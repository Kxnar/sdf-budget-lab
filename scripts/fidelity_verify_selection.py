"""Reconstruct d5 residual ranks independently from saved fitting data."""
import json
import struct
import numpy as np
from distortion_analysis import oracle
from tail_evidence import rows, close


def verify_selection(directory):
    config = json.loads((directory/'config.json').read_text())
    policy = json.loads((directory/'policy.json').read_text())
    selections = rows(directory/'selection.csv')
    allocations = rows(directory/'allocation.csv')
    costs = rows(directory/'construction.csv')
    # Reproduce only spatial node numbering, not the C++ residual implementation.
    ids = np.empty((32, 32, 32), dtype=np.int64)
    next_node = 1
    def visit(node, depth, xyz):
        nonlocal next_node
        if depth == 5:
            ids[tuple(xyz)] = node
            return
        first = next_node
        next_node += 8
        for c in range(8):
            side = np.array([bool(c & (1 << a)) for a in range(3)])
            visit(first+c, depth+1, xyz*2+side)
    visit(0, 0, np.zeros(3, dtype=np.int64))
    assert next_node == 37449
    for allocation in allocations:
        case, method = allocation['case'], allocation['method']
        k = policy['k'] if method == 'risk_q8' else 32
        shape, feature_text, _ = case.rsplit('_', 2)
        feature = float(feature_text)
        fit = np.loadtxt(directory/'samples'/f'{case}_fit.csv', delimiter=',', skiprows=1)
        coords = np.minimum(((fit[:, :3]+1)*16).astype(int), 31)
        lo = -1+coords/16
        t = (fit[:, :3]-lo)*16
        values = np.column_stack([oracle(lo+np.array([bool(c & (1 << a)) for a in range(3)])/16, shape, feature) for c in range(8)]).astype(np.float32)
        offset = values.min(axis=1)
        scale = ((values.max(axis=1).astype(float)-offset)/255).astype(np.float32)
        assert np.all(scale > 0)
        quant = np.clip(np.floor((values.astype(float)-offset[:, None])/scale[:, None]+.5), 0, 255)
        decoded = (offset[:, None].astype(float)+scale[:, None].astype(float)*quant).astype(np.float32)
        prediction = np.zeros(len(fit))
        for c in range(8):
            side = np.array([bool(c & (1 << a)) for a in range(3)])
            prediction += decoded[:, c]*np.prod(np.where(side, t, 1-t), axis=1)
        nodes = ids[coords[:, 0], coords[:, 1], coords[:, 2]]
        losses = np.bincount(nodes, weights=((prediction-fit[:, 3])/2)**2, minlength=37449)
        counts = np.bincount(nodes, minlength=37449)
        eligible = np.flatnonzero((counts > 0) & (losses > 0))
        ranked = sorted(eligible, key=lambda node: (-losses[node], node))[:k]
        chosen = sorted([r for r in selections if r['case']==case and r['method']==method], key=lambda r: int(r['rank']))
        assert [int(r['node']) for r in chosen] == ranked, case
        data = (directory/'fields'/f'{case}_6_{method}.bin').read_bytes()
        retained = 0
        for r in chosen:
            node = int(r['node'])
            assert int(r['fit_count']) == counts[node]
            close(losses[node], r['q8_loss'], 'ranked residual SSE')
            xyz = np.argwhere(ids == node)[0]
            cell_lo = -1+xyz/16
            np.testing.assert_array_equal(cell_lo, [float(r['lo_'+a]) for a in 'xyz'])
            np.testing.assert_array_equal(cell_lo+1/16, [float(r['hi_'+a]) for a in 'xyz'])
            # A selected split is retained iff its centre reaches a branch at d5.
            centre = cell_lo+1/32
            pos, lower, upper, depth = 0, np.full(3, -1.), np.full(3, 1.), 0
            while depth < 5:
                desc, = struct.unpack_from('<I', data, pos)
                if desc & 3:
                    break
                mid = (lower+upper)/2
                side = centre >= mid
                pos = (desc & ~3)+4*int(side @ np.array([1, 2, 4]))
                lower, upper = np.where(side, mid, lower), np.where(side, upper, mid)
                depth += 1
            desc, = struct.unpack_from('<I', data, pos)
            keep = depth == 5 and not desc & 3
            assert keep == bool(int(r['retained']))
            retained += keep
        assert int(allocation['retained_parents']) == retained
        assert int(allocation['selected_parents']) == len(chosen)
        assert int(allocation['candidate_nodes']) == 37449+8*len(chosen)
        cost = next(r for r in costs if r['case']==case and r['method']==method)
        q=policy.get('quadrature',0) if method=='risk_q8' else 0
        assert int(cost['oracle_calls']) == config['fit']+8*(37449+8*len(chosen))+q**3*len(chosen)
        if q:
            expected=[]
            for r in chosen:
                lo=np.array([float(r['lo_'+a]) for a in 'xyz'])
                for z in range(q):
                    for y in range(q):
                        for x in range(q):
                            expected.append(lo+(np.array([x,y,z])+.5)/q/16)
            extra=np.loadtxt(directory/'samples'/f'{case}_extra.csv',delimiter=',',skiprows=1)
            np.testing.assert_array_equal(extra[:,:3],expected)
    for cost in costs:
        total = sum(float(cost[key]) for key in ('label_seconds', 'candidate_seconds', 'refinement_seconds', 'allocation_seconds', 'pack_seconds'))
        close(total, cost['total_seconds'], 'construction sum')
    print('Independent residual rankings, selection, retention and costs verified', flush=True)
