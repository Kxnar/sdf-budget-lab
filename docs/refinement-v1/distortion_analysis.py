"""Independent packed verification, error decomposition, widths and fit coverage."""
import argparse
import json
import math
from pathlib import Path
import struct
import numpy as np
from tail_evidence import rows, write_csv, hashes, sha, close

SHIFT = np.array([.017, -.023, .011])


def decode(data, points, metadata=False):
    points = np.asarray(points, dtype=float)
    raw = np.frombuffer(data, dtype=np.uint8)
    words = np.frombuffer(data, dtype='<u4')
    floats = np.frombuffer(data, dtype='<f4')
    lo, hi = np.full_like(points, -1.), np.full_like(points, 1.)
    node = np.zeros(len(points), dtype=np.int64)
    depth = np.zeros(len(points), dtype=np.int64)
    for _ in range(7):
        desc = words[node // 4]
        active = (desc & 3) == 0
        if not active.any():
            break
        mid = (lo[active] + hi[active]) / 2
        side = points[active] >= mid
        node[active] = (desc[active] & ~np.uint32(3)).astype(np.int64) + 4 * (side @ np.array([1, 2, 4]))
        lo[active] = np.where(side, mid, lo[active])
        hi[active] = np.where(side, hi[active], mid)
        depth[active] += 1
    assert not active.any(), 'Unexpected depth/cycle'
    desc = words[node // 4]
    enc = desc & 3
    assert np.all((enc == 1) | (enc == 3))
    offset = (desc & ~np.uint32(3)).astype(np.int64)
    t = (points - lo) / (hi - lo)
    result = np.zeros(len(points))
    for c in range(8):
        value = np.zeros(len(points))
        f = enc == 1
        value[f] = floats[offset[f] // 4 + c].astype(float)
        q = enc == 3
        value[q] = (floats[offset[q] // 4].astype(float) + floats[offset[q] // 4 + 1].astype(float) * raw[offset[q] + 8 + c]).astype(np.float32).astype(float)
        weight = np.ones(len(points))
        for axis in range(3):
            weight *= t[:, axis] if c & (1 << axis) else 1 - t[:, axis]
        result += value * weight
    return (result, node, depth) if metadata else result


def oracle(points, shape, feature):
    p = np.asarray(points) - SHIFT
    if shape == 'torus':
        return np.hypot(np.hypot(p[:, 0], p[:, 1]) - .55, p[:, 2]) - feature
    c = np.array([.3 + feature/2, 0, 0])
    return np.minimum(np.linalg.norm(p-c, axis=1), np.linalg.norm(p+c, axis=1)) - .3


def double_interpolation(points, shape, feature, depth):
    h = 2 / 2**depth
    lo = -1 + np.minimum(np.floor((points + 1) / h), 2**depth - 1) * h
    t = (points-lo) / h
    result = np.zeros(len(points))
    for c in range(8):
        side = np.array([bool(c & (1 << i)) for i in range(3)])
        result += oracle(lo + side*h, shape, feature) * np.prod(np.where(side, t, 1-t), axis=1)
    return result


def layout(data, method, depth):
    used = np.zeros(len(data)//4, dtype=bool)
    counts = [0, 0]
    def claim(pos, n):
        assert pos % 4 == 0 and 0 <= pos <= len(data)-n
        assert not used[pos//4:(pos+n)//4].any()
        used[pos//4:(pos+n)//4] = True
    def visit(node, d):
        claim(node, 4)
        desc, = struct.unpack_from('<I', data, node)
        e, pos = desc & 3, desc & ~3
        if e:
            assert e == (1 if method == 'full_f32' else 3)
            claim(pos, 32 if e == 1 else 16)
            assert d <= depth
            if method.startswith('full'):
                assert d == depth
            counts[1] += 1
        else:
            assert d < depth
            counts[0] += 1
            for c in range(8):
                visit(pos+4*c, d+1)
    visit(0, 0)
    assert used.all()
    if method.startswith('full'):
        assert counts == [(8**depth-1)//7, 8**depth]


def intervals(data, centers, directions, feature, positive):
    """First sign transition on each side of the centre; zero if centre fails."""
    query = data if callable(data) else lambda points: decode(data, points)
    n = len(centers)
    boundaries = []
    valid_center = query(centers) * (1 if positive else -1) > 0
    for sign in (-1, 1):
        found = np.zeros(n, dtype=bool)
        lower, upper = np.zeros(n), np.full(n, np.nan)
        # Chunk the dense scan to avoid a multi-million-query temporary array.
        for start in range(0, n, 16):
            end = min(start+16, n)
            offsets = np.linspace(0, 2*feature, 4097)
            points = centers[start:end, None, :] + sign * offsets[None, :, None] * directions[start:end, None, :]
            inside = (query(points.reshape(-1, 3)).reshape(end-start, -1) * (1 if positive else -1)) > 0
            bad = ~inside
            has = bad.any(axis=1)
            idx = np.argmax(bad, axis=1)
            found[start:end] = has
            lower[start:end] = offsets[np.maximum(idx-1, 0)]
            upper[start:end] = np.where(has, offsets[idx], np.nan)
        good = valid_center & found
        for _ in range(40):
            mid = (lower[good]+upper[good])/2
            inside = query(centers[good]+sign*mid[:, None]*directions[good]) * (1 if positive else -1) > 0
            lower[good] = np.where(inside, mid, lower[good])
            upper[good] = np.where(inside, upper[good], mid)
        boundaries.append(np.where(valid_center, sign*(lower+upper)/2, 0))
    return boundaries[0], boundaries[1]


def analyze(directory, out):
    if out.exists():
        raise ValueError('Analysis output must be new')
    assert (directory/'COMPLETE').exists() and (directory/'EXIT_CODE').read_text() == '0'
    hashes(directory)
    hashes(directory.with_name(directory.name+'-source'))
    out.mkdir(parents=True)
    config = json.loads((directory/'config.json').read_text())
    metrics = rows(directory/'metrics.csv')
    fields = [r for r in metrics if r['dataset'] == 'probes']
    assert len(fields) == (10 if config['smoke'] else 50)
    summaries, coverage, profiles, widths, decompositions = [], [], [], [], []
    sample_cache = {}
    geometry_cache = {}
    verified = 0
    for row in fields:
        case, shape = row['case'], row['shape']
        feature, depth = float(row['feature']), int(row['depth'])
        method = row['method']
        name = f'{case}_{depth}_{method}'
        if case not in sample_cache:
            sample_cache.clear()
            sample_cache[case] = {s: np.loadtxt(directory/'samples'/f'{case}_{s}.csv', delimiter=',', skiprows=1, ndmin=2) for s in ('fit', 'uniform', 'near', 'probes')}
        samples = sample_cache[case]
        data = (directory/'fields'/f'{name}.bin').read_bytes()
        layout(data, method, depth)
        assert len(data)+config['field_object_bytes'] == int(row['bytes'])
        if method == 'mse_q8':
            assert int(row['bytes']) <= 32768
        for ds in ('uniform', 'near', 'probes'):
            sample = samples[ds]
            pred = decode(data, sample[:, :3])
            cpp = np.fromfile(directory/'predictions'/f'{name}_{ds}.bin', dtype='<f8')
            np.testing.assert_allclose(pred, cpp, rtol=2e-12, atol=2e-14)
            np.testing.assert_allclose(sample[:, 3], oracle(sample[:, :3], shape, feature), rtol=2e-12, atol=2e-14)
            errors = np.abs(pred-sample[:, 3])
            m = next(r for r in metrics if r['case'] == case and int(r['depth']) == depth and r['method'] == method and r['dataset'] == ds)
            assert len(sample) == int(m['count'])
            close(errors.mean(), m['mae'], 'MAE')
            close(np.sort(errors)[-math.ceil(.01*len(errors)):].mean(), m['cvar99'], 'CVaR99')
            for sign, key in ((-1, 'inside_recall'), (1, 'outside_recall')):
                mask = sample[:, 3]*sign > 2e-9
                close(np.mean(pred[mask]*sign > 0), m[key], key)
            valid = np.abs(sample[:, 3]) > 2e-9
            close(np.mean(pred[valid]*np.sign(sample[valid, 3]) <= 0), m['sign_error'], 'sign error')
            verified += len(pred)
        # Exact frozen probes, routed to the selected leaf, compared with fitting pool.
        probes, fit = samples['probes'], samples['fit']
        pred, nodes, depths = decode(data, probes[:, :3], True)
        _, fit_nodes, _ = decode(data, fit[:, :3], True)
        failed = pred*np.sign(probes[:, 3]) <= 0
        counts, same_counts = [], []
        for i in range(len(probes)):
            same_leaf = fit_nodes == nodes[i]
            count = int(same_leaf.sum())
            correct_sign = int(np.sum(same_leaf & (fit[:, 3]*np.sign(probes[i, 3]) > 2e-9)))
            counts.append(count); same_counts.append(correct_sign)
            coverage.append(dict(case=case, depth=depth, method=method, probe=i, x=probes[i, 0], y=probes[i, 1], z=probes[i, 2], truth=probes[i, 3], prediction=pred[i], failed=int(failed[i]), leaf_depth=depths[i], fit_count=count, same_sign_fit_count=correct_sign))
        summaries.append(dict(case=case, shape=shape, feature=feature, seed=row['seed'], depth=depth, method=method, failed_probes=int(failed.sum()), failing_with_fit=int(np.sum(failed & (np.array(counts)>0))), failing_with_same_sign_fit=int(np.sum(failed & (np.array(same_counts)>0))), failing_at_max_depth=int(np.sum(failed & (depths==depth)))))
        # Full fields are seed-independent; retain profiles/widths once for them.
        gkey = (shape, feature, depth, method)
        if method.startswith('full') and gkey in geometry_cache:
            assert sha(directory/'fields'/f'{name}.bin') == geometry_cache[gkey]
            continue
        geometry_cache[gkey] = sha(directory/'fields'/f'{name}.bin')
        if shape == 'torus':
            u = np.arange(1024)*2*np.pi/1024
            points = SHIFT + np.column_stack([.55*np.cos(u), .55*np.sin(u), np.zeros_like(u)])
        else:
            u = np.linspace(-feature, feature, 1025)
            points = SHIFT + np.column_stack([u, np.zeros_like(u), np.zeros_like(u)])
        truth, pred = oracle(points, shape, feature), decode(data, points)
        if method.startswith('full'):
            f32data = (directory/'fields'/f'{case}_{depth}_full_f32.bin').read_bytes()
            f32 = decode(f32data, points)
            exact = double_interpolation(points, shape, feature, depth)
            quant, rounding, interpolation = pred-f32, f32-exact, exact-truth
            np.testing.assert_allclose(pred-truth, quant+rounding+interpolation, rtol=1e-12, atol=1e-14)
            decompositions.append(dict(case=case, shape=shape, feature=feature, depth=depth, method=method, max_abs_quantization=float(np.max(np.abs(quant))), max_abs_f32_rounding=float(np.max(np.abs(rounding))), max_abs_interpolation=float(np.max(np.abs(interpolation))), mean_signed_error=float(np.mean(pred-truth))))
        else:
            quant = rounding = interpolation = np.full(len(points), np.nan)
        for i in range(len(points)):
            profiles.append(dict(case=case, shape=shape, feature=feature, depth=depth, method=method, coordinate=u[i], truth=truth[i], prediction=pred[i], signed_error=pred[i]-truth[i], quantization=quant[i], f32_rounding=rounding[i], interpolation=interpolation[i]))
        if shape == 'torus':
            a = np.arange(256)*2*np.pi/256
            radial = np.column_stack([np.cos(a), np.sin(a), np.zeros_like(a)])
            centers = SHIFT + .55*radial
            axes = [('radial', radial), ('vertical', np.tile([0., 0., 1.], (256, 1)))]
        else:
            a, centers = np.array([0.]), SHIFT[None, :]
            axes = [('gap', np.array([[1., 0., 0.]]))]
        for axis, directions in axes:
            if depth == 5 and method == 'full_q8':
                ol, ore = intervals(lambda p: oracle(p, shape, feature), centers, directions, feature, shape != 'torus')
                np.testing.assert_allclose(ore-ol, feature if shape != 'torus' else 2*feature, atol=1e-12)
            left, right = intervals(data, centers, directions, feature, shape != 'torus')
            for i in range(len(a)):
                widths.append(dict(case=case, shape=shape, feature=feature, depth=depth, method=method, axis=axis, angle=a[i], left=left[i], right=right[i], width=right[i]-left[i], expected=feature if shape != 'torus' else 2*feature))
        print(name, 'verified and localized', flush=True)
    for filename, values in [('coverage.csv', coverage), ('coverage_summary.csv', summaries), ('profiles.csv', profiles), ('widths.csv', widths), ('decomposition.csv', decompositions)]:
        write_csv(out/filename, values)
    # Frozen recovery across both seeds, final resolutions, oracle and larger features.
    surface = rows(directory/'surfaces.csv')
    recovery = []
    for shape in ('torus', 'double_sphere'):
        features = sorted({float(r['feature']) for r in fields if r['shape'] == shape})
        for method, depth in [('full_q8', 5), ('full_q8', 6), ('full_f32', 5), ('full_f32', 6), ('mse_q8', 5)]:
            passes = {}
            for f in features:
                selected = [r for r in surface if r['shape'] == shape and float(r['feature']) == f and int(r['cells']) in ((16, 24) if config['smoke'] else (64, 96)) and (r['method'] == 'oracle' or (r['method'] == method and int(r['depth']) == depth))]
                assert len(selected) == (4 if config['smoke'] else 8)
                passes[f] = all(int(r['success']) for r in selected)
            threshold = next((f for f in features if all(passes[g] for g in features if g>=f)), None)
            recovery.append(dict(shape=shape, method=method, depth=depth, minimum_tested_recovered_feature=threshold))
    write_csv(out/'recovery.csv', recovery)
    (out/'VERIFIED.json').write_text(json.dumps(dict(fields=len(fields), cpp_queries=verified, raw_manifest_sha256=sha(directory/'SHA256SUMS.json'), source_manifest_sha256=sha(directory.with_name(directory.name+'-source')/'SHA256SUMS.json'), verifier_sha256=sha(Path(__file__))), indent=2))
    print(f'Verified {len(fields)} fields and {verified} C++ queries', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('directory', type=Path)
    p.add_argument('out', type=Path)
    a = p.parse_args()
    analyze(a.directory, a.out)
