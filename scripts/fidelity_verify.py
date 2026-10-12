"""Independent NumPy decoding, objectives, accounting and geometry verification."""
import argparse
import json
import math
import struct
from pathlib import Path
import numpy as np
from distortion_analysis import decode, oracle, layout, intervals, SHIFT
from tail_evidence import rows, hashes, sha, close, write_csv


def mesh_check(path, shape):
    vertices, faces = [], []
    for line in path.read_text().splitlines():
        if line.startswith('v '):
            vertices.append([float(v) for v in line.split()[1:]])
        elif line.startswith('f '):
            faces.append([int(v)-1 for v in line.split()[1:]])
    parent = list(range(len(vertices)))
    def root(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v
    edges, links = {}, {}
    for f in faces:
        assert len(f) == len(set(f)) == 3
        a,b,c = f
        assert np.linalg.norm(np.cross(np.array(vertices[b])-vertices[a],np.array(vertices[c])-vertices[a])) > 0
        for u,v in ((a,b),(b,c),(c,a)):
            key = tuple(sorted((u,v)))
            edges.setdefault(key, []).append(1 if u<v else -1)
            parent[root(u)] = root(v)
        for u,v,w in ((a,b,c),(b,c,a),(c,a,b)):
            links.setdefault(u, []).append((v,w))
    # A closed triangle manifold has oppositely oriented edge pairs and a
    # single cyclic link at every vertex; disconnected vertex fans fail.
    closed = all(len(v)==2 and sum(v)==0 for v in edges.values())
    for pairs in links.values():
        adjacency = {}
        for a,b in pairs:
            adjacency.setdefault(a, []).append(b)
            adjacency.setdefault(b, []).append(a)
        seen, todo = set(), [next(iter(adjacency))]
        while todo:
            v = todo.pop()
            if v not in seen:
                seen.add(v)
                todo.extend(adjacency[v])
        closed &= len(seen)==len(adjacency) and all(len(v)==2 for v in adjacency.values())
    components = {}
    for v in links:
        components.setdefault(root(v), [0,0,0])[0] += 1
    for a,b in edges:
        components[root(a)][1] += 1
    for f in faces:
        components[root(f[0])][2] += 1
    euler = sorted(v-e+f for v,e,f in components.values())
    return closed and euler == ([0] if shape=='torus' else [2,2])


def verify(directory, out, widths=True):
    assert not out.exists(), 'Use new analysis directory'
    assert (directory/'COMPLETE').exists() and (directory/'EXIT_CODE').read_text()=='0'
    hashes(directory)
    hashes(directory.with_name(directory.name+'-source'))
    out.mkdir(parents=True)
    config=json.loads((directory/'config.json').read_text())
    policy=json.loads((directory/'policy.json').read_text())
    metrics, allocations, costs = [rows(directory/f) for f in ('metrics.csv','allocation.csv','construction.csv')]
    fields=[r for r in metrics if r['dataset']=='probes']
    assert len(fields)==20
    expected={(s,f,str(seed),m) for s,fs in [('torus',[.035,.06,.1]),('double_sphere',[.06,.1])] for f in fs for seed in [policy['seed_start'],policy['seed_start']+1] for m in ['baseline_q8','risk_q8']}
    assert {(r['shape'],float(r['feature']),r['seed'],r['method']) for r in fields}==expected
    summaries, width_rows = [], []
    verified=0
    for r in fields:
        case,method,shape=r['case'],r['method'],r['shape']
        feature=float(r['feature'])
        name=f'{case}_6_{method}'
        data=(directory/'fields'/f'{name}.bin').read_bytes()
        layout(data,method,6)
        assert len(data)+config['field_object_bytes']==int(r['bytes'])<=32768
        a=next(a for a in allocations if a['case']==case and a['method']==method)
        for ds in ('fit','uniform','near','probes'):
            sample=np.loadtxt(directory/'samples'/f'{case}_{ds}.csv',delimiter=',',skiprows=1,ndmin=2)
            pred=decode(data,sample[:,:3])
            cpp=np.fromfile(directory/'predictions'/f'{name}_{ds}.bin',dtype='<f8')
            np.testing.assert_allclose(pred,cpp,rtol=2e-12,atol=2e-14)
            np.testing.assert_allclose(sample[:,3],oracle(sample[:,:3],shape,feature),rtol=2e-12,atol=2e-14)
            m=next(m for m in metrics if m['case']==case and m['method']==method and m['dataset']==ds)
            error=np.abs(pred-sample[:,3])
            assert len(sample)==int(m['count'])
            close(error.mean(),m['mae'],'MAE')
            close(np.sort(error)[-math.ceil(.01*len(error)):].mean(),m['cvar99'],'CVaR99')
            for sign,key in ((-1,'inside_recall'),(1,'outside_recall')):
                mask=sample[:,3]*sign>2e-9
                close(np.mean(pred[mask]*sign>0),m[key],key)
            valid=np.abs(sample[:,3])>2e-9
            wrong=(pred*np.sign(sample[:,3])<=0)&valid
            close(wrong.sum()/valid.sum(),m['sign_error'],'sign_error')
            if ds=='fit':
                sse=np.sum(((pred-sample[:,3])/2)**2)
                close(sse,a['fit_sse'],'fitting SSE')
                objective=sse+(policy['lambda'] if method=='risk_q8' else 0)*wrong.sum()
                if method=='risk_q8' and policy.get('quadrature',0):
                    extra=np.loadtxt(directory/'samples'/f'{case}_extra.csv',delimiter=',',skiprows=1,ndmin=2)
                    np.testing.assert_allclose(extra[:,3],oracle(extra[:,:3],shape,feature),rtol=2e-12,atol=2e-14)
                    ep=decode(data,extra[:,:3])
                    ew=(ep*np.sign(extra[:,3])<=0)&(np.abs(extra[:,3])>2e-9)
                    objective+=np.sum(((ep-extra[:,3])/2)**2)+policy['lambda']*ew.sum()
                close(objective,a['objective'],'allocated objective')
            verified+=len(sample)
        cost=next(c for c in costs if c['case']==case and c['method']==method)
        k=policy['k'] if method=='risk_q8' else 32
        assert int(a['candidate_nodes'])==37449+8*int(a['selected_parents'])
        assert int(a['selected_parents'])<=k
        extra_count=policy.get('quadrature',0)**3*int(a['selected_parents']) if method=='risk_q8' else 0
        assert int(cost['oracle_calls'])==config['fit']+8*int(a['candidate_nodes'])+extra_count
        close(sum(float(cost[x]) for x in ('label_seconds','candidate_seconds','refinement_seconds','allocation_seconds','pack_seconds')),cost['total_seconds'],'construction sum')
        summaries.append(dict(case=case,shape=shape,feature=feature,seed=r['seed'],method=method,bytes=r['bytes'],inside_recall=r['inside_recall'],outside_recall=r['outside_recall'],retained=a['retained_parents'],oracle_calls=cost['oracle_calls'],build_ms=1000*float(cost['total_seconds'])))
        if widths and (not policy['quick'] or (shape=='torus' and feature==.035) or (shape=='double_sphere' and feature==.06)):
            if shape=='torus':
                angles=np.arange(256)*2*np.pi/256
                radial=np.column_stack([np.cos(angles),np.sin(angles),np.zeros(256)])
                centers=SHIFT+.55*radial
                axes=[('radial',radial),('vertical',np.tile([0.,0.,1.],(256,1)))]
            else:
                angles=np.array([0.]); centers=SHIFT[None,:]; axes=[('gap',np.array([[1.,0.,0.]]))]
            truth=2*feature if shape=='torus' else feature
            for axis,directions in axes:
                left,right=intervals(data,centers,directions,feature,shape!='torus')
                ol,orr=intervals(lambda p:oracle(p,shape,feature),centers,directions,feature,shape!='torus')
                np.testing.assert_allclose(orr-ol,truth,atol=1e-12)
                for angle,l,h in zip(angles,left,right):
                    width_rows.append(dict(case=case,shape=shape,feature=feature,method=method,axis=axis,angle=angle,left=l,right=h,width=h-l,expected=truth,relative_error=(h-l)/truth-1))
        print(name,'independently verified',flush=True)
    # Verify residual ranks and retained branches independently for both policies.
    from fidelity_verify_selection import verify_selection
    verify_selection(directory)
    if not policy['quick']:
        surfaces=rows(directory/'surfaces.csv')
        assert len(surfaces)==90
        for r in surfaces:
            if r['method']=='oracle':
                probe_pass=True
            else:
                m=next(m for m in fields if m['case']==r['case'] and m['method']==r['method'])
                probe_pass=min(float(m['inside_recall']),float(m['outside_recall']))>=.95
            assert bool(int(r['probe_pass']))==probe_pass
            mesh=directory/'meshes'/f"{r['case']}_{r['depth']}_{r['method']}_{r['cells']}.obj"
            passed=mesh_check(mesh,r['shape'])
            assert passed==bool(int(r['topology_pass']))
            assert bool(int(r['success']))==(passed and probe_pass)
        recovery=[]
        for shape in ('torus','double_sphere'):
            for method in ('baseline_q8','risk_q8'):
                features=sorted({float(r['feature']) for r in fields if r['shape']==shape})
                passes={}
                for f in features:
                    checks=[r for r in surfaces if r['shape']==shape and float(r['feature'])==f and r['method'] in ('oracle',method) and int(r['cells']) in (64,96)]
                    assert len(checks)==8
                    passes[f]=all(int(r['success']) for r in checks)
                for f in features:
                    recovery.append(dict(shape=shape,method=method,feature=f,individual_pass=passes[f],recovered=all(passes[g] for g in features if g>=f)))
        write_csv(out/'recovery.csv',recovery)
    write_csv(out/'summary.csv',summaries)
    if width_rows:
        write_csv(out/'widths.csv',width_rows)
    verifier_files=[Path(__file__),Path(__file__).with_name('fidelity_verify_selection.py'),Path(__file__).with_name('distortion_analysis.py'),Path(__file__).with_name('tail_evidence.py')]
    (out/'VERIFIED.json').write_text(json.dumps(dict(fields=len(fields),queries=verified,raw_manifest=sha(directory/'SHA256SUMS.json'),source_manifest=sha(directory.with_name(directory.name+'-source')/'SHA256SUMS.json'),verifier_sha256={f.name:sha(f) for f in verifier_files}),indent=2))
    (out/'SHA256SUMS.json').write_text(json.dumps({f.name:sha(f) for f in out.iterdir() if f.is_file()},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('directory',type=Path)
    p.add_argument('out',type=Path)
    p.add_argument('--no-widths',action='store_true')
    a=p.parse_args()
    verify(a.directory,a.out,not a.no_widths)
