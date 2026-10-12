"""Independent buffer, scalar/feature/ray/topology and decoded-oracle verification."""
from generalization_common import *
from compact_mesh import decode,source_check,ray_events,rational,exact_sign

def brute_distance(vertices,faces,points):
    xyz=vertices[faces];a,b,c=xyz[:,0],xyz[:,1],xyz[:,2];ab=b-a;ac=c-a;n=np.cross(ab,ac);nn=(n*n).sum(axis=1);result=[]
    for p in points:
        ap=p-a;height=(ap*n).sum(axis=1);projection=p-height[:,None]*n/nn[:,None]
        sides=np.stack([(np.cross(b-a,p-a)*n).sum(axis=1),(np.cross(c-b,p-b)*n).sum(axis=1),(np.cross(a-c,p-c)*n).sum(axis=1)])
        best=np.where(np.all(sides>=0,axis=0),height*height/nn,np.inf)
        for x,y in ((a,b),(b,c),(c,a)):
            e=y-x;t=((p-x)*e).sum(axis=1)/(e*e).sum(axis=1);delta=p-x-np.clip(t,0,1)[:,None]*e;best=np.minimum(best,(delta*delta).sum(axis=1))
        result.append(np.sqrt(best.min()))
    return np.array(result)

def verify(folder,ref,g,pilot=False):
    if (folder/'STATUS').exists() and (folder/'STATUS').read_text().splitlines()[0] in ('decoded_invalid','resource_unresolved','incomplete'):
        status=(folder/'STATUS').read_text().splitlines()[0]
        construction=json.loads((folder/'construction.json').read_text()) if (folder/'construction.json').exists() else {}
        layout_verified=False;displacement=None
        if (folder/'field.bin').exists():
            try:
                decoded=decode((folder/'field.bin').read_bytes());v,f=obj(g['normalized_mesh']);displacement=source_check(decoded,v,f);layout_verified=True
            except AssertionError:
                if status!='decoded_invalid':raise
        assessment=dict(geometry=g['id'],status=status,joint_pass=False,actual_bytes=construction.get('actual_bytes'),
          layout_independently_verified=layout_verified,maximum_vertex_displacement=displacement,
          native_evidence={k:v for k,v in manifest(folder).items() if k not in ('VERIFIED.json','SHA256SUMS.json')},evaluator_sha256=sha(Path(__file__)),pilot=pilot)
        json_write(folder/'VERIFIED.json',assessment);return assessment
    construction=json.loads((folder/'construction.json').read_text());data=(folder/'field.bin').read_bytes();decoded=decode(data);vertices,faces=decoded['vertices'],decoded['triangles'];source,source_faces=obj(g['normalized_mesh']);displacement=source_check(decoded,source,source_faces)
    assert construction['actual_bytes']==construction['object_bytes']+len(data)
    assert construction['buffer_bytes']==len(data) and construction['vertices']==len(vertices) and construction['faces']==len(faces) and construction['nodes']==len(decoded['nodes'])
    closed,euler=mesh_topology(vertices,faces);expected=[int(x) for x in (ref/'expected_euler.txt').read_text().split()];assert closed and euler==expected
    scalar=True;unresolved=False;metrics=[];independent_points=[];independent_native=[];probe_pass=True
    for name in ('uniform','near','probes','audit','boundary'):
        samples=np.fromfile(folder/'boundary.bin' if name=='boundary' else ref/(name+'.bin'),dtype='<f8').reshape(-1,4)
        pred=np.fromfile(folder/(name+'_pred.bin'),dtype='<f8');flags=np.fromfile(folder/(name+'_flags.bin'),dtype='u1');assert len(pred)==len(flags)
        if pilot and name in ('uniform','near'):samples=samples[:len(pred)]
        assert len(pred)==len(samples);bad=int(np.sum((flags&1)==0));nonfinite=int(np.sum(~np.isfinite(pred)));unresolved |= bad>0 or nonfinite>0
        if bad or nonfinite:record=dict(set=name,count=len(pred),unresolved=bad,nonfinite=nonfinite,maximum=None,cvar99=None);scalar=False
        else:
            error=np.abs(pred-samples[:,3]);tail=np.sort(error)[len(error)-max(1,int(np.ceil(.01*len(error)))):];maximum=float(error.max());cvar=float(tail.mean());record=dict(set=name,count=len(pred),unresolved=0,nonfinite=0,maximum=maximum,cvar99=cvar)
            scalar &= maximum<=.020 and (name!='uniform' or cvar<=.015) and (name!='near' or cvar<=.010)
            # Independent unpruned decoded-triangle distance on a stratified set,
            # complete probe witnesses, and every ray sign witness below.
            indices=np.arange(len(pred)) if name=='probes' else np.linspace(0,len(pred)-1,min(128,len(pred)),dtype=int)
            independent_points.extend(samples[indices,:3]);independent_native.extend(pred[indices])
        if name=='probes':
            groups=np.loadtxt(ref/'probe_groups.txt',dtype=int);recalls=[]
            for group in g['feature_groups']:
                for sign in (-1,1):
                    mask=(groups==group)&(samples[:,3]*sign>2e-9);assert mask.any()
                    recall=float(np.mean(pred[mask]*sign>0));probe_pass &= recall>=.95
                    recalls.append(dict(group=group,sign=sign,recall=recall))
        metrics.append(record)
    points=np.array(independent_points);native=np.array(independent_native);dist=brute_distance(vertices,faces,points)
    np.testing.assert_allclose(dist,np.abs(native),atol=3e-12,rtol=2e-8)
    tris=[tuple(rational(p) for p in vertices[t]) for t in faces];xyz=vertices[faces]
    signs=[]
    for p in points:signs.append(exact_sign(rational(p),tris,xyz))
    assert all(s is not None for s in signs),'Independent parity unresolved'
    signed=dist*np.array(signs);np.testing.assert_allclose(signed,native,atol=3e-12,rtol=2e-8)
    np.savez(folder/'independent_oracle.npz',points=points,native=native,independent=signed)
    rays=np.loadtxt(ref/'rays.txt');reference=[np.fromstring(r['roots'],sep=';') for r in rows(ref/'reference_roots.csv')];cpp=rows(folder/'roots.csv');assert len(cpp)==len(rays)
    ray_results=[];width=True;count_failures=0;sign_failures=0;numerical=0;errors={int(k):[] for k in np.unique(rays[:,0])};witnesses=rows(folder/'ray_witnesses.csv')
    for i,(ray,truth,record) in enumerate(zip(rays,reference,cpp)):
        result=ray_events(vertices,faces,ray);native_unresolved=bool(int(record['numerical_unresolved']))
        if not result['unresolved'] and not native_unresolved:
            np.testing.assert_array_equal(result['roots'],np.fromstring(record['roots'],sep=';'));np.testing.assert_array_equal(result['tangencies'],np.fromstring(record['tangencies'],sep=';'))
        numerical+=result['unresolved'] or native_unresolved;width &= not(result['unresolved'] or native_unresolved);roots=np.array(result['roots']) if not native_unresolved else np.fromstring(record['roots'],sep=';');correct_count=len(roots)==len(truth) and len(roots)>0 and len(roots)%2==0
        if not correct_count:count_failures+=1;width=False
        native_w=[r for r in witnesses if int(r['ray'])==i];times=[];expected_sign=[]
        for rr in (truth,roots):
            cuts=np.r_[ray[7],rr,ray[8]];times.extend((cuts[:-1]+cuts[1:])*.5);expected_sign.extend(np.where(np.arange(len(cuts)-1)%2,-1,1))
        times.extend(ray[[7,8]]);expected_sign.extend([1,1]);assert len(times)==len(native_w)
        sign_ok=True
        for t,sign,r in zip(times,expected_sign,native_w):
            assert float(r['t'])==t
            point=ray[1:4]+ray[4:7]*t;actual=exact_sign(rational(point),tris,xyz)
            native_sign=int(np.sign(float(r['prediction']))) if int(r['resolved']) and np.isfinite(float(r['prediction'])) else None
            if native_sign is not None:assert actual==native_sign,'Independent interval sign disagrees with standalone native query'
            sign_ok &= actual==sign and native_sign==sign
        sign_failures+=not sign_ok;width &= sign_ok
        if correct_count:errors[int(ray[0])].extend(abs(np.diff(roots)/np.diff(truth)-1))
        ray_results.append(dict(ray=i,roots=result['roots'],tangencies=result['tangencies'],numerical_unresolved=result['unresolved'],root_count_pass=bool(correct_count),signed_interval_pass=bool(sign_ok)))
    widths=[]
    for group,e in errors.items():
        if not e:width=False;continue
        mean=float(np.mean(e));p95=float(np.quantile(e,.95,method='inverted_cdf'));width &= mean<=.10 and p95<=.15;widths.append(dict(group=group,mean=mean,p95=p95))
    levels=[int(x) for x in (ref/'qualified_levels.txt').read_text().split()];topology_pass=len(levels)==2;topology=[]
    for level in levels:
        path=folder/f'surface_{level}.obj'
        if not path.exists():topology_pass=False;topology.append(dict(level=level,status='not_evaluated'));continue
        vv,ff=obj(path);closed,chi=mesh_topology(vv,ff);okay=closed and chi==expected;topology_pass &= okay;topology.append(dict(level=level,status='passed' if okay else 'observed_failure',closed=closed,component_euler=chi))
        nr=next(r for r in rows(folder/'topology.csv') if int(r['cells'])==level);assert closed==bool(int(nr['closed'])) and chi==[int(x) for x in nr['component_euler'].split(';')]
    joint=scalar and probe_pass and width and topology_pass and not unresolved
    status='reference_unresolved' if len(levels)!=2 else 'numerical_unresolved' if unresolved or numerical else 'passed' if joint else 'observed_gate_failure'
    assessment=dict(geometry=g['id'],codec=decoded['codec'],actual_bytes=construction['actual_bytes'],status=status,joint_pass=bool(joint),scalar_pass=bool(scalar),feature_pass=bool(probe_pass),width_pass=bool(width),topology_pass=bool(topology_pass),metrics=metrics,feature_recalls=recalls,width_groups=widths,incorrect_root_count_rays=count_failures,incorrect_signed_interval_rays=sign_failures,numerical_unresolved_rays=numerical,topology=topology,maximum_vertex_displacement=displacement,sections={k:decoded[k] for k in ('vertex_section_bytes','triangle_section_bytes','bvh_section_bytes','padding_bytes')},independent_decoded_queries=len(points),field_sha256=sha(folder/'field.bin'),source_sha256=sha(Path(g['normalized_mesh'])),reference_hashes=manifest(ref),evaluator_hashes={p.name:sha(p) for p in (Path(__file__),Path(__file__).with_name('compact_mesh.py'),Path(__file__).with_name('generalization_common.py'))},pilot=pilot)
    json_write(folder/'VERIFIED.json',assessment);json_write(folder/'VERIFY-RAYS.json',ray_results);return assessment
