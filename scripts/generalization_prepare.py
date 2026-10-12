"""Reference-only input qualification, development pilot and immutable freeze."""
import argparse, subprocess, urllib.request, platform, time, shutil
from generalization_common import *

REV='835213277571d625722d25880c9c3010c6f45a5a'
def mesh_intersections(v,f,o,d,lo,hi):
    a,b,c=v[f[:,0]],v[f[:,1]],v[f[:,2]];e1=b-a;e2=c-a;h=np.cross(np.broadcast_to(d,e2.shape),e2);det=np.einsum('ij,ij->i',e1,h);good=np.abs(det)>1e-12;inv=np.zeros(len(det));inv[good]=1/det[good];s=o-a;u=inv*np.einsum('ij,ij->i',s,h);q=np.cross(s,e1);w=inv*(q@d);t=inv*np.einsum('ij,ij->i',e2,q);hit=good&(u>=-1e-10)&(w>=-1e-10)&(u+w<=1+1e-10)&(t>lo)&(t<hi)
    ts=np.sort(t[hit]);return ts[np.r_[True,np.diff(ts)>1e-8]] if len(ts) else ts
def ray_bounds(o,d):
    mask=np.abs(d)>1e-12;t=np.stack([(-.999-o[mask])/d[mask],(.999-o[mask])/d[mask]]);return max(t.min(axis=0)),min(t.max(axis=0))
def analytic_rays(g):
    s=g['shape'];a,b,c=g['parameters'];local=[]
    def add(group,o,d,extent=.85):local.append([group,*o,*d,-extent,extent])
    for i in range(48):
        angle=2*np.pi*(i+.37)/48;rad=np.array([np.cos(angle),np.sin(angle),0]);q=(i+.5)/48
        if s=='thin_box':
            x=-.48+.96*((i*17)%48+.5)/48;y=-.48+.96*q;add(0,[x,y,0],[0,0,1],.4);add(1,[0,y,0],[1,0,0])
        elif s=='torus':add(0,.55*rad,rad,.20);add(1,.55*rad,[0,0,1],.20)
        elif s=='double_sphere':
            y=.08*np.cos(angle);z=.08*np.sin(angle);add(0,[0,y,z],[1,0,0]);add(1,[.3+a/2,0,0],[0,np.cos(angle),np.sin(angle)],.45)
        elif s=='capsule':add(0,[0,0,-b+2*b*q],rad,.4);add(1,[0,0,0],np.array([.3*np.cos(angle),.3*np.sin(angle),1])/np.sqrt(1.09))
        elif s=='cylinder':add(0,[0,0,-.8*b+1.6*b*q],rad,.6);add(1,[.65*a*np.cos(angle),.65*a*np.sin(angle),0],[0,0,1],.6)
        elif s=='rounded_box':add(0,[.7*b*np.cos(angle),.7*c*np.sin(angle),0],[0,0,1],.4);add(1,[0,0,0],np.array([b*np.cos(angle),c*np.sin(angle),.12])/np.linalg.norm([b*np.cos(angle),c*np.sin(angle),.12]))
        else:
            d=np.array([np.cos(angle)*np.sqrt(1-(2*q-1)**2),np.sin(angle)*np.sqrt(1-(2*q-1)**2),2*q-1]);add(0,[0,0,0],d);add(1,[0,0,0],np.roll(d,1))
    result=[]
    for group,x,y,z,dx,dy,dz,lo,hi in local:
        o=forward(g,[x,y,z]);d=np.array([dx,dy,dz])@np.array(g['rotation']).T;lower,upper=ray_bounds(o,d);result.append([group,*o,*d,max(lo,lower),min(hi,upper)])
    return np.array(result)
def mesh_rays(g,v,f):
    rng=np.random.default_rng(713007);candidates=[];centroids=v[f].mean(axis=1);_,_,basis=np.linalg.svd(v-v.mean(axis=0),full_matrices=False)
    for i in range(768):
        o=centroids[int(rng.integers(len(f)))];d=basis[i%3] if i<384 else rng.normal(size=3);d=d/np.linalg.norm(d);lo,hi=ray_bounds(o,d);ts=mesh_intersections(v,f,o,d,lo,hi)
        if len(ts)>=2 and len(ts)%2==0 and np.min(np.diff(ts))>1e-5:
            candidates.append((np.min(np.diff(ts)),[0,*o,*d,lo,hi],ts))
    assert len(candidates)>=96,'Insufficient unambiguous intersecting rays'
    candidates.sort(key=lambda x:x[0]);selected=candidates[:48]+[candidates[i] for i in np.linspace(48,len(candidates)-1,48,dtype=int)]
    rr=[];tt=[]
    for i,(_,r,t) in enumerate(selected):r[0]=0 if i<48 else 1 if i<72 else 2;rr.append(r);tt.append(t)
    return np.array(rr),tt
def write_spec(g,p):
    a,b,c=g['parameters'];mesh=g.get('normalized_mesh','');v=[g['id'],g['shape'],a,b,c,g['noise'],*np.array(g['rotation']).ravel(),*g['translation']];p.write_text(' '.join(str(x) for x in v)+' '+json.dumps(mesh)+'\n')
def create_reference(g,folder,exe,count=100000,seed=880001,levels=(64,96,128,192)):
    folder.mkdir(parents=True,exist_ok=False);spec=folder/'geometry.txt';write_spec(g,spec)
    if g['shape']=='mesh':
        v,f=obj(g['normalized_mesh']);ray,rr=mesh_rays(g,v,f)
    else:
        ray=analytic_rays(g);rr=[r for r,c in independent_roots(lambda p:exact(g,p),ray)]
    assert all(len(t)>=2 and len(t)%2==0 for t in rr)
    if g['shape']=='mesh':csv_write(folder/'expected_ray_roots.csv',[dict(ray=i,roots=';'.join(format(float(x),'.17g') for x in t)) for i,t in enumerate(rr)])
    np.savetxt(folder/'rays.txt',ray,fmt='%.17g');probes=[];groups=[]
    for r,t in zip(ray,rr):
        o=r[1:4];d=r[4:7];values=np.r_[r[7],t,r[8]];centers=(values[:-1]+values[1:])/2
        for center in centers:probes.append(o+center*d);groups.append(int(r[0]))
    np.savetxt(folder/'probes.points',probes,fmt='%.17g');np.savetxt(folder/'probe_groups.txt',groups,fmt='%d')
    # Fixed lattice plus all potential depth-six interface planes; same set for all methods.
    grid=np.stack(np.meshgrid(*([np.linspace(-.999,.999,25)]*3),indexing='ij'),axis=-1).reshape(-1,3);audit=[grid]
    for axis in range(3):
        p=np.stack(np.meshgrid(np.arange(1,64)/32-1,np.linspace(-.85,.85,7),np.linspace(-.85,.85,7),indexing='ij'),axis=-1).reshape(-1,3)
        p=np.roll(p,axis,axis=1)
        for e in (-1e-8,1e-8):q=p.copy();q[:,axis]+=e;audit.append(q)
    np.savetxt(folder/'audit.points',np.concatenate(audit),fmt='%.17g')
    (folder/'levels.txt').write_text(' '.join(map(str,levels)));(folder/'expected_euler.txt').write_text(' '.join(map(str,g['expected_euler'])))
    start=time.time();result=subprocess.run([str(exe),'qualify',str(spec),str(folder),str(count),str(seed)],stdout=(folder/'qualification.log').open('w'),stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError((folder/'qualification.log').read_text())
    ref=rows(folder/'reference_roots.csv')
    for i,t in enumerate(rr):
        native=np.fromstring(ref[i]['roots'],sep=';');np.testing.assert_allclose(native,t,atol=2e-8,rtol=1e-10);assert int(ref[i]['continuous'])==1
    samples=np.fromfile(folder/'probes.bin',dtype='<f8').reshape(-1,4)
    assert np.all(np.abs(samples[:,3])>1e-8),'Witness on surface'
    for group in sorted(set(groups)):
        vals=samples[np.array(groups)==group,3];assert np.any(vals<0) and np.any(vals>0),'Each group needs both signs'
    top=rows(folder/'reference_topology.csv');valid=[bool(int(r['closed'])) and list(map(int,r['component_euler'].split(';')))==g['expected_euler'] for r in top];qualified=[]
    for i in range(1,len(valid)):
        if valid[i-1] and valid[i]:qualified=[int(top[i-1]['cells']),int(top[i]['cells'])];break
    for r in top:
        closed,euler=mesh_topology(*obj(folder/f"oracle_0_{r['cells']}.obj"));assert closed==bool(int(r['closed']));assert euler==list(map(int,r['component_euler'].split(';')))
    (folder/'qualified_levels.txt').write_text(' '.join(map(str,qualified)))
    g['qualified_levels']=qualified;g['reference_status']='qualified' if qualified else 'oracle_resolution_unresolved';g['feature_groups']=sorted(set(groups));g['ray_count']=len(ray);g['qualification_seconds']=time.time()-start
    json_write(folder/'independent_reference.json',{'ray_reference':'actual triangle intersections' if g['shape']=='mesh' else 'independent analytic SDF roots','ray_count':len(ray),'reference_topology_verified':True})
    return g
def source_files():
    files=[ROOT/'CMakeLists.txt',ROOT/'docs/generalization-protocol.md']
    for folder in ('src','include','tests','scripts'):files.extend(p for p in (ROOT/folder).rglob('*') if p.suffix in ('.cpp','.hpp','.py'))
    files.extend(p for p in (ROOT/'fixtures').rglob('*') if p.is_file())
    return files
def main():
    ap=argparse.ArgumentParser();ap.add_argument('out',type=Path);ap.add_argument('--exe',type=Path,default=ROOT/'build-generalization-cgal/Release/sdf_generalization_bench.exe');ap.add_argument('--prepare-only',action='store_true');ap.add_argument('--freeze-existing',action='store_true');ap.add_argument('--selected',type=Path);a=ap.parse_args();out=a.out.resolve();exe=a.exe.resolve()
    if not a.freeze_existing:
        assert not out.exists();out.mkdir(parents=True);(out/'inputs').mkdir();geometries=[]
        definitions=[('thin_box',.042,0,0),('thin_box',.081,0,0),('torus',.047,0,0),('torus',.073,0,0),('double_sphere',.052,0,0),('double_sphere',.087,0,0),('capsule',.075,.46,0),('capsule',.14,.31,0),('cylinder',.13,.37,0),('cylinder',.22,.055,0),('rounded_box',.047,.38,.23),('shell',.042,.46,0)]
        rng=np.random.default_rng(99173021)
        for i,(shape,x,y,z) in enumerate(definitions):
            g=dict(id=f'{i:02d}_{shape}',shape=shape,stratum='familiar_stress' if i<6 else 'novel_family',parameters=[x,y,z],rotation=rotation(rng.uniform(-np.pi,np.pi,3)).tolist(),translation=rng.uniform(-.055,.055,3).tolist(),noise=min(x/2,.035),expected_euler=[0] if shape=='torus' else [2,2] if shape in ('double_sphere','shell') else [2]);geometries.append(g)
        licenses={};reject=[]
        for i,name in enumerate(('goathead','penguin','wingnut','brucewick')):
            source=f'https://raw.githubusercontent.com/odedstein/meshes/{REV}/objects/{name}/';raw=out/'inputs'/f'{name}.obj';readme=out/'inputs'/f'{name}-README.md'
            urllib.request.urlretrieve(source+name+'.obj',raw);urllib.request.urlretrieve(source+'README.md',readme);v,f=obj(raw);closed,euler=mesh_topology(v,f)
            assert closed,'Mesh edge/vertex manifold validation failed; record rejection before choosing replacement'
            center=(v.min(axis=0)+v.max(axis=0))/2;radius=np.linalg.norm(v-center,axis=1).max();scale=.75/radius;rotation_matrix=rotation(rng.uniform(-np.pi,np.pi,3));translation=rng.uniform(-.035,.035,3);normalized=(v-center)*scale@rotation_matrix.T+translation;normpath=out/'inputs'/f'{name}-normalized.obj'
            with normpath.open('w') as stream:
                for p in normalized:stream.write('v '+' '.join(format(x,'.17g') for x in p)+'\n')
                for p in f:stream.write('f '+' '.join(str(x+1) for x in p)+'\n')
            g=dict(id=f'{i+12:02d}_{name}',shape='mesh',stratum='real_mesh',parameters=[.04,0,0],rotation=rotation_matrix.tolist(),translation=translation.tolist(),noise=.02,expected_euler=euler,normalized_mesh=str(normpath),source_mesh=str(raw),source_revision=REV,source_url=source+name+'.obj',license_text=readme.read_text(),source_sha256=sha(raw),normalized_sha256=sha(normpath),normalization=dict(original_center=center.tolist(),original_radius=float(radius),scale_to_domain=float(scale),original_units='unspecified in upstream OBJ; distances convert by dividing scale_to_domain',normalized_bounding_sphere_radius=.75),vertices=len(v),triangles=len(f));geometries.append(g)
        json_write(out/'inputs/rejections.json',reject)
        for i,g in enumerate(geometries):
            print('Qualifying',g['id'],flush=True);geometries[i]=create_reference(g,out/'reference'/g['id'],exe,seed=880001+i*1000);json_write(out/'qualification.json',geometries)
    else:geometries=json.loads((out/'qualification.json').read_text())
    if a.prepare_only:return
    assert not (out/'FROZEN.json').exists()
    assert a.selected and (out/'REFERENCES_COMPLETE.json').exists(),'Final freeze requires completed bounded development selection and all reference inputs'
    selected=json.loads(a.selected.read_text());assert selected['bounded_intervention_complete'] and selected['no_further_policy_expansion'];candidate=selected['policy'];assert candidate['probe_max_gate'];candidate['lambda']=candidate.pop('lambda_');candidate['maximum_depth']=7 if candidate.get('k7') else 6
    (out/'FINAL_POLICY.txt').write_text(f"{candidate['k']} {candidate['lambda']} {candidate['global_quadrature']} {candidate['covered_ranking']} {candidate.get('k7',0)}\n")
    frozen=dict(protocol='generalization-v1',status='frozen_before_final_fitting',domain=[[-1]*3,[1]*3],length_scale=2,fit_seeds=SEEDS,query_namespace=880001,fit_count=10000,uniform_count=100000,near_count=100000,caps=CAPS,methods=METHODS,candidate=candidate,candidate_original=dict(depth=5,k=128,lambda_= .001,quadrature=4,global_quadrature=0,k7=0),selection_sha256=sha(a.selected),selection_declaration_sha256=selected['declaration_sha256'],final_policy_sha256=sha(out/'FINAL_POLICY.txt'),sign_baseline_penalty=.0025,cvar_grid_points=9,gates=dict(feature_inside_recall=.95,feature_outside_recall=.95,sampled_max=.020,sampled_max_datasets=['uniform','near','audit','probes','boundary'],uniform_cvar99=.015,near_cvar99=.010,width_mean_relative_error=.10,width_p95_relative_error=.15,topology='closed edge and vertex manifold and expected component Euler at both qualified levels'),geometries=geometries,source_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in source_files()},executable_sha256=sha(exe),reference_sha256=manifest(out/'reference'),reference_fast_sha256=manifest(out/'reference-fast'),oracle_audit_sha256=manifest(out/'oracle-audit'),input_sha256=manifest(out/'inputs'),platform=platform.platform(),machine=platform.processor(),numpy_version=np.__version__,created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    json_write(out/'FROZEN.json',frozen)
    snapshot=out/'source';snapshot.mkdir()
    for p in source_files():d=snapshot/p.relative_to(ROOT);d.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,d)
    for p in (exe,*exe.parent.glob('*.dll')):shutil.copy2(p,snapshot/p.name)
    shutil.copy2(out/'FROZEN.json',snapshot/'FROZEN.json');json_write(snapshot/'SHA256SUMS.json',manifest(snapshot));print('FROZEN',out/'FROZEN.json',flush=True)
if __name__=='__main__':main()
