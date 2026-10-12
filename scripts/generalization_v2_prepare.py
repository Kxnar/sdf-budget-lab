"""Reference-only reserved mesh acquisition/qualification, before any fresh fields."""
import argparse,subprocess,urllib.request,time,shutil
from generalization_common import *
from generalization_prepare import REV,mesh_rays,write_spec,create_reference
from compact_mesh import ray_events

FILES={'fish':'fish_low_resolution.obj','mushroom':'mushroom.obj','springer':'springer.obj','hammer':'hammer.obj','stuffedtoy':'stuffedtoy.obj','boat':'boat.obj'}
LICENSES={'fish':('Public domain','Keenan Crane'),'mushroom':('CC BY-NC 4.0','Josh Holinaty'),'springer':('GPLv2','TimEdwards'),'hammer':('CC BY 4.0','Andrea Spognetta'),'stuffedtoy':('Individual README follows libigl license; verify linked scope before use','Ilya Baran and Cosmic blobs'),'boat':('CC BY 4.0','Oded Stein, Eitan Grinspun and Keenan Crane')}
PRIMARY=('fish','mushroom','springer','hammer');FALLBACK=('stuffedtoy','boat')
class InputRejected(ValueError):pass

def execution_identity(folder,g,exe,phase):
    expected=dict(phase=phase,normalized_sha256=g['normalized_sha256'],source_sha256=g['source_sha256'],spec_sha256=sha(folder/'geometry.txt'),executable_sha256=sha(exe),dll_sha256={p.name:sha(p) for p in sorted(exe.parent.glob('*.dll'))},helper_sha256={name:sha(Path(__file__).with_name(name)) for name in ('generalization_v2_prepare.py','generalization_prepare.py','generalization_common.py','compact_mesh.py')})
    path=folder/'EXECUTION_IDENTITY.json'
    if path.exists():assert json.loads(path.read_text())==expected,'Reference/audit identity changed: preserve completed evidence and use its frozen executable/helpers'
    else:json_write(path,expected)

def declaration(out):
    rng=np.random.default_rng(39173021);slots=[]
    for i,name in enumerate(PRIMARY):
        angles=rng.uniform(-np.pi,np.pi,3);translation=rng.uniform(-.035,.035,3)
        slots.append(dict(slot=i,primary=name,rotation_angles=angles.tolist(),rotation=rotation(angles).tolist(),translation=translation.tolist(),query_seed=1880001+i*1000,kernel_audit_seed=591007+i))
    record=dict(protocol='generalization-v2-compact',revision=REV,slots=slots,fallback_order=list(FALLBACK),files=FILES,licenses={k:list(v) for k,v in LICENSES.items()},fit_seeds=list(range(340,345)),fit_count=10000,uniform_count=100000,near_count=100000,noise=.02,normalized_radius=.75,reference_ray_seed=713007,reference_min_interval_eligibility=1e-5,reference_hit_proximity_dedup=1e-8,reference_probe_min_absolute_label=1e-8,reference_ladder=[64,96,128,192],kernel_comparison_stratified_queries=8192,kernel_distance_threshold=1e-10,kernel_sign_band=1e-9,no_fresh_representation_fitting=True)
    path=out/'PREPARATION-DECLARATION.json'
    if path.exists():assert json.loads(path.read_text())==record
    else:json_write(path,record)
    return record

def acquire(out,name,slot):
    folder=out/'inputs';folder.mkdir(exist_ok=True);base=f'https://raw.githubusercontent.com/odedstein/meshes/{REV}/objects/{name}/'
    raw=folder/f'{name}.obj';readme=folder/f'{name}-README.md'
    for url,path in ((base+FILES[name],raw),(base+'README.md',readme)):
        if not path.exists():
            request=urllib.request.Request(url,headers={'User-Agent':'SDFGeneralization-v2'})
            path.write_bytes(urllib.request.urlopen(request).read())
    if name=='stuffedtoy':raise InputRejected('Fallback license scope must be resolved from its linked libigl license before use; no representation fitting')
    v,f=obj(raw);closed,euler=mesh_topology(v,f)
    if not closed:raise InputRejected('Source oriented edge/cyclic vertex-link/nondegeneracy validation failed')
    center=(v.min(axis=0)+v.max(axis=0))/2;radius=float(np.linalg.norm(v-center,axis=1).max());scale=.75/radius
    matrix=np.array(slot['rotation']);translation=np.array(slot['translation']);normalized=((v-center)*scale)@matrix.T+translation
    path=folder/f'{name}-normalized.obj'
    if not path.exists():
        with path.open('w',encoding='utf-8') as stream:
            for p in normalized:stream.write('v '+' '.join(format(float(x),'.17g') for x in p)+'\n')
            for t in f:stream.write('f '+' '.join(str(int(x)+1) for x in t)+'\n')
    saved_vertices,saved_faces=obj(path);np.testing.assert_array_equal(saved_vertices,normalized);np.testing.assert_array_equal(saved_faces,f)
    source_identity=folder/f'{name}-SOURCE.json';expected_source=dict(source_url=base+FILES[name],readme_url=base+'README.md',source_revision=REV,source_sha256=sha(raw),readme_sha256=sha(readme),normalized_sha256=sha(path),slot=slot)
    if source_identity.exists():assert json.loads(source_identity.read_text())==expected_source,'Previously acquired source/normalization changed; preserve and investigate'
    else:json_write(source_identity,expected_source)
    return dict(id=f"{slot['slot']:02d}_{name}",slot=slot['slot'],shape='mesh',stratum='fresh_real_mesh',parameters=[.04,0,0],rotation=slot['rotation'],translation=slot['translation'],noise=.02,expected_euler=euler,feature_groups=[0,1,2],normalized_mesh=str(path.resolve()),source_mesh=str(raw.resolve()),source_revision=REV,source_url=base+FILES[name],readme_url=base+'README.md',license_text=readme.read_text(encoding='utf-8'),license=LICENSES[name][0],author=LICENSES[name][1],source_sha256=sha(raw),readme_sha256=sha(readme),normalized_sha256=sha(path),vertices=len(v),triangles=len(f),normalization=dict(original_center=center.tolist(),original_radius=radius,scale_to_domain=scale,original_units='Unspecified upstream OBJ units; divide normalized distances by scale_to_domain',normalized_bounding_sphere_radius=.75),query_seed=slot['query_seed'],kernel_audit_seed=slot['kernel_audit_seed'])

def audit(g,out,exe):
    folder=out/'oracle-audit'/g['id'];folder.mkdir(parents=True,exist_ok=True);write_spec(g,folder/'geometry.txt');v,f=obj(g['normalized_mesh']);centers=v[f].mean(axis=1)
    execution_identity(folder,g,exe,'EPECK_input_validation_and_paired_kernel_audit')
    edges=sorted({tuple(sorted((int(a),int(b)))) for t in f for a,b in zip(t,np.roll(t,-1))});edge_points=np.array([(v[a]+v[b])*.5 for a,b in edges]);normal=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]);normal/=np.linalg.norm(normal,axis=1)[:,None]
    rays,roots=mesh_rays(g,v,f);witnesses=[]
    for ray,rr in zip(rays,roots):
        cuts=np.r_[ray[7],rr,ray[8]];witnesses.extend(ray[1:4]+((cuts[:-1]+cuts[1:])*.5)[:,None]*ray[4:7])
    points=np.concatenate([v,centers,edge_points,centers+1e-7*normal,centers-1e-7*normal,np.array(witnesses)]);np.savetxt(folder/'comparison.points',points,fmt='%.17g')
    if not (folder/'COMPLETE').exists():
        with (folder/'comparison.log').open('w',encoding='utf-8') as log:
            result=subprocess.run([str(exe),'compare',str(folder/'geometry.txt'),str(folder),'4096',str(g['kernel_audit_seed'])],stdout=log,stderr=subprocess.STDOUT)
        if result.returncode:
            message=(folder/'comparison.log').read_text()
            input_errors=('Mesh must be nonempty, closed, and triangulated','Normalize mesh vertices into [-1,1]^3 before loading','Degenerate triangle in mesh','Invalid triangle area','Self-intersecting mesh')
            if any(error in message for error in input_errors):raise InputRejected('Demonstrated EPECK input validation failure: '+message)
            json_write(folder/'UNRESOLVED.json',dict(status='runtime_resource_or_oracle_unresolved',native_exit=result.returncode,diagnostic=message,primary_preserved=True));raise RuntimeError('Native audit failed without demonstrated invalid input; investigation required before replacement')
    exact_samples=np.fromfile(folder/'epeck.bin',dtype='<f8').reshape(-1,4);fast=np.fromfile(folder/'epick.bin',dtype='<f8');assert len(fast)==len(exact_samples)==8192+len(points)
    difference=abs(fast-exact_samples[:,3]);away=abs(exact_samples[:,3])>1e-9;sign_mismatch=away&(np.sign(fast)!=np.sign(exact_samples[:,3]));native=json.loads((folder/'comparison.json').read_text());eligible=bool(difference.max()<=1e-10 and not sign_mismatch.any() and native['sampler_max_delta']<=1e-12)
    indices=np.r_[np.linspace(0,8191,128,dtype=int),np.arange(len(exact_samples)-min(96,len(witnesses)),len(exact_samples))];independent=brute_mesh(g,exact_samples[indices,:3]);np.testing.assert_allclose(independent,exact_samples[indices,3],atol=2e-10,rtol=2e-9)
    np.savez(folder/'independent.npz',points=exact_samples[indices,:3],epeck=exact_samples[indices,3],independent=independent)
    record=dict(geometry=g['id'],epeck_input_validated=True,epick_eligible=eligible,stratified_queries=8192,total_queries=len(fast),maximum_distance_difference=float(difference.max()),unexplained_sign_mismatches=int(sign_mismatch.sum()),tiny_sign_band=1e-9,near_band_sign_mismatches=int(np.sum((~away)&(np.sign(fast)!=np.sign(exact_samples[:,3])))),sampler_max_delta=native['sampler_max_delta'],independent_triangle_queries=len(indices),epeck_phase_seconds=native['epeck_seconds'],epick_distance_only_seconds=native['epick_seconds'],timing_scope='Scopes differ; no matched query-speed ratio',normalized_sha256=g['normalized_sha256'])
    json_write(folder/'RESULT.json',record);json_write(folder/'SHA256SUMS.json',manifest(folder));return eligible

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='runs/generalization-v2/fresh');ap.add_argument('--exe',default='build-generalization-v2-cgal/Release/sdf_generalization_v2_baseline.exe');ap.add_argument('--declaration-only',action='store_true');args=ap.parse_args();out=ROOT/args.out;out.mkdir(parents=True,exist_ok=True);decl=declaration(out)
    if args.declaration_only:return
    helpers=out/'preparation-source';identity=dict(executable_sha256=sha((ROOT/args.exe).resolve()),helpers={name:sha(Path(__file__).with_name(name)) for name in ('generalization_v2_prepare.py','generalization_prepare.py','generalization_common.py','compact_mesh.py')})
    if (out/'PREPARATION-SOURCE.json').exists():assert json.loads((out/'PREPARATION-SOURCE.json').read_text())==identity,'Preparation source changed; resume with its frozen modules/executable'
    else:
        helpers.mkdir(exist_ok=False)
        for name in identity['helpers']:shutil.copy2(Path(__file__).with_name(name),helpers/name)
        shutil.copy2((ROOT/args.exe).resolve(),helpers/Path(args.exe).name)
        for dll in (ROOT/args.exe).resolve().parent.glob('*.dll'):shutil.copy2(dll,helpers/dll.name)
        json_write(out/'PREPARATION-SOURCE.json',identity)
    assert not (out/'FROZEN.json').exists(),'Fresh fitting freeze already exists; reference evidence is immutable'
    exe=(ROOT/args.exe).resolve();geometries=json.loads((out/'qualification.json').read_text()) if (out/'qualification.json').exists() else [];reject=json.loads((out/'rejections.json').read_text()) if (out/'rejections.json').exists() else [];used={g['id'].split('_',1)[1] for g in geometries}
    for slot in decl['slots']:
        if any(g['slot']==slot['slot'] for g in geometries):continue
        accepted=None
        for name in (slot['primary'],*FALLBACK):
            if name in used:continue
            try:g=acquire(out,name,slot);eligible=audit(g,out,exe);accepted=g;break
            except (InputRejected,urllib.error.URLError) as e:
                reject.append(dict(slot=slot['slot'],name=name,reason=str(e),stage='acquisition_license_or_independent_input_validation',before_any_fresh_fields=True));json_write(out/'rejections.json',reject)
        if accepted is None:raise RuntimeError('Declared input/fallback list exhausted; protocol review required')
        g=accepted;g['parameters'][1]=1 if eligible else 0;ref=out/'reference'/g['id'];print('reference qualification',g['id'],flush=True)
        if not (ref/'QUALIFIED').exists():
            if ref.exists():raise RuntimeError('Incomplete reference attempt preserved; investigate before resuming')
            # create_reference creates the folder itself. Bind immutable input,
            # executable and helper hashes in an adjacent identity beforehand.
            identity_folder=out/'reference-identities'/g['id'];identity_folder.mkdir(parents=True,exist_ok=True);write_spec(g,identity_folder/'geometry.txt');execution_identity(identity_folder,g,exe,'reference_qualification')
            g=create_reference(g,ref,exe,seed=g['query_seed']);shutil.copy2(identity_folder/'EXECUTION_IDENTITY.json',ref/'EXECUTION_IDENTITY.json')
        else:
            identity_folder=out/'reference-identities'/g['id'];execution_identity(identity_folder,g,exe,'reference_qualification')
            if not (ref/'EXECUTION_IDENTITY.json').exists():shutil.copy2(identity_folder/'EXECUTION_IDENTITY.json',ref/'EXECUTION_IDENTITY.json')
            execution_identity(ref,g,exe,'reference_qualification')
            g['qualified_levels']=[int(x) for x in (ref/'qualified_levels.txt').read_text().split()];g['reference_status']='qualified' if g['qualified_levels'] else 'oracle_resolution_unresolved';g['ray_count']=96
        v,f=obj(g['normalized_mesh']);rs=np.loadtxt(ref/'rays.txt');native=rows(ref/'reference_roots.csv');records=[]
        for i,r in enumerate(rs):
            result=ray_events(v,f,r);expected=np.fromstring(native[i]['roots'],sep=';')
            if result['unresolved'] or len(result['roots'])!=len(expected):
                json_write(ref/'EXACT-REFERENCE-CHECKPOINT.json',dict(ray=i,reference=expected.tolist(),exact=result,requires_reference_consistency_replan=True));raise RuntimeError('Reference-only exact event inconsistency; preserve case and review reference consistency before fitting fields')
            np.testing.assert_allclose(result['roots'],expected,atol=2e-8,rtol=1e-10);records.append(result)
        json_write(ref/'EXACT-REFERENCE-EVENTS.json',records);json_write(ref/'SHA256SUMS.json',manifest(ref));used.add(g['id'].split('_',1)[1]);geometries.append(g);json_write(out/'qualification.json',geometries);json_write(out/'rejections.json',reject)
    assert len(geometries)==4;json_write(out/'REFERENCES_COMPLETE.json',dict(geometries=geometries,reference_hashes=manifest(out/'reference'),input_hashes=manifest(out/'inputs'),oracle_audit_hashes=manifest(out/'oracle-audit'),declaration_sha256=sha(out/'PREPARATION-DECLARATION.json'),no_fresh_fields=True))
if __name__=='__main__':main()
