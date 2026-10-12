"""Versioned complete-ray reassessment; never change any v1 artifact."""
import subprocess,time,shutil
from generalization_v2_rays import *

def ray_gate(query,rays,reference,results):
    errors={int(g):[] for g in np.unique(rays[:,0])};passed=True;counts=0;sign_failures=0;ambiguous=0
    for ray,truth,result in zip(rays,reference,results):
        ambiguous+=bool(result['ambiguous'])
        if result['ambiguous']:passed=False
        roots=result['roots'];ok=len(roots)==len(truth) and len(roots)>0 and len(roots)%2==0
        if not ok:counts+=1;passed=False;continue
        signs=signed_intervals(query,ray,truth,roots);sign_failures+=not signs;passed &= signs
        errors[int(ray[0])].extend(abs(np.diff(roots)/np.diff(truth)-1))
    means=[];p95=[]
    for group,e in errors.items():
        if not e:passed=False;continue
        means.append(float(np.mean(e)));p95.append(float(np.quantile(e,.95,method='inverted_cdf')));passed &= means[-1]<=.10 and p95[-1]<=.15
    return dict(width_pass=bool(passed),width_status='numerical_unresolved' if ambiguous else 'passed' if passed else 'observed_geometric_failure',incorrect_root_count_rays=counts,incorrect_signed_interval_rays=sign_failures,ambiguous_rays=ambiguous,width_worst_group_mean=max(means,default=None),width_worst_group_p95=max(p95,default=None))

def main():
    old=ROOT/'runs/generalization-v1';new=ROOT/'runs/generalization-v2/reassessment';new.mkdir(exist_ok=True);frozen=json.loads((old/'FROZEN.json').read_text());geometries={g['id']:g for g in frozen['geometries']};matrix=rows(old/'verification/matrix.csv');reproduced=json.loads((old/'REPRODUCED.json').read_text())['fields'];chosen={(r['geometry'],int(r['seed']),r['method'],int(r['cap'])) for r in matrix if r['pass']=='True'}|{(r['geometry'],r['seed'],r['method'],r['cap']) for r in reproduced}
    chosen|={(g,240,'candidate',CAPS[-1]) for g in ('12_goathead','14_spot','15_brucewick')};exe=ROOT/'build-generalization-cgal/Release/sdf_generalization_v2_diagnostic.exe';identity=dict(v1_frozen_sha256=sha(old/'FROZEN.json'),native_sha256=sha(exe),python_sha256=sha(Path(__file__).with_name('generalization_v2_rays.py')),script_sha256=sha(Path(__file__)));cache={};output=[]
    for i,key in enumerate(sorted(chosen)):
        name,seed,method,cap=key;g=geometries[name];folder=old/'jobs'/f'{name}_{seed}_{method}';ref=reference_dir(old,g);record=next(r for r in matrix if (r['geometry'],int(r['seed']),r['method'],int(r['cap']))==key);native=next(r for r in rows(folder/'fields.csv') if int(r['cap'])==cap);data=(folder/f'{cap}.bin').read_bytes();n=int(native['grid_n']);hashkey=(name,hashlib.sha256(data).hexdigest(),n)
        rays=np.loadtxt(ref/'rays.txt');reference=[np.fromstring(r['roots'],sep=';') for r in rows(ref/'reference_roots.csv')]
        query=lambda p:decode_grid(data,n,p) if n else decode_packed(data,p)
        if hashkey not in cache:
            target=new/f'{name}_{hashkey[1][:16]}';target.mkdir(exist_ok=True);out=target/'native'
            expected_identity=dict(identity=identity,field_sha256=hashkey[1],rays_sha256=sha(ref/'rays.txt'),reference_roots_sha256=sha(ref/'reference_roots.csv'))
            if (target/'NATIVE_IDENTITY.json').exists():assert json.loads((target/'NATIVE_IDENTITY.json').read_text())==expected_identity,'Assessment input/evaluator changed; preserve and use a new version'
            else:json_write(target/'NATIVE_IDENTITY.json',expected_identity)
            assert hashkey[1]==json.loads((folder/'SHA256SUMS.json').read_text())[f'{cap}.bin'],'V1 field changed'
            prior=ROOT/'runs/generalization-v2/reassessment-initial-interrupted'/target.name
            if not out.exists() and (prior/'native/COMPLETE').exists():
                provenance=json.loads((prior/'NATIVE_IDENTITY.json').read_text());assert provenance['identity']['native_sha256']==identity['native_sha256'] and all(provenance[k]==expected_identity[k] for k in ('field_sha256','rays_sha256','reference_roots_sha256'))
                shutil.copytree(prior/'native',out)
            if not (out/'COMPLETE').exists():
                assert not out.exists(),'Incomplete assessment preserved';subprocess.run([str(exe),'field',str(folder/f'{cap}.bin'),str(ref),str(out),str(n)],check=True,stdout=subprocess.DEVNULL)
            cpp=rows(out/'roots.csv');results=[]
            for j,ray in enumerate(rays):
                if n:
                    times=[ray[7],ray[8]]
                    for a in range(3):
                        if ray[4+a]:times.extend(t for t in ((np.arange(1,n-1)*2/(n-1)-1-ray[1+a])/ray[4+a]) if ray[7]<t<ray[8])
                    times=np.unique(times);segments=[];values=np.frombuffer(data,dtype='<f4').reshape(n,n,n)
                    for lo,hi in zip(times[:-1],times[1:]):
                        point=ray[1:4]+ray[4:7]*(lo+hi)/2;index=np.clip(np.floor((point+1)*(n-1)/2).astype(int),0,n-2);low=-1+index*2/(n-1);high=-1+(index+1)*2/(n-1);corners=[values[index[2]+bool(c&4),index[1]+bool(c&2),index[0]+bool(c&1)] for c in range(8)];segments.append([lo,hi,*low,*high,*corners])
                    segments=np.array(segments)
                else:segments=packed_segments(data,ray)
                result=transitions(segments,ray);native_roots=np.fromstring(cpp[j]['roots'],sep=';');native_ambiguous=bool(int(cpp[j]['ambiguous']))
                if not native_ambiguous and not result['ambiguous']:np.testing.assert_allclose(result['roots'],native_roots,rtol=2e-8,atol=3e-9)
                result['ambiguous'] |= native_ambiguous;results.append(result)
            assessment=ray_gate(query,rays,reference,results);cache[hashkey]=assessment;json_write(target/'ASSESSMENT.json',dict(identity=identity,field_sha256=hashkey[1],assessment=assessment))
        assessment=cache[hashkey];passed=assessment['width_pass'] and all(record[k]=='True' for k in ('scalar_pass','probe_pass','topology_pass','reference_qualified'));status='passed' if passed else 'observed_failure' if record['scalar_pass']!='True' or record['probe_pass']!='True' or assessment['width_status']=='observed_geometric_failure' else 'reference_unresolved' if record['reference_qualified']!='True' else 'numerical_unresolved' if assessment['ambiguous_rays'] else 'topology_not_evaluated_or_failed';output.append(dict(geometry=name,seed=seed,method=method,cap=cap,bytes=int(record['bytes']),v1_joint_pass=record['pass']=='True',v2_complete_ray_joint_pass=bool(passed),v2_joint_status=status,v1_width_pass=record['width_pass']=='True',topology_status=record['topology_status'],**assessment));print(f'{i+1}/{len(chosen)} reassessed {name} {seed} {method} {cap}',flush=True)
        if (i+1)%25==0:csv_write(new/'partial_matrix.csv',output)
    csv_write(new/'matrix.csv',output);json_write(new/'COMPLETE.json',dict(identity=identity,fields=len(output),unique_fields=len(cache),old_passes_rejected=sum(r['v1_joint_pass'] and not r['v2_complete_ray_joint_pass'] for r in output),complete_rays=True,prior_artifacts_unchanged=True));json_write(new/'SHA256SUMS.json',manifest(new))
if __name__=='__main__':main()
