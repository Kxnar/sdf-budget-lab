"""Serial query-local full-depth controls and independent complete-ray checks."""
import argparse,subprocess,time
from generalization_v2_rays import *

def virtual_query(segments,ray,encoding,times):
    values=decoded_corners(segments[:,8:],encoding);indices=np.searchsorted(segments[:,1],times,side='right').clip(0,len(segments)-1);s=segments[indices];v=values[indices];u=(ray[1:4]+np.asarray(times)[:,None]*ray[4:7]-s[:,2:5])/(s[:,5:8]-s[:,2:5]);result=np.zeros(len(times))
    for c in range(8):result+=v[:,c]*np.prod(np.where(np.array([bool(c&(1<<a)) for a in range(3)]),u,1-u),axis=1)
    return result

def verify_geometry(g,ref,folder):
    rays=np.loadtxt(ref/'rays.txt');reference=[np.fromstring(r['roots'],sep=';') for r in rows(ref/'reference_roots.csv')];probes=np.fromfile(ref/'probes.bin',dtype='<f8').reshape(-1,4);groups=np.loadtxt(ref/'probe_groups.txt',dtype=int);native=rows(folder/'roots.csv');predictions=rows(folder/'probes.csv');output=[];ray_details=[]
    for depth in (7,9,12):
        probe_cells=np.fromfile(folder/f'probe_cells_{depth}.bin',dtype='<f8').reshape(-1,16);assert len(probe_cells)==len(probes);width=probe_cells[:,5:8]-probe_cells[:,2:5];np.testing.assert_array_equal(width,np.full_like(width,2/(1<<depth)))
        for encoding in ('q8','f32','f64'):
            values=decoded_corners(probe_cells[:,8:],encoding);u=(probes[:,:3]-probe_cells[:,2:5])/width;pred=np.zeros(len(probes))
            for c in range(8):pred+=values[:,c]*np.prod(np.where(np.array([bool(c&(1<<a)) for a in range(3)]),u,1-u),axis=1)
            cpp=np.array([float(r['prediction']) for r in predictions if int(r['depth'])==depth and r['encoding']==encoding]);np.testing.assert_allclose(pred,cpp,atol=2e-13,rtol=2e-10);feature=True;recalls=[]
            for group in g['feature_groups']:
                mask=groups==group;inside=mask&(probes[:,3]<-2e-9);outside=mask&(probes[:,3]>2e-9);ir=float(np.mean(pred[inside]<0));orr=float(np.mean(pred[outside]>0));recalls.append(dict(group=group,inside_recall=ir,outside_recall=orr));feature &= ir>=.95 and orr>=.95
            results=[];all_segments=[]
            for i,ray in enumerate(rays):
                segments=np.fromfile(folder/f'segments_{depth}_{i}.bin',dtype='<f8').reshape(-1,16);assert segments[0,0]==ray[7] and segments[-1,1]==ray[8];np.testing.assert_array_equal(segments[:-1,1],segments[1:,0]);np.testing.assert_array_equal(segments[:,5:8]-segments[:,2:5],np.full((len(segments),3),2/(1<<depth)));result=transitions(segments,ray,encoding);record=next(r for r in native if int(r['depth'])==depth and r['encoding']==encoding and int(r['ray'])==i);cpp=np.fromstring(record['roots'],sep=';');native_ambiguous=bool(int(record['ambiguous']))
                if not native_ambiguous and not result['ambiguous']:np.testing.assert_allclose(result['roots'],cpp,atol=3e-9,rtol=2e-8)
                result['ambiguous'] |= native_ambiguous;results.append(result);all_segments.append(segments)
                ray_details.append(dict(depth=depth,encoding=encoding,ray=i,reference_count=len(reference[i]),predicted_count=len(result['roots']),ambiguous=bool(result['ambiguous']),roots=';'.join(format(float(v),'.17g') for v in result['roots'])))
            # The query callback is invoked for each ray in declared order, twice
            # for interval signs plus endpoints; coordinates identify its ray.
            errors={group:[] for group in g['feature_groups']};counts=0;signs=0;ambiguous=0
            for i,(ray,truth,result) in enumerate(zip(rays,reference,results)):
                ambiguous+=bool(result['ambiguous']);roots=result['roots'];ok=len(roots)==len(truth) and len(roots)>0 and len(roots)%2==0
                if not ok:counts+=1;continue
                endpoint=virtual_query(all_segments[i],ray,encoding,ray[[7,8]]);signed=bool(np.all(endpoint>0))
                for rr in (truth,roots):
                    edges=np.r_[ray[7],rr,ray[8]];times=(edges[:-1]+edges[1:])/2;expected=np.where(np.arange(len(times))%2,-1,1);signed &= bool(np.all(virtual_query(all_segments[i],ray,encoding,times)*expected>0))
                signs+=not signed;errors[int(ray[0])].extend(abs(np.diff(roots)/np.diff(truth)-1))
            means=[float(np.mean(v)) for v in errors.values() if len(v)];p95=[float(np.quantile(v,.95,method='inverted_cdf')) for v in errors.values() if len(v)];width_pass=not(counts or signs or ambiguous) and len(means)==len(errors) and max(means)<=.10 and max(p95)<=.15
            output.append(dict(geometry=g['id'],depth=depth,encoding=encoding,feature_pass=bool(feature),wrong_probes=int(np.sum(pred*np.sign(probes[:,3])<=0)),probe_sampled_max=float(np.max(abs(pred-probes[:,3]))),width_pass=bool(width_pass),width_status='numerical_unresolved' if ambiguous else 'passed' if width_pass else 'observed_geometric_failure',incorrect_root_count_rays=counts,incorrect_signed_interval_rays=signs,ambiguous_rays=ambiguous,width_worst_group_mean=max(means,default=None),width_worst_group_p95=max(p95,default=None),feature_recalls=json.dumps(recalls)))
            print(g['id'],depth,encoding,output[-1]['width_status'],'root mismatches',counts,'width',max(means,default=None),max(p95,default=None),flush=True)
    csv_write(folder/'independent-summary.csv',output);csv_write(folder/'independent-rays.csv',ray_details);json_write(folder/'VERIFIED.json',dict(geometry=g['id'],complete_96_rays=True,all_probes=True,rows=output,native_sha256=sha(ROOT/'build-generalization-cgal/Release/sdf_generalization_v2_diagnostic.exe'),python_sha256=sha(Path(__file__).with_name('generalization_v2_rays.py')),verification_sha256=sha(Path(__file__)),source_frozen_sha256=sha(ROOT/'runs/generalization-v1/FROZEN.json'),storage_contenders=False));json_write(folder/'SHA256SUMS.json',manifest(folder));return output

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--geometry',nargs='*');a=ap.parse_args();old=ROOT/'runs/generalization-v1';out=ROOT/'runs/generalization-v2/feasibility';out.mkdir(exist_ok=True);frozen=json.loads((old/'FROZEN.json').read_text());exe=ROOT/'build-generalization-cgal/Release/sdf_generalization_v2_diagnostic.exe';all_rows=[]
    for g in frozen['geometries']:
        if g['id'] not in ('12_goathead','14_spot','15_brucewick') or (a.geometry and g['id'] not in a.geometry):continue
        ref=reference_dir(old,g);folder=out/g['id']
        if not (folder/'COMPLETE').exists():
            assert not folder.exists(),'Incomplete diagnostic preserved; use a fresh version'
            with (out/(g['id']+'.log')).open('w') as log:subprocess.run([str(exe),'virtual',str(ref/'geometry.txt'),str(ref),str(folder)],check=True,stdout=log,stderr=subprocess.STDOUT)
        all_rows.extend(verify_geometry(g,ref,folder));csv_write(out/'summary.csv',all_rows)
    if not a.geometry:
        assert len(all_rows)==27;deep=[r for r in all_rows if r['depth']==12];json_write(out/'CHECKPOINT.json',dict(deep_controls=deep,requires_representation_replan=any(not r['feature_pass'] or r['width_status']!='passed' for r in deep if r['encoding'] in ('f32','f64')),q8_only_failures=[r for r in deep if r['encoding']=='q8' and not r['width_pass']],not_full_gate_certification=True));json_write(out/'SHA256SUMS.json',manifest(out))
if __name__=='__main__':main()
