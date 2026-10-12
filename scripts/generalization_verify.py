"""Independent full packed/grid decode, strict gates and Cartesian completeness."""
import argparse,time,sys,shutil
from generalization_common import *
from generalization_run import validate_freeze,expected_jobs
from distortion_analysis import layout
from generalization_common import decode_packed as decode

def metrics(pred,truth):
    e=np.abs(pred-truth);sorted_e=np.sort(e);inside=truth<-2e-9;outside=truth>2e-9;valid=inside|outside
    return dict(count=len(e),mae=float(e.mean()),p99=float(sorted_e[math.ceil(.99*len(e))-1]),cvar99=float(sorted_e[-math.ceil(.01*len(e)):].mean()),max_error=float(e.max()),inside_recall=float(np.mean(pred[inside]<0)) if inside.any() else float('nan'),outside_recall=float(np.mean(pred[outside]>0)) if outside.any() else float('nan'),sign_error=float(np.mean(pred[valid]*np.sign(truth[valid])<=0)) if valid.any() else float('nan'))
def close(actual,expected):np.testing.assert_allclose(actual,float(expected),rtol=5e-10,atol=2e-11,equal_nan=True)
def strict_gate(scalar,probe,width,topology,qualified):return scalar and probe and width and topology and qualified
def ray_sign_pattern(query,ray,reference_roots,recovered_roots):
    origin=ray[1:4];direction=ray[4:7];endpoint_values=query(origin+ray[[7,8]][:,None]*direction)
    if not np.all(endpoint_values>0):return False
    for crossing_roots in (reference_roots,recovered_roots):
        edges=np.r_[ray[7],crossing_roots,ray[8]];mid=(edges[:-1]+edges[1:])/2;expected_sign=np.where(np.arange(len(mid))%2,-1,1);values=query(origin+mid[:,None]*direction)
        if not np.all(values*expected_sign>0):return False
    return True
def evaluator_identity():
    names=['generalization_verify.py','generalization_common.py','generalization_run.py','generalization_prepare.py','distortion_analysis.py'];return dict(source_sha256={n:sha(Path(__file__).with_name(n)) for n in names},numpy_version=np.__version__,python_version=sys.version)
def fitting_evidence(folder,g,method):
    if method=='grid':return None,None
    fit=np.fromfile(folder/'fit.bin',dtype='<f8').reshape(-1,4);assert len(fit)==10000;policy=json.loads((folder/'policy.json').read_text());global_samples=np.fromfile(folder/'global.bin',dtype='<f8').reshape(-1,4) if (folder/'global.bin').exists() else np.empty((0,4));extra=np.fromfile(folder/'extra.bin',dtype='<f8').reshape(-1,4) if (folder/'extra.bin').exists() else np.empty((0,4))
    if len(global_samples):
        q=policy['global_quadrature'];assert len(global_samples)==32768*q**3
        expected=[]
        for iz in range(32):
            for iy in range(32):
                for ix in range(32):
                    for z in range(q):
                        for y in range(q):
                            for x in range(q):expected.append([-1+(ix+(x+.5)/q)/16,-1+(iy+(y+.5)/q)/16,-1+(iz+(z+.5)/q)/16])
        np.testing.assert_array_equal(global_samples[:,:3],np.array(expected))
    selection=rows(folder/'selection.csv') if (folder/'selection.csv').exists() else[]
    if len(extra):
        expected=[]
        for row in selection:
            low=np.array([float(row['lo_'+x]) for x in 'xyz']);high=np.array([float(row['hi_'+x]) for x in 'xyz'])
            for z in range(4):
                for y in range(4):
                    for x in range(4):expected.append(low+(high-low)*(np.array([x,y,z])+.5)/4)
        np.testing.assert_array_equal(extra[:,:3],np.array(expected));assert len(extra)==64*len(selection)
    for sample in (fit,global_samples,extra):
        if g['shape']!='mesh' and len(sample):np.testing.assert_allclose(sample[:,3],exact(g,sample[:,:3]),atol=2e-13,rtol=2e-11)
    rank_pool=np.concatenate([fit,global_samples]) if policy['covered_ranking'] else fit
    if selection:
        if (folder/'depth5_corners.bin').exists():corners=np.fromfile(folder/'depth5_corners.bin',dtype='<f4').reshape(-1,8);children=np.fromfile(folder/'depth5_children.bin',dtype='<u4')
        else:
            # Old development raws predate the diagnostic export. Reconstruct
            # this analytic-only tree independently; never rewrite raw jobs.
            assert g['shape']!='mesh';boxes=[(np.full(3,-1.),np.full(3,1.))];first=[0]
            def build(v,depth):
                if depth==5:return
                low,high=boxes[v];mid=(low+high)/2;start=len(boxes);first[v]=start
                for c in range(8):side=np.array([bool(c&(1<<i)) for i in range(3)]);boxes.append((np.where(side,mid,low),np.where(side,high,mid)));first.append(0)
                for c in range(8):build(start+c,depth+1)
            build(0,0);children=np.array(first);low=np.array([b[0] for b in boxes]);high=np.array([b[1] for b in boxes]);corners=np.stack([exact(g,np.where(np.array([bool(c&(1<<i)) for i in range(3)]),high,low)).astype(np.float32) for c in range(8)],axis=1)
        assert len(corners)==len(children)==37449
        lo=np.full_like(rank_pool[:,:3],-1.);hi=np.full_like(lo,1.);nodes=np.zeros(len(lo),dtype=int)
        for depth in range(5):mid=(lo+hi)/2;side=rank_pool[:,:3]>=mid;nodes=children[nodes]+side@np.array([1,2,4]);lo=np.where(side,mid,lo);hi=np.where(side,hi,mid)
        vv=corners[nodes].astype(float);offset=vv.min(axis=1);scale=((vv.max(axis=1)-offset)/255).astype(np.float32).astype(float);safe=np.where(scale>0,scale,1);quantized=np.clip(np.floor((vv-offset[:,None])/safe[:,None]+.5),0,255);decoded=(offset[:,None]+scale[:,None]*quantized).astype(np.float32).astype(float);t=(rank_pool[:,:3]-lo)/(hi-lo);pred=np.zeros(len(lo))
        for c in range(8):side=np.array([bool(c&(1<<i)) for i in range(3)]);pred+=decoded[:,c]*np.prod(np.where(side,t,1-t),axis=1)
        loss=np.bincount(nodes,weights=((pred-rank_pool[:,3])/2)**2,minlength=37449);counts=np.bincount(nodes,minlength=37449);ranked=sorted((int(v) for v in np.flatnonzero((loss>0)&(counts>0))),key=lambda v:(-loss[v],v));assert len(selection)==min(policy['k'] if method in ('candidate','candidate_original','coverage_q8') else 32,len(ranked))
        for i,row in enumerate(selection):
            node=int(row['node']);assert node==ranked[i],f'Independent residual rank mismatch at {i}';close(loss[node],row['q8_loss']);assert counts[node]==int(row['fit_count'])
    allocation=np.concatenate([rank_pool,extra]) if policy['covered_ranking'] else np.concatenate([fit,extra,global_samples])
    if policy.get('k7',0):
        selection7=rows(folder/'selection7.csv');extra7=np.fromfile(folder/'extra7.bin',dtype='<f8').reshape(-1,4);expected=[]
        for row in selection7:
            low=np.array([float(row['lo_'+x]) for x in 'xyz']);high=np.array([float(row['hi_'+x]) for x in 'xyz']);np.testing.assert_allclose(high-low,1/32,atol=0)
            for z in range(4):
                for y in range(4):
                    for x in range(4):expected.append(low+(high-low)*(np.array([x,y,z])+.5)/4)
        np.testing.assert_array_equal(extra7[:,:3],np.array(expected));assert len(extra7)==64*len(selection7)
        if g['shape']!='mesh':np.testing.assert_allclose(extra7[:,3],exact(g,extra7[:,:3]),atol=2e-13,rtol=2e-11)
        corners6=np.fromfile(folder/'depth6_corners.bin',dtype='<f4').reshape(-1,8);children6=np.fromfile(folder/'depth6_children.bin',dtype='<u4');assert len(corners6)==len(children6)==37449+8*len(selection)
        lo=np.full_like(allocation[:,:3],-1.);hi=np.full_like(lo,1.);node=np.zeros(len(lo),dtype=int);depth=np.zeros(len(lo),dtype=int)
        for _ in range(6):
            active=children6[node]>0
            if not active.any():break
            mid=(lo[active]+hi[active])/2;side=allocation[active,:3]>=mid;node[active]=children6[node[active]]+side@np.array([1,2,4]);lo[active]=np.where(side,mid,lo[active]);hi[active]=np.where(side,hi[active],mid);depth[active]+=1
        vv=corners6[node].astype(float);offset=vv.min(axis=1);scale=((vv.max(axis=1)-offset)/255).astype(np.float32).astype(float);q=np.clip(np.floor((vv-offset[:,None])/np.where(scale>0,scale,1)[:,None]+.5),0,255);decoded=(offset[:,None]+scale[:,None]*q).astype(np.float32).astype(float);t=(allocation[:,:3]-lo)/(hi-lo);pred=np.zeros(len(lo))
        for c in range(8):side=np.array([bool(c&(1<<i)) for i in range(3)]);pred+=decoded[:,c]*np.prod(np.where(side,t,1-t),axis=1)
        active=depth==6;loss=np.bincount(node[active],weights=((pred[active]-allocation[active,3])/2)**2,minlength=len(corners6));counts=np.bincount(node[active],minlength=len(corners6));ranked=sorted((int(v) for v in np.flatnonzero((loss>0)&(counts>0))),key=lambda v:(-loss[v],v));assert len(selection7)==min(policy['k7'],len(ranked))
        for i,row in enumerate(selection7):v=int(row['node']);assert v==ranked[i],'Independent depth-six residual rank mismatch';close(loss[v],row['q8_loss']);assert counts[v]==int(row['fit_count'])
        allocation=np.concatenate([allocation,extra7])
    return allocation,policy
def verify_job(directory,g,seed,method,verify_widths=True):
    folder=directory/'jobs'/f"{g['id']}_{seed}_{method}";assert (folder/'COMPLETE').exists() and (folder/'EXIT_CODE').read_text()=='0';check_hashes(folder,json.loads((folder/'SHA256SUMS.json').read_text()));assert json.loads((folder/'job.json').read_text())['frozen_sha256']==sha(directory/'FROZEN.json')
    ff=rows(folder/'fields.csv');job_caps=json.loads((folder/'job.json').read_text())['caps'];assert len(ff)==len(job_caps) and [int(r['cap']) for r in ff]==job_caps;abi=json.loads((folder/'abi.json').read_text());assert abi['configuration']=='Release' and abi['sample_bytes']==32
    ref=reference_dir(directory,g);samples={s:np.fromfile(ref/(s+'.bin'),dtype='<f8').reshape(-1,4) for s in ('uniform','near','audit','probes')};groups=np.loadtxt(ref/'probe_groups.txt',dtype=int);rays=np.loadtxt(ref/'rays.txt');reference=[np.fromstring(r['roots'],sep=';') for r in rows(ref/'reference_roots.csv')]
    native_metrics=rows(folder/'metrics.csv');native_groups=rows(folder/'groups.csv');native_widths=rows(folder/'widths.csv');native_topology=rows(folder/'topology.csv');summaries=[];queries=0;prediction_cache={};root_cache={};fit_prediction_cache={}
    allocation_fit,policy=fitting_evidence(folder,g,method)
    probe_max_gate=json.loads((folder/'policy.json').read_text()).get('probe_max_gate',False)
    for r in ff:
        cap=int(r['cap']);data=(folder/f'{cap}.bin').read_bytes();fieldhash=hashlib.sha256(data).hexdigest();n=int(r['grid_n']);assert int(r['bytes'])==len(data)+abi['grid_object_bytes' if method=='grid' else 'packed_object_bytes']<=cap
        if method=='grid':assert len(data)==n**3*4
        else:layout(data,method,7 if policy and policy.get('k7',0) else 6 if method in ('selective32','candidate','candidate_original','coverage_q8') else 5)
        query=lambda points:decode_grid(data,n,points) if method=='grid' else decode(data,points)
        if allocation_fit is not None:
            if fieldhash not in fit_prediction_cache:fit_prediction_cache[fieldhash]=query(allocation_fit[:,:3])
            prediction=fit_prediction_cache[fieldhash];error=np.abs(prediction-allocation_fit[:,3])
            if method=='cvar99':thresholds=rows(folder/f'{cap}_thresholds.csv');chosen=min(thresholds,key=lambda t:(float(t['objective']),int(t['bytes'])));objective=np.maximum(error-float(chosen['threshold']),0).sum()
            else:
                penalty=policy['lambda'] if method in ('candidate','candidate_original') else .0025 if method=='sign' else 0;valid=np.abs(allocation_fit[:,3])>2e-9;wrong=valid&(prediction*np.sign(allocation_fit[:,3])<=0);objective=np.sum((error/2)**2)+penalty*wrong.sum()
            close(objective,r['objective']);queries+=len(prediction)
        samples['boundary']=np.fromfile(folder/f'{cap}_boundary.bin',dtype='<f8').reshape(-1,4);assert len(samples['boundary'])==int(r['boundary_oracle_calls'])<=4096
        mm={};preds={};scalar=True
        for ds,s in samples.items():
            key=(fieldhash,ds)
            if key not in prediction_cache:prediction_cache[key]=(s.copy() if ds=='boundary' else None,query(s[:,:3]))
            original,pred=prediction_cache[key]
            if original is not None:np.testing.assert_array_equal(original,s)
            cpp=np.fromfile(folder/f'{cap}_{ds}.pred',dtype='<f8');np.testing.assert_allclose(pred,cpp,rtol=2e-12,atol=2e-14);m=metrics(pred,s[:,3]);native=next(x for x in native_metrics if int(x['cap'])==cap and x['dataset']==ds)
            for key,value in m.items():close(value,native[key])
            if g['shape']!='mesh':np.testing.assert_allclose(s[:,3],exact(g,s[:,:3]),rtol=2e-11,atol=2e-13)
            if ds!='probes' or probe_max_gate:scalar &= m['max_error']<=.020 and (ds!='uniform' or m['cvar99']<=.015) and (ds!='near' or m['cvar99']<=.010)
            mm[ds]=m;preds[ds]=pred;queries+=len(pred)
        probe=True
        for group in g['feature_groups']:
            mask=groups==group;m=metrics(preds['probes'][mask],samples['probes'][mask,3]);native=next(x for x in native_groups if int(x['cap'])==cap and int(x['group'])==group);passed=m['inside_recall']>=.95 and m['outside_recall']>=.95
            for key in ('inside_recall','outside_recall'):close(m[key],native[key])
            assert passed==bool(int(native['pass']));probe &= passed
        if fieldhash not in root_cache:root_cache[fieldhash]=independent_roots(query,rays)
        width=True;ray_signs=True;errors={group:[] for group in g['feature_groups']};root_results=root_cache[fieldhash]
        for i,(v,cont) in enumerate(root_results):
            native=next(x for x in native_widths if int(x['cap'])==cap and int(x['ray'])==i);cpp=np.fromstring(native['predicted_roots'],sep=';');np.testing.assert_allclose(v,cpp,rtol=2e-9,atol=2e-9);assert cont==bool(int(native['continuous']));truth=reference[i];ok=len(v)==len(truth) and len(v)>0 and len(v)%2==0
            assert ok==bool(int(native['interval_count_ok']));width &= ok
            if ok:
                e=np.abs(np.diff(v)/np.diff(truth)-1);errors[int(rays[i,0])].extend(e);close(e.mean(),native['mean_relative_error']);close(np.quantile(e,.95,method='inverted_cdf'),native['p95_relative_error'])
                ray_signs &= ray_sign_pattern(query,rays[i],truth,v)
            else:ray_signs=False
        for error in errors.values():width &= bool(error) and np.mean(error)<=.10 and np.quantile(error,.95,method='inverted_cdf')<=.15
        assert scalar==bool(int(r['scalar_pass']));assert probe==bool(int(r['probe_pass']));assert width==bool(int(r['width_pass']))
        native_width_pass=width;width &= ray_signs
        tops=[x for x in native_topology if int(x['cap'])==cap];top_pass=False
        if scalar and probe and native_width_pass:
            if len(g['qualified_levels'])==2:
                assert r['topology_status']=='evaluated' and [int(x['cells']) for x in tops]==g['qualified_levels'];top_pass=True
                for t in tops:
                    closed,euler=mesh_topology(*obj(folder/f"{method}_{cap}_{t['cells']}.obj"));assert closed==bool(int(t['closed']));assert euler==list(map(int,t['component_euler'].split(';')));top_pass &= closed and euler==g['expected_euler']
            else:assert r['topology_status']=='oracle_resolution_unresolved' and not tops
        else:assert r['topology_status']=='not_evaluated_due_to_prior_failure' and not tops
        qualified=len(g['qualified_levels'])==2;passed=strict_gate(scalar,probe,width,top_pass,qualified)
        failures=[name for name,ok in [('scalar',scalar),('features',probe),('width',width),('reference_qualified',qualified)] if not ok]
        if r['topology_status']=='evaluated' and not top_pass:failures.append('observed_topology_failure')
        elif r['topology_status']=='not_evaluated_due_to_prior_failure':failures.append('topology_not_evaluated_due_to_prior_failure')
        wrmean=max((float(np.mean(e)) for e in errors.values() if e),default=None);wrp95=max((float(np.quantile(e,.95,method='inverted_cdf')) for e in errors.values() if e),default=None)
        summary=dict(geometry=g['id'],stratum=g['stratum'],seed=seed,method=method,cap=cap,bytes=int(r['bytes']),scalar_pass=bool(scalar),probe_pass=bool(probe),width_pass=bool(width),topology_pass=bool(top_pass),reference_qualified=qualified,failures=';'.join(failures),uniform_mae=mm['uniform']['mae'],uniform_cvar99=mm['uniform']['cvar99'],uniform_max=mm['uniform']['max_error'],near_cvar99=mm['near']['cvar99'],near_max=mm['near']['max_error'],audit_max=mm['audit']['max_error'],boundary_max=mm['boundary']['max_error'],width_worst_group_mean=wrmean,width_worst_group_p95=wrp95,build_tree_seconds=float(r['build_tree_seconds']),loss_seconds=float(r['loss_seconds']),dp_seconds=float(r['dp_seconds']),pack_seconds=float(r['pack_seconds']),oracle_calls=int(r['oracle_calls']),boundary_oracle_calls=int(r['boundary_oracle_calls']))
        summary['pass']=bool(passed)
        summary['topology_status']=r['topology_status'];summary['topology_observed_failure']=r['topology_status']=='evaluated' and not top_pass
        summary['native_interval_width_pass']=bool(native_width_pass);summary['ray_signed_interval_pass']=bool(ray_signs)
        summaries.append(summary)
    # Fit separation and packed construction accounting, including optional quadrature.
    if method!='grid':
        fit=np.fromfile(folder/'fit.bin',dtype='<f8').reshape(-1,4);assert len(fit)==10000
        if g['shape']!='mesh':np.testing.assert_allclose(fit[:,3],exact(g,fit[:,:3]),atol=2e-13,rtol=2e-11)
        for s in (samples['uniform'],samples['near']):assert not np.array_equal(fit[:10,:3],s[:10,:3])
        for r in ff:
            policy=json.loads((folder/'policy.json').read_text());candidate=method in ('candidate','candidate_original','coverage_q8');selected=int(r['selected_parents']);selected7=int(r.get('selected_d6_parents',0));nodes=int(r['nodes']);extra=64*(selected+selected7) if candidate else 0;global_count=32768*policy['global_quadrature']**3 if candidate else 0;assert nodes==37449+8*(selected+selected7);assert selected<=(policy['k'] if candidate else 32 if method=='selective32' else 0);assert selected7<=policy.get('k7',0);assert int(r['oracle_calls'])==8*nodes+10000+extra+global_count;assert int(r['fit_count'])==10000+extra+global_count
    else:
        for r in ff:assert int(r['oracle_calls'])==int(r['grid_n'])**3
    return summaries,queries
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);ap.add_argument('--out',type=Path);ap.add_argument('--only',nargs='*');ap.add_argument('--geometry',nargs='*');ap.add_argument('--wait-for-jobs',action='store_true');ap.add_argument('--merge-shards',type=Path);a=ap.parse_args();directory=a.directory.resolve();frozen=validate_freeze(directory,check_source=False);out=(a.out or directory/'verification').resolve();out.mkdir(parents=True,exist_ok=True);all_rows=[];total_queries=0;jobs=expected_jobs(frozen);evaluator=evaluator_identity();frozen_hash=sha(directory/'FROZEN.json');partial=bool(a.only or a.geometry)
    if not partial:assert (directory/'RUN_COMPLETE.json').exists()
    if a.merge_shards:
        for cache in a.merge_shards.rglob('*.json'):
            if not cache.name.endswith(('_grid.json','_mse.json','_sign.json','_cvar99.json','_selective32.json','_candidate_original.json','_coverage_q8.json','_candidate.json')):continue
            value=json.loads(cache.read_text());identity=value['identity'];job=directory/'jobs'/cache.stem
            assert identity==dict(frozen_sha256=frozen_hash,job_manifest_sha256=sha(job/'SHA256SUMS.json'),evaluator=evaluator),'Shard evaluator/input identity mismatch';shutil.copy2(cache,out/cache.name)
    for g in frozen['geometries']:
        if a.geometry and g['id'] not in a.geometry:continue
        if a.only and not any(name.startswith(g['id']+'_') for name in a.only):continue
        ref=reference_dir(directory,g)
        if g['shape']=='mesh':
            sample=np.concatenate([np.fromfile(ref/(s+'.bin'),dtype='<f8').reshape(-1,4)[np.linspace(0,99999,64,dtype=int)] for s in ('uniform','near')]);np.testing.assert_allclose(brute_mesh(g,sample[:,:3]),sample[:,3],atol=2e-10,rtol=2e-9);json_write(out/f"{g['id']}_independent_mesh.json",dict(points=128,method='brute-force point triangle distances and independent ray parity sign',passed=True))
    for i,(g,s,m) in enumerate(jobs):
        name=f"{g['id']}_{s}_{m}"
        if a.only and name not in a.only:continue
        if a.geometry and g['id'] not in a.geometry:continue
        job=directory/'jobs'/name
        while a.wait_for_jobs and not (job/'SHA256SUMS.json').exists():time.sleep(1)
        identity=dict(frozen_sha256=frozen_hash,job_manifest_sha256=sha(job/'SHA256SUMS.json'),evaluator=evaluator)
        cache=out/(name+'.json');
        if cache.exists():
            v=json.loads(cache.read_text());assert v['identity']==identity,'Cached evaluation identity mismatch';rr,q=v['rows'],v['queries']
        else:
            start=time.time();rr,q=verify_job(directory,g,s,m);json_write(cache,dict(rows=rr,queries=q,identity=identity));print(f'{i+1}/{len(jobs)} independently verified {name} {time.time()-start:.1f}s',flush=True)
        all_rows.extend(rr);total_queries+=q
    csv_write(out/'matrix.csv',all_rows)
    if not partial:
        assert len(all_rows)==len(jobs)*len(CAPS);assert {(r['geometry'],r['seed'],r['method'],r['cap']) for r in all_rows}=={(g['id'],s,m,c) for g,s,m in jobs for c in CAPS}
        json_write(out/'VERIFIED.json',dict(fields=len(all_rows),jobs=len(jobs),queries=total_queries,frozen_sha256=frozen_hash,evaluator=evaluator,verifier_sha256=sha(Path(__file__)),all_gates_independently_recomputed=True));json_write(out/'SHA256SUMS.json',manifest(out))
if __name__=='__main__':main()
