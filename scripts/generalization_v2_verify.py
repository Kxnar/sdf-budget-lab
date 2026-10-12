"""Independent full packed/grid decode, strict gates and Cartesian completeness."""
import argparse,time,sys,shutil
from generalization_common import *
from generalization_v2_run import validate_freeze,expected_jobs,COMPACT
from generalization_v2_rays import packed_segments,transitions
from compact_verify import verify as verify_compact
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
    names=['generalization_v2_verify.py','generalization_common.py','generalization_v2_run.py','generalization_v2_rays.py','generalization_prepare.py','distortion_analysis.py','compact_mesh.py','compact_verify.py'];return dict(source_sha256={n:sha(Path(__file__).with_name(n)) for n in names},numpy_version=np.__version__,python_version=sys.version)
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
        if fieldhash not in root_cache:
            results=[]
            for ray in rays:
                if method=='grid':
                    times=[ray[7],ray[8]]
                    for axis in range(3):
                        if ray[4+axis]:times.extend(t for t in ((np.arange(1,n-1)*2/(n-1)-1-ray[1+axis])/ray[4+axis]) if ray[7]<t<ray[8])
                    times=np.unique(times);segments=[];values=np.frombuffer(data,dtype='<f4').reshape(n,n,n)
                    for lo,hi in zip(times[:-1],times[1:]):
                        point=ray[1:4]+ray[4:7]*(lo+hi)/2;index=np.clip(np.floor((point+1)*(n-1)/2).astype(int),0,n-2);low=-1+index*2/(n-1);high=-1+(index+1)*2/(n-1);corners=[values[index[2]+bool(c&4),index[1]+bool(c&2),index[0]+bool(c&1)] for c in range(8)];segments.append([lo,hi,*low,*high,*corners])
                    segments=np.array(segments)
                else:segments=packed_segments(data,ray)
                results.append(transitions(segments,ray))
            root_cache[fieldhash]=results
        width=True;ray_signs=True;numerical=0;errors={group:[] for group in g['feature_groups']}
        for i,result in enumerate(root_cache[fieldhash]):
            native=next(x for x in native_widths if int(x['cap'])==cap and int(x['ray'])==i);cpp=np.fromstring(native['predicted_roots'],sep=';');native_unresolved=bool(int(native['numerical_unresolved']));ambiguous=result['ambiguous'] or native_unresolved;numerical+=ambiguous
            v=cpp if native_unresolved else np.array(result['roots'])
            if not ambiguous:np.testing.assert_allclose(v,cpp,rtol=2e-8,atol=3e-9)
            truth=reference[i];ok=len(v)==len(truth) and len(v)>0 and len(v)%2==0
            signs=ray_sign_pattern(query,rays[i],truth,v)
            if not ambiguous:
                assert ok==bool(int(native['interval_count_ok']));assert signs==bool(int(native['signed_intervals_ok']))
            width &= ok and signs and not ambiguous;ray_signs &= signs
            if ok:
                e=np.abs(np.diff(v)/np.diff(truth)-1);errors[int(rays[i,0])].extend(e)
                if not ambiguous:
                    native_e=np.abs(np.diff(cpp)/np.diff(truth)-1);close(native_e.mean(),native['mean_relative_error']);close(np.quantile(native_e,.95,method='inverted_cdf'),native['p95_relative_error'])
        for error in errors.values():width &= bool(error) and np.mean(error)<=.10 and np.quantile(error,.95,method='inverted_cdf')<=.15
        assert scalar==bool(int(r['scalar_pass']));assert probe==bool(int(r['probe_pass']))
        if not numerical:assert width==bool(int(r['width_pass']))
        native_width_pass=bool(int(r['width_pass']))
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
        summary['native_interval_width_pass']=bool(native_width_pass);summary['ray_signed_interval_pass']=bool(ray_signs);summary['numerical_unresolved_rays']=int(numerical);summary['width_status']='numerical_unresolved' if numerical else 'passed' if width else 'observed_geometric_failure';summary['status']='reference_unresolved' if not qualified else 'numerical_unresolved' if numerical else 'passed' if passed else 'observed_gate_failure'
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
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);ap.add_argument('--only',nargs='*');a=ap.parse_args();d=a.directory.resolve();f=validate_freeze(d,check_source=False);out=d/'verification';out.mkdir(exist_ok=True);all_rows=[];total_queries=0;identity_evaluator=evaluator_identity();jobs=expected_jobs(f)
    for i,(g,s,m) in enumerate(jobs):
        name=f"{g['id']}_{s}_{m}"
        if a.only and name not in a.only:continue
        job=d/'jobs'/name;assert (job/'SHA256SUMS.json').exists();check_hashes(job,json.loads((job/'SHA256SUMS.json').read_text()));identity=dict(frozen_sha256=sha(d/'FROZEN.json'),job_manifest_sha256=sha(job/'SHA256SUMS.json'),evaluator=identity_evaluator);cache=out/(name+'.json')
        if cache.exists():
            value=json.loads(cache.read_text());assert value['identity']==identity;rr=value['rows'];q=value['queries']
        elif m in COMPACT:
            assessment=verify_compact(job,reference_dir(d,g),g);rr=[];q=assessment.get('independent_decoded_queries',0)
            for cap in CAPS:
                size=assessment['actual_bytes'];feasible=size is not None and size<=cap
                rr.append(dict(geometry=g['id'],stratum=g['stratum'],seed=s,method=m,cap=cap,bytes=size,storage_feasible=feasible,status=assessment['status'] if feasible else 'storage_infeasible' if size is not None else assessment['status'],representation_status=assessment['status'],scalar_pass=assessment.get('scalar_pass',False),probe_pass=assessment.get('feature_pass',False),width_pass=assessment.get('width_pass',False),topology_pass=assessment.get('topology_pass',False),reference_qualified=len(g['qualified_levels'])==2,**{'pass':bool(feasible and assessment['joint_pass'])}))
            json_write(cache,dict(identity=identity,rows=rr,queries=q,assessment=assessment))
        else:
            rr,q=verify_job(d,g,s,m);json_write(cache,dict(identity=identity,rows=rr,queries=q))
        all_rows.extend(rr);total_queries+=q;print(f'{i+1}/{len(jobs)} verified {name}',flush=True)
        # Independent evidence is additive; retain the native immutable manifest.
    keys=list(dict.fromkeys(k for r in all_rows for k in r));csv_write(out/'matrix.partial.csv',[{k:r.get(k,'') for k in keys} for r in all_rows])
    if not a.only:
        assert len(all_rows)==2028 and {(r['geometry'],r['seed'],r['method'],r['cap']) for r in all_rows}=={(g['id'],s,m,c) for g,s,m in jobs for c in CAPS}
        csv_write(out/'matrix.csv',[{k:r.get(k,'') for k in keys} for r in all_rows]);json_write(out/'VERIFIED.json',dict(jobs=156,logical_cap_rows=2028,baseline_decoded_queries=total_queries,frozen_sha256=sha(d/'FROZEN.json'),evaluator=identity_evaluator,complete=True));json_write(out/'SHA256SUMS.json',manifest(out))
if __name__=='__main__':main()
