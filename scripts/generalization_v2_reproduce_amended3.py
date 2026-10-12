"""Fresh-process frozen-source reproduction; no study verifier runs concurrently."""
import argparse,subprocess,shutil
from generalization_common import *
from generalization_v2_run import validate_freeze,COMPACT,owned_process

def deterministic_rows(records):
    return [{k:v for k,v in r.items() if k not in ('seconds','full_query_seconds','accelerated_seconds')} for r in records]

def same_acceptance(actual,expected):
    assert deterministic_rows(actual)==deterministic_rows(expected),'Deterministic acceptance evidence mismatch'

def verified_identity(d):
    value=json.loads((d/'verification/VERIFIED.json').read_text());assert value['frozen_sha256']==sha(d/'FROZEN.json');check_hashes(d/'verification',json.loads((d/'verification/SHA256SUMS.json').read_text()))
    from generalization_v2_verify import evaluator_identity
    assert value['evaluator']==evaluator_identity();return value

def selection(d,f):
    v=verified_identity(d);matrix=rows(d/'verification/matrix.csv');result=[]
    for g in f['geometries']:
        for m in f['methods']:
            s=0 if m=='grid' or m in COMPACT else 340
            rr=[r for r in matrix if r['geometry']==g['id'] and r['method']==m and int(r['seed'])==s]
            passing=[r for r in rr if r['pass']=='True']
            caps={int(min(passing,key=lambda r:(int(r['bytes']),int(r['cap'])))['cap'])} if passing else {CAPS[-1]}
            if passing:
                cap=next(iter(caps));i=CAPS.index(cap)
                if i:caps.add(CAPS[i-1])
            result.append((g,s,m,sorted(caps)))
    # Preserve the44 declared jobs, then include every contributing adaptive fit.
    for g in f['geometries']:
        for m in METHODS:
            if m=='grid':continue
            for s in f['fit_seeds'][1:]:
                rr=[r for r in matrix if r['geometry']==g['id'] and r['method']==m and int(r['seed'])==s and r['pass']=='True']
                if not rr:continue
                chosen=min(rr,key=lambda r:(int(r['bytes']),int(r['cap'])));cap=int(chosen['cap']);caps={cap};index=CAPS.index(cap)
                if index:caps.add(CAPS[index-1])
                result.append((g,s,m,sorted(caps)))
    return result

def baseline_ray_evidence(folder,ref,cap,method):
    from generalization_v2_rays import packed_segments,transitions
    from generalization_common import decode_packed,decode_grid
    data=(folder/f'{cap}.bin').read_bytes();rays=np.loadtxt(ref/'rays.txt');native=next(r for r in rows(folder/'fields.csv') if int(r['cap'])==cap);n=int(native['grid_n'])
    query=lambda points:decode_grid(data,n,points) if method=='grid' else decode_packed(data,points)
    references=[np.fromstring(r['roots'],sep=';') for r in rows(ref/'reference_roots.csv')];evidence=[]
    for ray,truth in zip(rays,references):
        if method=='grid':
            times=[ray[7],ray[8]]
            for axis in range(3):
                if ray[4+axis]:times.extend(t for t in ((np.arange(1,n-1)*2/(n-1)-1-ray[1+axis])/ray[4+axis]) if ray[7]<t<ray[8])
            times=np.unique(times);segments=[];values=np.frombuffer(data,dtype='<f4').reshape(n,n,n)
            for lo,hi in zip(times[:-1],times[1:]):
                point=ray[1:4]+ray[4:7]*(lo+hi)/2;index=np.clip(np.floor((point+1)*(n-1)/2).astype(int),0,n-2);low=-1+index*2/(n-1);high=-1+(index+1)*2/(n-1);corners=[values[index[2]+bool(c&4),index[1]+bool(c&2),index[0]+bool(c&1)] for c in range(8)];segments.append([lo,hi,*low,*high,*corners])
            segments=np.array(segments)
        else:segments=packed_segments(data,ray)
        result=transitions(segments,ray);record={k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in result.items()};witness=[]
        for rr in (truth,result['roots']):
            edges=np.r_[ray[7],rr,ray[8]];times=np.r_[ray[7],ray[8],(edges[:-1]+edges[1:])/2];witness.append(dict(times=times.tolist(),values=query(ray[1:4]+times[:,None]*ray[4:7]).tolist()))
        record['signed_witnesses']=witness;evidence.append(record)
    return evidence

def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);a=ap.parse_args();d=a.directory.resolve();f=validate_freeze(d,False);out=d/'reproduction-amended3';out.mkdir(exist_ok=True);records=[]
    for i,(g,s,m,caps) in enumerate(selection(d,f)):
        ref=reference_dir(d,g);original=d/'jobs'/f"{g['id']}_{s}_{m}";target=out/f"{g['id']}_{s}_{m}";exe=d/'rebuilt/Release'/('sdf_compact_bench.exe' if m in COMPACT else 'sdf_generalization_v2_baseline.exe')
        assert exe.exists(),'Rebuild frozen source separately before reproduction'
        if not (target/'COMPLETE').exists():
            assert not target.exists();command=[str(exe),str(ref/'geometry.txt'),str(ref),str(target),str(COMPACT[m])] if m in COMPACT else [str(exe),'fit',str(ref/'geometry.txt'),str(ref),str(target),str(s),m,*map(str,caps),'--probe-max']
            if m=='candidate':command+=['--policy',str(d/'FINAL_POLICY.txt')]
            process=owned_process(command,out/(target.name+'.log'));assert process['exit_code']==0;json_write(target/'PROCESS.json',process)
        if m in COMPACT:
            assert (target/'STATUS').read_text()==(original/'STATUS').read_text(),'Decoded rejection/resolution status mismatch'
            candidates=['field.bin','decoded.obj','roots.csv','ray_witnesses.csv']+[ds+suffix for ds in ('uniform','near','audit','probes','boundary') for suffix in ('_pred.bin','_flags.bin')]+[p.name for p in original.glob('surface_*.obj')]
            for csv_name in ('metrics.csv','topology.csv','extraction.csv'):
                if (original/csv_name).exists():same_acceptance(rows(target/csv_name),rows(original/csv_name))
            files=[name for name in candidates if (original/name).exists()]
            assert all((target/name).exists() for name in files)
        else:
            files=[name for cap in caps for name in [f'{cap}.bin',f'{cap}_boundary.bin',*[f'{cap}_{ds}.pred' for ds in ('uniform','near','audit','probes','boundary')]]]+[p.name for p in original.glob(f'{m}_*.obj') if any(f'_{cap}_' in p.name for cap in caps)]
            for csv_name in ('widths.csv','metrics.csv','groups.csv','topology.csv'):
                expected=[r for r in rows(original/csv_name) if int(r['cap']) in caps];actual=rows(target/csv_name);same_acceptance(actual,expected)
        for name in files:assert sha(target/name)==sha(original/name),f'Reproduction mismatch {target/name}'
        if m not in COMPACT:
            original_events={str(cap):baseline_ray_evidence(original,ref,cap,m) for cap in caps};replayed_events={str(cap):baseline_ray_evidence(target,ref,cap,m) for cap in caps};assert original_events==replayed_events,'Complete independent baseline events/witnesses mismatch';json_write(target/'INDEPENDENT-RAY-REPLAY.json',replayed_events);files.append('INDEPENDENT-RAY-REPLAY.json')
        records.append(dict(geometry=g['id'],seed=s,method=m,caps=caps,exact_field_and_predictions=True,executable_sha256=sha(exe),timing_scope='fresh process, study verification stopped; background load uncontrolled',hashes={name:sha(target/name) for name in files}));print(f'{i+1}/{len(selection(d,f))} reproduced {target.name}',flush=True)
    assert len(records)>=44;json_write(d/'REPRODUCED-AMENDED.json',dict(frozen_sha256=sha(d/'FROZEN.json'),jobs=len(records),records=records,helper_sha256=sha(Path(__file__)),verification_manifest_sha256=sha(d/'verification/SHA256SUMS.json')));json_write(out/'SHA256SUMS.json',manifest(out))
if __name__=='__main__':main()
