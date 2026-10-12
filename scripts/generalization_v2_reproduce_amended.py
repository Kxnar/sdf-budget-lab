"""Fresh-process frozen-source reproduction; no study verifier runs concurrently."""
import argparse,subprocess,shutil
from generalization_common import *
from generalization_v2_run import validate_freeze,COMPACT,owned_process

def selection(d,f):
    v=json.loads((d/'verification/VERIFIED.json').read_text());assert v['frozen_sha256']==sha(d/'FROZEN.json');check_hashes(d/'verification',json.loads((d/'verification/SHA256SUMS.json').read_text()));matrix=rows(d/'verification/matrix.csv');result=[]
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

def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);a=ap.parse_args();d=a.directory.resolve();f=validate_freeze(d,False);out=d/'reproduction-amended';out.mkdir(exist_ok=True);records=[]
    for i,(g,s,m,caps) in enumerate(selection(d,f)):
        ref=reference_dir(d,g);original=d/'jobs'/f"{g['id']}_{s}_{m}";target=out/f"{g['id']}_{s}_{m}";exe=d/'rebuilt/Release'/('sdf_compact_bench.exe' if m in COMPACT else 'sdf_generalization_v2_baseline.exe')
        assert exe.exists(),'Rebuild frozen source separately before reproduction'
        if not (target/'COMPLETE').exists():
            assert not target.exists();command=[str(exe),str(ref/'geometry.txt'),str(ref),str(target),str(COMPACT[m])] if m in COMPACT else [str(exe),'fit',str(ref/'geometry.txt'),str(ref),str(target),str(s),m,*map(str,caps),'--probe-max']
            if m=='candidate':command+=['--policy',str(d/'FINAL_POLICY.txt')]
            process=owned_process(command,out/(target.name+'.log'));assert process['exit_code']==0;json_write(target/'PROCESS.json',process)
        if m in COMPACT:
            assert (target/'STATUS').read_text()==(original/'STATUS').read_text(),'Decoded rejection/resolution status mismatch'
            candidates=['field.bin','decoded.obj','roots.csv','ray_witnesses.csv','metrics.csv','topology.csv']+[ds+suffix for ds in ('uniform','near','audit','probes','boundary') for suffix in ('_pred.bin','_flags.bin')]+[p.name for p in original.glob('surface_*.obj')]
            files=[name for name in candidates if (original/name).exists()]
            assert all((target/name).exists() for name in files)
        else:
            files=[name for cap in caps for name in [f'{cap}.bin',f'{cap}_boundary.bin',*[f'{cap}_{ds}.pred' for ds in ('uniform','near','audit','probes','boundary')]]]+[p.name for p in original.glob(f'{m}_*.obj') if any(f'_{cap}_' in p.name for cap in caps)]
            for csv_name in ('widths.csv','metrics.csv','groups.csv','topology.csv'):
                expected=[r for r in rows(original/csv_name) if int(r['cap']) in caps];actual=rows(target/csv_name);assert actual==expected,f'Gate evidence mismatch {csv_name}'
        for name in files:assert sha(target/name)==sha(original/name),f'Reproduction mismatch {target/name}'
        records.append(dict(geometry=g['id'],seed=s,method=m,caps=caps,exact_field_and_predictions=True,executable_sha256=sha(exe),timing_scope='fresh process, study verification stopped; background load uncontrolled',hashes={name:sha(target/name) for name in files}));print(f'{i+1}/{len(selection(d,f))} reproduced {target.name}',flush=True)
    assert len(records)>=44;json_write(d/'REPRODUCED-AMENDED.json',dict(frozen_sha256=sha(d/'FROZEN.json'),jobs=len(records),records=records,helper_sha256=sha(Path(__file__)),verification_manifest_sha256=sha(d/'verification/SHA256SUMS.json')));json_write(out/'SHA256SUMS.json',manifest(out))
if __name__=='__main__':main()
