"""Fresh-process frozen-source reproduction; no study verifier runs concurrently."""
import argparse,subprocess,shutil
from generalization_common import *
from generalization_v2_run import validate_freeze,COMPACT,owned_process

def selection(d,f):
    matrix=rows(d/'verification/matrix.csv');result=[]
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
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);a=ap.parse_args();d=a.directory.resolve();f=validate_freeze(d,False);out=d/'reproduction';out.mkdir(exist_ok=True);records=[]
    for i,(g,s,m,caps) in enumerate(selection(d,f)):
        ref=reference_dir(d,g);original=d/'jobs'/f"{g['id']}_{s}_{m}";target=out/f"{g['id']}_{s}_{m}";exe=d/'rebuilt/Release'/('sdf_compact_bench.exe' if m in COMPACT else 'sdf_generalization_v2_baseline.exe')
        assert exe.exists(),'Rebuild frozen source separately before reproduction'
        if not (target/'COMPLETE').exists():
            assert not target.exists();command=[str(exe),str(ref/'geometry.txt'),str(ref),str(target),str(COMPACT[m])] if m in COMPACT else [str(exe),'fit',str(ref/'geometry.txt'),str(ref),str(target),str(s),m,*map(str,caps),'--probe-max']
            if m=='candidate':command+=['--policy',str(d/'FINAL_POLICY.txt')]
            process=owned_process(command,out/(target.name+'.log'));assert process['exit_code']==0;json_write(target/'PROCESS.json',process)
        files=['field.bin']+[x+'_pred.bin' for x in ('uniform','near','audit','probes','boundary')] if m in COMPACT else [name for cap in caps for name in [f'{cap}.bin',*[f'{cap}_{s}.pred' for s in ('uniform','near','audit','probes','boundary')]]]
        for name in files:assert sha(target/name)==sha(original/name),f'Reproduction mismatch {target/name}'
        records.append(dict(geometry=g['id'],seed=s,method=m,caps=caps,exact_field_and_predictions=True,executable_sha256=sha(exe),timing_scope='fresh process, study verification stopped; background load uncontrolled',hashes={name:sha(target/name) for name in files}));print(f'{i+1}/44 reproduced {target.name}',flush=True)
    assert len(records)==44;json_write(d/'REPRODUCED.json',dict(frozen_sha256=sha(d/'FROZEN.json'),jobs=44,records=records));json_write(out/'SHA256SUMS.json',manifest(out))
if __name__=='__main__':main()
