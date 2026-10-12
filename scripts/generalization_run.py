"""Serial immutable/resumable final jobs; accepts no method tuning arguments."""
import argparse, subprocess, time
from generalization_common import *
def expected_jobs(frozen):
    return [(g,0 if m=='grid' else s,m) for g in frozen['geometries'] for m in frozen.get('methods',METHODS) for s in ([0] if m=='grid' else SEEDS)]
def validate_freeze(directory,check_source=True,exe=None):
    frozen=json.loads((directory/'FROZEN.json').read_text());assert frozen['protocol']=='generalization-v1'
    assert frozen['fit_seeds']==SEEDS and frozen['caps']==CAPS and frozen['methods']==METHODS
    if check_source:check_hashes(ROOT,frozen['source_sha256'])
    if exe:assert sha(exe)==frozen['executable_sha256'],'Executable changed after freeze'
    check_hashes(directory/'reference',frozen['reference_sha256']);check_hashes(directory/'reference-fast',frozen.get('reference_fast_sha256',{}));check_hashes(directory/'oracle-audit',frozen.get('oracle_audit_sha256',{}));check_hashes(directory/'inputs',frozen['input_sha256']);assert sha(directory/'FINAL_POLICY.txt')==frozen['final_policy_sha256'];return frozen
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);ap.add_argument('--exe',type=Path);ap.add_argument('--only',nargs='*');a=ap.parse_args();directory=a.directory.resolve();exe=(a.exe or directory/'source/sdf_generalization_bench.exe').resolve();frozen=validate_freeze(directory,exe=exe);freeze_hash=sha(directory/'FROZEN.json');jobs=expected_jobs(frozen);start=time.time()
    for i,(g,seed,method) in enumerate(jobs):
        name=f"{g['id']}_{seed}_{method}"
        if a.only and name not in a.only:continue
        folder=directory/'jobs'/name;identity=dict(frozen_sha256=freeze_hash,geometry=g['id'],seed=seed,method=method,caps=CAPS,executable_sha256=sha(exe))
        if folder.exists():
            if (folder/'COMPLETE').exists() and (folder/'SHA256SUMS.json').exists() and (folder/'EXIT_CODE').read_text()=='0':
                assert json.loads((folder/'job.json').read_text())==identity,'Resume identity mismatch';check_hashes(folder,json.loads((folder/'SHA256SUMS.json').read_text()));continue
            folder.rename(folder.with_name(name+f'.incomplete-{time.time_ns()}'))
        ref=reference_dir(directory,g);folder.parent.mkdir(exist_ok=True);log=directory/'jobs'/f'{name}.log';cache=directory/'shared-cache'/(g['id']+'.bin');cachekey=cache_identity(g,ref/'geometry.txt');command=[str(exe),'fit',str(ref/'geometry.txt'),str(ref),str(folder),str(seed),method,*map(str,CAPS),'--probe-max','--cache',str(cache),'--cache-key',cachekey]
        if cache.exists():
            size=cache.stat().st_size
            with cache.open('rb') as stream:header=stream.read(65)
            if size<65 or (size-65)%32 or header!=cachekey.encode()+b'\n':
                quarantine=cache.with_name(cache.name+f'.invalid-{time.time_ns()}');cache.rename(quarantine)
                with (directory/'CACHE_RECOVERY.jsonl').open('a') as recovery:recovery.write(json.dumps(dict(geometry=g['id'],preserved_cache=str(quarantine),reason='truncated record or spec/mesh/kernel identity mismatch; rebuilt from frozen oracle'))+'\n')
        if method=='candidate':command+=['--policy',str(directory/'FINAL_POLICY.txt')]
        t=time.time()
        with log.open('w') as output:result=subprocess.run(command,stdout=output,stderr=subprocess.STDOUT)
        if not folder.exists():raise RuntimeError(log.read_text())
        json_write(folder/'job.json',identity);(folder/'EXIT_CODE').write_text(str(result.returncode));json_write(folder/'SHA256SUMS.json',manifest(folder))
        if result.returncode:raise RuntimeError(log.read_text())
        assert (folder/'COMPLETE').exists();print(f'{i+1}/{len(jobs)} {name} {time.time()-t:.1f}s elapsed {time.time()-start:.1f}s',flush=True)
    if not a.only:
        assert all((directory/'jobs'/f"{g['id']}_{s}_{m}"/'COMPLETE').exists() for g,s,m in jobs)
        json_write(directory/'RUN_COMPLETE.tmp',dict(jobs=len(jobs),frozen_sha256=freeze_hash,elapsed_seconds=time.time()-start));(directory/'RUN_COMPLETE.tmp').replace(directory/'RUN_COMPLETE.json')
if __name__=='__main__':main()
