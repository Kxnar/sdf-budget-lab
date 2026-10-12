"""Procedural v2 amendment: strict completion, immutable freeze, job-boundary stop."""
import argparse,sys,time
from generalization_common import *
import generalization_v2_run as original

def validate(d):
    f=original.validate_freeze(d);a=json.loads((d/'PROCEDURAL-AMENDMENT-1.json').read_text())
    assert a['frozen_sha256']==sha(d/'FROZEN.json') and a['runner_sha256']==sha(Path(__file__))
    assert sha(d/'PREPARATION-DECLARATION.json')==f['preparation_declaration_sha256']
    assert sha(ROOT/'runs/generalization-v2/compact-development-v2/SELECTION.json')==f['selected_development_sha256']
    for g,s,m in original.expected_jobs(f):
        folder=d/'jobs'/f"{g['id']}_{s}_{m}"
        if not (folder/'SHA256SUMS.json').exists():continue
        assert (folder/'COMPLETE').exists() and (folder/'EXIT_CODE').read_text()=='0',f'Explicit nonpass/incomplete checkpoint requires investigation: {folder}'
        identity=json.loads((folder/'job.json').read_text());assert identity['frozen_sha256']==sha(d/'FROZEN.json');check_hashes(folder,json.loads((folder/'SHA256SUMS.json').read_text()))
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);args=ap.parse_args();d=args.directory.resolve();f=validate(d);start=time.time()
    for g,s,m in original.expected_jobs(f):
        name=f"{g['id']}_{s}_{m}";folder=d/'jobs'/name
        if (folder/'SHA256SUMS.json').exists():continue
        if (d/'STOP_AFTER_CURRENT_JOB').exists():print('Stopped at requested native-job boundary',flush=True);return
        sys.argv=[str(Path(original.__file__)),str(d),'--only',name];original.main()
        assert (folder/'COMPLETE').exists() and (folder/'EXIT_CODE').read_text()=='0';check_hashes(folder,json.loads((folder/'SHA256SUMS.json').read_text()))
    validate(d);json_write(d/'RUN_COMPLETE.tmp',dict(jobs=len(original.expected_jobs(f)),frozen_sha256=sha(d/'FROZEN.json'),procedural_amendment_sha256=sha(d/'PROCEDURAL-AMENDMENT-1.json'),elapsed_resumed_seconds=time.time()-start,all_complete_zero_exit=True));(d/'RUN_COMPLETE.tmp').replace(d/'RUN_COMPLETE.json')
if __name__=='__main__':main()
