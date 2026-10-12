"""Frozen held-out jobs: one owned heavy process, no persistent oracle cache."""
import argparse, ctypes, shutil, subprocess, time
from ctypes import wintypes
from generalization_common import *

FIT_SEEDS=list(range(340,345))
COMPACT={'compact_f64':1,'compact_f32':2,'compact_Q24':3}
GUARD=1610612736

def expected_jobs(frozen):
    return [(g,s,m) for g in frozen['geometries'] for m in frozen['methods'] for s in ([0] if m=='grid' or m in COMPACT else frozen['fit_seeds'])]

def validate_freeze(directory,check_source=True):
    f=json.loads((directory/'FROZEN.json').read_text())
    assert f['protocol']=='generalization-v2-compact' and f['fit_seeds']==FIT_SEEDS and f['caps']==CAPS
    assert f['methods']==METHODS+list(COMPACT) and f['selected_codec']=='Q24' and f['process_guard_bytes']==GUARD
    check_hashes(directory/'source',f['source_sha256'])
    if check_source:check_hashes(ROOT,f['source_sha256'])
    for folder,key in [('reference','reference_sha256'),('inputs','input_sha256'),('oracle-audit','oracle_audit_sha256')]:check_hashes(directory/folder,f[key])
    for name,h in f['executable_sha256'].items():assert sha(directory/'source'/name)==h
    assert sha(directory/'FINAL_POLICY.txt')==f['final_policy_sha256']
    return f

def freeze(directory):
    assert (directory/'REFERENCES_COMPLETE.json').exists() and not (directory/'FROZEN.json').exists()
    assert not (directory/'jobs').exists(),'Freeze must precede all held-out fields'
    selection=ROOT/'runs/generalization-v2/compact-development-v2/SELECTION.json';assert json.loads(selection.read_text())['winner']=='Q24'
    target=directory/'source';target.mkdir()
    files=[ROOT/'CMakeLists.txt',ROOT/'docs/generalization-v2-protocol.md',ROOT/'docs/generalization-v2-compact-plan.md',ROOT/'docs/generalization-v2-final-protocol.md']
    for folder in ('src','include','scripts','tests'):
        files.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and p.suffix in ('.cpp','.hpp','.py'))
    identity={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)}
    for name in identity:
        dest=target/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
    executables={}
    release=ROOT/'build-generalization-v2-cgal/Release'
    for name in ('sdf_generalization_v2_baseline.exe','sdf_compact_bench.exe'):
        shutil.copy2(release/name,target/name);executables[name]=sha(target/name)
    for p in release.glob('*.dll'):shutil.copy2(p,target/p.name);executables[p.name]=sha(target/p.name)
    shutil.copy2(ROOT/'runs/generalization-v1/FINAL_POLICY.txt',directory/'FINAL_POLICY.txt')
    f=dict(protocol='generalization-v2-compact',geometries=json.loads((directory/'qualification.json').read_text()),methods=METHODS+list(COMPACT),caps=CAPS,fit_seeds=FIT_SEEDS,selected_codec='Q24',selected_development_sha256=sha(selection),process_guard_bytes=GUARD,source_sha256=identity,executable_sha256=executables,reference_sha256=manifest(directory/'reference'),input_sha256=manifest(directory/'inputs'),oracle_audit_sha256=manifest(directory/'oracle-audit'),final_policy_sha256=sha(directory/'FINAL_POLICY.txt'),preparation_declaration_sha256=sha(directory/'PREPARATION-DECLARATION.json'),reuse='grid once per geometry/cap; each compact codec once per geometry; seven adaptive methods five fits',oracle_cache='no persistent cache; baseline direct oracle; compact owns decoded geometry only',expected_physical_jobs=156,expected_logical_cap_rows=2028,expected_unique_encoded_fields=1884,timing_scope='serial native processes; nine compact batches of 4096 versus baseline batches of 100000, no unmatched speed ratio')
    assert len(f['geometries'])==4 and len(expected_jobs(f))==156
    json_write(directory/'FROZEN.json',f);validate_freeze(directory)
    print('FROZEN',directory/'FROZEN.json',sha(directory/'FROZEN.json'),flush=True)

class Counters(ctypes.Structure):
    _fields_=[('cb',wintypes.DWORD),('PageFaultCount',wintypes.DWORD)]+[(n,ctypes.c_size_t) for n in ('PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage')]

def owned_process(command,log):
    start=time.time();peak=0;resource=False
    with log.open('w',encoding='utf-8') as stream:
        p=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT)
        while p.poll() is None:
            if hasattr(p,'_handle'):
                counter=Counters();counter.cb=ctypes.sizeof(counter)
                if ctypes.windll.psapi.GetProcessMemoryInfo(wintypes.HANDLE(int(p._handle)),ctypes.byref(counter),counter.cb):
                    peak=max(peak,int(counter.PeakWorkingSetSize),int(counter.PeakPagefileUsage))
                    if peak>GUARD:resource=True;p.terminate();p.wait();break
            if shutil.disk_usage(ROOT).free<30*1024**3:resource=True;p.terminate();p.wait();break
            time.sleep(.5)
    return dict(exit_code=p.returncode,resource_unresolved=resource,observed_peak_committed_or_working_bytes=peak,seconds=time.time()-start,guard_bytes=GUARD)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);ap.add_argument('--freeze',action='store_true');ap.add_argument('--only',nargs='*');a=ap.parse_args();d=a.directory.resolve()
    if a.freeze:freeze(d);return
    f=validate_freeze(d);fh=sha(d/'FROZEN.json');jobs=expected_jobs(f);(d/'jobs').mkdir(exist_ok=True)
    for i,(g,s,m) in enumerate(jobs):
        name=f"{g['id']}_{s}_{m}"
        if a.only and name not in a.only:continue
        folder=d/'jobs'/name;exe=d/'source'/('sdf_compact_bench.exe' if m in COMPACT else 'sdf_generalization_v2_baseline.exe');identity=dict(frozen_sha256=fh,geometry=g['id'],seed=s,method=m,caps=CAPS,executable_sha256=sha(exe))
        if (folder/'SHA256SUMS.json').exists():
            assert json.loads((folder/'job.json').read_text())==identity;check_hashes(folder,json.loads((folder/'SHA256SUMS.json').read_text()));continue
        assert not folder.exists(),'Incomplete attempt retained; investigate before resume'
        ref=reference_dir(d,g)
        command=[str(exe),str(ref/'geometry.txt'),str(ref),str(folder),str(COMPACT[m])] if m in COMPACT else [str(exe),'fit',str(ref/'geometry.txt'),str(ref),str(folder),str(s),m,*map(str,CAPS),'--probe-max']
        if m=='candidate':command+=['--policy',str(d/'FINAL_POLICY.txt')]
        result=owned_process(command,d/'jobs'/(name+'.log'));folder.mkdir(exist_ok=True);json_write(folder/'job.json',identity);json_write(folder/'PROCESS.json',result);(folder/'EXIT_CODE').write_text(str(result['exit_code']),encoding='utf-8')
        if result['exit_code']:
            status='resource_unresolved' if result['resource_unresolved'] else 'incomplete';(folder/'STATUS').write_text(status+'\n',encoding='utf-8')
            json_write(folder/'SHA256SUMS.json',manifest(folder));raise RuntimeError(f'{name}: {status}; preserve and checkpoint implementation before continuing')
        assert (folder/'COMPLETE').exists();json_write(folder/'SHA256SUMS.json',manifest(folder));print(f'{i+1}/{len(jobs)} {name} {result["seconds"]:.1f}s',flush=True)
    if not a.only:
        assert all((d/'jobs'/f"{g['id']}_{s}_{m}"/'SHA256SUMS.json').exists() for g,s,m in jobs)
        json_write(d/'RUN_COMPLETE.tmp',dict(jobs=len(jobs),frozen_sha256=fh));(d/'RUN_COMPLETE.tmp').replace(d/'RUN_COMPLETE.json')
if __name__=='__main__':main()
