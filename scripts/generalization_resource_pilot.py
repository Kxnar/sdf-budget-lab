"""Development-only1MiB frontier pilot with OS process peak memory evidence."""
import argparse,subprocess,time,ctypes
from ctypes import wintypes
from generalization_common import *
class Counters(ctypes.Structure):
    _fields_=[('cb',wintypes.DWORD),('PageFaultCount',wintypes.DWORD),*[(x,ctypes.c_size_t) for x in ('PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage','QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage','PrivateUsage')]]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('out',type=Path);a=ap.parse_args();out=a.out.resolve();assert not out.exists();out.mkdir();g=json.loads((ROOT/'runs/generalization-development-screen/reference/development_0_torus/metadata.json').read_text());ref=ROOT/'runs/generalization-development-screen/reference'/g['id'];exe=ROOT/'build-generalization-core/Release/sdf_generalization_bench.exe';policy=ROOT/'runs/generalization-development-d7/confirmation/d7_2048/policy.txt';get=ctypes.WinDLL('psapi').GetProcessMemoryInfo;get.argtypes=[wintypes.HANDLE,ctypes.POINTER(Counters),wintypes.DWORD];get.restype=wintypes.BOOL;results=[]
    for method in ('coverage_q8','candidate','cvar99'):
        folder=out/method;command=[str(exe),'fit',str(ref/'geometry.txt'),str(ref),str(folder),'40',method,*map(str,CAPS),'--probe-max','--cache',str(out/'cache.bin'),'--cache-key',cache_identity(g,ref/'geometry.txt')]
        if method=='candidate':command+=['--policy',str(policy)]
        start=time.time();peak=0
        with (out/(method+'.log')).open('w') as log:
            process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
            while process.poll() is None:
                counters=Counters();counters.cb=ctypes.sizeof(counters)
                if get(int(process._handle),ctypes.byref(counters),ctypes.sizeof(counters)):peak=max(peak,counters.PeakWorkingSetSize)
                time.sleep(.1)
        assert process.returncode==0;fields=rows(folder/'fields.csv');last=fields[-1];result=dict(method=method,elapsed_seconds=time.time()-start,peak_process_working_set_bytes=peak,maximum_cap=int(last['cap']),candidate_nodes=int(last['nodes']),dp_states=int(last['dp_states']),dp_transitions=int(last['dp_transitions']),dp_seconds=float(last['dp_seconds']),tree_seconds=float(last['build_tree_seconds']),loss_seconds=float(last['loss_seconds']),cache=json.loads((folder/'cache.json').read_text()),scope='OS process peak includes oracle, cache, labels, candidate tree, loss and DP scratch; excluded from retained field bytes');results.append(result);print(method,result['elapsed_seconds'],'s',peak/1024**2,'MiB peak',flush=True)
    json_write(out/'RESULT.json',dict(development_geometry=g,seed=40,caps=CAPS,executable_sha256=sha(exe),results=results));json_write(out/'SHA256SUMS.json',manifest(out))
if __name__=='__main__':main()
