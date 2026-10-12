"""Rebuild frozen source separately and reproduce crossings and adjacent failures."""
import argparse, os, subprocess
from generalization_common import *
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path)
    ap.add_argument('--toolchain',type=Path,default=os.environ.get('CMAKE_TOOLCHAIN_FILE'),help='CMake toolchain file; defaults to CMAKE_TOOLCHAIN_FILE, otherwise installed CGAL')
    ap.add_argument('--vcpkg-classic',action='store_true',help='Disable vcpkg manifest mode for a classic-mode toolchain')
    a=ap.parse_args();directory=a.directory.resolve();frozen=json.loads((directory/'FROZEN.json').read_text());source=directory/'source';check_hashes(source,json.loads((source/'SHA256SUMS.json').read_text()));build=directory/'rebuild'
    configure=['cmake','-S',str(source),'-B',str(build),'-DSDF_WITH_CGAL=ON','-DCMAKE_BUILD_TYPE=Release']
    if a.toolchain is not None:configure+=['-DCMAKE_TOOLCHAIN_FILE='+str(a.toolchain.resolve())]
    if a.vcpkg_classic:configure+=['-DVCPKG_MANIFEST_MODE=OFF']
    commands=[configure,['cmake','--build',str(build),'--config','Release','-j','8'],['ctest','--test-dir',str(build),'-C','Release','--output-on-failure']]
    with (directory/'rebuild.log').open('w') as log:
        for cmd in commands:subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT)
    exe=build/('Release/sdf_generalization_bench.exe' if os.name=='nt' else 'sdf_generalization_bench');matrix=rows(directory/'verification/matrix.csv');reproduced=[]
    for g in frozen['geometries']:
        for method in (METHODS if g['shape']=='mesh' or g['id'] in ('00_thin_box','06_capsule') else ['candidate']):
            seed=0 if method=='grid' else 240;name=f"{g['id']}_{seed}_{method}";original=directory/'jobs'/name;target=directory/'reproduction'/name;ref=reference_dir(directory,g);target.parent.mkdir(exist_ok=True)
            if not (target/'COMPLETE').exists():
                assert not target.exists(),'Preserve incomplete reproduction and rerun with a fresh target';command=[str(exe),'fit',str(ref/'geometry.txt'),str(ref),str(target),str(seed),method,*map(str,CAPS),'--probe-max']
                if method=='candidate':command+=['--policy',str(directory/'FINAL_POLICY.txt')]
                subprocess.run(command,check=True,stdout=(target.parent/(name+'.log')).open('w'),stderr=subprocess.STDOUT)
            rr=[r for r in matrix if r['geometry']==g['id'] and int(r['seed'])==seed and r['method']==method];passing=[r for r in rr if r['pass']=='True'];cap=int(min(passing,key=lambda r:(int(r['bytes']),int(r['cap'])))['cap']) if passing else CAPS[-1];index=CAPS.index(cap)
            failed=[int(r['cap']) for r in rr if r['pass']!='True' and int(r['cap'])!=cap]
            adjacent=min(failed,key=lambda c:(abs(CAPS.index(c)-index),c)) if failed else None
            selected=([adjacent] if adjacent is not None else[])+[cap]
            for c in selected:
                timing_target=target
                if method=='grid':
                    timing_target=directory/'reproduction-cold-grid'/f'{name}_{c}'
                    timing_target.parent.mkdir(exist_ok=True)
                    if not (timing_target/'COMPLETE').exists():
                        assert not timing_target.exists(),'Preserve incomplete cold grid reproduction'
                        with (timing_target.parent/f'{name}_{c}.log').open('w') as log:
                            subprocess.run([str(exe),'fit',str(ref/'geometry.txt'),str(ref),str(timing_target),str(seed),method,str(c),'--probe-max'],check=True,stdout=log,stderr=subprocess.STDOUT)
                names=[f'{c}.bin',*(f'{c}_{ds}.pred' for ds in ('uniform','near','audit','probes','boundary')),f'{c}_boundary.bin']
                for file in names:
                    assert sha(original/file)==sha(target/file),f'Reproduction mismatch: {name}/{file}'
                    assert sha(original/file)==sha(timing_target/file),f'Cold timing reproduction mismatch: {name}/{file}'
                reproduced.append(dict(geometry=g['id'],method=method,seed=seed,cap=c,field_sha256=sha(target/f'{c}.bin'),predictions_exact=True,timing_folder=str(timing_target.relative_to(directory)),status='crossing' if any(int(r['cap'])==c and r['pass']=='True' for r in rr) else 'adjacent_or_censored_failure'))
            print(f'Frozen-source job {name}: exact fields/predictions checked; {len(reproduced)} selected fields total',flush=True)
    json_write(directory/'REPRODUCED.json',dict(fields=reproduced,source_manifest_sha256=sha(source/'SHA256SUMS.json'),frozen_sha256=sha(directory/'FROZEN.json'),rebuild_executable_sha256=sha(exe),fresh_tests_passed=True,timing='remeasured raw timing.csv in reproduction jobs'))
    print(f'Exact frozen-source reproduction verified for {len(reproduced)} selected fields',flush=True)
if __name__=='__main__':main()
