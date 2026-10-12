"""Final bounded d7 DEVELOPMENT intervention; no further policy expansion."""
import argparse,subprocess,shutil,time
from generalization_common import *
from generalization_screen import rank_policy,link_reference
from generalization_verify import verify_job
CAPS=[1024*x for x in (8,12,16,24,32,48,64,96,128,192,256)]
def job(base,g,seed,exe,out):
    ref=base/'reference'/g['id'];folder=base/'jobs'/f"{g['id']}_{seed}_candidate";folder.parent.mkdir(exist_ok=True)
    if (folder/'COMPLETE').exists():check_hashes(folder,json.loads((folder/'SHA256SUMS.json').read_text()));return
    assert not folder.exists();cache=out/'shared-cache'/(g['id']+'.bin');identity=cache_identity(g,ref/'geometry.txt');command=[str(exe),'fit',str(ref/'geometry.txt'),str(ref),str(folder),str(seed),'candidate',*map(str,CAPS),'--policy',str(base/'policy.txt'),'--probe-max','--cache',str(cache),'--cache-key',identity]
    with (folder.parent/(folder.name+'.log')).open('w') as log:result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
    (folder/'EXIT_CODE').write_text(str(result.returncode));json_write(folder/'job.json',dict(frozen_sha256=sha(base/'FROZEN.json'),geometry=g['id'],seed=seed,method='candidate',caps=CAPS,cache_identity=identity));json_write(folder/'SHA256SUMS.json',manifest(folder));assert result.returncode==0,(folder.parent/(folder.name+'.log')).read_text()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('out',type=Path);ap.add_argument('--previous',type=Path,default=ROOT/'runs/generalization-development-screen');ap.add_argument('--exe',type=Path,default=ROOT/'build-generalization-core/Release/sdf_generalization_bench.exe');a=ap.parse_args();out=a.out.resolve();out.mkdir(exist_ok=True);previous=a.previous.resolve();exe=a.exe.resolve();policies=[dict(k=4096,lambda_=0,global_quadrature=1,covered_ranking=1,k7=k7) for k7 in (0,128,512,2048)];declaration=out/'SCREEN_FROZEN.json'
    if not declaration.exists():json_write(declaration,dict(protocol='bounded-development-d7-v1',policies=policies,fit_seeds=[40,41],caps=CAPS,maximum_screen_fields=440,probe_max_gate=True,original_screen_sha256=sha(previous/'SCREEN_FROZEN.json'),executable_sha256=sha(exe),source_sha256={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'src/generalization_benchmark.cpp',ROOT/'src/precision.cpp',ROOT/'include/sdf/precision.hpp',Path(__file__)]},selection='count common-capgeometry crossings, actualretainedbytes+512KiBuncrossedselectionpenalty, physical/logicalcalls, buildtime, fixedorder; confirmrankedwinner at100kqueries then freeze one candidate regardlessall5'))
    else:assert json.loads(declaration.read_text())['policies']==policies and json.loads(declaration.read_text())['executable_sha256']==sha(exe)
    source=out/'source'
    if not source.exists():
        source.mkdir()
        for folder in ('src','include','tests','scripts'):
            for p in (ROOT/folder).rglob('*'):
                if p.suffix in ('.cpp','.hpp','.py'):target=source/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
        shutil.copy2(exe,source/exe.name);shutil.copy2(ROOT/'CMakeLists.txt',source/'CMakeLists.txt');json_write(source/'SHA256SUMS.json',manifest(source))
    geometries=[json.loads((ref/'metadata.json').read_text()) for ref in sorted((previous/'reference').iterdir())];ranking=[]
    for i,policy in enumerate(policies):
        base=out/'policies'/f'd7_{policy["k7"]}';base.mkdir(parents=True,exist_ok=True);(base/'policy.txt').write_text(f'4096 0 1 1 {policy["k7"]}\n')
        if not (base/'FROZEN.json').exists():
            for g in geometries:link_reference(previous/'reference'/g['id'],base/'reference'/g['id'],screen=True)
            json_write(base/'FROZEN.json',dict(declaration_sha256=sha(declaration),policy=policy))
        for g in geometries:
            for seed in (40,41):job(base,g,seed,exe,out)
        ranked=rank_policy(base,geometries,i);ranking.append(ranked);json_write(out/'ranking-progress.json',ranking);print(base.name,ranked['passing_geometry_count'],'/5 strict common-capcrossings',flush=True)
    ranking.sort(key=lambda r:(-r['passing_geometry_count'],r['selection_penalized_bytes'],r['oracle_calls'],r['build_time'],r['order']));json_write(out/'RANKING.json',ranking);winner=ranking[0];policy=policies[winner['order']];base=out/'confirmation'/winner['policy'];base.mkdir(parents=True,exist_ok=True);shutil.copy2(out/'policies'/winner['policy']/'policy.txt',base/'policy.txt')
    if not (base/'FROZEN.json').exists():
        for g in geometries:link_reference(previous/'reference'/g['id'],base/'reference'/g['id'])
        json_write(base/'FROZEN.json',dict(declaration_sha256=sha(declaration),policy=policy,confirmation=True))
    independent=[]
    for g in geometries:
        for seed in (40,41):
            job(base,g,seed,exe,out);rr,q=verify_job(base,g,seed,'candidate');independent.extend(rr);print('Independently confirmed',winner['policy'],g['id'],seed,flush=True)
    csv_write(base/'independent_matrix.csv',independent);confirmed=rank_policy(base,geometries,winner['order']);selected=dict(depth=5,k=4096,k7=policy['k7'],lambda_=0,quadrature=4,global_quadrature=1,covered_ranking=1,probe_max_gate=True);json_write(out/'SELECTED.json',dict(policy=selected,screen_winner=winner,confirmation=confirmed,declaration_sha256=sha(declaration),bounded_intervention_complete=True,no_further_policy_expansion=True));print('FIXED FINAL CANDIDATE',selected,flush=True)
if __name__=='__main__':main()
