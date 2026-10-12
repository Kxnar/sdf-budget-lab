"""Predeclared bounded DEVELOPMENT selection; never reads final method outputs."""
import argparse,subprocess,time,shutil
from generalization_common import *
from generalization_prepare import create_reference
from generalization_verify import verify_job
CAPS=[1024*x for x in (8,12,16,24,32,48,64,96,128,192,256)]
POLICIES=[dict(k=128,lambda_= .001,global_quadrature=0,covered_ranking=0)]
POLICIES += [dict(k=k,lambda_=lam,global_quadrature=0,covered_ranking=0) for k in (512,1024,4096) for lam in (0,.001)]
POLICIES += [dict(k=4096,lambda_=.001,global_quadrature=1,covered_ranking=0)]
POLICIES += [dict(k=4096,lambda_=lam,global_quadrature=q,covered_ranking=1) for q in (1,2) for lam in (0,.001)]
def run_job(base,g,seed,exe,policy):
    ref=base/'reference'/g['id'];folder=base/'jobs'/f"{g['id']}_{seed}_candidate";folder.parent.mkdir(exist_ok=True)
    if (folder/'COMPLETE').exists():check_hashes(folder,json.loads((folder/'SHA256SUMS.json').read_text()));return
    assert not folder.exists(),'Incomplete development evidence retained; use separate named screen attempt'
    with (folder.parent/(folder.name+'.log')).open('w') as log:
        result=subprocess.run([str(exe),'fit',str(ref/'geometry.txt'),str(ref),str(folder),str(seed),'candidate',*map(str,CAPS),'--policy',str(base/'policy.txt')],stdout=log,stderr=subprocess.STDOUT)
    (folder/'EXIT_CODE').write_text(str(result.returncode));json_write(folder/'job.json',dict(frozen_sha256=sha(base/'FROZEN.json'),geometry=g['id'],seed=seed,method='candidate',caps=CAPS));json_write(folder/'SHA256SUMS.json',manifest(folder));assert result.returncode==0
def native_pass(base,g,seed):
    folder=base/'jobs'/f"{g['id']}_{seed}_candidate";top=rows(folder/'topology.csv');result=[]
    for r in rows(folder/'fields.csv'):
        passed=all(int(r[k]) for k in ('scalar_pass','probe_pass','width_pass')) and r['topology_status']=='evaluated';tt=[t for t in top if t['cap']==r['cap']];passed &= len(tt)==2 and all(int(t['closed']) and list(map(int,t['component_euler'].split(';')))==g['expected_euler'] for t in tt);result.append(dict(cap=int(r['cap']),bytes=int(r['bytes']),pass_=passed,oracle_calls=int(r['oracle_calls']),build_time=sum(float(r[k]) for k in ('build_tree_seconds','loss_seconds','dp_seconds','pack_seconds'))))
    return result
def rank_policy(base,geometries,order):
    count=0;penalized_bytes=0;details=[];calls=0;timecost=0
    for g in geometries:
        seedrows=[native_pass(base,g,s) for s in (40,41)];common=[]
        for i,c in enumerate(CAPS):
            r=[rr[i] for rr in seedrows]
            if all(x['pass_'] for x in r):common.append((max(x['bytes'] for x in r),c))
        best=min(common) if common else None;count+=bool(best);penalized_bytes+=best[0] if best else 512*1024
        calls+=sum(rr[-1]['oracle_calls'] for rr in seedrows);timecost+=sum(rr[-1]['build_time'] for rr in seedrows);details.append(dict(geometry=g['id'],common_passing_actual_bytes=best[0] if best else None,cap=best[1] if best else None,censored=best is None))
    return dict(policy=base.name,order=order,passing_geometry_count=count,selection_penalized_bytes=penalized_bytes,oracle_calls=calls,build_time=timecost,details=details)
def link_reference(source,target,screen=False):
    target.mkdir(parents=True)
    for p in source.iterdir():
        if p.is_file():
            if screen and p.name in ('uniform.bin','near.bin'):
                data=np.fromfile(p,dtype='<f8').reshape(-1,4)[:10000];data.tofile(target/p.name)
            else:shutil.copy2(p,target/p.name)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('out',type=Path);ap.add_argument('--exe',type=Path,default=ROOT/'build-generalization-cgal/Release/sdf_generalization_bench.exe');a=ap.parse_args();out=a.out.resolve();exe=a.exe.resolve();out.mkdir(exist_ok=True)
    declaration=out/'SCREEN_FROZEN.json'
    if not declaration.exists():json_write(declaration,dict(protocol='bounded-development-screen-v1',policies=POLICIES,caps=CAPS,seeds=[40,41],fit_count=10000,screen_uniform_count=10000,screen_near_count=10000,confirm_uniform_count=100000,confirm_near_count=100000,maximum_screen_fields=1320,ranking='all5geometry crossing bothseeds samecap; then sum actual crossing bytes plus512KiB selection-only penalty per censoredcase; then feweroraclecalls, buildtime, fixedorder',gates='generalization-v1 unchanged',source_at_declaration_sha256={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'src/generalization_benchmark.cpp',ROOT/'include/sdf/generalization.hpp',Path(__file__)]}))
    else:assert json.loads(declaration.read_text())['policies']==POLICIES
    geometries=[]
    for i,(shape,feature) in enumerate([('torus',.035),('torus',.060),('torus',.100),('double_sphere',.060),('double_sphere',.100)]):
        g=dict(id=f'development_{i}_{shape}',shape=shape,stratum='development',parameters=[feature,0,0],rotation=np.eye(3).tolist(),translation=[.017,-.023,.011],noise=feature/2,expected_euler=[0] if shape=='torus' else[2,2]);ref=out/'reference'/g['id']
        if (ref/'metadata.json').exists():g=json.loads((ref/'metadata.json').read_text())
        else:g=create_reference(g,ref,exe,count=100000,seed=778001+i*1000,levels=(64,96));json_write(ref/'metadata.json',g)
        geometries.append(g)
    ranking=[]
    for i,policy in enumerate(POLICIES):
        base=out/'policies'/f'p{i:02d}';base.mkdir(parents=True,exist_ok=True);(base/'policy.txt').write_text(f"{policy['k']} {policy['lambda_']} {policy['global_quadrature']} {policy['covered_ranking']}\n")
        if not (base/'FROZEN.json').exists():
            for g in geometries:link_reference(out/'reference'/g['id'],base/'reference'/g['id'],screen=True)
            json_write(base/'FROZEN.json',dict(declaration_sha256=sha(declaration),policy=policy))
        for g in geometries:
            for seed in (40,41):run_job(base,g,seed,exe,policy)
        result=rank_policy(base,geometries,i);ranking.append(result);json_write(out/'ranking-progress.json',ranking);print(f"Screen {base.name}: {result['passing_geometry_count']}/5 geometry crossings",flush=True)
    ranking.sort(key=lambda r:(-r['passing_geometry_count'],r['selection_penalized_bytes'],r['oracle_calls'],r['build_time'],r['order']));json_write(out/'RANKING.json',ranking)
    for rank in ranking:
        if rank['passing_geometry_count']<5:break
        i=rank['order'];policy=POLICIES[i];base=out/'confirmation'/f'p{i:02d}';base.mkdir(parents=True,exist_ok=True);shutil.copy2(out/'policies'/f'p{i:02d}'/'policy.txt',base/'policy.txt')
        if not (base/'FROZEN.json').exists():
            for g in geometries:link_reference(out/'reference'/g['id'],base/'reference'/g['id'])
            json_write(base/'FROZEN.json',dict(declaration_sha256=sha(declaration),policy=policy,confirmation=True))
        independent=[]
        for g in geometries:
            for seed in (40,41):
                run_job(base,g,seed,exe,policy);rr,q=verify_job(base,g,seed,'candidate');independent.extend(rr);print('Confirmed independently',base.name,g['id'],seed,flush=True)
        csv_write(base/'independent_matrix.csv',independent);result=rank_policy(base,geometries,i)
        if result['passing_geometry_count']==5:
            selected={'depth':5,'k':policy['k'],'lambda':policy['lambda_'],'quadrature':4,'global_quadrature':policy['global_quadrature'],'covered_ranking':policy['covered_ranking']};json_write(out/'SELECTED.json',dict(policy=selected,confirmation=result,declaration_sha256=sha(declaration),rankings=ranking));print('DEVELOPMENT READY',selected,flush=True);return
    json_write(out/'NO_POLICY_PASSES.json',dict(reason='No predeclared policy independently confirmed all5development cases at unchanged strict joint gates',rankings=ranking));print('NO POLICY PASSES: bounded screen complete',flush=True)
if __name__=='__main__':main()
