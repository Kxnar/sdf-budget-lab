"""Immutable resumable deterministic development jobs; no persistent oracle cache."""
import argparse,subprocess,shutil,os,time
from generalization_common import *
from compact_verify import verify
CODECS={'f64':1,'f32':2,'Q24':3,'Q28':4,'XOR32':5}
DEVELOPMENT=('12_goathead','14_spot','15_brucewick')

def sources():
    files=[ROOT/'CMakeLists.txt',ROOT/'docs/generalization-v2-compact-plan.md']
    for folder in ('src','include','tests'):
        files.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and p.suffix in ('.cpp','.hpp','.py'))
    files.extend(ROOT/'scripts'/name for name in ('compact_run.py','compact_mesh.py','compact_verify.py','generalization_common.py'))
    return {p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)}

def snapshot(out,exe):
    target=out/'source';identity=dict(sources=sources(),executable_sha256=sha(exe),declaration_sha256=sha(ROOT/'runs/generalization-v2/COMPACT-DECLARATION.json'))
    if (out/'SOURCE.json').exists():
        assert json.loads((out/'SOURCE.json').read_text())==identity,'Changed native/evaluator/source; preserve old run and use frozen source to resume'
    else:
        assert not target.exists();target.mkdir(parents=True)
        for name in identity['sources']:
            dest=target/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/name,dest)
        for p in exe.parent.glob('*.dll'):shutil.copy2(p,target/p.name)
        shutil.copy2(exe,target/exe.name);json_write(out/'SOURCE.json',identity)
    return target/exe.name,identity

def selection(out,assessments):
    rr=[]
    for name,codec in CODECS.items():
        if name in ('f64','f32'):continue
        mesh=[];cost=0.;query=0.
        for geometry in DEVELOPMENT:
            a=next(a for a in assessments if a['geometry']==geometry and a['codec']==codec)
            okay=a['joint_pass'] and a['actual_bytes']<=CAPS[-1]
            mesh.append(dict(geometry=geometry,passing=okay,selection_bytes=a['actual_bytes'] if okay else 2*1024*1024))
            if (out/'jobs'/f'{geometry}_{name}'/'construction.json').exists():
                c=json.loads((out/'jobs'/f'{geometry}_{name}'/'construction.json').read_text());cost+=sum(c.get(k,0.) for k in ('source_seconds','encoding_bvh_packing_seconds','decoded_validation_seconds'))
                tt=rows(out/'jobs'/f'{geometry}_{name}'/'timing.csv') if (out/'jobs'/f'{geometry}_{name}'/'timing.csv').exists() else []
                query+=float(np.median([float(r['ns']) for r in tt])) if tt else float('inf')
        rr.append(dict(name=name,passing_meshes=sum(m['passing'] for m in mesh),selection_bytes=sum(m['selection_bytes'] for m in mesh),construction_seconds=cost,warm_query_ns_sum=query,mesh_scores=mesh,order=list(CODECS).index(name)))
    ranked=sorted(rr,key=lambda r:(-r['passing_meshes'],r['selection_bytes'],r['construction_seconds'],r['warm_query_ns_sum'],r['order']))
    json_write(out/'SELECTION.json',dict(ranked=ranked,winner=ranked[0]['name'],selection_penalty_only=2*1024*1024,no_fresh_fitting=True))
    return ranked[0]['name']

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='runs/generalization-v2/compact-development');ap.add_argument('--exe',default='build-generalization-v2-cgal/Release/sdf_compact_bench.exe');ap.add_argument('--geometry',nargs='*');ap.add_argument('--codec',nargs='*');ap.add_argument('--surface-check',action='store_true');args=ap.parse_args()
    out=ROOT/args.out;out.mkdir(parents=True,exist_ok=True);exe,identity=snapshot(out,ROOT/args.exe)
    old=ROOT/'runs/generalization-v1';frozen=json.loads((old/'FROZEN.json').read_text());assessments=[];matrix=[]
    for geometry in DEVELOPMENT:
        if args.geometry and geometry not in args.geometry:continue
        g=next(g for g in frozen['geometries'] if g['id']==geometry);ref=reference_dir(old,g)
        for name,codec in CODECS.items():
            if args.codec and name not in args.codec:continue
            job=out/'jobs'/f'{geometry}_{name}';job.parent.mkdir(exist_ok=True)
            check_surface=args.surface_check and name=='f64'
            job_identity=dict(source_identity=identity,geometry_sha256=sha(Path(g['normalized_mesh'])),reference_hashes=manifest(ref),codec=codec,counts={'uniform':100000,'near':100000},surface_check=check_surface)
            saved=out/'jobs'/f'{geometry}_{name}.json'
            if saved.exists():assert json.loads(saved.read_text())==job_identity
            else:json_write(saved,job_identity)
            if not (job/'COMPLETE').exists():
                assert not job.exists(),'Incomplete attempt preserved; explicit investigation required'
                command=[str(exe),str(ref/'geometry.txt'),str(ref),str(job),str(codec)]
                if check_surface:command.append('--surface-check')
                with (out/'jobs'/f'{geometry}_{name}.log').open('w',encoding='utf-8') as log:
                    result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
                if result.returncode:
                    text=(out/'jobs'/f'{geometry}_{name}.log').read_text();status='resource_unresolved' if 'resource_unresolved' in text else 'incomplete'
                    (job/'STATUS').write_text(status+'\n',encoding='utf-8');print(geometry,name,status,flush=True);raise RuntimeError('Native failure retained; investigate before continuation')
            print('independent verification',geometry,name,flush=True);a=verify(job,ref,g);a['codec']=codec;assessments.append(a)
            json_write(job/'SHA256SUMS.json',manifest(job))
            if name=='f64' and not a['joint_pass']:raise RuntimeError('f64 control failed: implementation/reference checkpoint required')
            for cap in CAPS:
                feasible=a['actual_bytes'] is not None and a['actual_bytes']<=cap
                matrix.append(dict(geometry=geometry,method=name,cap=cap,actual_bytes=a['actual_bytes'],status=a['status'] if feasible else 'storage_infeasible' if a['actual_bytes'] is not None else a['status'],joint_pass=bool(feasible and a['joint_pass']),unique_builds=1,seed='deterministic'))
            csv_write(out/'matrix.partial.csv',matrix);print('verified',geometry,name,a['status'],a['actual_bytes'],flush=True)
    if not args.geometry and not args.codec:
        assert len(assessments)==15 and len(matrix)==195
        csv_write(out/'matrix.csv',matrix);winner=selection(out,assessments);json_write(out/'COMPLETE.json',dict(unique_fields=15,logical_cap_rows=195,winner=winner,source_identity=identity,complete=True));json_write(out/'SHA256SUMS.json',manifest(out))
if __name__=='__main__':main()
