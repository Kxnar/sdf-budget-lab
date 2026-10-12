"""Explicit pre-freeze EPECK versus EPICK construction audit on unchanged meshes."""
import argparse,subprocess,time
from generalization_common import *
from generalization_prepare import write_spec,mesh_rays
def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);ap.add_argument('--exe',type=Path,default=ROOT/'build-generalization-cgal/Release/sdf_generalization_bench.exe');a=ap.parse_args();directory=a.directory.resolve();exe=a.exe.resolve();geometries=json.loads((directory/'qualification.json').read_text());audit=directory/'oracle-audit';audit.mkdir(exist_ok=True);results=[]
    for g in geometries:
        if g['shape']!='mesh':continue
        folder=audit/g['id'];folder.mkdir(exist_ok=True);spec=folder/'geometry.txt';write_spec(g,spec);v,f=obj(g['normalized_mesh']);centroids=v[f].mean(axis=1);edges=sorted({tuple(sorted((int(a),int(b)))) for a,b,c in f for a,b in ((a,b),(b,c),(c,a))});edge_points=np.array([(v[a]+v[b])/2 for a,b in edges]);normals=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]);normals/=np.linalg.norm(normals,axis=1)[:,None];ray,roots=mesh_rays(g,v,f);witness=[]
        for r,t in zip(ray,roots):
            times=np.r_[r[7],t,r[8]];witness.extend(r[1:4]+((times[:-1]+times[1:])/2)[:,None]*r[4:7])
        points=np.concatenate([v,centroids,edge_points,centroids+1e-7*normals,centroids-1e-7*normals,np.array(witness)]);np.savetxt(folder/'comparison.points',points,fmt='%.17g');json_write(folder/'groups.json',dict(uniform=4096,near=4096,vertices=len(v),face_centers=len(centroids),edge_centers=len(edge_points),two_sided_faces=2*len(centroids),witnesses=len(witness)))
        if not (folder/'COMPLETE').exists():
            print('Comparing exact and explicit fast oracle',g['id'],len(points)+8192,'queries',flush=True)
            with (folder/'comparison.log').open('w') as log:subprocess.run([str(exe),'compare',str(spec),str(folder),'4096',str(991007+int(g['id'][:2]))],check=True,stdout=log,stderr=subprocess.STDOUT)
        exact_samples=np.fromfile(folder/'epeck.bin',dtype='<f8').reshape(-1,4);fast=np.fromfile(folder/'epick.bin',dtype='<f8');assert len(exact_samples)==len(fast)==len(points)+8192;delta=np.abs(exact_samples[:,3]-fast);band=1e-9;away=np.abs(exact_samples[:,3])>band;sign_mismatch=away&(np.sign(exact_samples[:,3])!=np.sign(fast));near_mismatch=(~away)&(np.sign(exact_samples[:,3])!=np.sign(fast));native=json.loads((folder/'comparison.json').read_text());passed=bool(delta.max()<=1e-10 and not sign_mismatch.any() and native['sampler_max_delta']<=1e-12)
        # Independent geometry audit of representative queries, distinct from kernel comparison.
        ix=np.r_[np.linspace(0,8191,128,dtype=int),np.arange(len(exact_samples)-min(32,len(witness)),len(exact_samples))];independent=brute_mesh(g,exact_samples[ix,:3]);np.testing.assert_allclose(independent,exact_samples[ix,3],atol=2e-10,rtol=2e-9)
        result=dict(geometry=g['id'],passed=passed,queries=len(fast),stratified_queries=8192,maximum_distance_difference=float(delta.max()),distance_threshold=1e-10,unexplained_sign_mismatches=int(sign_mismatch.sum()),near_band_sign_mismatches=int(near_mismatch.sum()),tiny_sign_band=band,epeck_seconds=native['epeck_seconds'],epick_seconds=native['epick_seconds'],sampler_max_delta=native['sampler_max_delta'],epeck_geometric_validation=True,independent_brute_force_queries=len(ix),source_sha256=g['source_sha256'],normalized_sha256=g['normalized_sha256'],executable_sha256=sha(exe));json_write(folder/'RESULT.json',result);json_write(folder/'SHA256SUMS.json',manifest(folder));results.append(result);print(g['id'],'pass',passed,'maxdiff',delta.max(),'speed ratio',native['epeck_seconds']/native['epick_seconds'],flush=True)
    assert len(results)==4;json_write(audit/'SUMMARY.json',dict(all_passed=all(r['passed'] for r in results),chosen_oracle='EPICK-audited' if all(r['passed'] for r in results) else'EPECK',results=results,documentation='EPICK exact predicates; constructions generally inexact. This audit bounds observed kernel differences only, not all possible mesh distance queries.'));assert all(r['passed'] for r in results),'Fast oracle not eligible'
if __name__=='__main__':main()
