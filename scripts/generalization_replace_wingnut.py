"""Record pre-freeze geometric rejection and attributed public-domain replacement."""
import urllib.request,shutil
from generalization_common import *
from generalization_prepare import REV
def main():
    directory=ROOT/'runs/generalization-v1';assert not (directory/'FROZEN.json').exists();geometries=json.loads((directory/'qualification.json').read_text());old=next(g for g in geometries if g['id']=='14_wingnut');log=(directory/'oracle-audit/14_wingnut/comparison.log').read_text();assert 'Self-intersecting mesh' in log
    json_write(directory/'inputs/REJECTED-wingnut.json',dict(geometry=old,reason='CGAL EPECK exact-predicate does_self_intersect rejected unchanged normalized mesh before any final field fitting',log=log))
    source=f'https://raw.githubusercontent.com/odedstein/meshes/{REV}/objects/spot/';raw=directory/'inputs/spot_low_resolution.obj';readme=directory/'inputs/spot-README.md';urllib.request.urlretrieve(source+'spot_low_resolution.obj',raw);urllib.request.urlretrieve(source+'README.md',readme);v,f=obj(raw);closed,euler=mesh_topology(v,f);assert closed;center=(v.min(axis=0)+v.max(axis=0))/2;radius=np.linalg.norm(v-center,axis=1).max();scale=.75/radius;r=np.array(old['rotation']);t=np.array(old['translation']);normalized=(v-center)*scale@r.T+t;normpath=directory/'inputs/spot-normalized.obj'
    with normpath.open('w') as stream:
        for p in normalized:stream.write('v '+' '.join(format(x,'.17g') for x in p)+'\n')
        for p in f:stream.write('f '+' '.join(str(x+1) for x in p)+'\n')
    g=dict(old);g.update(id='14_spot',normalized_mesh=str(normpath),source_mesh=str(raw),source_url=source+'spot_low_resolution.obj',license_text=readme.read_text(),source_sha256=sha(raw),normalized_sha256=sha(normpath),expected_euler=euler,normalization=dict(original_center=center.tolist(),original_radius=float(radius),scale_to_domain=float(scale),original_units='unspecified in upstream OBJ; distances convert by dividing scale_to_domain',normalized_bounding_sphere_radius=.75),vertices=len(v),triangles=len(f));geometries[geometries.index(old)]=g;shutil.copy2(directory/'qualification.json',directory/'qualification-before-wingnut-replacement.json');json_write(directory/'qualification.json',geometries);rejections=json.loads((directory/'inputs/rejections.json').read_text());rejections.append(dict(id='14_wingnut',reason='EPECK self-intersection validation',source_sha256=old['source_sha256'],replacement='14_spot'));json_write(directory/'inputs/rejections.json',rejections);print(g['id'],len(v),len(f),euler,flush=True)
if __name__=='__main__':main()
