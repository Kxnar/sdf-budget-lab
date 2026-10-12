"""Independent protocol utilities. No selection on representation results."""
import csv, hashlib, json, math, struct
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
METHODS=['grid','mse','sign','cvar99','selective32','candidate_original','coverage_q8','candidate']
CAPS=[1024*x for x in (8,12,16,24,32,48,64,96,128,192,256,512,1024)]
SEEDS=list(range(240,245))
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()
def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f))
def manifest(folder):return {p.relative_to(folder).as_posix():sha(p) for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.json'}
def check_hashes(folder,expected):
    for name,h in expected.items():assert sha(folder/name)==h, f'Hash mismatch: {folder/name}'
def json_write(p,value):Path(p).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def cache_identity(g,spec):
    record=dict(spec_sha256=sha(spec),normalized_mesh_sha256=g.get('normalized_sha256'),oracle='EPICK-audited' if g['shape']=='mesh' and g['parameters'][1]==1 else'EPECK' if g['shape']=='mesh' else'analytic-generalization-v1')
    return hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()
def reference_dir(directory,g):return directory/g.get('reference_subdir','reference/'+g['id'])
def csv_write(p,rr):
    if not rr:return
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,rr[0].keys());w.writeheader();w.writerows(rr)
def rotation(angles):
    x,y,z=angles;cx,sx,cy,sy,cz,sz=np.cos(x),np.sin(x),np.cos(y),np.sin(y),np.cos(z),np.sin(z)
    return np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]])@np.array([[cy,0,sy],[0,1,0],[-sy,0,cy]])@np.array([[1,0,0],[0,cx,-sx],[0,sx,cx]])
def forward(g,p):return np.asarray(p)@np.array(g['rotation']).T+g['translation']
def exact(g,points):
    p=(np.asarray(points)-g['translation'])@np.array(g['rotation']);a,b,c=g['parameters'];s=g['shape']
    if s=='torus':return np.hypot(np.hypot(p[...,0],p[...,1])-.55,p[...,2])-a
    if s=='double_sphere':return np.minimum(np.linalg.norm(p-[.3+a/2,0,0],axis=-1),np.linalg.norm(p+[.3+a/2,0,0],axis=-1))-.3
    if s=='capsule':
        q=p.copy();q[...,2]-=np.clip(q[...,2],-b,b);return np.linalg.norm(q,axis=-1)-a
    if s=='cylinder':
        q=np.stack([np.linalg.norm(p[...,:2],axis=-1)-a,np.abs(p[...,2])-b],axis=-1);return np.linalg.norm(np.maximum(q,0),axis=-1)+np.minimum(q.max(axis=-1),0)
    if s=='shell':return np.abs(np.linalg.norm(p,axis=-1)-b)-a
    h=np.array([.6,.6,a]) if s=='thin_box' else np.array([b,c,.12]);q=np.abs(p)-h
    return np.linalg.norm(np.maximum(q,0),axis=-1)+np.minimum(q.max(axis=-1),0)-(a if s=='rounded_box' else 0)
def obj(path):
    vv,ff=[],[]
    for l in Path(path).read_text().splitlines():
        if l.startswith('v '):vv.append([float(x) for x in l.split()[1:4]])
        if l.startswith('f '):
            ids=[int(x.split('/')[0])-1 for x in l.split()[1:]]
            assert len(ids)==3,'Nontriangle input; no triangulation allowed'
            ff.append(ids)
    return np.array(vv),np.array(ff,dtype=int)
def mesh_topology(vv,ff):
    parent=list(range(len(vv)))
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    edges,links={},{}
    closed=len(ff)>0
    for a,b,c in ff:
        closed &= len({int(a),int(b),int(c)})==3 and np.linalg.norm(np.cross(vv[b]-vv[a],vv[c]-vv[a]))>0
        for u,v,w in ((a,b,c),(b,c,a),(c,a,b)):
            edges.setdefault(tuple(sorted((u,v))),[]).append(1 if u<v else -1);parent[root(u)]=root(v);links.setdefault(u,[]).append((v,w))
    closed &= all(len(v)==2 and sum(v)==0 for v in edges.values())
    for pairs in links.values():
        adj={}
        for a,b in pairs:adj.setdefault(a,[]).append(b);adj.setdefault(b,[]).append(a)
        seen,todo=set(),[next(iter(adj))]
        while todo:
            u=todo.pop()
            if u not in seen:seen.add(u);todo.extend(adj[u])
        closed &= len(seen)==len(adj) and all(len(v)==2 for v in adj.values())
    chi={}
    for u in links:chi[root(u)]=chi.get(root(u),0)+1
    for a,b in edges:chi[root(a)]-=1
    for a,b,c in ff:chi[root(a)]+=1
    return bool(closed),sorted(chi.values())
def decode_grid(data,n,points):
    values=np.frombuffer(data,dtype='<f4');p=np.asarray(points);u=(p+1)*(n-1)/2;base=np.minimum(np.floor(u).astype(int),n-2);t=u-base;result=np.zeros(len(p))
    for c in range(8):
        side=np.array([c&1,(c>>1)&1,(c>>2)&1]);v=base+side
        result+=values[v[:,0]+n*(v[:,1]+n*v[:,2])]*np.prod(np.where(side,t,1-t),axis=1)
    return result
def decode_packed(data,points,metadata=False):
    """Independent decoder with explicit depth-seven terminal validation."""
    p=np.asarray(points);raw=np.frombuffer(data,dtype=np.uint8);words=np.frombuffer(data,dtype='<u4');floats=np.frombuffer(data,dtype='<f4');lo=np.full_like(p,-1.);hi=np.full_like(p,1.);node=np.zeros(len(p),dtype=np.int64);depth=np.zeros(len(p),dtype=np.int64)
    for _ in range(8):
        desc=words[node//4];active=(desc&3)==0
        if not active.any():break
        mid=(lo[active]+hi[active])/2;side=p[active]>=mid;node[active]=(desc[active]&~np.uint32(3)).astype(np.int64)+4*(side@np.array([1,2,4]));lo[active]=np.where(side,mid,lo[active]);hi[active]=np.where(side,hi[active],mid);depth[active]+=1
    desc=words[node//4];assert np.all(desc&3),'Unexpected cycle/depth beyond seven';enc=desc&3;assert np.all((enc==1)|(enc==3));offset=(desc&~np.uint32(3)).astype(np.int64);t=(p-lo)/(hi-lo);result=np.zeros(len(p))
    for c in range(8):
        v=np.zeros(len(p));f=enc==1;v[f]=floats[offset[f]//4+c];q=enc==3;v[q]=(floats[offset[q]//4].astype(float)+floats[offset[q]//4+1].astype(float)*raw[offset[q]+8+c]).astype(np.float32).astype(float);side=np.array([bool(c&(1<<i)) for i in range(3)]);result+=v*np.prod(np.where(side,t,1-t),axis=1)
    return (result,node,depth) if metadata else result
def independent_roots(query,rays):
    rays=np.asarray(rays);o=rays[:,1:4];d=rays[:,4:7];t=rays[:,7,None]+(rays[:,8]-rays[:,7])[:,None]*np.linspace(0,1,1025)
    val=query((o[:,None,:]+t[...,None]*d[:,None,:]).reshape(-1,3)).reshape(len(rays),1025);ri,ix=np.nonzero((val[:,:-1]<0)!=(val[:,1:]<0));a,b=t[ri,ix].copy(),t[ri,ix+1].copy();sign=val[ri,ix]<0
    for _ in range(36):
        m=(a+b)/2;s=query(o[ri]+m[:,None]*d[ri])<0;a=np.where(s==sign,m,a);b=np.where(s==sign,b,m)
    v=(a+b)/2;values=query(o[ri]+v[:,None]*d[ri]);return [(v[ri==i],bool(np.all(np.abs(values[ri==i])<=1e-7))) for i in range(len(rays))]
def brute_mesh(g,points):
    """Independent brute-force Euclidean triangle distance and ray parity sign."""
    vv,ff=obj(g['normalized_mesh']);a,b,c=vv[ff[:,0]],vv[ff[:,1]],vv[ff[:,2]];ab=b-a;ac=c-a;n=np.cross(ab,ac);nn=np.einsum('ij,ij->i',n,n);values=[]
    direction=np.array([.738179,.521087,.429771]);direction/=np.linalg.norm(direction)
    # Deliberately no CGAL, BVH or triangulation dependency.
    for p in points:
        ap=p-a;s=np.einsum('ij,ij->i',ap,n)/nn;projection=p-s[:,None]*n
        v0=ab;v1=ac;v2=projection-a;d00=np.einsum('ij,ij->i',v0,v0);d01=np.einsum('ij,ij->i',v0,v1);d11=np.einsum('ij,ij->i',v1,v1);d20=np.einsum('ij,ij->i',v2,v0);d21=np.einsum('ij,ij->i',v2,v1);den=d00*d11-d01*d01;u=(d11*d20-d01*d21)/den;v=(d00*d21-d01*d20)/den;inside=(u>=0)&(v>=0)&(u+v<=1);best=np.where(inside,s*s*nn,np.inf)
        for x,y in ((a,b),(b,c),(c,a)):
            edge=y-x;t=np.clip(np.einsum('ij,ij->i',p-x,edge)/np.einsum('ij,ij->i',edge,edge),0,1);delta=p-x-t[:,None]*edge;best=np.minimum(best,np.einsum('ij,ij->i',delta,delta))
        from generalization_prepare import mesh_intersections
        hits=mesh_intersections(vv,ff,p,direction,1e-10,10);values.append(np.sqrt(best.min())*(-1 if len(hits)%2 else 1))
    return np.array(values)
def crossing(rr):
    passing=[r for r in rr if r['pass']]
    return min((r['bytes'] for r in passing),default=None)
