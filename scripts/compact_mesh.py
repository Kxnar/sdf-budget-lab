"""Independent complete serialized-mesh decoder and exact triangle-ray events."""
from fractions import Fraction as F
import struct
import numpy as np
from generalization_common import mesh_topology

def decode(data):
    assert len(data)>=80 and data[:4]==b'SDFM'
    version,codec,iw=struct.unpack_from('<HBB',data,4)
    assert version==1 and codec in range(1,6)
    v,f,n,vp,tp,bp=struct.unpack_from('<6I',data,8)
    assert v and f and n and iw==(2 if v<=65536 else 4)
    assert struct.unpack_from('<6d',data,32)==(-1.,-1.,-1.,1.,1.,1.)
    assert vp==80 and vp<=tp<=bp and tp%4==bp%4==0
    assert bp+20*n==len(data) and tp+3*f*iw<=bp
    def integers(raw,count,width):
        value=int.from_bytes(raw,'little');mask=(1<<width)-1
        assert value>>(count*width)==0,'Nonzero packed tail'
        return np.array([(value>>(i*width))&mask for i in range(count)],dtype=np.uint64)
    blocks=[];params=None
    if codec in (1,2):
        end=vp+v*(24 if codec==1 else 12)
        assert end<=tp
        vertices=np.frombuffer(data[vp:end],dtype='<f8' if codec==1 else '<f4').astype(float).reshape(v,3)
    elif codec in (3,4):
        width=24 if codec==3 else 28;end=vp+48+(3*v*width+7)//8;assert end<=tp
        params=np.frombuffer(data[vp:vp+48],dtype='<f8').reshape(2,3).copy()
        assert np.isfinite(params).all() and np.all(params[1]>=0)
        q=integers(data[vp+48:end],3*v,width).reshape(v,3)
        assert np.all(q[:,params[1]==0]==0),'Nonzero constant-axis code'
        vertices=params[0]+q.astype(float)*params[1]
    else:
        count=(v+31)//32;directory=struct.unpack_from('<'+'I'*(count+1),data,vp)
        assert directory[0]==4*(count+1) and vp+directory[-1]<=tp
        words=[]
        for i in range(count):
            begin,endblock=directory[i:i+2];raw=data[vp+begin:vp+endblock]
            assert 2<=len(raw) and raw[0] in (0,1) and raw[1]==min(32,v-32*i)
            size=raw[1];blocks.append(dict(mode=raw[0],count=size,length=len(raw)))
            if raw[0]==0:
                assert len(raw)==2+24*size
                words.extend(struct.unpack_from('<'+'Q'*(3*size),raw,2))
            else:
                assert len(raw)>=29
                base=struct.unpack_from('<3Q',raw,2);widths=raw[26:29];assert max(widths)<=64
                bitcount=(size-1)*sum(widths);assert len(raw)==29+(bitcount+7)//8
                residual=int.from_bytes(raw[29:],'little');assert residual>>bitcount==0
                words.extend(base);cursor=0
                for j in range(1,size):
                    for a,w in enumerate(widths):
                        words.append(base[a]^((residual>>cursor)&((1<<w)-1)));cursor+=w
        end=vp+directory[-1]
        vertices=np.array(words,dtype='<u8').view('<f8').reshape(v,3)
    assert 0<=tp-end<=3 and not any(data[end:tp])
    triangles=np.frombuffer(data[tp:tp+3*f*iw],dtype='<u2' if iw==2 else '<u4').astype(np.int64).reshape(f,3)
    assert np.all(triangles<v) and np.isfinite(vertices).all() and np.all(abs(vertices)<=1)
    assert not any(data[tp+3*f*iw:bp]) and bp-(tp+3*f*iw)<=3
    nodes=[]
    for i in range(n):
        raw=struct.unpack_from('<6H2I',data,bp+20*i)
        assert all(raw[a]<=raw[a+3] for a in range(3));nodes.append(raw)
    seen=set();covered=[]
    def visit(i):
        assert i<n and i not in seen;seen.add(i);node=nodes[i];first,count=node[6:]
        if count:
            assert 1<=count<=16 and first+count<=f
            covered.extend(range(first,first+count))
            for p in vertices[np.unique(triangles[first:first+count])]:
                for a in range(3):
                    coordinate=(F(float(p[a]))+1)*65535
                    assert 2*node[a]<=coordinate<=2*node[a+3],'Nonconservative exact BVH bound'
        else:
            assert i+1<first<n
            for child in (i+1,first):
                assert all(node[a]<=nodes[child][a] and node[a+3]>=nodes[child][a+3] for a in range(3))
                visit(child)
    visit(0);assert len(seen)==n and sorted(covered)==list(range(f))
    return dict(vertices=vertices,triangles=triangles,nodes=nodes,codec=codec,index_width=iw,
                vertex_section_bytes=tp-vp,triangle_section_bytes=bp-tp,bvh_section_bytes=20*n,
                padding_bytes=(tp-end)+(bp-(tp+3*f*iw)),blocks=blocks,parameters=params)

def morton_order(vertices):
    lo=vertices.min(axis=0);w=vertices.max(axis=0)-lo
    q=np.floor(np.divide(vertices-lo,w,out=np.zeros_like(vertices),where=w!=0).clip(0,1)*2097151).astype(np.uint64)
    code=np.zeros(len(vertices),dtype=np.uint64)
    for b in range(21):
        for a in range(3):code|=((q[:,a]>>b)&1)<<(3*b+a)
    return np.lexsort((np.arange(len(vertices)),code))

def source_check(decoded,vertices,triangles):
    order=morton_order(vertices);inverse=np.empty(len(order),dtype=int);inverse[order]=np.arange(len(order))
    source=vertices[order];codec=decoded['codec']
    if codec in (1,5):np.testing.assert_array_equal(decoded['vertices'].view('u8'),source.view('u8'))
    elif codec==2:np.testing.assert_array_equal(decoded['vertices'],source.astype('f4').astype(float))
    else:
        p=decoded['parameters'];maximum=(1<<(24 if codec==3 else 28))-1
        np.testing.assert_array_equal(p[0],source.min(axis=0));np.testing.assert_array_equal(p[1],(source.max(axis=0)-source.min(axis=0))/maximum)
        codes=np.floor(np.divide(source-p[0],p[1],out=np.zeros_like(source),where=p[1]!=0)+.5).clip(0,maximum)
        np.testing.assert_array_equal(decoded['vertices'],p[0]+codes*p[1])
    def cycle(t):return min(tuple(t),tuple(np.roll(t,1)),tuple(np.roll(t,2)))
    assert sorted(cycle(t) for t in decoded['triangles'])==sorted(cycle(t) for t in inverse[triangles]),'Changed connectivity/winding'
    return float(np.linalg.norm(decoded['vertices']-source,axis=1).max())

def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def rational(p):return tuple(F(float(x)) for x in p)
def exact_hit(origin,direction,triangle):
    a,b,c=triangle;e1=sub(b,a);e2=sub(c,a);s=sub(origin,a);h=cross(direction,e2);det=dot(e1,h)
    if not det:return (None,False,dot(s,cross(e1,e2))==0)
    u=dot(s,h)/det;q=cross(s,e1);v=dot(direction,q)/det
    if u<0 or v<0 or u+v>1:return None,False,False
    return dot(e2,q)/det,u==0 or v==0 or u+v==1,False

def coplanar_clip(o,d,triangle,low,high):
    a,b,c=triangle;n=cross(sub(b,a),sub(c,a))
    for x,y in ((a,b),(b,c),(c,a)):
        e=sub(y,x);value=dot(cross(e,sub(o,x)),n);slope=dot(cross(e,d),n)
        if not slope:
            if value<0:return None
        elif slope>0:low=max(low,-value/slope)
        else:high=-value/slope if high is None else min(high,-value/slope)
        if high is not None and high<low:return None
    return low,high

def candidates(point,direction,triangles,low=0.,high=float('inf'),xyz=None):
    # Independent per-triangle slabs (no serialized BVH), rounded outwards at
    # each subtraction/division. Account for exact rational midpoint rounding.
    p=np.array([float(x) for x in point]);error=[]
    for x,y in zip(point,p):
        difference=abs(x-F(float(y)));e=np.nextafter(float(difference),np.inf)
        while F(float(e))<difference:e=np.nextafter(e,np.inf)
        error.append(e)
    if xyz is None:xyz=np.array(triangles,dtype=float)
    lower=np.nextafter(xyz.min(axis=1)-error,-np.inf);upper=np.nextafter(xyz.max(axis=1)+error,np.inf)
    enter=np.full(len(triangles),low);exit=np.full(len(triangles),high);mask=np.ones(len(triangles),dtype=bool)
    for a,dd in enumerate(direction):
        d=float(dd)
        if not d:mask &= (p[a]>=lower[:,a])&(p[a]<=upper[:,a]);continue
        l=np.nextafter(lower[:,a]-p[a],-np.inf);h=np.nextafter(upper[:,a]-p[a],np.inf)
        if d>0:l=np.nextafter(l/d,-np.inf);h=np.nextafter(h/d,np.inf)
        else:l,h=np.nextafter(h/d,-np.inf),np.nextafter(l/d,np.inf)
        enter=np.maximum(enter,l);exit=np.minimum(exit,h)
    return np.flatnonzero(mask&(enter<=exit))

def exact_sign(point,triangles,xyz=None):
    # Independent directions, not native retry sequence.
    for direction in ((.917137,.231719,.419731),(-.313177,.851719,.423173),(.731117,-.519731,.449177)):
        d=rational(direction);hits=0;bad=False
        for i in candidates(point,d,triangles,xyz=xyz):
            tri=triangles[i]
            t,edge,plane=exact_hit(point,d,tri)
            if plane:
                if coplanar_clip(point,d,tri,F(0),None) is not None:bad=True;break
            elif t is not None and t==0:return 0
            elif t is not None and t>0:
                if edge:bad=True;break
                hits+=1
        if not bad:return -1 if hits%2 else 1
    return None

def ray_events(vertices,faces,ray):
    o,d=rational(ray[1:4]),rational(ray[4:7]);lo,hi=F(float(ray[7])),F(float(ray[8]));events=set();coplanar=False
    triangles=[tuple(rational(p) for p in vertices[t]) for t in faces]
    xyz=vertices[faces]
    for i in candidates(o,d,triangles,np.nextafter(float(lo),-np.inf),np.nextafter(float(hi),np.inf),xyz):
        tri=triangles[i]
        t,edge,plane=exact_hit(o,d,tri)
        if plane:
            span=coplanar_clip(o,d,tri,lo,hi)
            if span is not None:
                coplanar=True
                if span[0]!=span[1]:return dict(unresolved=True,roots=[],tangencies=[],coplanar=True)
                if lo<span[0]<hi:events.add(span[0])
        elif t is not None and lo<t<hi:events.add(t)
    cuts=[lo]+sorted(events)+[hi];signs=[]
    for a,b in zip(cuts[:-1],cuts[1:]):
        t=(a+b)/2;p=tuple(x+t*y for x,y in zip(o,d));signs.append(exact_sign(p,triangles,xyz))
    if any(s is None or s==0 for s in signs):return dict(unresolved=True,roots=[],tangencies=[],coplanar=coplanar)
    roots=[];tangencies=[]
    for i,t in enumerate(cuts[1:-1]):(roots if signs[i]!=signs[i+1] else tangencies).append(float(t))
    if len(set(float(t) for t in cuts))!=len(cuts):return dict(unresolved=True,roots=[],tangencies=[],coplanar=coplanar)
    return dict(unresolved=False,roots=roots,tangencies=tangencies,coplanar=coplanar)
