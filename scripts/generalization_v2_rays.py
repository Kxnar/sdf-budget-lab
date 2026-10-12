"""Independent segmented field rays: algebraic polynomial roots, no sample scan."""
import struct
from fractions import Fraction
from generalization_common import *

def decoded_corners(values,encoding):
    if encoding=='f64':return values.astype(float)
    v=values.astype(np.float32).astype(float)
    if encoding=='f32':return v
    offset=v.min(axis=1);scale=((v.max(axis=1)-offset)/255).astype(np.float32).astype(float)
    q=np.clip(np.floor((v-offset[:,None])/np.where(scale>0,scale,1)[:,None]+.5),0,255)
    return (offset[:,None]+scale[:,None]*q).astype(np.float32).astype(float)

def polynomials(segments,ray,values=None):
    """Explicit multilinear-to-power-basis substitution, batched over cells."""
    v=segments[:,8:] if values is None else values
    # f(x,y,z)=a0+a1*x+a2*y+a3*z+a12*x*y+a13*x*z+a23*y*z+a123*x*y*z.
    a=np.stack([v[:,0],v[:,1]-v[:,0],v[:,2]-v[:,0],v[:,4]-v[:,0],v[:,3]-v[:,1]-v[:,2]+v[:,0],v[:,5]-v[:,1]-v[:,4]+v[:,0],v[:,6]-v[:,2]-v[:,4]+v[:,0],v[:,7]-v[:,3]-v[:,5]-v[:,6]+v[:,1]+v[:,2]+v[:,4]-v[:,0]],axis=1)
    width=segments[:,5:8]-segments[:,2:5];u=(ray[1:4]+segments[:,0,None]*ray[4:7]-segments[:,2:5])/width;b=(segments[:,1]-segments[:,0])[:,None]*ray[4:7]/width
    p=np.zeros((len(v),4));p[:,0]=a[:,0]
    for k,axes in enumerate(((0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)),1):
        for mask in range(1<<len(axes)):
            term=a[:,k].copy();power=0
            for j,axis in enumerate(axes):
                chosen=bool(mask&(1<<j));term*=b[:,axis] if chosen else u[:,axis];power+=chosen
            p[:,power]+=term
    return p

def trim(p):
    while len(p)>1 and p[-1]==0:p.pop()
    return p
def divide(a,b):
    a=a.copy();q=[Fraction(0)]*max(1,len(a)-len(b)+1)
    while len(a)>=len(b) and any(a):
        i=len(a)-len(b);v=a[-1]/b[-1];q[i]+=v
        for j in range(len(b)):a[i+j]-=v*b[j]
        trim(a)
    return trim(q),trim(a)
def gcd(a,b):
    while any(b):a,b=b,divide(a,b)[1]
    return [x/a[-1] for x in a]
def value(p,x):
    v=Fraction(0)
    for c in p[::-1]:v=v*x+c
    return v
def sturm(p):
    result=[p,[i*p[i] for i in range(1,len(p))]]
    while any(result[-1]):
        remainder=divide(result[-2],result[-1])[1]
        if not any(remainder):break
        result.append([-x for x in remainder])
    return result
def variation(seq,x):
    signs=[1 if v>0 else -1 for p in seq if (v:=value(p,x))!=0]
    return sum(a!=b for a,b in zip(signs,signs[1:]))
def isolated_factor_roots(p):
    # Strip exact endpoints; ownership belongs to the interface pass.
    for endpoint in (Fraction(0),Fraction(1)):
        if len(p)>1 and value(p,endpoint)==0:p=divide(p,[-endpoint,Fraction(1)])[0]
    if len(p)<=1:return []
    seq=sturm(p);count=variation(seq,Fraction(0))-variation(seq,Fraction(1))
    eigen=np.polynomial.polynomial.polyroots([float(x) for x in p]);rr=sorted(float(x.real) for x in eigen if abs(x.imag)<1e-12 and 0<x.real<1)
    if len(rr)==count and all(abs(value(p,Fraction(x)))<Fraction(1,10**10) for x in rr):return rr
    result=[];pending=[(Fraction(0),Fraction(1),count)]
    while pending:
        lo,hi,n=pending.pop()
        if not n:continue
        if n==1 and float(hi-lo)<1e-15:result.append(float((lo+hi)/2));continue
        mid=(lo+hi)/2
        if value(p,mid)==0:
            quotient=divide(p,[-mid,Fraction(1)])[0]
            return sorted([float(mid),*isolated_factor_roots(quotient)])
        left=variation(seq,lo)-variation(seq,mid);pending.extend([(lo,mid,left),(mid,hi,n-left)])
    return sorted(result)
def polynomial_roots(c):
    # Square-free factors have exact multiplicities for these independently
    # formed binary64 coefficients. Companion roots are Sturm-count checked;
    # rational interval isolation resolves a disagreement, without root merging.
    p=trim([Fraction(float(x)) for x in c])
    if len(p)==1:return [],[],p[0]==0
    common=gcd(p,[i*p[i] for i in range(1,len(p))]);w=divide(p,common)[0];multiplicity=1;crossings=[];tangencies=[]
    while len(w)>1:
        y=gcd(w,common);factor=divide(w,y)[0];roots=isolated_factor_roots(factor)
        (crossings if multiplicity%2 else tangencies).extend(roots)
        w=y;common=divide(common,y)[0];multiplicity+=1
    return sorted(crossings),sorted(tangencies),False

def endpoint_side(coefficients,endpoint,direction):
    p=[Fraction(float(c)) for c in coefficients];x=Fraction(endpoint)
    for order in range(4):
        v=value(p,x)
        if v:return (1 if v>0 else -1)*(direction**order)
        p=[i*p[i] for i in range(1,len(p))]
        if not p:break
    return 0

def transitions(segments,ray,encoding='f64'):
    values=decoded_corners(segments[:,8:],encoding);p=polynomials(segments,ray,values)
    crossings=[];tangencies=[];jumps=[];zero=int(np.sum(np.all(p==0,axis=1)));ambiguous=zero>0
    # A same-strict-sign cube cannot have an interior root while its ray stays in the cube.
    mixed=~(np.all(values>0,axis=1)|np.all(values<0,axis=1))
    for i in np.flatnonzero(mixed):
        rr,tt,aa=polynomial_roots(p[i]);crossings.extend(segments[i,0]+np.array(rr)*(segments[i,1]-segments[i,0]));tangencies.extend(segments[i,0]+np.array(tt)*(segments[i,1]-segments[i,0]));ambiguous|=aa
    left=p.sum(axis=1);right=p[:,0];scale=np.sum(np.abs(p),axis=1);uncertain=128*np.finfo(float).eps*scale
    for i in np.flatnonzero((np.sign(left[:-1])!=np.sign(right[1:])) | (left[:-1]==0) | (right[1:]==0)):
        if segments[i,1]!=segments[i+1,0]:ambiguous=True;continue
        a,b=endpoint_side(p[i],1,-1),endpoint_side(p[i+1],0,1)
        if a*b<0:
            crossings.append(segments[i,1])
            if left[i]!=0 or right[i+1]!=0:jumps.append(segments[i,1])
        elif a*b==0:ambiguous=True
        elif left[i]==0 and right[i+1]==0:tangencies.append(segments[i,1])
    # No proximity-based root coalescing: two close eigenvalue roots remain two.
    return dict(roots=np.sort(crossings),tangencies=np.sort(tangencies),jumps=np.sort(jumps),zero_segments=zero,ambiguous=ambiguous)

def packed_segments(data,ray):
    result=[]
    def visit(offset,lo,hi):
        low,high=ray[7:9]
        for axis in range(3):
            if ray[4+axis]==0:
                if ray[1+axis]<lo[axis] or ray[1+axis]>hi[axis] or (ray[1+axis]==hi[axis] and hi[axis]!=1):return
            else:
                times=(np.array([lo[axis],hi[axis]])-ray[1+axis])/ray[4+axis];low=max(low,times.min());high=min(high,times.max())
        if high<=low:return
        descriptor=struct.unpack_from('<I',data,offset)[0];tag=descriptor&3;payload=descriptor&~3
        if not tag:
            mid=(lo+hi)/2
            for c in range(8):side=np.array([bool(c&(1<<a)) for a in range(3)]);visit(payload+4*c,np.where(side,mid,lo),np.where(side,hi,mid))
        else:
            if tag==1:values=np.frombuffer(data,dtype='<f4',count=8,offset=payload).astype(float)
            else:
                offset_value,scale=struct.unpack_from('<ff',data,payload);q=np.frombuffer(data,dtype='<u2' if tag==2 else'u1',count=8,offset=payload+8);values=(offset_value+float(scale)*q.astype(float)).astype(np.float32).astype(float)
            result.append([low,high,*lo,*hi,*values])
    visit(0,np.full(3,-1.),np.full(3,1.));return np.array(sorted(result,key=lambda s:s[0]))

def signed_intervals(query,ray,truth,roots):
    if len(roots)!=len(truth) or not len(roots) or len(roots)%2:return False
    if not np.all(query(ray[1:4]+ray[[7,8]][:,None]*ray[4:7])>0):return False
    for rr in (truth,roots):
        edges=np.r_[ray[7],rr,ray[8]];mid=(edges[:-1]+edges[1:])/2;expected=np.where(np.arange(len(mid))%2,-1,1)
        if not np.all(query(ray[1:4]+mid[:,None]*ray[4:7])*expected>0):return False
    return True
