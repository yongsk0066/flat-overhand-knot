# ref15_farkas.py — REFEREE: is each straggler a BLOCKED point?  Farkas test, oracle-free.
# blocked  <=>  grad(ell) in cone{ grad g_k : k active }.  If not, a first-order feasible
# descent direction exists and the point is a mere stall of the rejection-sampling descent.
import sys, math, json, itertools
import numpy as np
sys.path.insert(0,'/Users/yongseok/dev/pentagon/verify')
from n4_tighten2 import config, pieces, ell as ELL
ELLS=4*math.sqrt(15)/3
t0=math.acos(-0.25); s0=2*math.pi-2*t0; L1=6/math.sqrt(15); L2=8/math.sqrt(15)
STAR=[t0,s0,-s0,-t0,L1,L2,L1]; MIRR=[-t0,-s0,s0,t0,L1,L2,L1]

def gvec(th, order):
    c=config(tuple(th))
    segs=pieces(c); rank={p:i for i,p in enumerate(order)}; out=[]
    for i in range(4):
        ra,rb=rank[i],rank[i+1]; lo,hi=min(ra,rb),max(ra,rb)
        for k in [order[r] for r in range(lo+1,hi)]:
            P,dd,ln=segs[k]; n=(-dd[1],dd[0])
            for e in (0,1):
                p=c['CH'][i][e]
                out.append(abs((p[0]-P[0])*n[0]+(p[1]-P[1])*n[1])-0.5)
    return out

def analyse(th, order, tol=1e-4):
    g0=gvec(th,order); act=[k for k,v in enumerate(g0) if abs(v)<tol]
    h=1e-7; J=np.zeros((len(act),7))
    for j in range(7):
        tp=list(th); tm=list(th); tp[j]+=h; tm[j]-=h
        gp=gvec(tp,order); gm=gvec(tm,order)
        for c_,k in enumerate(act): J[c_,j]=(gp[k]-gm[k])/(2*h)
    b=np.array([0,0,0,0,1.,1.,1.])
    A=J.T                                    # 7 x m
    best=(np.inf,None)
    m=len(act)
    for r in range(1, m+1):
        for S in itertools.combinations(range(m), r):
            As=A[:,list(S)]
            lam,_,_,_=np.linalg.lstsq(As,b,rcond=None)
            if np.all(lam>=-1e-9):
                res=np.linalg.norm(As@lam-b)
                if res<best[0]: best=(res,(S,lam))
    # steepest feasible descent: minimise grad(ell).v s.t. Jv>=0, |v|<=1 (projected iteration)
    v=-b/np.linalg.norm(b)
    for _ in range(60000):
        viol=J@v
        bad=np.where(viol<0)[0]
        if len(bad)==0: break
        k=bad[np.argmin(viol[bad])]
        n=J[k]/np.linalg.norm(J[k])
        v=v-(v@n)*n
        nv=np.linalg.norm(v)
        if nv<1e-12: v=np.zeros(7); break
        v=v/nv
    dell = float(b@v) if np.linalg.norm(v)>1e-9 else 0.0
    feas = float(np.min(J@v)) if np.linalg.norm(v)>1e-9 else 0.0
    return act, best[0], dell, feas

print("Farkas / blockedness test  (residual ~0 and lambda>=0  <=>  BLOCKED)")
print(f"  {'point':>10} {'ell-ell*':>11} {'#active':>8} {'Farkas resid':>13} {'best d(ell) on cone':>21} {'min Jv':>9}")
for nm,th,o in (("theta*",STAR,[1,4,2,0,3]),("mirror",MIRR,[3,0,2,4,1])):
    act,res,dell,feas=analyse(th,o)
    print(f"  {nm:>10} {ELL(th)-ELLS:>+11.3e} {len(act):>8} {res:>13.3e} {dell:>21.4e} {feas:>9.1e}")
ends=json.load(open('/private/tmp/claude-501/-Users-yongseok-dev-pentagon/0c3f9975-7617-44d6-be99-ab855225d015/scratchpad/ref10_endpoints.json'))
print()
nb=0; ns=0
for k,e in enumerate(ends):
    d=e['ell']-ELLS
    if d<=1e-3: continue
    act,res,dell,feas=analyse(e['theta'], e['order'])
    blocked = res<1e-6
    ns+=1; nb+= blocked
    print(f"  {('stall#%d'%k):>10} {d:>+11.3e} {len(act):>8} {res:>13.3e} {dell:>21.4e} {feas:>9.1e}"
          f"   {'BLOCKED' if blocked else 'NOT blocked -> descent direction exists'}")
print(f"\n  stragglers that are genuinely blocked: {nb}/{ns}")
print("  'best d(ell) on cone' < 0 means a first-order feasible direction strictly decreases ell,")
print("  i.e. the point is a STALL of the rejection sampler and not a rival tight state.")
