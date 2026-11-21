"""Compare boundary tolerances and knot-type oracles at recorded states."""
from probe_files import VERIFY, read_probe

import sys, math, random, json, collections
sys.path.insert(0, str(VERIFY))
import n4_tighten2 as A
from n4_tighten2 import config, pieces, crossings, code_str, canon, ev_over_of, fox5_events, ell as ELL
ELLS=4*math.sqrt(15)/3
t0=math.acos(-0.25); s0=2*math.pi-2*t0
STAR=(t0,s0,-s0,-t0,6/math.sqrt(15),8/math.sqrt(15),6/math.sqrt(15)); O=[1,4,2,0,3]

DELTA=1e-12
def axiom_fixed(c, order, tol=1e-11, delta=DELTA):
    """same interval logic, but the slab is shrunk by delta so that exact boundary contact
    (|lat| == 1/2) is never mistaken for penetration by round-off."""
    segs=pieces(c); rank={p:q for q,p in enumerate(order)}
    for i in range(4):
        ra,rb=rank[i],rank[i+1]; lo,hi=min(ra,rb),max(ra,rb)
        between=[order[r] for r in range(lo+1,hi)]
        if not between: continue
        (ax,ay),(bx,by)=c['CH'][i]
        for k in between:
            P,dd,l=segs[k]; n=(-dd[1],dd[0])
            latA=(ax-P[0])*n[0]+(ay-P[1])*n[1]; latB=(bx-P[0])*n[0]+(by-P[1])*n[1]
            lonA=(ax-P[0])*dd[0]+(ay-P[1])*dd[1]; lonB=(bx-P[0])*dd[0]+(by-P[1])*dd[1]
            def clip(lo_,hi_,c0,c1,f0,f1):
                d=c1-c0
                if abs(d)<1e-15:
                    if not (lo_<c0<hi_): return 1.0,0.0
                    return f0,f1
                fa=(lo_-c0)/d; fb=(hi_-c0)/d
                if fa>fb: fa,fb=fb,fa
                return max(f0,fa),min(f1,fb)
            g0,g1=clip(-0.5+delta,0.5-delta,latA,latB,0.0,1.0)
            if g0>=g1: continue
            g0,g1=clip(0.0,l,lonA,lonB,g0,g1)
            if g1-g0>tol: return False
    return True

def valid_fixed3(th, order):
    c=config(th)
    if c is None: return None
    X=crossings(c); cd=code_str(X)
    if cd is None or canon(cd)!=A.TARGET_CANON: return None
    if fox5_events(ev_over_of(X,order))!=25: return None
    if not axiom_fixed(c,order): return None
    return dict(c=c,X=X)

print("Actual oblique panels and a rectangular slab approximation at the reference state")
print(f"  axiom_exact_ok(theta*) = {A.axiom_exact_ok(config(STAR),O)}   (actual panel geometry)")
print(f"  axiom_fixed   (theta*) = {axiom_fixed(config(STAR),O)}   (slab shrunk by delta={DELTA:.0e})")

rnd=random.Random(4242)
def rate(oracle, radius, n=400):
    ok=0
    for _ in range(n):
        v=[rnd.gauss(0,1) for _ in range(7)]; nn=math.sqrt(sum(x*x for x in v))
        th2=tuple(STAR[j]+radius*v[j]/nn for j in range(7))
        if oracle(th2,O) is not None: ok+=1
    return ok
for r in (1e-3,1e-4,1e-5):
    print(f"  neighbourhood validity rate at radius {r:.0e}: actual panels {rate(A.valid_fixed2,r)}/400"
          f"   rectangular slabs {rate(valid_fixed3,r)}/400")

# Descend using the rectangular slab approximation.
A.valid_fixed2 = valid_fixed3           # tighten_fixed2 calls valid_fixed2 by name
from n4_tighten2 import tighten_fixed2
independent_hits=read_probe("ref_hits.json")
starts=[(tuple(t),o) for t,o,cd in independent_hits]
import multiprocessing as mp
def dw(a):
    seed,th0,order=a
    rnd=random.Random(seed)
    thf,ef=tighten_fixed2(th0,order,rnd)
    for x in range(4):
        t2,e2=tighten_fixed2(list(thf),order,random.Random(seed*13+x),rounds=10)
        if e2<ef: thf,ef=t2,e2
    return ef
if __name__=="__main__":
    with mp.get_context("fork").Pool(8) as pool:
        res=sorted(pool.map(dw,[(990000+i,th,o) for i,(th,o) in enumerate(starts)]))
    print(f"\n[descent with the rectangular slab approximation] {len(res)} independent starts")
    print(f"  best ell = {res[0]:.9f}   ell* = {ELLS:.9f}   gap = {res[0]-ELLS:.3e}")
    print(f"  within 1e-4 of ell*: {sum(1 for e in res if e-ELLS<1e-4)}/{len(res)};"
          f"  within 1e-5: {sum(1 for e in res if e-ELLS<1e-5)}/{len(res)};"
          f"  within 0.01: {sum(1 for e in res if e-ELLS<0.01)}/{len(res)}")
    print("  gaps:", [f"{e-ELLS:.2e}" for e in res])
