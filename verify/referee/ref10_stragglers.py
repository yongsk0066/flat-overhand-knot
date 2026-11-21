"""Recount sampled states and descend from the independent harvest."""
from probe_files import VERIFY, read_probe, write_probe

import sys, math, random, json, collections
sys.path.insert(0, str(VERIFY))
from n4_tighten2 import (config, crossings, code_str, canon, orders_making_41,
                         axiom_exact_ok, tighten_fixed2, valid_fixed2, ell, TARGET_CANON)
import multiprocessing as mp
ELL_STAR=4*math.sqrt(15)/3

# ---------- (i) does the _hworker bookkeeping double count? ----------
def hw(args):
    seed,n=args; rnd=random.Random(seed); hits=[]
    for _ in range(n):
        if rnd.random()<0.5:
            rho=rnd.choice([1,-1]); sg=[rho]*4
        else: sg=[rnd.choice([1,-1]) for _ in range(4)]
        th=tuple([sg[k]*rnd.uniform(0.3,3.1) for k in range(4)]+[rnd.uniform(0.3,4.0) for _ in range(3)])
        c=config(th)
        if c is None: continue
        X=crossings(c); cd=code_str(X)
        if cd is None or canon(cd)!=TARGET_CANON: continue
        for m,exts in orders_making_41(X):        # <-- the original loop, no outer break
            for o in exts:
                if axiom_exact_ok(c,o): hits.append((list(th),o)); break
    return hits

def dw(args):
    seed,th0,order=args
    rnd=random.Random(seed)
    thf,ef=tighten_fixed2(th0,order,rnd)
    for extra in range(4):                                   # heavy re-polish
        t2,e2=tighten_fixed2(list(thf),order,random.Random(seed*7+extra+1),rounds=10)
        if e2<ef: thf,ef=t2,e2
    return ef,list(thf),order

if __name__=="__main__":
    P=8
    with mp.get_context("fork").Pool(P) as pool:
        parts=pool.map(hw,[(7000+i,2_000_000//P) for i in range(P)])   # seeds 7000 onward
    hits=[h for p in parts for h in p]
    key=lambda t: tuple(round(v,10) for v in t)
    dis=collections.Counter(key(t) for t,o in hits)
    ords=collections.Counter(tuple(o) for t,o in hits)
    print("="*74)
    print("(i) harvest bookkeeping, reproduced on 2M samples (seeds 7000 onward)")
    print("="*74)
    print(f"  hit records: {len(hits)}   distinct theta: {len(dis)}")
    print(f"  records-per-theta histogram: {dict(collections.Counter(dis.values()))}")
    print(f"  layer orders: {dict(ords)}")
    print("  => each geometric state is emitted once per over/under mask (the mirror pair),")
    print("     i.e. twice, once per top-bottom flip of the SAME configuration.")

    independent_hits=read_probe("ref_hits.json")
    starts=[(tuple(t),o) for t,o,cd in independent_hits]
    print()
    print("="*74)
    print(f"(ii) descent from {len(starts)} independent starts (seed 424242), heavy re-polish")
    print("="*74)
    with mp.get_context("fork").Pool(P) as pool:
        res=pool.map(dw,[(556677+i,th,o) for i,(th,o) in enumerate(starts)])
    res.sort(key=lambda r:r[0])
    cl=collections.Counter(round(r[0],2) for r in res)
    print(f"  endpoint ell distribution (0.01): {dict(sorted(cl.items()))}")
    print(f"  ell* = {ELL_STAR:.9f}")
    for k,(ef,thf,o) in enumerate(res):
        d=ef-ELL_STAR
        sym=f"t4+t1={thf[3]+thf[0]:+.2e} t3+t2={thf[2]+thf[1]:+.2e} L3-L1={thf[6]-thf[4]:+.2e}"
        print(f"   {k:2}: ell={ef:.7f}  ell-ell*={d:+.3e}  {sym}")
    hit=sum(1 for r in res if abs(r[0]-ELL_STAR)<0.01)
    hit3=sum(1 for r in res if abs(r[0]-ELL_STAR)<1e-3)
    print(f"  within 0.01 of ell*: {hit}/{len(res)};  within 1e-3: {hit3}/{len(res)}")
    t0=math.acos(-0.25); s0=2*math.pi-2*t0
    TH0=[t0,s0,-s0,-t0,6/math.sqrt(15),8/math.sqrt(15),6/math.sqrt(15)]
    mirror = [-v for v in TH0[:4]] + TH0[4:]
    write_probe("ref10_endpoints.json", [{"ell": ef, "theta": th, "order": o} for ef, th, o in res])
    print("\n  distance to the reference for the endpoint's layer order (max-norm):")
    for k,(ef,thf,o) in enumerate(res):
        reference = {(1,4,2,0,3): TH0, (3,0,2,4,1): mirror}[tuple(o)]
        dd=max(abs(thf[j]-reference[j]) for j in range(7))
        print(f"   {k:2}: ell-ell*={ef-ELL_STAR:+.3e}   |theta-theta*|_inf = {dd:.4e}")
