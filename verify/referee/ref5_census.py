"""Independent four-crossing harvest with a finite-tail comparison."""
from probe_files import VERIFY, write_probe

import sys, math, random, collections, json
sys.path.insert(0, str(VERIFY))
import n4_tighten2 as A
from n4_tighten2 import (config, code_str, canon, events, linear_extensions,
                         fox5_events, axiom_exact_ok, ell, TARGET_CANON)
import multiprocessing as mp

def crossings_R(c, RAY):
    """Intersections after truncating both tail rays to the specified length."""
    d,M = c['d'], c['M']
    segs={0:(M[0],(-d[0][0],-d[0][1]),RAY),1:(M[0],d[1],c['L'][0]),2:(M[1],d[2],c['L'][1]),
          3:(M[2],d[3],c['L'][2]),4:(M[3],d[4],RAY)}
    X=[]
    for i in range(5):
        for j in range(i+2,5):
            P1,d1,l1=segs[i]; P2,d2,l2=segs[j]
            det=d1[0]*d2[1]-d1[1]*d2[0]
            if abs(det)<1e-12: continue
            rx,ry=P2[0]-P1[0],P2[1]-P1[1]
            s=(rx*d2[1]-ry*d2[0])/det; u=(rx*d1[1]-ry*d1[0])/det
            if 1e-9<s<l1-1e-9 and 1e-9<u<l2-1e-9: X.append(((i,j),s,u))
    return X

def orders41_general(X):
    """all over/under masks giving fox5=25 AND an acyclic layer graph — no canon assumption."""
    res=[]
    n=len(X)
    for mask in range(1<<n):
        pairs=[]
        for cid,((i,j),_,_) in enumerate(X):
            pairs.append((i,j) if (mask>>cid)&1 else (j,i))
        exts=linear_extensions(pairs)
        if not exts: continue
        evo=[(cid,(side=='i')==((mask>>cid)&1==1)) for _,_,cid,side in events(X)]
        if fox5_events(evo)==25: res.append((mask,exts))
    return res

def worker(args):
    seed,n=args
    rnd=random.Random(seed)
    shadow=collections.Counter()          # canon -> count (4 crossings, any shadow)
    valid_by_canon=collections.Counter()  # canon -> count of fox5=25 + exact axiom
    cells=collections.Counter()
    ray_disagree=0; ray_examples=[]
    hits=[]
    for _ in range(n):
        if rnd.random()<0.5:
            rho=rnd.choice([1,-1]); sg=[rho]*4
        else: sg=[rnd.choice([1,-1]) for _ in range(4)]
        th=tuple([sg[k]*rnd.uniform(0.3,3.1) for k in range(4)]+[rnd.uniform(0.3,4.0) for _ in range(3)])
        c=config(th)
        if c is None: continue
        X=crossings_R(c,60.0)
        XL=crossings_R(c,1e7)
        if len(X)!=len(XL):
            ray_disagree+=1
            if len(ray_examples)<5: ray_examples.append((list(th),len(X),len(XL)))
        if len(X)!=4: continue
        cd=canon(code_str(X)); shadow[cd]+=1
        for m,exts in orders41_general(X):
            done=False
            for o in exts:
                if axiom_exact_ok(c,o):
                    valid_by_canon[cd]+=1
                    cells[(tuple(sorted(p for p,_,_ in X)),tuple(o))]+=1
                    hits.append((list(th),o,cd)); done=True; break
            if done: break
    return shadow, valid_by_canon, cells, ray_disagree, ray_examples, hits

if __name__=="__main__":
    TOTAL=8_000_000; P=8
    with mp.get_context("fork").Pool(P) as pool:
        parts=pool.map(worker,[(424242+i,TOTAL//P) for i in range(P)])
    S=collections.Counter(); V=collections.Counter(); C=collections.Counter()
    RD=0; RE=[]; H=[]
    for s,v,c,rd,re,h in parts:
        S.update(s); V.update(v); C.update(c); RD+=rd; RE+=re[:2]; H+=h
    print(f"[independent harvest] {TOTAL:,} samples, seed family 424242 (the recorded sweep starts at 7000)")
    print(f"\n(a) shadow-canon distribution over ALL 4-crossing configurations (no prefilter):")
    for cd,n in S.most_common(20): print(f"    {cd}: {n}")
    print(f"\n(b) sampled shadows admitting a Fox-5 layering and passing the panel predicate")
    for cd,n in V.most_common(20):
        print(f"    {cd}: {n}   {'<= canonical target' if cd==TARGET_CANON else '*** NON-TARGET SHADOW WITH A VALID 4_1 STATE ***'}")
    if not V: print("    (none)")
    print(f"\n(c) cells (pass-pairs, layer order) among valid states:")
    for k,n in C.most_common(20): print(f"    {k}: {n}")
    print(f"\n(d) RAY-truncation attack: RAY=60 vs RAY=1e7 gave a different crossing count in "
          f"{RD} / {TOTAL} sampled configurations")
    for e in RE[:6]: print(f"    example: n60={e[1]} n1e7={e[2]}  theta={[round(v,3) for v in e[0]]}")
    write_probe("ref_hits.json", H[:300])
    print(f"\n  valid states harvested: {len(H)}")
