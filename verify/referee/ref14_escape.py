"""Probe descent from stalled endpoints with perturbed searches."""
from probe_files import VERIFY, read_probe

import sys, math, json, random
sys.path.insert(0, str(VERIFY))
from n4_tighten2 import config, pieces, valid_fixed2, ell as ELL
ELLS=4*math.sqrt(15)/3
t0=math.acos(-0.25); s0=2*math.pi-2*t0; L1=6/math.sqrt(15); L2=8/math.sqrt(15)
STAR=[t0,s0,-s0,-t0,L1,L2,L1]; MIRR=[-t0,-s0,s0,t0,L1,L2,L1]

def gvec(th, order):
    c=config(tuple(th))
    if c is None: return None
    segs=pieces(c); rank={p:i for i,p in enumerate(order)}; out=[]
    for i in range(4):
        ra,rb=rank[i],rank[i+1]; lo,hi=min(ra,rb),max(ra,rb)
        for k in [order[r] for r in range(lo+1,hi)]:
            P,dd,ln=segs[k]; n=(-dd[1],dd[0])
            for e in (0,1):
                p=c['CH'][i][e]
                out.append(abs((p[0]-P[0])*n[0]+(p[1]-P[1])*n[1])-0.5)
    return out

def escape(th0, order, seed, tries=40000):
    """search for a feasible point with strictly smaller ell near th0"""
    rnd=random.Random(seed); e0=ELL(th0); best=None
    g0=gvec(th0,order); act=[k for k,v in enumerate(g0) if abs(v)<1e-4]
    # numeric gradients of the active contacts
    h=1e-7; Gr=[]
    for j in range(7):
        tp=list(th0); tm=list(th0); tp[j]+=h; tm[j]-=h
        gp=gvec(tp,order); gm=gvec(tm,order)
        Gr.append([(gp[k]-gm[k])/(2*h) for k in act])
    for _ in range(tries):
        v=[rnd.gauss(0,1) for _ in range(7)]
        # bias toward decreasing ell
        for j in (4,5,6): v[j]-=abs(rnd.gauss(0,1))*0.7
        n=math.sqrt(sum(x*x for x in v)); v=[x/n for x in v]
        if sum(v[4:])>=0: continue                      # need d(ell)<0
        # first-order feasibility wrt the active contacts
        if any(sum(Gr[j][ai]*v[j] for j in range(7)) < -1e-9 for ai in range(len(act))): continue
        for r in (3e-3,1e-3,3e-4,1e-4,3e-5,1e-5,3e-6):
            th2=[th0[j]+r*v[j] for j in range(7)]
            if any(th2[j]<0.05 for j in (4,5,6)): continue
            if valid_fixed2(tuple(th2),order) is not None:
                de=ELL(th2)-e0
                if de<-1e-12 and (best is None or de<best[0]): best=(de,r)
                break
    return e0, act, best

ends=read_probe("ref10_endpoints.json")
print("Reference: does the SAME search escape from the certified tight state theta*?")
for nm,th,o in (("theta*",STAR,[1,4,2,0,3]),("mirror",MIRR,[3,0,2,4,1])):
    e0,act,best=escape(th,o,11,tries=40000)
    print(f"  {nm:8}: ell-ell*={e0-ELLS:+.3e}  active={len(act)}  escape found: "
          f"{'YES d(ell)=%.3e'%best[0] if best else 'NO'}")
print("\nStragglers (endpoints above ell*+1e-3):")
print(f"  {'#':>3} {'ell-ell*':>11} {'dist to theta*/mirror':>22} {'#active':>8}  escape?")
nesc=0; nstall=0
for k,e in enumerate(ends):
    d=e['ell']-ELLS
    if d<=1e-3: continue
    th=e['theta']; o=e['order']
    dist=min(max(abs(th[j]-R[j]) for j in range(7)) for R in (STAR,MIRR))
    e0,act,best=escape(th,o,1000+k,tries=40000)
    nstall+=1
    if best: nesc+=1
    print(f"  {k:>3} {d:>+11.3e} {dist:>22.4e} {len(act):>8}  "
          f"{'ESCAPED d(ell)=%.2e at r=%.0e'%(best[0],best[1]) if best else 'no escape found'}")
print(f"\n  A lower-length sample passed the numerical panel tests for {nesc}/{nstall} endpoints.")
