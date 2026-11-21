"""Compare truncation rules on the trefoil and classify mirror-related endpoints."""
from probe_files import VERIFY, read_probe

import sys, math, json, collections
sys.path.insert(0, str(VERIFY))

print("="*76); print("(A) trefoil pentagon: outline rule vs outermost-crossing rule"); print("="*76)
import check_edge as CE
x0=math.sqrt(1-2/math.sqrt(5)); th=CE.slice_theta(x0); c=CE.config(th)
def sub(a,b): return (a[0]-b[0],a[1]-b[1])
def add(a,b): return (a[0]+b[0],a[1]+b[1])
def mul(a,s): return (a[0]*s,a[1]*s)
def meet(P1,d1,P2,d2):
    det=d1[0]*(-d2[1])-d1[1]*(-d2[0])
    if abs(det)<1e-14: return None,None
    rhs=sub(P2,P1)
    s=(rhs[0]*(-d2[1])-rhs[1]*(-d2[0]))/det
    u=(d1[0]*rhs[1]-d1[1]*rhs[0])/det
    return s,u
M1,M2,M3=c['M1'],c['M2'],c['M3']
_,_=0,0
s_in,_=meet(M1,mul(c['d2'],-1.0), c['G'], sub(c['Q1'],c['G']))
t_out,_=meet(M3,c['v3'],(0.0,0.0),(1.0,0.0))
ell=math.hypot(*sub(M2,M1))+math.hypot(*sub(M3,M2))
print(f"  outline rule (as implemented for n=3): s_in={s_in:.9f}  core={ell:.9f}  t_out={t_out:.9f}")
print(f"     L_trunc = {s_in+ell+t_out:.10f}   Kauffman 4/sqrt(5-2sqrt5) = {4/math.sqrt(5-2*math.sqrt(5)):.10f}")
# the trefoil passes: P0 = ray back from M1 along -d2, P1=[M1,M2], P2=[M2,M3], P3 = ray from M3 along v3
segs={0:(M1,mul(c['d2'],-1.0),60.0), 1:(M1,mul(sub(M2,M1),1/math.hypot(*sub(M2,M1))),math.hypot(*sub(M2,M1))),
      2:(M2,mul(sub(M3,M2),1/math.hypot(*sub(M3,M2))),math.hypot(*sub(M3,M2))),
      3:(M3,mul(c['v3'],1/math.hypot(*c['v3'])),60.0)}
outer_in=0.0; outer_out=0.0
for j in (2,3):
    P1_,d1_,l1=segs[0]; P2_,d2_,l2=segs[j]
    s,u=meet(P1_,d1_,P2_,d2_)
    if s is not None and 1e-9<s<l1 and 1e-9<u<l2: outer_in=max(outer_in,s)
for j in (0,1):
    P1_,d1_,l1=segs[3]; P2_,d2_,l2=segs[j]
    s,u=meet(P1_,d1_,P2_,d2_)
    if s is not None and 1e-9<s<l1 and 1e-9<u<l2: outer_out=max(outer_out,s)
print(f"  outermost centreline CROSSING rule: s_in={outer_in:.9f}  t_out={outer_out:.9f}")
print(f"     L_trunc = {outer_in+ell+outer_out:.10f}")
print(f"  the two rules coincide at the pentagon: {abs(outer_in-s_in)<1e-9 and abs(outer_out-t_out)<1e-9}")
print(f"  difference per stub: in {s_in-outer_in:+.9f}, out {t_out-outer_out:+.9f}")

print(); print("="*76); print("(B) straggler diagnosis"); print("="*76)
from n4_tighten2 import config, pieces, ell as ELL
t0=math.acos(-0.25); s0=2*math.pi-2*t0; L1=6/math.sqrt(15); L2=8/math.sqrt(15)
STAR={(1,4,2,0,3):[t0,s0,-s0,-t0,L1,L2,L1], (3,0,2,4,1):[-t0,-s0,s0,t0,L1,L2,L1]}
ELLS=4*math.sqrt(15)/3
def contacts(th, order):
    cc=config(tuple(th))
    if cc is None: return None
    segs=pieces(cc); rank={p:i for i,p in enumerate(order)}; out=[]
    for i in range(4):
        ra,rb=rank[i],rank[i+1]; lo,hi=min(ra,rb),max(ra,rb)
        for k in [order[r] for r in range(lo+1,hi)]:
            P,dd,ln=segs[k]; n=(-dd[1],dd[0])
            for e in (0,1):
                p=cc['CH'][i][e]
                out.append(abs((p[0]-P[0])*n[0]+(p[1]-P[1])*n[1])-0.5)
    return out
for nm,o in (("theta*",(1,4,2,0,3)),("mirror",(3,0,2,4,1))):
    cs=contacts(STAR[o],list(o))
    print(f"  reference {nm:8}: ell-ell*={ELL(STAR[o])-ELLS:+.2e}  active(<1e-9)={sum(1 for x in cs if abs(x)<1e-9)}/12")
ends=read_probe("ref10_endpoints.json")
print(f"\n  {'#':>3} {'ell-ell*':>11} {'|dtheta| to own-branch theta*':>30} {'active<1e-4':>12} {'active<1e-2':>12}")
worst=[]
for k,e in enumerate(ends):
    o=tuple(e['order']); th=e['theta']
    ref=STAR[o]
    dd=max(abs(th[j]-ref[j]) for j in range(7))
    cs=contacts(th,list(o))
    a4=sum(1 for x in cs if abs(x)<1e-4); a2=sum(1 for x in cs if abs(x)<1e-2)
    print(f"  {k:>3} {e['ell']-ELLS:>+11.3e} {dd:>30.4e} {a4:>12} {a2:>12}")
    worst.append((e['ell']-ELLS, dd, a4, a2))
far=[w for w in worst if w[0]>1e-3]
print(f"\n  endpoints stalled above ell*+1e-3: {len(far)}/{len(ends)}")
print(f"  max |dtheta| among them: {max((w[1] for w in far), default=0):.4e}"
      f"   (distance to the reference with the same layer order)")
print(f"  active contacts (<1e-4) among the stalls: {[w[2] for w in far]}")
print(f"  none of the stalls reaches 8 active contacts: {all(w[2]<8 for w in far)}")
