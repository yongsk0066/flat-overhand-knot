"""Historical line-intersection diagnostic. Its segment-parameter sign is incorrect; use ref12_trunc2.py for the clipping result."""
from probe_files import VERIFY, read_probe

import sympy as sp, sys, math
r15=sp.sqrt(15)
print("="*76); print("(A) L_trunc convention audit at theta*_{4_1}"); print("="*76)
# exact hexagon (scaled x = 30 x / sqrt15)
HEXs=[(-13,sp.Rational(1,2)),(-7,sp.Rational(3,2)),(1,sp.Rational(3,2)),
      (3,sp.Rational(1,2)),(-3,sp.Rational(-1,2)),(-11,sp.Rational(-1,2))]
HEX=[sp.Matrix([sp.Integer(x)*r15/30, y]) for x,y in HEXs]
M0=sp.Matrix([0,0]); M3=sp.Matrix([-r15/3,1])
d0=sp.Matrix([1,0]); d4=sp.Matrix([1,0])
def exit_param(P, dirv):
    """Historical intersection formula with a reversed edge parameter."""
    ts=[]
    for k in range(6):
        A=HEX[k]; B=HEX[(k+1)%6]; e=B-A
        den=dirv[0]*(-e[1])+dirv[1]*e[0]
        if sp.simplify(den)==0: continue
        rhs=A-P
        u=sp.simplify((rhs[0]*(-e[1])+rhs[1]*e[0])/den)          # P + u*dir  on the line AB
        if u is sp.nan: continue
        w=sp.simplify(((rhs[0]*(-dirv[1])+rhs[1]*dirv[0]))/(e[0]*(-dirv[1])+e[1]*dirv[0]))
        if sp.simplify(w)>=0 and sp.simplify(w)<=1 and sp.simplify(u)>0: ts.append(sp.simplify(u))
    return sp.simplify(sp.Max(*ts)) if ts else None
s_in=exit_param(M0, -d0)     # in-tail: ray from M0 along -d0
t_out=exit_param(M3, d4)     # out-tail: ray from M3 along d4
core=4*r15/3
computed = sp.simplify(s_in + core + t_out)
outline = 32*r15/15
print(f"  Historical signed-segment calculation: incoming {s_in}, outgoing {t_out}, length {computed}")
print(f"  Convex clipping (ref12_trunc2.py): {outline}")
print(f"  Historical result minus clipping result: {sp.simplify(computed - outline)}")
print(f"  Outermost centreline crossing convention: {28*r15/15}")

print("\nContact counts at the two mirror-related reference states")
sys.path.insert(0, str(VERIFY))
from n4_tighten2 import config, pieces, ell
t0=math.acos(-0.25); s0=2*math.pi-2*t0; L1=6/math.sqrt(15); L2=8/math.sqrt(15)
STAR=[t0,s0,-s0,-t0,L1,L2,L1]
MIRROR=[-t0,-s0,s0,t0,L1,L2,L1]
def contacts(th, order):
    """all 12 (crease, intervening pass, endpoint) slack values; >=0 required."""
    c=config(th)
    if c is None: return None
    segs=pieces(c); rank={p:i for i,p in enumerate(order)}; out=[]
    for i in range(4):
        ra,rb=rank[i],rank[i+1]; lo,hi=min(ra,rb),max(ra,rb)
        for k in [order[r] for r in range(lo+1,hi)]:
            P,dd,ln=segs[k]; n=(-dd[1],dd[0])
            for e in (0,1):
                p=c['CH'][i][e]
                la=(p[0]-P[0])*n[0]+(p[1]-P[1])*n[1]
                out.append(abs(la)-0.5)
    return out
ends = [row["ell"] for row in read_probe("ref10_endpoints.json")]
print(f"  (endpoint ell values recovered from the census run: {len(ends)})")
ELLS=4*math.sqrt(15)/3
print(f"\n  ell* = {ELLS:.9f};  mirror of theta* has the same ell by construction: "
      f"{abs(ell(MIRROR)-ELLS)<1e-12}")
print(f"\n  {'endpoint':>10} {'ell-ell*':>12} {'#contacts < 1e-4':>18} {'#contacts < 1e-2':>18}")
for nm,th,o in [("theta*",STAR,[1,4,2,0,3]),("mirror",MIRROR,[3,0,2,4,1])]:
    cs=contacts(tuple(th),o)
    print(f"  {nm:>10} {ell(th)-ELLS:>12.3e} {sum(1 for x in cs if abs(x)<1e-4):>18} {sum(1 for x in cs if abs(x)<1e-2):>18}")
print("  (theta* shows 8 active of 12 -- the certified active set)")
