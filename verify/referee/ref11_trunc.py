# SUPERSEDED by ref12_trunc2.py — the edge-parameter test in exit_param() has a sign error
# and accepts intersections with edge LINES outside the segment.  Kept only for the record.
# ref11_trunc.py — REFEREE:
#  (A) L_trunc at theta*_{4_1} under the paper's OWN convention (tail cut where the centreline
#      LEAVES THE CORE OUTLINE, doc 11 SS2 / Def 7.1 / Remark rem:trunclines) vs. the rule
#      actually coded in n4_exact_kkt.py ("outermost crossing with another pass").
#  (B) straggler diagnosis: mirror-aware distance to theta*, and active-set / Farkas test
#      (a genuine second attractor must be a BLOCKED point, not just a stalled one).
import sympy as sp, sys, math, json, itertools
r15=sp.sqrt(15)
print("="*76); print("(A) L_trunc convention audit at theta*_{4_1}"); print("="*76)
# exact hexagon (scaled x = 30 x / sqrt15)
HEXs=[(-13,sp.Rational(1,2)),(-7,sp.Rational(3,2)),(1,sp.Rational(3,2)),
      (3,sp.Rational(1,2)),(-3,sp.Rational(-1,2)),(-11,sp.Rational(-1,2))]
HEX=[sp.Matrix([sp.Integer(x)*r15/30, y]) for x,y in HEXs]
M0=sp.Matrix([0,0]); M3=sp.Matrix([-r15/3,1])
d0=sp.Matrix([1,0]); d4=sp.Matrix([1,0])
def exit_param(P, dirv):
    """largest t>0 with P + t*dirv on the hexagon boundary (convex => the exit point)"""
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
print(f"  in-tail  exits the hexagon at material distance s_in  = {sp.radsimp(s_in)} = {sp.N(s_in,10)}")
print(f"  out-tail exits the hexagon at material distance t_out = {sp.radsimp(t_out)} = {sp.N(t_out,10)}")
print(f"  s_in == t_out == L1 == L3 == 2*sqrt(15)/5 : {sp.simplify(s_in-2*r15/5)==0 and sp.simplify(t_out-2*r15/5)==0}")
Ltr_outline=sp.radsimp(sp.simplify(s_in+core+t_out))
print(f"\n  L_trunc (paper's outline convention) = {Ltr_outline} = {sp.N(Ltr_outline,15)}")
print(f"  L_trunc (n4_exact_kkt.py, outermost CROSSING) = {sp.radsimp(28*r15/15)} = {sp.N(28*r15/15,15)}")
print(f"  ratio L_trunc/core (outline) = {sp.simplify(Ltr_outline/core)}")
print("  the two rules DISAGREE:", sp.simplify(Ltr_outline-28*r15/15)!=0)
print("  reason: along P_0 the outermost centreline CROSSING (with P_3) is at 4*sqrt15/15,")
print("  but the tail leaves the hexagon later, through the BAND EDGE E3+ -> E4- , at 6*sqrt15/15.")
# does the confusion also occur for the trefoil?  (check where the pentagon's rule lands)
print("\n  trefoil cross-check (why the conflation was not caught at n=3):")
sys.path.insert(0,'/Users/yongseok/dev/pentagon')
print("   doc 11 SS2 defines the stub as reaching 'the point where the centreline leaves the")
print("   core outline'; exact_kauffman.py implements it as the meet with the OUTLINE LINE")
print("   l2=[G,Q1] and with {y=0}=[E,H] -- outline edges, not centreline crossings. So the")
print("   n=3 code follows the outline rule and the n=4 code does not.")

print()
print("="*76); print("(B) straggler diagnosis: stalls or genuine second attractors?"); print("="*76)
sys.path.insert(0,'/Users/yongseok/dev/pentagon/verify')
from n4_tighten2 import config, pieces, valid_fixed2, ell
t0=math.acos(-0.25); s0=2*math.pi-2*t0; L1=6/math.sqrt(15); L2=8/math.sqrt(15)
STAR=[t0,s0,-s0,-t0,L1,L2,L1]
MIRROR=[-t0,-s0,s0,t0,L1,L2,L1]
ORD={ (1,4,2,0,3):[1,4,2,0,3], (3,0,2,4,1):[3,0,2,4,1] }
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
res=json.load(open('/private/tmp/claude-501/-Users-yongseok-dev-pentagon/0c3f9975-7617-44d6-be99-ab855225d015/scratchpad/ref10_endpoints.json')) \
    if False else None
# re-read the endpoints printed by ref10 from its output file
import re
ends=[]
for line in open('/private/tmp/claude-501/-Users-yongseok-dev-pentagon/0c3f9975-7617-44d6-be99-ab855225d015/scratchpad/ref10_out.txt'):
    m=re.match(r"\s+(\d+): ell=([\d.]+)\s+ell-ell\*=([-+\d.e]+)\s+t4",line)
    if m: ends.append((int(m.group(1)), float(m.group(2))))
print(f"  (endpoint ell values recovered from the census run: {len(ends)})")
print("  NOTE: ref10 printed |theta-theta*| without accounting for the MIRROR branch;")
print("  states harvested with layer order [3,0,2,4,1] converge to the mirror of theta*,")
print("  whose tau signs are flipped -- hence the spurious |dtheta|~5.27 entries.")
ELLS=4*math.sqrt(15)/3
print(f"\n  ell* = {ELLS:.9f};  mirror of theta* has the same ell by construction: "
      f"{abs(ell(MIRROR)-ELLS)<1e-12}")
print("  Any endpoint whose ell exceeds ell* must, to be a SECOND ATTRACTOR, be a blocked")
print("  (KKT) point: i.e. carry enough ACTIVE contacts to span grad(ell).  A stall does not.")
print(f"\n  {'endpoint':>10} {'ell-ell*':>12} {'#contacts < 1e-4':>18} {'#contacts < 1e-2':>18}")
for nm,th,o in [("theta*",STAR,[1,4,2,0,3]),("mirror",MIRROR,[3,0,2,4,1])]:
    cs=contacts(tuple(th),o)
    print(f"  {nm:>10} {ell(th)-ELLS:>12.3e} {sum(1 for x in cs if abs(x)<1e-4):>18} {sum(1 for x in cs if abs(x)<1e-2):>18}")
print("  (theta* shows 8 active of 12 -- the certified active set)")
