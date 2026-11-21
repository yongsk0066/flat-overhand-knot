"""Reconstruct the figure-eight contact gradients and KKT multipliers."""
from probe_files import VERIFY

import sympy as sp, itertools

r15 = sp.sqrt(15)
T = sp.symbols('T1 T2 T3 T4', real=True)
LL = sp.symbols('L1 L2 L3', positive=True)
TH = list(T)+list(LL)

phi = [sp.Integer(0)]
for t in T: phi.append(phi[-1]+t)
d = [sp.Matrix([sp.cos(p), sp.sin(p)]) for p in phi]
M = [sp.Matrix([0,0])]
for i in range(3): M.append(M[-1]+LL[i]*d[i+1])

# fixed signs of (d_i+d_{i+1}).n_i in a neighbourhood of theta* (verified numerically below)
SIG = [1,1,-1,-1]
CH = []
for i in range(4):
    v = d[i]+d[i+1]; n = sp.Matrix([-d[i][1], d[i][0]])
    hv = v/(2*SIG[i]*v.dot(n))
    CH.append((M[i]-hv, M[i]+hv))

tstar = sp.acos(sp.Rational(-1,4)); sstar = 2*sp.pi-2*tstar
STAR = {T[0]:tstar, T[1]:sstar, T[2]:-sstar, T[3]:-tstar,
        LL[0]:6/r15, LL[1]:8/r15, LL[2]:6/r15}
SUB = {sp.cos(T[0]):sp.Rational(-1,4), sp.sin(T[0]):r15/4,
       sp.cos(T[1]):sp.Rational(-7,8), sp.sin(T[1]):r15/8,
       sp.cos(T[2]):sp.Rational(-7,8), sp.sin(T[2]):-r15/8,
       sp.cos(T[3]):sp.Rational(-1,4), sp.sin(T[3]):-r15/4,
       LL[0]:6/r15, LL[1]:8/r15, LL[2]:6/r15}
def EV(e): return sp.radsimp(sp.simplify(sp.expand_trig(sp.expand(e)).subs(SUB)))
def NUM(e): return complex(sp.N(e.subs(STAR),30)).real

print("[sign check] (d_i+d_{i+1}).n_i at theta*:",
      [sp.nsimplify(EV((d[i]+d[i+1]).dot(sp.Matrix([-d[i][1],d[i][0]])))) for i in range(4)],
      " -> assumed signs", SIG)

# passes as (basepoint, unit direction, length); P0 runs backwards along -d0
PASS = {0:(M[0], -d[0], sp.Integer(60)), 1:(M[0], d[1], LL[0]), 2:(M[1], d[2], LL[1]),
        3:(M[2], d[3], LL[2]), 4:(M[3], d[4], sp.Integer(60))}
lat = lambda p,P,dd: -(p[0]-P[0])*dd[1]+(p[1]-P[1])*dd[0]
lon = lambda p,P,dd:  (p[0]-P[0])*dd[0]+(p[1]-P[1])*dd[1]

ORDER = [1,4,2,0,3]; rank = {p:i for i,p in enumerate(ORDER)}
print("\n[axiom scan] every (crease, intervening pass, endpoint) triple, exact lat/lon at theta*")
G = []
for i in range(4):
    ra, rb = rank[i], rank[i+1]; lo, hi = min(ra,rb), max(ra,rb)
    between = [ORDER[r] for r in range(lo+1,hi)]
    for k in between:
        P, dd, ln = PASS[k]
        la = [EV(lat(CH[i][e],P,dd)) for e in (0,1)]
        lo_ = [EV(lon(CH[i][e],P,dd)) for e in (0,1)]
        sgn = 1 if min(NUM(la[0]),NUM(la[1])) > 0 else -1     # which side of the slab
        for e in (0,1):
            g = sgn*lat(CH[i][e],P,dd) - sp.Rational(1,2)
            gv = EV(g)
            inwin = 0 < NUM(lo_[e]) < NUM(ln)
            tag = "ACTIVE" if gv==0 else f"slack {gv}"
            print(f"  crease{i+1} vs P{k} end{e}: lat={la[e]}  lon={sp.nsimplify(lo_[e])} in(0,{sp.nsimplify(ln)})={inwin}  {tag}")
            if gv==0:
                assert inwin, "ACTIVE contact with endpoint OUTSIDE the longitudinal window!"
                G.append((f"c{i+1}/P{k}/e{e}", g))
print(f"\n  -> {len(G)} active constraints found by exhaustive scan (eight at the reference state)")

# Exact Jacobian of the active contact functions.
J = sp.Matrix([[EV(sp.diff(g,v)) for v in TH] for _,g in G])
print("  Jacobian rank (exact):", J.rank(), " of", J.shape)
grad_ell = sp.Matrix([0,0,0,0,1,1,1])

print("\n[Farkas] solve grad(ell) = sum lambda_k grad(g_k),  lambda >= 0")
found=[]
for sub in itertools.combinations(range(len(G)), 7):
    A = J[list(sub),:].T
    if A.det()==0: continue
    lam = sp.simplify(A.solve(grad_ell))
    if all(sp.simplify(x)>0 for x in lam):
        found.append((sub, lam))
for sub,lam in found:
    print("  subset:", [G[k][0] for k in sub])
    print("   lambda =", [sp.radsimp(sp.nsimplify(x,[r15])) for x in lam])
    print("   lambda*45/sqrt(15) =", [sp.nsimplify(sp.radsimp(x*45/r15)) for x in lam])
    print("   residual:", list(sp.simplify(J[list(sub),:].T*lam - grad_ell)))
if not found: print("  !!! NO positive Farkas certificate on any 7-subset")

print("\n[strict blocking] is  {v : J v >= 0}  cap  {grad_ell . v <= 0}  =  {0} ?")
# with a strictly positive lambda on a rank-7 subset: Jv>=0 => grad_ell.v = sum lam_k (Jv)_k >= 0,
# equality forces (Jv)_k = 0 for all k in the subset, and the subset matrix is invertible => v=0.
sub,lam = found[0]
print("   rank of chosen 7-subset:", J[list(sub),:].rank(), "(=7 => LICQ + strict first-order min)")

# ---- CROSS-CHECK: is g >= 0 really the feasible side?  brute numeric probe ---
print("\n[feasible-side probe] tiny steps along +/- grad g_k, evaluated by the actual-panel validity predicate")
import sys, math, random
sys.path.insert(0, str(VERIFY))
from n4_tighten2 import valid_fixed2, config, crossings, ell as ELL, code_str, canon, ev_over_of, fox5_events
t0 = math.acos(-0.25); s0 = 2*math.pi-2*t0
TH0 = [t0, s0, -s0, -t0, 6/math.sqrt(15), 8/math.sqrt(15), 6/math.sqrt(15)]
Jn = [[float(sp.N(J[k,j],30)) for j in range(7)] for k in range(len(G))]
for k,(nm,_) in enumerate(G):
    ok_plus=ok_minus=None
    gk = Jn[k]; nn = math.sqrt(sum(x*x for x in gk))
    for sgn,label in ((1,'plus'),(-1,'minus')):
        th2 = [TH0[j]+sgn*1e-6*gk[j]/nn for j in range(7)]
        r = valid_fixed2(tuple(th2), ORDER) is not None
        if sgn==1: ok_plus=r
        else: ok_minus=r
    print(f"   {nm:14} step along +grad g: valid={ok_plus}   along -grad g: valid={ok_minus}"
          f"   {'OK (g>=0 is feasible side)' if ok_plus and not ok_minus else 'CHECK'}")

print("\n[blockedness probe] 4000 random directions, exact-axiom oracle, radius sweep")
rnd = random.Random(20260804)
worst = None; nfeas=0
for _ in range(4000):
    v = [rnd.gauss(0,1) for _ in range(7)]; n = math.sqrt(sum(x*x for x in v)); v=[x/n for x in v]
    for rad in (1e-3, 1e-4, 1e-5, 1e-6):
        th2 = [TH0[j]+rad*v[j] for j in range(7)]
        if valid_fixed2(tuple(th2), ORDER) is not None:
            nfeas += 1
            de = ELL(th2)-ELL(TH0)
            if worst is None or de < worst[0]: worst = (de, rad, list(v))
            break
print(f"   feasible probes: {nfeas}/4000;  most negative d(ell) found = {worst[0]:.3e} at radius {worst[1]:.0e}")
print("   (a genuinely blocked point must show d(ell) >= 0 for every feasible probe)")
