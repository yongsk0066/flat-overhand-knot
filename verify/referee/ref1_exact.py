"""Exact geometry of the historical figure-eight candidate."""
import sympy as sp
from fractions import Fraction as F

print("="*78)
print("Independent construction of the four-crease figure-eight state")
print("="*78)

# Folding geometry in material coordinates.
# centerline: phi_0 = 0, phi_i = tau_1+..+tau_i ; d_i = (cos phi_i, sin phi_i)
# M_0 = 0, M_i = M_{i-1} + L_i d_i (i=1,2,3)
# pass P_0 = ray from M_0 along -d_0 ; P_i = [M_{i-1},M_i] ; P_4 = ray from M_3 along d_4
# crease at M_i (i=0..3) is the fold line: reflection across it maps d_i -> d_{i+1},
#   hence its unit direction u_i is parallel to d_i + d_{i+1}.
# crease CHORD = crease line clipped to band of width 1 about the centerline:
#   half-length h_i with h_i*|u_i . n_i| = 1/2,  n_i = rot90(d_i).

def build(taus, Ls):
    phi = [sp.Integer(0)]
    for t in taus: phi.append(sp.simplify(phi[-1] + t))
    d = [sp.Matrix([sp.cos(p), sp.sin(p)]) for p in phi]
    M = [sp.Matrix([0, 0])]
    for i in range(3): M.append(sp.simplify(M[-1] + Ls[i]*d[i+1]))
    CH = []
    for i in range(4):
        v = d[i] + d[i+1]
        n = sp.Matrix([-d[i][1], d[i][0]])
        # u = v/|v| ; h = 1/(2|u.n|) ; chord halfvector = h*u = v/(2|v.n|)  (|v| cancels)
        hv = sp.simplify(v / (2*sp.Abs(v.dot(n))))
        CH.append((sp.simplify(M[i]-hv), sp.simplify(M[i]+hv)))
    return dict(phi=phi, d=d, M=M, CH=CH, L=Ls)

t = sp.acos(sp.Rational(-1,4))
s = 2*sp.pi - 2*t
r15 = sp.sqrt(15)
L1 = 6/r15; L2 = 8/r15
G = build([t, s, -s, -t], [L1, L2, L1])
sim = lambda e: sp.radsimp(sp.simplify(sp.expand_trig(sp.simplify(e))))

print("\n--- fold points from reflected directions ---")
for i, m in enumerate(G['M']):
    print(f"  M{i} = ({sim(m[0])}, {sim(m[1])})")
print("  claim: M3 = (-sqrt(15)/3, 1)  ->",
      "MATCH" if sim(G['M'][3][0]+r15/3)==0 and sim(G['M'][3][1]-1)==0 else "MISMATCH")

print("\n--- crease chord endpoints, in units x = xs*sqrt(15)/30 ---")
def scaled(p):
    xs = sim(p[0]*30/r15); ys = sim(p[1])
    return (sp.nsimplify(xs), sp.nsimplify(ys))
lbl = [("E1-","E1+"),("E2-","E2+"),("E3-","E3+"),("E4-","E4+")]
SC = {}
for i,(a,b) in enumerate(G['CH']):
    for nm,p in zip(lbl[i],(a,b)):
        SC[nm] = scaled(p); print(f"  {nm} = {SC[nm]}  (scaled)")

claim = {"E1+":(3,sp.Rational(1,2)), "E1-":(-3,sp.Rational(-1,2)),
         "E2-":(1,sp.Rational(3,2)), "E2+":(-7,sp.Rational(3,2)),
         "E3-":(-3,sp.Rational(-1,2)),"E3+":(-11,sp.Rational(-1,2)),
         "E4-":(-13,sp.Rational(1,2)),"E4+":(-7,sp.Rational(3,2))}
bad = [k for k in claim if (sp.simplify(SC[k][0]-claim[k][0])!=0 or sp.simplify(SC[k][1]-claim[k][1])!=0)]
print("  comparison with the paper's endpoints:", "ALL MATCH" if not bad else f"MISMATCH {bad}")
print("  double points: E2+=E4+ ?", SC["E2+"]==SC["E4+"], "  E1-=E3- ?", SC["E1-"]==SC["E3-"])

# ---- exact tiling / no-notch, in EXACT RATIONAL scaled coordinates ----------
print("\n--- no-notch: union of core panels vs hexagon, EXACT rational arithmetic ---")
def lat(p, P, dd):  # signed lateral coordinate wrt line (P,dd), dd unit
    return -(p[0]-P[0])*dd[1] + (p[1]-P[1])*dd[0]

# assemble each panel from the two bounding crease chords, sorted by lat sign
panels = []
for i in (1,2,3):
    P, dd = G['M'][i-1], G['d'][i]
    ends = list(G['CH'][i-1]) + list(G['CH'][i])
    pos = [e for e in ends if sim(lat(e,P,dd)) > 0]
    neg = [e for e in ends if sim(lat(e,P,dd)) < 0]
    assert len(pos)==2 and len(neg)==2, (i,[sim(lat(e,P,dd)) for e in ends])
    for e in ends:
        assert sp.Abs(sim(lat(e,P,dd)))==sp.Rational(1,2), "crease endpoint not on band edge!"
    lon = lambda e: sim((e[0]-P[0])*dd[0]+(e[1]-P[1])*dd[1])
    neg.sort(key=lambda e: float(lon(e))); pos.sort(key=lambda e: float(lon(e)), reverse=True)
    quad = [scaled(e) for e in neg+pos]
    panels.append([(F(sp.Rational(a).p, sp.Rational(a).q), F(sp.Rational(b).p, sp.Rational(b).q)) for a,b in quad])
    print(f"  panel {i} (scaled) = {panels[-1]}")

HEX = [(-13,F(1,2)),(-7,F(3,2)),(1,F(3,2)),(3,F(1,2)),(-3,F(-1,2)),(-11,F(-1,2))]
HEX = [(F(a),F(b)) for a,b in HEX]

def area(poly):
    A = F(0)
    for k in range(len(poly)):
        x1,y1 = poly[k]; x2,y2 = poly[(k+1)%len(poly)]
        A += x1*y2 - x2*y1
    return A/2
def orient(poly):
    return poly if area(poly) > 0 else poly[::-1]
def clip(subject, clipper):   # Sutherland-Hodgman, both convex, CCW
    out = list(subject)
    C = orient(list(clipper))
    for k in range(len(C)):
        A = C[k]; B = C[(k+1)%len(C)]
        inp, out = out, []
        if not inp: break
        side = lambda P: (B[0]-A[0])*(P[1]-A[1]) - (B[1]-A[1])*(P[0]-A[0])
        for j in range(len(inp)):
            P = inp[j]; Q = inp[(j+1)%len(inp)]
            sp_, sq = side(P), side(Q)
            if sp_ >= 0: out.append(P)
            if (sp_ > 0 and sq < 0) or (sp_ < 0 and sq > 0):
                tt = sp_/(sp_-sq)
                out.append((P[0]+tt*(Q[0]-P[0]), P[1]+tt*(Q[1]-P[1])))
    return out

PN = [orient(p) for p in panels]; H = orient(HEX)
print(f"  hexagon area (scaled) = {area(H)}   panel areas = {[area(p) for p in PN]}")
# (a) each panel subset of hexagon: every panel vertex is a hexagon vertex
allin = all(v in set(H) for p in PN for v in p)
print(f"  (a) every panel vertex is a hexagon vertex -> panels <= conv(hex) = hex : {allin}")
# (b) area of union by exact inclusion-exclusion
A1,A2,A3 = [area(p) for p in PN]
A12 = area(clip(PN[0],PN[1])) if clip(PN[0],PN[1]) else F(0)
A13 = area(clip(PN[0],PN[2])) if clip(PN[0],PN[2]) else F(0)
A23 = area(clip(PN[1],PN[2])) if clip(PN[1],PN[2]) else F(0)
t123 = clip(clip(PN[0],PN[1]),PN[2]); A123 = area(t123) if t123 else F(0)
U = A1+A2+A3-A12-A13-A23+A123
print(f"  pairwise areas  A12={A12} A13={A13} A23={A23}  triple A123={A123}")
print(f"  (b) area(union) = {U}   area(hex) = {area(H)}   EQUAL: {U==area(H)}")
print(f"  => union == hexagon EXACTLY (no notch, no sliver): {allin and U==area(H)}")
print(f"  average overlap multiplicity = {(A1+A2+A3)}/{area(H)} = {sp.Rational((A1+A2+A3).numerator,(A1+A2+A3).denominator)/sp.Rational(area(H).numerator,area(H).denominator)}")

# ---- hexagon metrics --------------------------------------------------------
print("\n--- hexagon metrics (from scaled coords, unscaled back) ---")
def real(pt): return sp.Matrix([sp.Rational(pt[0].numerator,pt[0].denominator)*r15/30,
                                sp.Rational(pt[1].numerator,pt[1].denominator)])
Hr = [real(p) for p in HEX]
per = sp.S(0); sides=[]
for k in range(6):
    e = sim((Hr[(k+1)%6]-Hr[k]).norm()); sides.append(e); per += e
print("  sides   =", sides)
print("  claimed  [2sqrt10/5, 4/sqrt15, 4/sqrt15, 2sqrt10/5, 4/sqrt15, 4/sqrt15] ->",
      all(sim(sides[k]-c)==0 for k,c in enumerate([2*sp.sqrt(10)/5,4/r15,4/r15,2*sp.sqrt(10)/5,4/r15,4/r15])))
print("  perimeter =", sim(per), " claimed 4sqrt10/5+16/sqrt15 ->", sim(per-(4*sp.sqrt(10)/5+16/r15))==0)
for k in range(6):
    a=Hr[(k-1)%6]-Hr[k]; b=Hr[(k+1)%6]-Hr[k]
    print(f"  interior angle at V{k}: cos = {sim(a.dot(b)/(a.norm()*b.norm()))}")
ar = sp.Rational(area(H).numerator, area(H).denominator)*r15/30
print("  area =", sim(ar), "=", sp.N(ar,12), " claimed 4sqrt15/5 ->", sim(ar-4*r15/5)==0)
ell = sim(2*L1+L2)
print("  ell* =", ell, "=", sp.N(ell,20), " claimed 4sqrt15/3 ->", sim(ell-4*r15/3)==0)
print("  area/ell* =", sim(ar/ell), " claimed 3/5 ->", sim(ar/ell-sp.Rational(3,5))==0)
