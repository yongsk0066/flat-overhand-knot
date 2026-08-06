# ref12_trunc2.py — REFEREE: L_trunc convention audit, done robustly (convex half-plane clip).
import sympy as sp
r15=sp.sqrt(15)
HEXs=[(-13,sp.Rational(1,2)),(-7,sp.Rational(3,2)),(1,sp.Rational(3,2)),
      (3,sp.Rational(1,2)),(-3,sp.Rational(-1,2)),(-11,sp.Rational(-1,2))]
HEX=[sp.Matrix([sp.Integer(x)*r15/30, sp.Rational(y)]) for x,y in HEXs]
# orient CCW
def area(P):
    return sum(P[k][0]*P[(k+1)%len(P)][1]-P[(k+1)%len(P)][0]*P[k][1] for k in range(len(P)))/2
if sp.simplify(area(HEX))<0: HEX=HEX[::-1]
print("hexagon (CCW) area =", sp.simplify(area(HEX)), "=", sp.N(area(HEX),10))

def ray_exit(P, v):
    """largest u>=0 with P+u v inside the convex hexagon (half-plane clip). Exact."""
    lo, hi = sp.S(0), sp.oo
    for k in range(6):
        A=HEX[k]; B=HEX[(k+1)%6]; e=B-A
        # inward normal for CCW polygon: n = (-e_y, e_x); inside <=> n.(X-A) >= 0
        n=sp.Matrix([-e[1], e[0]])
        a=sp.simplify(n.dot(P-A)); b=sp.simplify(n.dot(v))
        if b==0:
            if a<0: return None
            continue
        t=sp.simplify(-a/b)
        if b>0: lo=sp.Max(lo,t) if sp.simplify(t)>lo else lo
        else:   hi=sp.Min(hi,t)
    return sp.simplify(hi)

M0=sp.Matrix([sp.S(0),sp.S(0)]); M3=sp.Matrix([-r15/3,sp.S(1)])
d0=sp.Matrix([sp.S(1),sp.S(0)])
s_in = ray_exit(M0, -d0)      # in-tail: from the first fold point, outward along -d0
t_out= ray_exit(M3,  d0)      # out-tail: from the last fold point, outward along d4=d0
core=4*r15/3
print(f"\nM0 on the outline boundary: {any(sp.simplify((M0-HEX[k]).dot(sp.Matrix([-(HEX[(k+1)%6]-HEX[k])[1],(HEX[(k+1)%6]-HEX[k])[0]])))==0 for k in range(6))}")
print(f"s_in  (in-tail leaves the hexagon at) = {sp.radsimp(s_in)} = {sp.N(s_in,12)}")
print(f"t_out (out-tail leaves the hexagon at)= {sp.radsimp(t_out)} = {sp.N(t_out,12)}")
print(f"s_in == t_out == L1 == 2*sqrt(15)/5 ? {sp.simplify(s_in-2*r15/5)==0 and sp.simplify(t_out-2*r15/5)==0}")
Lo=sp.radsimp(sp.simplify(s_in+core+t_out))
Lc=sp.radsimp(28*r15/15)
print(f"\nL_trunc, paper's OUTLINE convention (doc 11 SS2 / Def 7.1 / Rem rem:trunclines)")
print(f"   = {sp.radsimp(s_in)} + {sp.radsimp(core)} + {sp.radsimp(t_out)} = {Lo} = {sp.N(Lo,15)}")
print(f"L_trunc as coded in n4_exact_kkt.py (outermost centreline CROSSING)")
print(f"   = 4sqrt15/15 + {sp.radsimp(core)} + 4sqrt15/15 = {Lc} = {sp.N(Lc,15)}")
print(f"\nDISAGREE: {sp.simplify(Lo-Lc)!=0};  difference = {sp.radsimp(sp.simplify(Lo-Lc))} = {sp.N(Lo-Lc,10)}")
print(f"ratio L_trunc/core under the outline rule = {sp.simplify(Lo/core)}")
# where exactly does the in-tail leave, and through which side?
X=M0-s_in*d0
print(f"\nin-tail exit point = ({sp.radsimp(X[0])}, {X[1]}) ; scaled x = {sp.simplify(X[0]*30/r15)}")
for k in range(6):
    A=HEX[k]; B=HEX[(k+1)%6]
    cr=sp.simplify((B[0]-A[0])*(X[1]-A[1])-(B[1]-A[1])*(X[0]-A[0]))
    if cr==0:
        print(f"   lies on side {k}: ({sp.simplify(A[0]*30/r15)},{A[1]}) -> ({sp.simplify(B[0]*30/r15)},{B[1]})")
print("   sides of the hexagon in scaled x: V0(-13,1/2) V1(-7,3/2) V2(1,3/2) V3(3,1/2) V4(-3,-1/2) V5(-11,-1/2)")
print("   crease chords = V3V4, V1V2, V4V5, V0V1 ; band edges = V2V3, V5V0")
print("   -> the in-tail exits through the BAND EDGE V5V0 (= E3+ -> E4-), at 2*sqrt(15)/5;")
print("      the outermost centreline crossing (P0 x P3) sits earlier, at 4*sqrt(15)/15.")
