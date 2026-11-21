# n4_krawczyk.py — θ*_{4_1}의 존재·유일성 구간 증서 (세잎 정리 10.1의 유사물) — interval certificate for existence and uniqueness of θ*_{4_1} (analogue of Theorem 10.1 for the trefoil)
# 활성 7벽(랭크7 부분집합: c1P2E- 제외) 다항계 W(u,L)=0 에 크라우칙 연산자 — Krawczyk operator applied to the polynomial system W(u,L)=0 of the 7 active walls (the rank-7 subset: c1P2E- excluded)
import sympy as sp, mpmath as mp
import sys
sys.path.insert(0,'verify')

# n4_chart 재구성 (다항 벽) — rebuild n4_chart (polynomial walls)
exec(open('verify/n4_chart.py').read().split('# 수치 대조용')[0])
S7=["a:c1P4E+","a:c2P4E+","a:c2P4E-","a:c4P0E-","a:c4P2E+","a:c3P0E+","a:c3P0E-"]
def swap_e(i,e):
    if i in (2,3): return '-' if e=='+' else '+'
    return e
POLYS=[]
for name in S7:
    i,e0,k=WALL_DEFS[name]
    e=swap_e(i,e0)
    base,dd=BASE[k]
    g=sp.together(lat(E[(i,e)],base,dd))
    num,den=sp.fraction(g)
    # 부호 s: θ* 근방 수치로 결정 (정배향 셀) — sign s: fixed numerically near θ* (positively oriented cell)
    W=sp.expand(2*num-den)  # s=+1 가정 후 θ* 대입으로 s 판정 — assume s=+1, then decide s by substituting θ*
    z=sp.simplify(W.subs(SUB))
    if z!=0:
        W=sp.expand(-2*num-den)
        z2=sp.simplify(W.subs(SUB))
        assert z2==0, f"{name}: 양쪽 부호 모두 실패 {z} {z2}"
    # 정수 계수화 — clear denominators to integer coefficients
    dens=[sp.denom(c) for c in sp.Poly(W,[u1,u2,u3,u4,l1,l2,l3]).coeffs()]
    lcm=sp.ilcm(*[int(d) for d in dens]) if dens else 1
    POLYS.append(sp.expand(W*lcm))
print(f"[다항계] {len(POLYS)}개, 차수 {[sp.total_degree(p) for p in POLYS]}")
VARS=[u1,u2,u3,u4,l1,l2,l3]
J=sp.Matrix([[sp.diff(p,v) for v in VARS] for p in POLYS])
Pf=[sp.lambdify(VARS,p,'mpmath') for p in POLYS]
Jf=[[sp.lambdify(VARS,J[r,c],'mpmath') for c in range(7)] for r in range(7)]

mp.mp.dps=60
r15=mp.sqrt(15)
ctr=[r15/3, r15, -r15, -r15/3, 2*r15/5, 8*r15/15, 2*r15/5]
print("[중심 잔차]", [mp.nstr(abs(f(*ctr)),3) for f in Pf])
RAD=mp.mpf(10)**-25
iv=mp.iv
iv.dps=60
X=[iv.mpf([str(c-RAD),str(c+RAD)]) for c in ctr]
root15 = iv.sqrt(15)
exact_point = [root15/3, root15, -root15, -root15/3,
               2*root15/5, 8*root15/15, 2*root15/5]
assert all(box.a < point.a and point.b < box.b
           for box, point in zip(X, exact_point))
# All quantities in the Krawczyk expression use one consistent, exactly
# represented mpf centre and preconditioner. In particular, F(m) is evaluated
# in outward-rounded interval arithmetic, not converted from a rounded mpf.
assert all(sp.simplify(p.subs(SUB)) == 0 for p in POLYS)
Yc = mp.matrix([[Jf[r][c](*ctr) for c in range(7)] for r in range(7)])
Y = Yc ** -1
m = [iv.mpf(v) for v in ctr]
Yi = [[iv.mpf(Y[r, c]) for c in range(7)] for r in range(7)]
Fm = [f(*m) for f in Pf]
JX = [[Jf[r][c](*X) for c in range(7)] for r in range(7)]
IY = [[(1 if r == c else 0) - sum(Yi[r][q] * JX[q][c] for q in range(7))
       for c in range(7)] for r in range(7)]
K = [m[r] - sum(Yi[r][q] * Fm[q] for q in range(7))
     + sum(IY[r][c] * (X[c] - m[c]) for c in range(7)) for r in range(7)]
inside = all(K[r].a > X[r].a and K[r].b < X[r].b for r in range(7))
# A strict contraction also establishes nonsingularity of the square Y.
row_bounds = [sum(abs(entry).b for entry in row) for row in IY]
contractive = all(bound.b < 1 for bound in row_bounds)
print(f"[Krawczyk] outward-rounded K(X) in int(X): {inside}")
print(f"[Krawczyk] interval infinity-norm bound below one: {contractive}")
assert inside and contractive
print("  radius 1e-25; widths:", [mp.nstr(mp.mpf(k.b)-mp.mpf(k.a), 4) for k in K])
print("PASS: existence and uniqueness in the stated box for the seven-wall polynomial system")
