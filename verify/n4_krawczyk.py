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
# 크라우칙: K = y - Y f(y) + (I - Y J(X))(X - y) — Krawczyk operator
Yc=mp.matrix([[mp.mpf(str(Jf[r][c](*ctr))) for c in range(7)] for r in range(7)])
Y=Yc**-1
fy=mp.matrix([mp.mpf(str(Pf[r](*ctr))) for r in range(7)])
Yfy=Y*fy
JX=[[iv.mpf(str(Jf[r][c](*ctr)))*0 for c in range(7)] for r in range(7)]
# 구간 야코비안: iv 평가 — interval Jacobian: evaluate with iv
def iv_eval(fn):
    return fn(*X)
JXiv=[[Jf[r][c](*X) for c in range(7)] for r in range(7)]
# I - Y*J(X) (구간) — (interval)
IY=[[ (1 if r==c else 0) - sum(iv.mpf(mp.nstr(Y[r,q],40))*JXiv[q][c] for q in range(7)) for c in range(7)] for r in range(7)]
Xy=[X[c]-iv.mpf(mp.nstr(ctr[c],50)) for c in range(7)]
K=[ iv.mpf(mp.nstr(ctr[r],50)) - iv.mpf(mp.nstr(Yfy[r],40)) + sum(IY[r][c]*Xy[c] for c in range(7)) for r in range(7)]
inside=all(K[r].a>X[r].a and K[r].b<X[r].b for r in range(7))
wid=[mp.nstr(mp.mpf(K[r].b)-mp.mpf(K[r].a),3) for r in range(7)]
print(f"[크라우칙] K(X) ⊂ int X: {inside}")
print("  상자 반경 1e-25, K 폭:", wid)
if inside:
    print("  ⟹ 상자 내 근 존재+유일 [구간 엄밀] — θ*(u-차트) 증서 완료")
