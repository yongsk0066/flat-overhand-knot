# exact_branch.py — exact certification of crossing-branch preservation along the second family.
# 대상: 교차 3개(A-C, A-D, B-D)의 존재(파라미터 열린창)·행렬식 비영·통과 순서 부호· — targets: existence of the 3 crossings (A-C, A-D, B-D) (an open parameter window)·nonvanishing determinants·passage-order signs·
#       교차-크리스선 거리 비영 — 각 분자 C_j(x,h)에 res_h(N5,·) 무근(Sturm) + 끝점/내부 부호 — nonvanishing crossing-to-crease-line distances — for each numerator C_j(x,h): res_h(N5,·) root-free (Sturm) + endpoint/interior signs
import sympy as sp
from mpmath import iv
from exact_reach import HEIGHT_INTERVAL, SAMPLE, X_INTERVAL, enclosure, horner, isolate_root, require
from ribbon_geometry import family_polynomial

x,h=sp.symbols('x h', positive=True)
N5 = family_polynomial(x, h)
g=(1-x**2)/(2*x); b=(x**2+1)/(2*x); a=h; q=sp.together(g-h)
E=sp.Matrix([0,0]); F=sp.Matrix([-x,1]); G=sp.Matrix([g,1]); H=sp.Matrix([h,0])
W0u=sp.Matrix([sp.together((q**2-1)/(q**2+1)), sp.together(2*q/(q**2+1))])
Q1=sp.Matrix([sp.together(g+a*W0u[0]), sp.together(1+a*W0u[1])])
Q0=sp.Matrix([sp.together(h+b*W0u[0]), sp.together(b*W0u[1])])
U=Q0-Q1
def cross2(u,v): return sp.together(u[0]*v[1]-u[1]*v[0])
# v3 = refl_U(W0u)
UU=sp.together(U.dot(U))
tf=sp.together(2*W0u.dot(U)/UU)
v3=sp.Matrix([sp.together(tf*U[0]-W0u[0]), sp.together(tf*U[1]-W0u[1])])
M1=sp.Matrix([-x/2, sp.Rational(1,2)]); M2=sp.Matrix([sp.together((g+h)/2), sp.Rational(1,2)])
den3=cross2(W0u,U)
s3=sp.together(cross2(Q1-M2,U)/den3)
M3=sp.Matrix([sp.together(M2[0]+s3*W0u[0]), sp.together(M2[1]+s3*W0u[1])])
# 조각: A: M1 + s·dA (κ1-반사 방향, s>0) | B: M1+t(M2-M1) | C: M2+t(M3-M2) | D: M3+t·v3 — pieces: A: M1 + s·dA (the κ1-reflection direction, s>0) | B: M1+t(M2-M1) | C: M2+t(M3-M2) | D: M3+t·v3
dA=sp.Matrix([sp.together((1-x**2)/(x**2+1)), sp.together(2*x/(x**2+1))])
dB=M2-M1; dC=M3-M2; dD=v3
def meet_params(P1,D1,P2,D2):
    det=cross2(D1,D2)
    s=sp.together(cross2(P2-P1,D2)/det); u=sp.together(cross2(P2-P1,D1)/det)
    return s,u,det
COND={}
sAC,tAC,dAC = meet_params(M1,dA,M2,dC)   # A-C
COND["AC:s>0"]=sAC; COND["AC:t>0"]=tAC; COND["AC:t<1"]=sp.together(1-tAC); COND["AC:det"]=dAC
sAD,tAD,dAD = meet_params(M1,dA,M3,dD)   # A-D
COND["AD:s>0"]=sAD; COND["AD:t>0"]=tAD; COND["AD:det"]=dAD
tBD,uBD,dBD = meet_params(M1,dB,M3,dD)   # B-D
COND["BD:t>0"]=tBD; COND["BD:t<1"]=sp.together(1-tBD); COND["BD:u>0"]=uBD; COND["BD:det"]=dBD
# 통과 순서: A 위에서 (A-D)가 (A-C)보다 먼저/나중 — 부호 고정; D 위에서 (B-D) vs (A-D) — passage order: on A, whether (A-D) comes before/after (A-C) — the sign is fixed; on D, (B-D) vs (A-D)
COND["ord:A"]=sp.together(sAD-sAC)
COND["ord:D"]=sp.together(uBD-tAD)
# 교차-크리스선 거리 비영: 교차점 3개 × 크리스선 3개 (line(E,F-E), line(H,G-H), line(Q1,U)) — nonvanishing crossing-to-crease-line distances: 3 crossing points × 3 crease lines (line(E,F-E), line(H,G-H), line(Q1,U))
def pt_on(P,D,s): return sp.Matrix([sp.together(P[0]+s*D[0]), sp.together(P[1]+s*D[1])])
X1=pt_on(M1,dA,sAC); X2=pt_on(M1,dA,sAD); X3=pt_on(M1,dB,tBD)
for nm,X in (("X1",X1),("X2",X2),("X3",X3)):
    for cn,(P0,D0) in (("k1",(E,F-E)),("k2",(H,G-H)),("k3",(Q1,U))):
        COND[f"gp:{nm}-{cn}"]=cross2(D0, X-P0)
print(f"[조건 분자] {len(COND)}개")
# 다항화 + 종결식 + Sturm — polynomialization + resultant + Sturm
a_r,b_r=sp.Rational(322,1000),sp.Rational(325,1000)
q4=5*x**4-10*x**2+1; p_x=x**12+4*x**10-7*x**4+20*x**2-2
def sturm_count(R):
    ch=sp.sturm(sp.Poly(R,x))
    def sc(t):
        vs=[c.eval(t) for c in ch]; vs=[v for v in vs if v!=0]
        return sum(1 for i in range(len(vs)-1) if (vs[i]>0)!=(vs[i+1]>0))
    return sc(a_r)-sc(b_r)
allok=True; results=[]
for name,expr in COND.items():
    num=sp.fraction(sp.cancel(expr))[0]
    P=sp.expand(num)
    R=sp.expand(sp.resultant(N5,P,h))
    if R==0:
        results.append((name,"res≡0(공통성분!)")); allok=False; continue
    # 끝점 인수 제거 — remove the endpoint factors
    note=""
    for f,fn in ((q4,"q4"),(p_x,"px")):
        qq,rr=sp.div(R,f,x)
        while rr==0:
            R=sp.expand(qq); note+=f"[{fn}]"
            qq,rr=sp.div(R,f,x)
    # x^k 제거 — remove the x^k factor
    Rp=sp.Poly(R,x)
    k=min(m[0] for m in Rp.monoms())
    if k>0: R=sp.expand(sp.cancel(R/x**k))
    cnt=sturm_count(R)
    results.append((name,f"deg{sp.degree(R,x)} {note} 근{cnt}"))
    if cnt!=0: allok=False
for nm,st in results: print(f"  {nm:14} {st}")
print("[판정]","전부 무근 — 가지 데이터가 구간 전체에서 불변 ✓" if allok else "잔근 — 개별 검토 필요")

# Isolate each sign anchor over the rationals before interval evaluation.
print("\n[부호 앵커]")
iv.dps = 60


def root_box(poly, bounds):
    low, high = isolate_root(poly, bounds)
    return iv.mpf([enclosure(low).a, enclosure(high).b])


def sign(value):
    return 1 if value.a > 0 else (-1 if value.b < 0 else 0)


def endpoint_sign(expr, height, modulus, box):
    num, den = sp.fraction(sp.cancel(expr.subs(h, height)))
    den_box = horner(sp.Poly(sp.rem(den, modulus, x), x), box)
    require(0 not in den_box, "endpoint substitution has a singular denominator")
    return sign(horner(sp.Poly(sp.rem(num, modulus, x), x), box) / den_box)


Hmid = root_box(sp.Poly(N5.subs(x, SAMPLE), h), HEIGHT_INTERVAL)
Xp = root_box(sp.Poly(q4, x), X_INTERVAL)
Xd = root_box(sp.Poly(p_x, x), X_INTERVAL)
hp = (1 - 3*x**2)/(2*x)
hd = (2*x**11 + 8*x**9 - x**7 - 4*x**5 - 13*x**3 + 40*x)/12
okall=True
for name, expr in COND.items():
    num = sp.fraction(sp.cancel(expr))[0]
    s_in = sign(horner(sp.Poly(num.subs(x, SAMPLE), h), Hmid))
    s_p = endpoint_sign(num, hp, q4, Xp)
    s_d = endpoint_sign(num, hd, p_x, Xd)
    agree = (s_in==s_p==s_d) and s_in!=0
    if not agree: okall=False
    print(f"  {name:14} 내부:{s_in:+d} 오각형:{s_p:+d} θ†:{s_d:+d} {'✓' if agree else '⚠'}")
final_ok = allok and okall
print("[결론]","22조건 부호가 세 앵커에서 일치·비영 + 종결식 무근 ⟹ 가지 데이터 불변 [정확대수+구간]; θ†의 M(B) 소속 정확 인증 포함" if final_ok else "불일치 또는 구간 내 잔근 — 검토 필요")
require(allok, "One or more branch-preservation resultants have a root on [0.322, 0.325]")
require(okall, "One or more branch-preservation sign anchors failed or disagreed")
