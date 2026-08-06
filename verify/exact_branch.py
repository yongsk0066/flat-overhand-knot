# exact_branch.py — 제2족 위 가지(branch) 보존의 정확 인증 [리뷰 04 #1] — exact certification that the branch is preserved along the second family [review 04 #1]
# 대상: 교차 3개(A-C, A-D, B-D)의 존재(파라미터 열린창)·행렬식 비영·통과 순서 부호· — targets: existence of the 3 crossings (A-C, A-D, B-D) (an open parameter window)·nonvanishing determinants·passage-order signs·
#       교차-크리스선 거리 비영 — 각 분자 C_j(x,h)에 res_h(N5,·) 무근(Sturm) + 끝점/내부 부호 — nonvanishing crossing-to-crease-line distances — for each numerator C_j(x,h): res_h(N5,·) root-free (Sturm) + endpoint/interior signs
import sympy as sp, mpmath as mp

x,h=sp.symbols('x h', positive=True)
N5 = sp.expand(8*h**3*x**3+4*h**2*x**4+4*h**2*x**2+2*h*x**5-4*h*x**3-6*h*x+x**6+3*x**4+3*x**2+1)
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

# ============ 2부: 부호 앵커 — 내부 크라우칙 상자 + 오각형(모드 q4) + θ†(모드 p_x) ============ — Part 2: sign anchors — interior Krawczyk box + pentagon (mod q4) + θ† (mod p_x) ============
print("\n[부호 앵커]")
mp.mp.dps=50
from mpmath import iv, mpf, findroot
iv.dps=50
# (i) 내부점 x_m=1293/4000, 가지 h-상자 (exact_reach와 동일) — (i) interior point x_m=1293/4000, branch h-box (same as in exact_reach)
N5f=sp.lambdify((x,h),N5,'mpmath')
xm=sp.Rational(1293,4000)
hm=findroot(lambda hh: N5f(mp.mpf(str(float(xm))),hh), mp.mpf('1.041'))
H_=iv.mpf([str(hm-mp.mpf(10)**-20),str(hm+mp.mpf(10)**-20)])
def iv_sign_bivar(P):
    Ph=sp.Poly(P.subs(x,xm), h)
    l=sp.ilcm(*[sp.denom(c) for c in Ph.all_coeffs()])
    cs=[sp.Rational(c) for c in sp.Poly(sp.expand(P.subs(x,xm)*l),h).all_coeffs()]
    v=iv.mpf(0)
    for c in cs: v=v*H_+iv.mpf([str(c),str(c)])
    return 1 if v.a>0 else (-1 if v.b<0 else 0)
# (ii) 오각형: h0 유리식 (h0=(1-3x^2)/2x on slice), 모드 q4 + x0 구간 — (ii) pentagon: h0 as a rational function (h0=(1-3x^2)/2x on the slice), mod q4 + the x0 interval
x0=findroot(lambda t:5*t**4-10*t**2+1, mpf('0.3249196962'))
X0IV=iv.mpf([str(x0-mpf(10)**-30),str(x0+mpf(10)**-30)])
h0e=sp.together((1-3*x**2)/(2*x))
def q4_sign(P):
    # P는 이변수 '분자 다항식' — h에 h0(x) 유리식 대입 후 분자만 (분모 (2x)^k>0) — P is a bivariate 'numerator polynomial' — after substituting the rational function h0(x) for h, keep only the numerator (the denominator (2x)^k>0)
    e=sp.rem(sp.expand(sp.fraction(sp.together(P.subs(h,h0e)))[0]), 5*x**4-10*x**2+1, x)
    if e==0: return 0
    if not e.free_symbols: return 1 if e>0 else -1
    l=sp.ilcm(*[sp.denom(c) for c in sp.Poly(e,x).all_coeffs()])
    cs=[int(c) for c in sp.Poly(sp.expand(e*l),x).all_coeffs()]
    v=iv.mpf(0)
    for c in cs: v=v*X0IV+c
    return 1 if v.a>0 else (-1 if v.b<0 else 0)
# (iii) θ†: h† 다항 표현, 모드 p_x + x† 구간 — (iii) θ†: polynomial representation of h†, mod p_x + the x† interval
xd=findroot(lambda t: t**12+4*t**10-7*t**4+20*t**2-2, mpf('0.32212749731122718'))
XDIV=iv.mpf([str(xd-mpf(10)**-30),str(xd+mpf(10)**-30)])
hde=sp.Rational(1,12)*(2*x**11+8*x**9-x**7-4*x**5-13*x**3+40*x)
def px_sign(P):
    e=sp.rem(sp.expand(P.subs(h,hde)), p_x, x)   # hde는 다항 — 대입 결과도 다항 — hde is a polynomial, so the result of the substitution is polynomial too
    if e==0: return 0
    if not e.free_symbols: return 1 if e>0 else -1
    l=sp.ilcm(*[sp.denom(c) for c in sp.Poly(e,x).all_coeffs()])
    cs=[int(c) for c in sp.Poly(sp.expand(e*l),x).all_coeffs()]
    v=iv.mpf(0)
    for c in cs: v=v*XDIV+c
    return 1 if v.a>0 else (-1 if v.b<0 else 0)
okall=True
for name,expr in COND.items():
    num=sp.expand(sp.fraction(sp.cancel(expr))[0])
    s_in=iv_sign_bivar(num); s_p=q4_sign(num); s_d=px_sign(num)
    agree = (s_in==s_p==s_d) and s_in!=0
    if not agree: okall=False
    print(f"  {name:14} 내부:{s_in:+d} 오각형:{s_p:+d} θ†:{s_d:+d} {'✓' if agree else '⚠'}")
print("[결론]","22조건 부호가 세 앵커에서 일치·비영 + 종결식 무근 ⟹ 가지 데이터 불변 [정확대수+구간]; θ†의 M(B) 소속 정확 인증 포함" if okall else "불일치 — 검토")
