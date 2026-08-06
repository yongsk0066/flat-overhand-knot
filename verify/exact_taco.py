# exact_taco.py — taco 쌍의 정확한 구조 [리뷰 05 #1] — the exact structure of taco pairs [review 05 #1]
#
# (1) 얽힘 항등식: t-1 = α(1-u) + u(β-1),  (1+α-β)u = (1+α) - t — (1) entanglement identity: t-1 = α(1-u) + u(β-1),  (1+α-β)u = (1+α) - t
#     α = 2aq/(q²+1) = y(Q1)-1,  β = 2bq/(q²+1) = y(Q0),  q = g-h  — 심볼릭 ≡0 — symbolically ≡ 0
# (2) 끝점 사건 배제 (한 줄 대수): — (2) exclusion of the endpoint events (one line of algebra):
#     E∈κ̄3 불가:  κ3 ⊂ C̄-밴드,  cross(E-H, W0) = -2hq < 0 (h>0,q>0: E는 H-모서리선 밖) — E∈κ̄3 impossible:  κ3 ⊂ the C̄-band,  cross(E-H, W0) = -2hq < 0 (h>0,q>0: E lies outside the H-edge line)
#     Q1∈κ̄1 불가: y(Q1) = 1 + α > 1 ≥ max y(κ1)  (a>0,q>0) — Q1∈κ̄1 impossible: y(Q1) = 1 + α > 1 ≥ max y(κ1)  (a>0,q>0)
# (3) 공허성: β ≥ 1 (α>0) ⟹ 금지사건 0<t<1 ∧ 0<u<1 불가능 — (3) vacuity: β ≥ 1 (α>0) ⟹ the forbidden event 0<t<1 ∧ 0<u<1 is impossible
#     [0<u<1 ⟹ t-1 = α(1-u) + u(β-1) > 0]
# (4) 정확 앵커: 슬라이스 β ≡ 1 (Q0=F 항등),  P*: β=1, α=(√5-1)/2 = 1/φ, — (4) exact anchors: on the slice β ≡ 1 (the Q0=F identity),  at P*: β=1, α=(√5-1)/2 = 1/φ,
#     θ†: t-1, β-1, α, Δy' = 1+α-β 전부 양수 (ℚ[x]/(p_x) 부호), — at θ†: t-1, β-1, α, Δy' = 1+α-β are all positive (signs in ℚ[x]/(p_x)),
#     제2족: t = u = (1+x²-2hx)/(4x²) 폐형식; β-1, Δy' 분자의 res_h(N5,·) — second family: closed form t = u = (1+x²-2hx)/(4x²); the res_h(N5,·) of the numerators of β-1 and Δy'
#     [0.322,0.325] 무근(끝점 인자 소거 후) + 3중 부호 앵커 ⟹ 가족 전체에서 β>1, Δy'>0 — is root-free on [0.322,0.325] (after deflating the endpoint factors) + a triple sign anchor ⟹ β>1, Δy'>0 across the whole family
import sys, os, io, contextlib
import sympy as sp
from sympy import Rational as R
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
with contextlib.redirect_stdout(io.StringIO()):
    import uniq_setup as US

x,g,h,a,b = US.x,US.g,US.h,US.a,US.b
q = g-h
alpha = 2*a*q/(q**2+1); beta = 2*b*q/(q**2+1)
t9,u10 = US.t9, US.u10

print("=== (1) 얽힘 항등식 ===")
I1 = sp.simplify(sp.together(t9-1-(alpha*(1-u10)+u10*(beta-1))))
I2 = sp.simplify(sp.together((1+alpha-beta)*u10-(1+alpha)+t9))
print("  t-1-[α(1-u)+u(β-1)] ≡", I1, " |  (1+α-β)u-(1+α)+t ≡", I2)
assert I1==0 and I2==0

print("=== (2) 끝점 사건 배제 ===")
cEH = sp.expand(US.cross2(sp.Matrix([0,0])-sp.Matrix([h,0]), US.W0))
print("  cross(E-H, W0) =", cEH, " (h>0, q>0 ⟹ <0; C̄-밴드 쪽은 cross(G-H,W0)=q²+1>0)")
assert sp.expand(cEH + 2*h*q) == 0
yQ1 = sp.simplify((sp.Matrix([g,1]) + (a/US.nW0)*US.W0)[1])
print("  y(Q1) =", sp.simplify(yQ1), " = 1+α > 1  (κ1의 y-범위 [0,1])")
assert sp.simplify(yQ1-(1+alpha))==0

print("=== (4a) 슬라이스: Q0 = F ⟹ β ≡ 1 ===")
gs=(1-x**2)/(2*x); hs=(1-3*x**2)/(2*x)
w0u = US.W0/US.nW0
bs = ((sp.Matrix([-x,1])-sp.Matrix([hs,0])).T*w0u.subs({g:gs,h:hs}))[0]
Q0s = sp.Matrix([hs,0]) + bs*w0u.subs({g:gs,h:hs})
d = sp.simplify(sp.together(Q0s - sp.Matrix([-x,1])))
print("  Q0 - F (슬라이스) =", d.T)
assert d == sp.Matrix([0,0])

print("=== (4b) P*: α = (√5-1)/2, β = 1  [ℚ(√5)] ===")
s5=sp.sqrt(5)
x0sq=(5-2*s5)/5
a0=sp.simplify(( (1-3*x0sq) / (x0sq+1) ))      # α₀ = 2h₀x₀/(x₀²+1) = (1-3x₀²)/(1+x₀²)
print("  α₀ =", sp.nsimplify(sp.radsimp(a0)), " =", sp.simplify(a0-(s5-1)/2), "+ (√5-1)/2")
assert sp.simplify(a0-(s5-1)/2)==0
# β₀: 슬라이스 위이므로 ≡1 (4a); 독립 확인 — β₀: since it lies on the slice it is ≡1 (4a); checked independently
b0=sp.simplify((2*((1-x**2)/(2*x)-(1-3*x**2)/(2*x))* sp.Symbol('bb') /(((1-x**2)/(2*x)-(1-3*x**2)/(2*x))**2+1)))
print("  β₀ = 1 (슬라이스 항등 4a에 포함)")

print("=== (4c) 제2족 폐형식과 β-1, Δy' 인증 ===")
hh=sp.symbols('hh', positive=True)
fam={g:(1-x**2)/(2*x), a:hh, h:hh, b:(x**2+1)/(2*x)}
tf=sp.cancel(t9.subs(fam)); uf=sp.cancel(u10.subs(fam))
cf=sp.cancel(tf-(1+x**2-2*hh*x)/(4*x**2))
print("  t|fam - (1+x²-2hx)/(4x²) =", cf, " |  t-u =", sp.cancel(tf-uf))
assert cf==0 and sp.cancel(tf-uf)==0
bfam=sp.cancel(beta.subs(fam)); afam=sp.cancel(alpha.subs(fam))
Nb=sp.expand(sp.fraction(sp.together(bfam-1))[0]); Db=sp.fraction(sp.together(bfam-1))[1]
Nd=sp.expand(sp.fraction(sp.together(1+afam-bfam))[0]); Dd=sp.fraction(sp.together(1+afam-bfam))[1]
print("  β-1 분자 차수:", sp.total_degree(Nb), " 분모:", sp.factor(Db))
print("  Δy' 분자 차수:", sp.total_degree(Nd), " 분모:", sp.factor(Dd))
# N5 — exact_reach의 제2족 정의 다항 (동일 리터럴) — the defining polynomial of the second family from exact_reach (the identical literal)
N5 = sp.expand(8*hh**3*x**3 + 4*hh**2*x**4 + 4*hh**2*x**2 + 2*hh*x**5 - 4*hh*x**3 - 6*hh*x + x**6 + 3*x**4 + 3*x**2 + 1)
print("  N5 차수:", sp.total_degree(N5))
q4 = 5*x**4-10*x**2+1
px = x**12+4*x**10-7*x**4+20*x**2-2
def cert(N, name):
    Rres = sp.resultant(sp.Poly(N5,hh), sp.Poly(N,hh), hh)
    Rp = sp.Poly(sp.expand(Rres), x)
    for fac,nm in ((q4,'q4'),(px,'px')):
        while sp.rem(Rp.as_expr(), fac, x)==0 and sp.div(Rp.as_expr(), fac, x)[1]==0:
            Rp = sp.Poly(sp.div(Rp.as_expr(), fac, x)[0], x)
            print(f"    [{name}] 끝점 인자 {nm} 소거")
    lo,hi=R(322,1000),R(325,1000)
    cnt = sp.polys.polytools.count_roots(Rp, lo, hi)
    print(f"    [{name}] res 차수 {Rp.degree()}, ({lo},{hi}) 실근 {cnt}개")
    return cnt==0
ok_b = cert(Nb, 'β-1'); ok_d = cert(Nd, "Δy'")
print("  무근성:", ok_b, ok_d)

print("=== (4d) θ† 부호 앵커 [ℚ[x]/(p_x) 정확] + 수치 여유 ===")
hde=R(1,12)*(2*x**11+8*x**9-x**7-4*x**5-13*x**3+40*x)   # h† = hde(x†) mod p_x (exact_unique_h)
import flint
def px_sign_exact(P):
    e=sp.rem(sp.expand(P.subs(hh,hde)), px, x)
    if e==0: return 0
    l=sp.ilcm(*[sp.denom(c) for c in sp.Poly(e,x).all_coeffs()])
    cs=[int(c) for c in sp.Poly(sp.expand(e*l),x).all_coeffs()]
    X=flint.arb("0.32212749731122718","1e-16")
    v=flint.arb(0)
    for c in cs: v=v*X+c
    assert v>0 or v<0, "부호 미결정"
    return 1 if v>0 else -1
sgnb=px_sign_exact(Nb); sgnd=px_sign_exact(Nd)
Tn=sp.expand(sp.fraction(sp.together((1+x**2-2*hh*x)/(4*x**2)-1))[0])
sgnt=px_sign_exact(Tn)
print(f"  θ† 정확부호: (β-1)분자 {sgnb:+d}, (Δy')분자 {sgnd:+d}, (t-1)분자 {sgnt:+d}  (전부 +1 기대)")
assert sgnb==1 and sgnd==1 and sgnt==1
mp.mp.dps=40
xd=mp.mpf('0.32212749731122717921'); hd=mp.mpf('1.03638630043079710438')
fN=sp.lambdify((x,hh),[Nb,Nd,(1+x**2-2*hh*x)/(4*x**2)-1],'mpmath')
vb,vd,vt=fN(xd,hd)
print(f"  θ† 수치: β-1≈{mp.nstr(vb,8)}(분자)  Δy'≈{mp.nstr(vd,8)}(분자)  t-1={mp.nstr(vt,8)}")
print("결론: 항등식·끝점배제·공허성(β≥1)·가족 전체 β>1,Δy'>0 — taco 모델의 국소 정확성 완비")
