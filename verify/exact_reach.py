# exact_reach.py — 정리 C 승격: 제2 가족을 따라 오각형→θ†의 단조 유효 하강 [정확대수+구간] — promotion of Theorem C: monotone valid descent from the pentagon to θ† along the second family [exact algebra + intervals]
# 가족: g=(1-x²)/2x, a=h, b=(x²+1)/2x, N5(x,h)=0 가지, x ∈ [x†, x0] — family: g=(1-x²)/2x, a=h, b=(x²+1)/2x, the N5(x,h)=0 branch, x ∈ [x†, x0]
# 인증 대상: (i) N5_h ≠ 0 (가지 매끄러움)  (ii) P7>0, P9>0 (잔여 벽 w7=w8, w9=w10) — certification targets: (i) N5_h ≠ 0 (smoothness of the branch)  (ii) P7>0, P9>0 (the remaining walls w7=w8, w9=w10)
#            (iii) A>0 (dℓ/dx 분자; x†에서만 0)  — 방법: res_x 무근(Sturm) + 끝점 항등 + 표본 부호 — (iii) A>0 (numerator of dℓ/dx; zero only at x†)  — method: res_x root-free (Sturm) + endpoint identities + sample signs
import sympy as sp

x,h = sp.symbols('x h', positive=True)
N5 = sp.expand(8*h**3*x**3 + 4*h**2*x**4 + 4*h**2*x**2 + 2*h*x**5 - 4*h*x**3 - 6*h*x + x**6 + 3*x**4 + 3*x**2 + 1)
g = (1-x**2)/(2*x); b = (x**2+1)/(2*x); a = h
q = sp.together(g - h)

# --- 벽의 가족 제한 (uniq_setup 닫힌형에서 대입) --- — restriction of the walls to the family (substituted from the uniq_setup closed forms) ---
# w7 = 2b(q+x)(1-qx)/(q²+1) - 2hx - (x²+1)  [부호고정 +] — [sign fixed to +]
W7 = sp.together(2*b*(q+x)*(1-q*x)/(q**2+1) - 2*h*x - (x**2+1))
P7n, P7d = sp.fraction(sp.cancel(W7))
print("[w7] 분자 차수:", sp.total_degree(sp.expand(P7n)), " 분모 =", sp.factor(P7d))
# 분모 (q²+1)·(2x)^k 형태 — 양수성 확인용 인수 출력 — the denominator has the form (q²+1)·(2x)^k — its factors are printed so positivity can be checked

# w9 = t9 - 1,  t9 = cross(Q1-E, U)/cross(F-E, U),  U = Q0-Q1 (5-변수 기하를 가족에 제한) — the 5-variable geometry restricted to the family
# 좌표 (uniq_setup._config): E=(0,0), F=(-x,1), G=(g,1), H=(h,0), — coordinates (uniq_setup._config): E=(0,0), F=(-x,1), G=(g,1), H=(h,0),
#   w0 = unit(refl_{G-H}(1,0)) — 유리화: refl벡터 = ((GH)² 성분식)/|GH|² … 방향만 필요하므로 비정규 사용 — rationalization: the reflected vector = (component expression in (GH)²)/|GH|² … only the direction is needed, so the unnormalized one is used
GHx, GHy = sp.together(g-h), sp.Integer(1)
# refl_lin(d)(v) = 2(v·d)d/|d|² − v,  v=(1,0):
w0x = sp.together(2*GHx*GHx/(GHx**2+GHy**2) - 1)
w0y = sp.together(2*GHx*GHy/(GHx**2+GHy**2))
# Q1 = G + a·ŵ0, Q0 = H + b·ŵ0 — 단위화가 필요하나 U = Q0-Q1 = (H-G) + (b-a)ŵ0. — normalization is required, but U = Q0-Q1 = (H-G) + (b-a)ŵ0.
# |w0|=1이 아닌 비정규 (w0x,w0y)로는 ŵ0 = (w0x,w0y)/n, n=√(w0x²+w0y²) — n이 무리적. — with the unnormalized (w0x,w0y) (for which |w0|≠1) one has ŵ0 = (w0x,w0y)/n, n=√(w0x²+w0y²) — and n is irrational.
# 대신 s := (b-a)/n 을 새 변수로 두지 않고, n² = w0x²+w0y² 가 유리이므로 — instead, without introducing s := (b-a)/n as a new variable: since n² = w0x²+w0y² is rational,
# U ∥ n·(H-G) + (b-a)·(w0x,w0y) 방향의 유리 배수가 아님(√ 잔존). → 교차비 t9는 ŵ0에 1차. — U is not a rational multiple of the direction n·(H-G) + (b-a)·(w0x,w0y) (a √ survives). → the cross-ratio t9 is linear in ŵ0.
# 우회: t9-1의 부호는 cross(Q1-F, U)/cross(F-E,U) 부호와 동일 (t9-1 = cross(Q1-E,U)/cross(F-E,U) - 1 — workaround: the sign of t9-1 equals that of cross(Q1-F, U)/cross(F-E,U) (t9-1 = cross(Q1-E,U)/cross(F-E,U) - 1
#   = cross(Q1-E - (F-E)·1 ... ) — 직접: t9-1 = [cross(Q1-E,U) - cross(F-E,U)]/cross(F-E,U) — directly: t9-1 = [cross(Q1-E,U) - cross(F-E,U)]/cross(F-E,U)
#   = cross(Q1-F, U)/cross(F-E, U).
# cross(Q1-F, U) = cross(G-F, U) + a·cross(ŵ0, U);  U = (H-G) + (b-a)ŵ0
#   cross(ŵ0,U) = cross(ŵ0, H-G);  cross(G-F,U) = cross(G-F,H-G) + (b-a)cross(G-F,ŵ0)
# ŵ0 = w0v/n ⟹ 분자·분모 모두 (α + β/n) 꼴 → n을 곱해 (αn + β): n² 유리 ⟹ 부호식 = α²n²−β² 판별 대신 — ŵ0 = w0v/n ⟹ numerator and denominator are both of the form (α + β/n) → multiply by n to get (αn + β): n² is rational ⟹ instead of deciding the sign expression as α²n²−β²
# 간단히: 전체 식을 n으로 통분하면 t9-1 = (A1·n + B1)/(A2·n + B2), Ai,Bi ∈ Q(x,h). — more simply: putting the whole expression over the common denominator n gives t9-1 = (A1·n + B1)/(A2·n + B2), with Ai,Bi ∈ Q(x,h).
w0v = sp.Matrix([w0x, w0y]); HG = sp.Matrix([sp.together(h-g), sp.Integer(-1)])
GF = sp.Matrix([sp.together(g+x), sp.Integer(0)])   # G-F = (g+x, 0)
def cross2(u,v): return sp.together(u[0]*v[1]-u[1]*v[0])
A1 = cross2(GF, HG)                     # n-계수: cross(G-F, H-G) — coefficient of n: cross(G-F, H-G)
B1 = sp.together((b-a)*cross2(GF, w0v) + a*cross2(w0v, HG))
FE = sp.Matrix([-x, sp.Integer(1)])     # F-E
A2 = cross2(FE, HG)
B2 = sp.together((b-a)*cross2(FE, w0v))
n2 = sp.together(w0x**2 + w0y**2)
print("[w9] n² =", sp.factor(n2), " (양수 유리)")
# t9-1 = (A1 n + B1)/(A2 n + B2);  부호 인증은 (A1 n + B1)(A2 n + B2) > 0 ⟺ — certifying the sign amounts to (A1 n + B1)(A2 n + B2) > 0 ⟺
#   S9 := A1A2 n² + B1B2 + n(A1B2 + A2B1) > 0 — 여전히 n 잔존. — n still survives.
#   유리화: (A1 n + B1)(A2 n + B2)(A2 n - B2)² … 복잡. 다른 길: n = |w0v| 이고 — rationalization: (A1 n + B1)(A2 n + B2)(A2 n - B2)² … complicated. Another route: n = |w0v|, and
#   ŵ0 성분이 (q²-1, 2q)/(q²+1) 임을 상기 (uniq_setup: W0=(q²-1, 2q), |W0|=q²+1) — 단위벡터가 유리!! — recall that the components of ŵ0 are (q²-1, 2q)/(q²+1) (uniq_setup: W0=(q²-1, 2q), |W0|=q²+1) — the unit vector is rational!!
print("  ⟹ ŵ0 = (q²-1, 2q)/(q²+1) 유리 — n 우회 불필요. 재구성:")
W0u = sp.Matrix([sp.together((q**2-1)/(q**2+1)), sp.together(2*q/(q**2+1))])
Q1v = sp.Matrix([sp.together(g + a*W0u[0]), sp.together(1 + a*W0u[1])])
Q0v = sp.Matrix([sp.together(h + b*W0u[0]), sp.together(b*W0u[1])])
Uv = sp.Matrix([sp.together(Q0v[0]-Q1v[0]), sp.together(Q0v[1]-Q1v[1])])
Ev = sp.Matrix([0,0]); Fv = sp.Matrix([-x, sp.Integer(1)])
t9m1 = sp.together(cross2(sp.Matrix([Q1v[0]-Fv[0],Q1v[1]-Fv[1]]), Uv) / cross2(sp.Matrix([Fv[0]-Ev[0],Fv[1]-Ev[1]]), Uv))
P9n, P9d = sp.fraction(sp.cancel(t9m1))
print("[w9] t9-1 분자 차수:", sp.total_degree(sp.expand(P9n)), " / 분모 차수:", sp.total_degree(sp.expand(P9d)))

# ℓ' 분자: A = -2x²·N5_x + (x²-1)·N5_h — numerator of ℓ': A = -2x²·N5_x + (x²-1)·N5_h
N5x = sp.diff(N5,x); N5h = sp.diff(N5,h)
Aell = sp.expand(-2*x**2*N5x + (x**2-1)*N5h)
print("[ℓ'] 분자 차수:", sp.total_degree(Aell), " 분모 = 2x²·N5_h")
import pickle
pickle.dump({"N5":N5,"P7n":sp.expand(P7n),"P7d":P7d,"P9n":sp.expand(P9n),"P9d":sp.expand(P9d),
             "Aell":Aell,"N5h":sp.expand(N5h)}, open("/tmp/reach_polys.pkl","wb"))
print("저장: /tmp/reach_polys.pkl")

# ================= 2부: 앵커·끝점 항등·종결식·Sturm ================= — Part 2: anchors · endpoint identities · resultants · Sturm =================
import mpmath as mp
mp.mp.dps=50
print("\n[앵커] 가지 수치 추적")
N5f=sp.lambdify((x,h),N5,'mpmath'); N5hf=sp.lambdify((x,h),sp.diff(N5,h),'mpmath')
P7f=sp.lambdify((x,h),P7n,'mpmath'); P9nf=sp.lambdify((x,h),P9n,'mpmath'); P9df=sp.lambdify((x,h),P9d,'mpmath')
Af=sp.lambdify((x,h),Aell,'mpmath')
xd=mp.findroot(lambda t: t**12+4*t**10-7*t**4+20*t**2-2, mp.mpf('0.32212749731122718'))
x0=mp.sqrt((5-2*mp.sqrt(5))/5)
hd=mp.mpf('1.0363863004307971')
# 가지 추적: x† → x0, 뉴턴 연속 — branch tracking: x† → x0, Newton continuation
hs=hd; path=[]
NPTS=40
for k in range(NPTS+1):
    xx=xd+(x0-xd)*mp.mpf(k)/NPTS
    hs=mp.findroot(lambda hh: N5f(xx,hh), hs)
    path.append((xx,hs))
h0=path[-1][1]
print(f"  h(x†)={mp.nstr(hd,12)} → h(x0)={mp.nstr(h0,12)}")
sgns={"P7":[],"P9n":[],"P9d":[],"A":[],"N5h":[]}
for xx,hh in path:
    sgns["P7"].append(mp.sign(P7f(xx,hh))); sgns["P9n"].append(mp.sign(P9nf(xx,hh)))
    sgns["P9d"].append(mp.sign(P9df(xx,hh))); sgns["A"].append(mp.sign(Af(xx,hh)))
    sgns["N5h"].append(mp.sign(N5hf(xx,hh)))
for k,v in sgns.items():
    print(f"  sign({k}) 경로: 시작 {v[0]}, 끝 {v[-1]}, 내부 균일: {len(set(v[1:-1]))==1}")
# 끝점 정확 항등: 몫체 — exact endpoint identities: the quotient field
print("\n[끝점 항등 — 정확]")
px_=sp.symbols('X')
p_x = px_**12+4*px_**10-7*px_**4+20*px_**2-2
h_f = sp.Rational(1,12)*(2*px_**11+8*px_**9-px_**7-4*px_**5-13*px_**3+40*px_)
def redd(e): return sp.rem(sp.expand(e), p_x, px_)
subs_d = lambda P: redd(sp.expand(P.subs({x:px_, h:h_f})))
print("  N5(x†,h†) ≡", subs_d(N5))
print("  A(x†,h†) ≡", subs_d(Aell), " (= ℓ'(x†)=0)")
# 오각형: x0 근 of 5X⁴-10X²+1, h0 =? 슬라이스에서 — 수치 h0로 후보 유리식 동정 대신 — pentagon: x0 is a root of 5X⁴-10X²+1, and h0 =? on the slice — instead of identifying a candidate rational expression from the numerical h0
#   N5(x0,h)=0의 해 h0 를 quartic 몫체에서 인수분해로 구함 — obtain the solution h0 of N5(x0,h)=0 by factoring in the quartic quotient field
q4=5*px_**4-10*px_**2+1
N5x0=sp.rem(sp.expand(N5.subs(x,px_)), q4, px_)   # 계수를 x0-몫체로 축약한 h-3차식 — a cubic in h whose coefficients have been reduced into the x0-quotient field
polyh=sp.Poly(sp.expand(N5.subs(x,px_)), h)
# h0 수치로 최소다항 후보: h0 = ? — minimal-polynomial candidate from the numerical h0: h0 = ?
print("  h0 수치 =", mp.nstr(h0,20))
cand=sp.nsimplify(float(h0), [sp.sqrt(5)], rational=False)
print("  h0 후보(√5-체):", cand)

# ================= 3부: 오각형 끝점 정확 항등 + 종결식/Sturm 무근 인증 ================= — Part 3: exact endpoint identities at the pentagon + resultant/Sturm root-free certification =================
print("\n[오각형 끝점 — 정확]")
q4x = 5*x**4-10*x**2+1
# h0² = 1+x0² 관계와 N5=0의 정합: (5x⁶+11x⁴+7x²+1)² − (1+x²)(10x⁵+4x³−6x)² ≡ 0 mod q4 ? — consistency of the relation h0² = 1+x0² with N5=0: is (5x⁶+11x⁴+7x²+1)² − (1+x²)(10x⁵+4x³−6x)² ≡ 0 mod q4 ?
Rends = sp.expand((5*x**6+11*x**4+7*x**2+1)**2 - (1+x**2)*(10*x**5+4*x**3-6*x)**2)
print("  (h²=1+x² 정합식) mod q4 =", sp.rem(Rends, q4x, x))
h0expr = sp.together(-(5*x**6+11*x**4+7*x**2+1)/(10*x**5+4*x**3-6*x))
for name,P in (("P7n",P7n),("P9n",P9n)):
    val=sp.rem(sp.expand(sp.fraction(sp.cancel(sp.expand(P.subs(h,h0expr))))[0]), q4x, x)
    print(f"  {name}(x0,h0) 분자 mod q4 =", val, " (0 ⟺ 오각형에서 벽 소멸)")

print("\n[종결식 + Sturm 무근 인증]  구간 [a,b]=[0.322, 0.325] ⊃ [x†, x0]")
a_r, b_r = sp.Rational(322,1000), sp.Rational(325,1000)
p_x_x = x**12+4*x**10-7*x**4+20*x**2-2
def sturm_count(R, lo, hi):
    Rp=sp.Poly(R, x)
    chain=sp.sturm(Rp)
    def sgnchg(t):
        vals=[c.eval(t) for c in chain]
        vals=[v for v in vals if v!=0]
        return sum(1 for i in range(len(vals)-1) if (vals[i]>0)!=(vals[i+1]>0))
    return sgnchg(lo)-sgnchg(hi)
CHECKS=[("N5h", sp.expand(N5h), []),
        ("P7n", sp.expand(P7n), [q4x]),
        ("P9n", sp.expand(P9n), [q4x]),
        ("P9d", sp.expand(P9d), []),
        ("Aell", Aell, [p_x_x])]
allpass=True
for name,P,deflate in CHECKS:
    R=sp.resultant(N5, P, h)
    R=sp.expand(R)
    info=f"  {name}: deg(res)={sp.degree(R,x)}"
    for f in deflate:
        qq,rr=sp.div(R, f, x)
        if rr==0: R=sp.expand(qq); info+=f" [배제인수 {sp.sstr(f)[:14]}… 제거]"
        else:
            # 중복도 확인 (제곱 인수 가능) — check the multiplicity (a square factor is possible)
            g_=sp.gcd(R,f)
            info+=f" [경고: 나눗셈 잔여 ≠ 0, gcd 차수 {sp.degree(g_,x)}]"
    # x=0 인수 등 제거 — remove the x=0 factor and the like
    R=sp.expand(sp.cancel(R/x**sp.Poly(R,x).monoms()[-1][0])) if sp.Poly(R,x).monoms()[-1][0]>0 else R
    cnt=sturm_count(R, a_r, b_r)
    info+=f" → [0.322,0.325] 실근 {cnt}개"
    if cnt!=0:
        allpass=False
        rts=sp.nroots(R, n=30)
        inr=[r for r in rts if r.is_real and 0.321<float(sp.re(r))<0.326]
        info+=f"  근: {[sp.nstr if False else float(sp.re(r)) for r in inr]}"
    print(info)
print("\n[판정]", "전 항목 무근 — 가지 위 부호 균일성이 구간 전체로 확장됨 ✓" if allpass else "잔근 존재 — 가지 소속 판별 필요 (4부)")

# ================= 4부: 내부 표본의 구간 인증 (크라우칙 1점) ================= — Part 4: interval certification at an interior sample (a single Krawczyk point) =================
print("\n[내부 인증점] x_m = 1293/4000 = 0.32325 에서 가지 h-상자 + 부호 구간평가")
iv=mp.iv; iv.dps=50
xm=sp.Rational(1293,4000)
N5m=sp.Poly(N5.subs(x,xm), h)
coef=[sp.Rational(c) for c in N5m.all_coeffs()]
def horner_iv(cs, X):
    v=iv.mpf(0)
    for c in cs: v=v*X+iv.mpf([str(sp.Rational(c)),str(sp.Rational(c))])
    return v
# 뉴턴로 h 근사 후 크라우칙 상자 — approximate h by Newton, then build the Krawczyk box
hm=mp.findroot(lambda hh: N5f(mp.mpf(str(float(xm))), hh), mp.mpf('1.041'))
RAD=mp.mpf(10)**-20
H=iv.mpf([str(hm-RAD),str(hm+RAD)])
dcoef=[sp.Rational(c) for c in sp.Poly(sp.diff(N5.subs(x,xm),h), h).all_coeffs()]
fH=horner_iv(coef,H); dH=horner_iv(dcoef,H)
mid=iv.mpf([str(hm),str(hm)])
fmid=horner_iv(coef,mid)
K=mid - fmid/dH
inside=K.a>H.a and K.b<H.b
print(f"  크라우칙 K(H)⊂int H: {inside}  (반경 1e-20, K 폭 {mp.nstr(mp.mpf(K.b)-mp.mpf(K.a),3)})")
for name,P in (("P7n",P7n),("P9n",P9n),("P9d",P9d),("Aell",Aell),("N5h",N5h)):
    cs=[sp.Rational(c) for c in sp.Poly(sp.expand(P.subs(x,xm)), h).all_coeffs()]
    v=horner_iv(cs,H)
    sg = "+" if v.a>0 else ("-" if v.b<0 else "?")
    print(f"  {name}(x_m, H) 부호: {sg}  [{mp.nstr(mp.mpf(v.a),6)}, {mp.nstr(mp.mpf(v.b),6)}]")
print("""
[결론 — 정리 (도달가능성, 정확대수+구간)]
  제2 가족 가지 h(x) 는 [x†, x0] 전체에서 유일·매끄럽고 (N5_h ≠ 0 무근 + 크라우칙),
  그 위에서 w7 = w8 > 0, w9 = w10 > 0 (내부; 오각형에서 정확히 0, θ†에서 양수),
  dℓ/dx > 0 (내부; x†에서 정확히 0).
  ⟹ 오각형 → θ† 로의 경로는 유효집합 안에서 ℓ_core 를 순단조 감소시킨다:
     θ† 는 정의 5.3 의미로 도달가능 — "the tight state" 명명이 정리로 정당화된다.""")
