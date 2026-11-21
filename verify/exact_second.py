# exact_second.py — exact-algebra and interval certification of the second-order
# sufficient condition at theta-dagger (Theorem 6.6 / Theorem D).
#
# On the second family ell(x) = h(x) + (x^2+1)/(2x) with N5(x,h) = 0, computes
# ell''(x-dagger) in Q(x-dagger) = Q[x]/(p_x) and checks:
#   (i) exact identities N5(x-dagger, h-dagger) = 0 and ell'(x-dagger) = 0 mod p_x;
#   (ii) exact degree-11 polynomial representative of ell''(x-dagger);
#   (iii) interval-Newton enclosure of x-dagger and positivity of ell''(x-dagger);
#   (iv) independent 80-digit numerical implicit-function second difference.
#
# With rank-4 LICQ and strictly positive KKT multipliers on {w1, w4, w5, w6},
# the critical cone is spanned by c'(x-dagger), so ell''(x-dagger) > 0 certifies
# strict local minimality (Nocedal-Wright, Theorem 12.6).

import sympy as sp
from mpmath import mp, iv, mpf, findroot, nstr as mstr

x, h = sp.symbols('x h')

# Exact defining polynomials of the second family and theta-dagger.
px = sp.expand(x**12 + 4*x**10 - 7*x**4 + 20*x**2 - 2)          # x†의 최소다항식 — minimal polynomial of x†
N5 = sp.expand(8*h**3*x**3 + 4*h**2*x**4 + 4*h**2*x**2
               + 2*h*x**5 - 4*h*x**3 - 6*h*x
               + x**6 + 3*x**4 + 3*x**2 + 1)                      # 제2족 평면곡선 — the plane curve of the second family
h_field = sp.Rational(1,12)*(2*x**11 + 8*x**9 - x**7 - 4*x**5 - 13*x**3 + 40*x)  # h† ∈ ℚ(x†)

def red(e):
    """ℚ[x]/(px) 감축 (작은 다항식만 다루므로 expand+rem으로 충분). — Reduction into ℚ[x]/(px) (only small polynomials occur, so expand+rem suffices)."""
    return sp.rem(sp.expand(e), px, x)

def inv(e):
    """몫체 역원 — px 기약이므로 0 아닌 원소는 전부 가역. 실패 자체가 0-증명. — Inverse in the quotient field — px is irreducible, so every nonzero element is invertible; a failure is itself a proof of being 0."""
    return sp.invert(red(e), px, x)

def sub_h(e):
    """h ↦ h† 대입 후 감축. — Substitute h ↦ h†, then reduce."""
    return red(e.subs(h, h_field))

# ---- (i) 정확 항등식 검사 ---- — (i) check the exact identities ----
chk_curve = sub_h(N5)
assert chk_curve == 0, f"N5(x†,h†) != 0: {chk_curve}"
print("[i] N5(x†, h†) ≡ 0  (mod p_x)   — 곡선 위 확인: OK")

Fx  = sub_h(sp.diff(N5, x));   Fh  = sub_h(sp.diff(N5, h))
Fxx = sub_h(sp.diff(N5, x, 2)); Fhh = sub_h(sp.diff(N5, h, 2))
Fxh = sub_h(sp.diff(N5, x, h))

invFh = inv(Fh)   # 성공 ⇒ ∂N5/∂h(x†,h†) ≠ 0 인증 (음함수 정리 적용 가능) — success ⇒ certifies ∂N5/∂h(x†,h†) ≠ 0 (so the implicit function theorem applies)
print("[i] ∂N5/∂h(x†,h†) ≠ 0            — 몫체 가역성으로 인증: OK")

hp = red(-Fx * invFh)                                  # h'(x†)
# b(x) = (x²+1)/(2x) ⇒ b' = 1/2 − 1/(2x²),  b'' = 1/x³
ellp = red(hp + sp.Rational(1,2) - inv(2*x**2))        # ℓ'(x†)
assert ellp == 0, f"ell'(x†) != 0: {ellp}"
print("[i] ℓ'(x†) ≡ 0  (mod p_x)        — 임계점 항등식: OK")

# ---- (ii) ℓ''(x†) 정확 계산 ---- — (ii) exact computation of ℓ''(x†) ----
hpp = red(-(Fxx + 2*Fxh*hp + Fhh*red(hp*hp)) * invFh)  # h''(x†)
ellpp = red(hpp + inv(x**3))                           # ℓ''(x†) = h'' + 1/x³
P = sp.Poly(ellpp, x)
print("\n[ii] ℓ''(x†) ∈ ℚ(x†) 정확 표현 (기저 1, x, …, x^11):")
print("     ℓ''(x†) =", ellpp)
coeffs = [sp.Rational(c) for c in P.all_coeffs()]      # 내림차순 — descending order

# ---- (iii) 구간 인증 ---- — (iii) interval certification ----
mp.dps = 80
pxf  = lambda t: t**12 + 4*t**10 - 7*t**4 + 20*t**2 - 2
xr = findroot(pxf, mpf('0.3221274973112271792397142057801235136978'))

iv.dps = 80
def iv_rat(q):
    return iv.mpf(int(q.p)) / iv.mpf(int(q.q))

def iv_poly(cs, X):
    """내림차순 유리 계수 cs 의 구간 호너 평가. — Interval Horner evaluation of the rational coefficients cs, given in descending order."""
    acc = iv.mpf(0)
    for c in cs:
        acc = acc * X + iv_rat(sp.Rational(c))
    return acc

px_desc  = [sp.Rational(c) for c in sp.Poly(px, x).all_coeffs()]
dpx_desc = [sp.Rational(c) for c in sp.Poly(sp.diff(px, x), x).all_coeffs()]

rad = mpf('1e-60')
X  = iv.mpf([xr - rad, xr + rad])
m  = iv.mpf(xr)
Nx = m - iv_poly(px_desc, m) / iv_poly(dpx_desc, X)    # 구간 뉴턴 연산자 — interval Newton operator

def endpoints(I):
    a, b = I.a, I.b
    a = a.a if hasattr(a, 'a') else a
    b = b.b if hasattr(b, 'b') else b
    return mpf(a), mpf(b)

Na, Nb = endpoints(Nx); Xa, Xb = endpoints(X)
assert Xa < Na and Nb < Xb, "구간 뉴턴 축약 실패"
print("\n[iii] 구간 뉴턴: N(X) ⊂ int(X)  — X 안에 p_x의 근이 유일 존재: OK")
print(f"      X = [{mstr(Xa, 25)}, {mstr(Xb, 25)}]")

E = iv_poly(coeffs, X)
Ea, Eb = endpoints(E)
assert Ea > 0, f"ℓ''(x†) 양수성 인증 실패: [{Ea}, {Eb}]"
print(f"[iii] ℓ''(x†) ∈ [{mstr(Ea, 30)}, {mstr(Eb, 30)}]  ⊂ (0,∞): OK")

# ---- (iv) 독립 크로스체크: 수치 음함수 + 중심 2차 차분 ---- — (iv) independent cross-check: numerical implicit function + centred second difference ----
N5f  = sp.lambdify((x, h), N5, 'mpmath')
href = mpf('1.036386300430797104445252761398636908419')
def h_of(t):
    return findroot(lambda y: N5f(t, y), href)
def ell_of(t):
    return h_of(t) + (t**2 + 1) / (2*t)
d = mpf('1e-20')
num2 = (ell_of(xr + d) - 2*ell_of(xr) + ell_of(xr - d)) / d**2
exact_val = iv_poly(coeffs, iv.mpf(xr))
ex_mid = (endpoints(exact_val)[0] + endpoints(exact_val)[1]) / 2
print("\n[iv] 크로스체크 (서로 독립한 두 경로):")
print("     정확대수 평가   :", mstr(ex_mid, 40))
print("     2차 차분(80자리):", mstr(num2, 40))
print("     |차|            :", mstr(abs(ex_mid - num2), 5))
assert abs(ex_mid - num2) < mpf('1e-30'), "크로스체크 불일치"   # 2차 차분 절단오차 O(d²·ℓ⁗) ≈ 1e-34 수준 — the truncation error of the second difference is on the order of O(d²·ℓ⁗) ≈ 1e-34

print("\nPASS: ell''(x-dagger) > 0 certified in Q[x]/(p_x) and interval arithmetic (strict local minimum)")
