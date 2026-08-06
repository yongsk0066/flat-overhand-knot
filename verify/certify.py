"""
10 검증: theta†의 형식 증서 3종 + Lean 증서 계산.
10 verification: the three formal certificates for theta† + computation of the Lean certificate.

C1 (가지 기약성): N5(x,h) 가 QQ 위 기약 다항식 => 항등식 a=h, b=g+x, ell=h+b 는
C1 (irreducibility of the branch): N5(x,h) is irreducible over QQ => the identities a=h, b=g+x, ell=h+b
    (다항식 나눗셈으로 증명됐으므로) 곡선 {N5=0} '전체'에서 성립 — 가지 기약성 간극 해소.
    hold on the 'whole' curve {N5=0} (since they were proved by polynomial division) — closing the branch-irreducibility gap.
C2 (정확 식별): res_h(N5, K) = -2^18 * x^12 * (x^2+1)^2 * p_x  [정확 대수 항등식]
C2 (exact identification): res_h(N5, K) = -2^18 * x^12 * (x^2+1)^2 * p_x  [an exact algebraic identity]
    => (N5,K)=0 의 임의의 해의 x좌표는 {0, ±i} 또는 p_x의 근.
    => the x-coordinate of any solution of (N5,K)=0 lies in {0, ±i} or is a root of p_x.
C3 (존재 증서, Krawczyk): 상자 X = [x†±r]x[h†±r] (r=1e-30, iv.dps=80)에서
C3 (existence certificate, Krawczyk): on the box X = [x†±r]x[h†±r] (r=1e-30, iv.dps=80),
    K(X) ⊂ int(X) => (N5,K)=0 의 정확해가 X 안에 유일하게 존재.
    K(X) ⊂ int(X) => an exact solution of (N5,K)=0 exists uniquely inside X.
    C2와 결합, 상자가 x=0을 제외하므로 그 해의 x좌표는 p_x의 근 [완전 증서].
    combined with C2, and since the box excludes x=0, the x-coordinate of that solution is a root of p_x [complete certificate].
C4 (Lean 증서): 오각형 비근성 — 125*f(x) - (9664 - 4320*s) 를
C4 (Lean certificate): the pentagon is not a root — the coefficients that express 125*f(x) - (9664 - 4320*s)
    A = s*x^2-(s-2), B = s^2-5 의 결합으로 표현하는 계수 (linear_combination 용).
    as a combination of A = s*x^2-(s-2) and B = s^2-5 (for use with linear_combination).
"""
from sympy import (symbols, Rational, S, cancel, together, fraction, diff,
                   resultant, factor_list, Poly, groebner, reduced, expand, sqrt, nsimplify)

x, h, b, s = symbols('x h b s')

# ---------- N5 구성 (03/09의 기하: g=(1-x^2)/2x 대입, sI - h 의 분자) ---------- — construction of N5 (geometry from 03/09: substitute g=(1-x^2)/2x, take the numerator of sI - h) ----------
def refl_lin(d):
    dx, dy = d
    n = dx**2 + dy**2
    return lambda v: (cancel(((dx**2-dy**2)*v[0] + 2*dx*dy*v[1])/n),
                      cancel((2*dx*dy*v[0] + (dy**2-dx**2)*v[1])/n))
def cross(u, v): return u[0]*v[1] - u[1]*v[0]
def sub2(u, v): return (u[0]-v[0], u[1]-v[1])

g = (1 - x**2)/(2*x)
E = (S(0), S(0)); F = (-x, S(1)); G = (g, S(1)); H = (h, S(0))
d2 = refl_lin(sub2(F, E))((S(1), S(0)))
w0 = refl_lin(sub2(G, H))((S(1), S(0)))
sI = cancel(cross(sub2(F, G), d2) / cross(w0, d2))
N5 = Poly(fraction(together(sI - h))[0], x, h)
# 부호/배율 정규화 (정수 계수, 내용 제거) — sign/scale normalization (integer coefficients, content removed)
N5 = N5.primitive()[1]
print("N5 =", N5.as_expr())

print("\n=== C1: N5 기약성 ===")
fl = factor_list(N5.as_expr(), x, h)
print("  factor_list:", [(f.as_expr() if hasattr(f,'as_expr') else f, m) for f, m in fl[1]])
nontrivial = [(f, m) for f, m in fl[1] if Poly(f, x, h).total_degree() > 0]
assert len(nontrivial) == 1 and nontrivial[0][1] == 1, "N5가 기약이 아님!"
print("  [PASS] N5는 QQ 위 기약 — 항등식들은 곡선 {N5=0} 전체에서 성립")

print("\n=== C2: res_h(N5, K)의 정확 인수분해 ===")
K = expand(2*x**2*diff(N5.as_expr(), x) + (1 - x**2)*diff(N5.as_expr(), h))
R = resultant(N5.as_expr(), K, h)
p_x = x**12 + 4*x**10 - 7*x**4 + 20*x**2 - 2
quot = cancel(R / (x**12 * (x**2+1)**2 * p_x))
print(f"  res / (x^12 (x^2+1)^2 p_x) = {quot}")
assert quot.is_number and quot != 0, "종결식 인수 구조 불일치!"
print(f"  [PASS] res_h(N5,K) = {quot} * x^12 * (x^2+1)^2 * p_x  [정확 항등식]")

print("\n=== C3: Krawczyk 존재 증서 ===")
from mpmath import iv, mp
iv.dps = 80
mp.dps = 80
# 계수 사전 추출 — extract the coefficients in advance
def poly_terms(e):
    return [(m, Rational(c)) for m, c in Poly(e, x, h).terms()]
T_N5 = poly_terms(N5.as_expr())
T_K = poly_terms(K)
T_N5x = poly_terms(diff(N5.as_expr(), x)); T_N5h = poly_terms(diff(N5.as_expr(), h))
T_Kx = poly_terms(diff(K, x)); T_Kh = poly_terms(diff(K, h))

def iv_eval(terms, X, H):
    tot = iv.mpf(0)
    for (i, j), c in terms:
        t = iv.mpf(int(c.p)) / iv.mpf(int(c.q))
        tot += t * (X**i) * (H**j)
    return tot

xd = mp.mpf('0.3221274973112271792397142057801235136978')
hd = mp.mpf('1.036386300430797104445252761398636908419')
# 중심 정련 (mp 뉴턴 2D, 80자리) — refine the centre (2D Newton in mp, 80 digits)
def mpev(terms, X, H):
    tot = mp.mpf(0)
    for (i, j), c in terms:
        tot += mp.mpf(int(c.p))/mp.mpf(int(c.q)) * (X**i) * (H**j)
    return tot
for _ in range(8):
    F1 = mpev(T_N5, xd, hd); F2 = mpev(T_K, xd, hd)
    J11 = mpev(T_N5x, xd, hd); J12 = mpev(T_N5h, xd, hd)
    J21 = mpev(T_Kx, xd, hd); J22 = mpev(T_Kh, xd, hd)
    det = J11*J22 - J12*J21
    xd -= ( J22*F1 - J12*F2)/det
    hd -= (-J21*F1 + J11*F2)/det
print(f"  정련 중심: x = {mp.nstr(xd, 40)}")
print(f"            h = {mp.nstr(hd, 40)}")

r = mp.mpf('1e-30')
Xiv = iv.mpf([xd - r, xd + r]); Hiv = iv.mpf([hd - r, hd + r])
# Y = 중점 야코비안의 역 (점 행렬) — Y = inverse of the Jacobian at the midpoint (a point matrix)
J11m = mpev(T_N5x, xd, hd); J12m = mpev(T_N5h, xd, hd)
J21m = mpev(T_Kx, xd, hd); J22m = mpev(T_Kh, xd, hd)
detm = J11m*J22m - J12m*J21m
Y = [[ J22m/detm, -J12m/detm], [-J21m/detm, J11m/detm]]
Yiv = [[iv.mpf(mp.nstr(Y[i][j], 78)) for j in range(2)] for i in range(2)]
# F(m) 구간 평가 (m = 퇴화 구간) — interval evaluation of F(m) (m = a degenerate interval)
xm = iv.mpf(mp.nstr(xd, 78)); hm = iv.mpf(mp.nstr(hd, 78))
Fm = [iv_eval(T_N5, xm, hm), iv_eval(T_K, xm, hm)]
# J(X) 구간 야코비안 — the interval Jacobian J(X)
JX = [[iv_eval(T_N5x, Xiv, Hiv), iv_eval(T_N5h, Xiv, Hiv)],
      [iv_eval(T_Kx, Xiv, Hiv), iv_eval(T_Kh, Xiv, Hiv)]]
# I - Y*J(X)
M = [[ (iv.mpf(1) if i == j else iv.mpf(0))
       - (Yiv[i][0]*JX[0][j] + Yiv[i][1]*JX[1][j]) for j in range(2)] for i in range(2)]
DX = [Xiv - xm, Hiv - hm]
Kop = [xm - (Yiv[0][0]*Fm[0] + Yiv[0][1]*Fm[1]) + (M[0][0]*DX[0] + M[0][1]*DX[1]),
       hm - (Yiv[1][0]*Fm[0] + Yiv[1][1]*Fm[1]) + (M[1][0]*DX[0] + M[1][1]*DX[1])]
inside = (Xiv.a < Kop[0].a and Kop[0].b < Xiv.b and Hiv.a < Kop[1].a and Kop[1].b < Hiv.b)
wid_x = mp.mpf(Kop[0].delta.a); wid_h = mp.mpf(Kop[1].delta.a)
print(f"  K(X) 폭: x측 {mp.nstr(wid_x, 5)}, h측 {mp.nstr(wid_h, 5)}  vs 상자 폭 {mp.nstr(2*r, 5)}  -> 포함 = {inside}")
assert inside, "Krawczyk 포함 실패 — r 조정 필요"
print("  [PASS] K(X) ⊂ int(X): (N5,K)=0의 정확해가 상자 안에 유일 존재")
print(f"  상자 x-구간 = [{mp.nstr(mp.mpf(Xiv.a.a),35)}, {mp.nstr(mp.mpf(Xiv.b.a),35)}]  (0 제외)")
print("  ==> C2와 결합: 그 해의 x좌표는 p_x의 근 [완전 증서]")

print("\n=== C4: Lean 증서 — 오각형 x0는 p_x의 근이 아님 ===")
A = s*x**2 - (s - 2)
B = s**2 - 5
P = expand(125*p_x - (9664 - 4320*s))
gb = groebner([A, B], x, s, order='lex')
q, rem = reduced(P, gb)
print(f"  reduced 나머지 = {rem}  (0이어야)")
assert rem == 0
# P = alpha*A + beta*B, alpha = a0(x)+a1(x)s, beta = b0(x): s-차수별 계수 비교의 삼각 구조 — the triangular structure obtained by comparing coefficients degree by degree in s
# s^2: a1(x^2-1) + b0 = 0;  s^1: a0(x^2-1) + 2a1 = 4320;  s^0: 2a0 - 5b0 = 125p_x - 9664
# a1 = 2160 + (x^2-1)*gamma, a0 = -2*gamma 로 두면 — setting a1 = 2160 + (x^2-1)*gamma and a0 = -2*gamma gives
# gamma * (5(x^2-1)^2 - 4) = 125p_x - 9664 - 10800(x^2-1),  그리고 5(x^2-1)^2-4 = 5x^4-10x^2+1 (!) — and 5(x^2-1)^2-4 = 5x^4-10x^2+1 (!)
from sympy import div
pent = 5*x**4 - 10*x**2 + 1
num = expand(125*p_x - 9664 - 10800*(x**2 - 1))
gamma, remg = div(num, pent, x)
print(f"  분모 5(x^2-1)^2-4 = {expand(5*(x**2-1)**2-4)}  (= 오각형 접촉 다항식!)")
print(f"  gamma = {gamma},  나눗셈 나머지 = {remg}")
assert remg == 0, "gamma가 다항식이 아님 — 안사츠 확장 필요"
a1 = expand(2160 + (x**2 - 1)*gamma)
a0 = expand(-2*gamma)
b0 = expand(-a1*(x**2 - 1))
alpha = a0 + a1*s
beta = b0
check = expand(P - (alpha*A + beta*B))
print(f"  alpha = {alpha}")
print(f"  beta  = {beta}")
print(f"  잔차  = {check}  (0이어야)")
assert check == 0
print("  [PASS] 125*p_x = (9664-4320s) + alpha*A + beta*B")
print("  수치 확인: 9664^2 =", 9664**2, ", 4320^2*5 =", 4320**2*5, " (다름 => f(x0)≠0)")
print("\n==== ALL CERTIFICATES PASS ====")
