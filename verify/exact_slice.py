"""
2단계 검증: 슬라이스의 전 항등식을 sympy 정확 대수로 유도·증명.
Stage-2 verification: derive and prove every identity of the slice with sympy exact algebra.

수치 검증(check_*.py)과의 차이:
difference from the numerical verification (check_*.py):
  - G, H, k, t 공식을 '가정 후 확인'하지 않고 접촉 등식 C1-C4에서 '유도'한다.
  - the formulas for G, H, k, t are 'derived' from the contact equalities C1-C4, not 'assumed and then checked'.
  - 모든 항등식을 Q(x) 위 유리식 등식 / Q(sqrt5) 위 대수 등식으로 증명한다 (표본점 없음).
  - every identity is proved as a rational-function equality over Q(x) / an algebraic equality over Q(sqrt5) (no sample points).
  - Lean 형식화에 쓸 linear_combination 증서를 사전 검증한다.
  - the linear_combination certificates intended for the Lean formalization are pre-verified.

실행: verify/.venv/bin/python verify/exact_slice.py
Run: verify/.venv/bin/python verify/exact_slice.py
"""
from sympy import (symbols, Rational, sqrt, simplify, cancel, expand, factor,
                   together, solve, diff, S, nsimplify, radsimp)

x = symbols('x', positive=True)
R = Rational

ok_all = True
def check(name, expr):
    global ok_all
    z = simplify(expr)
    good = (z == 0)
    print(f"[{'PASS' if good else 'FAIL'}] {name}" + ("" if good else f"   residual = {z}"))
    ok_all = ok_all and good

def report(name, val):
    print(f"       {name} = {val}")

# ---------- 기하 도구 (전부 유리식 — 정규화 불필요) ---------- — geometric tools (all rational functions — no normalization needed) ----------
def refl_lin(d):
    dx, dy = d
    n = dx**2 + dy**2
    return lambda v: (cancel(((dx**2-dy**2)*v[0] + 2*dx*dy*v[1])/n),
                      cancel((2*dx*dy*v[0] + (dy**2-dx**2)*v[1])/n))

def refl_pt(P, d):
    L = refl_lin(d)
    def rho(Q):
        w = L((Q[0]-P[0], Q[1]-P[1]))
        return (cancel(P[0]+w[0]), cancel(P[1]+w[1]))
    return rho

def cross(a, b): return a[0]*b[1] - a[1]*b[0]
def dot(a, b):  return a[0]*b[0] + a[1]*b[1]
def sub(a, b):  return (a[0]-b[0], a[1]-b[1])
def addm(a, d, s): return (a[0]+s*d[0], a[1]+s*d[1])

def line_meet_param(P1, d1, P2, d2):
    """P1 + s d1 = P2 + t d2 의 s. — The s in P1 + s d1 = P2 + t d2."""
    return cancel(cross(sub(P2, P1), d2) / cross(d1, d2))

print("======== A. 슬라이스의 유도 (C1-C4를 풀어서) ========")
E = (S(0), S(0)); F = (-x, S(1))
L1 = refl_lin(sub(F, E))
d2 = L1((S(1), S(0)))
check("A1: d2 = ((x^2-1)/(1+x^2), -2x/(1+x^2))",
      (d2[0] - (x**2-1)/(1+x**2)) + (d2[1] + 2*x/(1+x**2)))

# C1에서 g 유도: G=(g,1)가 e_A=(E,d2) 위 — derive g from C1: G=(g,1) lies on e_A=(E,d2)
gs = symbols('gs')
g_roots = solve(cross((gs, S(1)), d2), gs)
report("C1의 해 g", g_roots)
g = [r0 for r0 in g_roots if simplify(r0 - (1-x**2)/(2*x)) == 0]
assert g, "g = (1-x^2)/(2x)가 C1의 해가 아님"
g = g[0]
check("A2: C1 ⟹ g = (1-x^2)/(2x)  (해 유일)", S(0) if len(g_roots) == 1 else S(1))

# C2에서 h 유도: F가 rho_HG({y=0}) = (H, w0) 위 — derive h from C2: F lies on rho_HG({y=0}) = (H, w0)
hs_ = symbols('hs_')
cd2_h = (g - hs_, S(1))                       # 크리스2 방향 (미지 h) — direction of crease2 (h still unknown)
w0_h = refl_lin(cd2_h)((S(1), S(0)))
h_roots = solve(together(cross(sub(F, (hs_, S(0))), w0_h)), hs_)
report("C2의 해 h", h_roots)
h_expect = (1-3*x**2)/(2*x)
h_good = [r0 for r0 in h_roots if simplify(r0 - h_expect) == 0]
assert h_good, "h = (1-3x^2)/(2x)가 C2의 해가 아님"
h = h_good[0]
print(f"       (다른 근 = 다른 가지: {[simplify(r0) for r0 in h_roots if simplify(r0-h_expect)!=0]})")

H = (h, S(0)); G = (g, S(1))
w0 = refl_lin(sub(G, H))((S(1), S(0)))        # pi2 진행방향 r — r, the progression direction of pi2
check("A3: r = ((x^2-1)/(1+x^2), 2x/(1+x^2))",
      (w0[0] - (x**2-1)/(1+x**2)) + (w0[1] - 2*x/(1+x**2)))
check("A4: r는 단위벡터", dot(w0, w0) - 1)

# I와 lambda — I and lambda
lam = line_meet_param(G, w0, F, d2)
check("A5: lambda = (1+x^2)^2/(4x(1-x^2))", lam - (1+x**2)**2/(4*x*(1-x**2)))

# C3+C4에서 크리스3 유도: F를 지나고, rho_3(H)가 l2=(G,w0) 위 — derive crease3 from C3+C4: it passes through F, and rho_3(H) lies on l2=(G,w0)
kap = symbols('kap')
Q1 = addm(G, w0, kap)                         # 크리스3의 l2쪽 끝점 후보 — candidate endpoint of crease3 on the l2 side
rho3 = refl_pt(F, sub(Q1, F))
c4 = together(cross(sub(rho3(H), G), w0))
kap_roots = solve(c4, kap)
report("C4의 해 kappa", [simplify(r0) for r0 in kap_roots])
k_expect = (1-3*x**2)/(2*x)
k_good = [r0 for r0 in kap_roots if simplify(r0 - k_expect) == 0]
assert k_good, "k = (1-3x^2)/(2x)가 C4의 해가 아님"
k = k_good[0]
K = addm(G, w0, k)
rho3 = refl_pt(F, sub(K, F))
t_val = cancel(dot(sub(rho3(H), G), w0))
check("A6: t = (1-x^2)/x  [asset-1의 가지]", t_val - (1-x**2)/x)

print("\n======== B. 간극 함수와 단조성 ========")
Delta = cancel(k - lam)
check("B1: Delta = (5x^4-10x^2+1)/(4x(1-x^2))", Delta - (5*x**4-10*x**2+1)/(4*x*(1-x**2)))
dDelta = cancel(diff(Delta, x))
check("B2: Delta' = -(5x^6-5x^4+7x^2+1)/(4x^2(1-x^2)^2)",
      dDelta + (5*x**6-5*x**4+7*x**2+1)/(4*x**2*(1-x**2)**2))
check("B3: 분자 SOS 분해  5x^6-5x^4+7x^2+1 = 5x^2(x^2-1/2)^2 + (23/4)x^2 + 1",
      expand(5*x**6-5*x**4+7*x**2+1 - (5*x**2*(x**2-R(1,2))**2 + R(23,4)*x**2 + 1)))
print("       ⟹ 분자 > 0 (x>0에서 각 항 양수) ⟹ Delta' < 0  [정확 증명 완료]")

print("\n======== C. ell_core ========")
M1 = (-x/2, R(1,2))
M2 = ((h+g)/2, R(1,2))
s12 = cancel(M2[0] - M1[0])
check("C1': |M1M2| = (1-x^2)/(2x)", s12 - (1-x**2)/(2*x))
s23 = line_meet_param(M2, w0, F, sub(K, F))
check("C2': |M2M3| = (1-x^2)/(2x)  (전방: 계수 양수)", cancel(s23) - (1-x**2)/(2*x))
check("C3': ell_core = 1/x - x", cancel(s12 + s23) - (1/x - x))

print("\n======== D. 타이트 상태 (x0, 정확 대수: Q(sqrt5)) ========")
x0 = sqrt(1 - 2/sqrt(5))
def at0(e): return e.subs(x, x0)
P5 = [E, (at0(h), S(0)), (at0(g), S(1)),
      (at0(K[0]), at0(K[1])), (-x0, S(1))]     # E, H, G, K, F
names = ['E','H','G','K','F']
L2ref = simplify((P5[1][0]-P5[0][0])**2 + (P5[1][1]-P5[0][1])**2)  # |EH|^2
report("|EH|^2 (기준)", simplify(nsimplify(L2ref, [sqrt(5)])))
for i in range(5):
    a, b = P5[i], P5[(i+1) % 5]
    d2s = simplify((b[0]-a[0])**2 + (b[1]-a[1])**2 - L2ref)
    check(f"D1: |{names[i]}{names[(i+1)%5]}|^2 = |EH|^2", d2s)
cos108 = -(sqrt(5)-1)/4
for i in range(5):
    a, v, c = P5[(i-1) % 5], P5[i], P5[(i+1) % 5]
    dp = (a[0]-v[0])*(c[0]-v[0]) + (a[1]-v[1])*(c[1]-v[1])
    check(f"D2: 꼭짓점 {names[i]}의 내각 = 108도  (cos = -(sqrt5-1)/4)",
          simplify(dp - L2ref*cos108))
exitd = refl_lin(sub((at0(K[0]), at0(K[1])), (-x0, S(1))))((at0(w0[0]), at0(w0[1])))
d2_0 = (at0(d2[0]), at0(d2[1]))
check("D3: 꼬리 사이각 = 72도  (cos = (sqrt5-1)/4)",
      simplify(dot(exitd, d2_0) - (sqrt(5)-1)/4))

print("\n======== E. Lean 증서(linear_combination) 사전 검증 ========")
s_, X = symbols('s_ X')
Gp = 5*X**4 - 10*X**2 + 1
Ap = s_*X**2 - (s_ - 2)      # hx : s*x^2 = s-2 의 (좌-우) — hx : (left − right) of s*x^2 = s-2
Bp = s_**2 - 5               # hs : s^2 = 5 의 (좌-우) — hs : (left − right) of s^2 = 5
check("E1: x0_root 증서  G = (s x^2 - s - 2)·A - (x^2-1)^2·B",
      expand(Gp - ((s_*X**2 - s_ - 2)*Ap - (X**2-1)**2*Bp)))
check("E2: 변 동일성 증서  (1-3x^2)^2 - (1+x^2)(2x)^2 = P(x)",
      expand((1-3*X**2)**2 - (1+X**2)*(2*X)**2 - Gp))
z = symbols('z')
check("E3: 폐포 증서  (z-1)(1+z+z^2+z^3+z^4) = z^5-1",
      expand((z-1)*(1+z+z**2+z**3+z**4) - (z**5-1)))
check("E4: p(y) SOS  5y^3-5y^2+7y+1 = 5y(y-1/2)^2 + (23/4)y + 1",
      expand(5*X**3-5*X**2+7*X+1 - (5*X*(X-R(1,2))**2 + R(23,4)*X + 1)))

print("\n========", "ALL PASS — 슬라이스 전 항목 정확 대수 증명 완료" if ok_all else "실패 있음", "========")
