"""
07 검증 1: 오각형 차단의 Farkas 증서 (정확 대수).
07 verification 1: Farkas certificate for the blocking of the pentagon (exact algebra).

내용:
Contents:
  P1. 5변수 유리식으로 벽 함수 w1,w4,w5,w6,w7,w8과 ell_core를 구성
  P1. Build the wall functions w1,w4,w5,w6,w7,w8 and ell_core as rational functions of the 5 variables
      (단위화 생략 — 양수 배율은 벽 집합과 원뿔을 바꾸지 않음; 반사·방향은 전부 유리식).
      (normalization is omitted — a positive rescaling changes neither the wall set nor the cone; reflections and directions are all rational).
  P2. 오각형 점에서 grad ell = (1/2,...,1/2) 를 '정확히' 검증 (05의 미해명 관찰의 승격).
  P2. Verify 'exactly' that grad ell = (1/2,...,1/2) at the pentagon point (promotion of the unexplained observation in 05).
  P3. Farkas: grad ell = sum lambda_j grad w_j, lambda_j >= 0 의 정확해 탐색
  P3. Farkas: search for an exact solution of grad ell = sum lambda_j grad w_j with lambda_j >= 0
      (6벽 중 5-부분집합 선형계를 정확히 풀고 부호를 고정밀로 판정).
      (solve the linear system for each 5-subset of the 6 walls exactly, then decide the signs in high precision).
  P4. 슬라이스 위에서 d(w_j)/dx (편미분)의 부호를 x의 함수로 정확 분석
  P4. Exact analysis, as a function of x, of the sign of d(w_j)/dx (partial derivative) on the slice
      — 전역 잼 없음의 벡터장 논증(V1 = -e_x)의 적용 범위 확정.
      — fixing the range of applicability of the vector-field argument (V1 = -e_x) for the absence of a global jam.

실행: verify/.venv/bin/python verify/exact_farkas.py
Run: verify/.venv/bin/python verify/exact_farkas.py
"""
from sympy import (symbols, sqrt, Rational, simplify, cancel, together, diff,
                   Matrix, S, nsimplify, factor, solve)

x, g, h, a, b = symbols('x g h a b', real=True)
R = Rational

def refl_lin(d):
    dx, dy = d
    n = dx**2 + dy**2
    return lambda v: (cancel(((dx**2-dy**2)*v[0] + 2*dx*dy*v[1])/n),
                      cancel((2*dx*dy*v[0] + (dy**2-dx**2)*v[1])/n))

def cross(u, v): return u[0]*v[1] - u[1]*v[0]
def dot(u, v):  return u[0]*v[0] + u[1]*v[1]
def sub2(u, v): return (u[0]-v[0], u[1]-v[1])
def add2(u, v): return (u[0]+v[0], u[1]+v[1])

E = (S(0), S(0)); F = (-x, S(1)); G = (g, S(1)); H = (h, S(0))
d2 = refl_lin(sub2(F, E))((S(1), S(0)))          # 단위 (반사는 노름 보존) — unit (reflections preserve the norm)
w0 = refl_lin(sub2(G, H))((S(1), S(0)))          # 단위 — unit
Q1 = (g + a*w0[0], 1 + a*w0[1])
Q0 = (h + b*w0[0], b*w0[1])
rl3 = refl_lin(sub2(Q0, Q1))
rho3 = lambda P: add2(Q1, rl3(sub2(P, Q1)))
dD = rl3(w0)                                     # 단위 — unit
e1p = rho3(G)

# 방향 상수(부호)는 오각형 근방에서 고정 — 수치로 한 번 결정해 정확 상수로 사용 — the orientation constants (signs) are fixed near the pentagon — determined numerically once, then used as exact constants
import math
x0n = math.sqrt(1 - 2/math.sqrt(5))
g0n = (1-x0n**2)/(2*x0n); h0n = (1-3*x0n**2)/(2*x0n)
w0n = None
def _num(e, X):  return float(e.subs(X).evalf(30))
b0expr = dot(sub2(F, H), w0)
X0n = {x: x0n, g: g0n, h: h0n, a: h0n, b: float(b0expr.subs({x: x0n, g: g0n, h: h0n}).evalf(30))}
sA = 1 if _num(cross(sub2(F, E), d2), X0n) > 0 else -1
pd_probe = (S(0), S(0))
# D 내부 표본: M3 + v3 대신, 레이아웃 진행 = rl3(w0) 방향으로 Q1에서 진입한 점 — interior sample of D: instead of M3 + v3, the point entered from Q1 along the layout progression = the rl3(w0) direction
pD = add2(Q1, dD)
sD = 1 if _num(cross(sub2(pD, e1p), dD), X0n) > 0 else -1

# ---- 벽 함수 (실현측 >= 0, 배율 양수 자유) ---- — wall functions (realizable side >= 0, positive rescaling is free) ----
w1 = -sA * cross(sub2(G, E), d2)
w7 =  sA * cross(sub2(Q0, E), d2) - 1
w4 = -sD * cross(sub2(H, e1p), dD)
w8 =  sD * cross(sub2(F, e1p), dD) - 1
sI = cancel(cross(sub2(F, G), d2) / cross(w0, d2))     # I = G + sI*w0
w5 = a - sI
w6raw = cross(sub2(E, rho3(H)), dD)
s6 = 1 if _num(w6raw, {x:0.25, g:(1-0.0625)/0.5, h:(1-3*0.0625)/0.5,
                       a:(1-3*0.0625)/0.5,
                       b: float(b0expr.subs({x:0.25, g:(1-0.0625)/0.5, h:(1-3*0.0625)/0.5}).evalf(30))}) > 0 else -1
w6 = s6 * w6raw

M2 = ((g+h)/2, R(1,2))
s23 = cancel(cross(sub2(Q1, M2), sub2(Q0, Q1)) / cross(w0, sub2(Q0, Q1)))
ell = (g + h + x)/2 + s23

# ---- 오각형 점 (정확) ---- — the pentagon point (exact) ----
s5 = sqrt(5)
x0 = sqrt(1 - 2/s5)
g0 = (1 - x0**2)/(2*x0); h0 = (1 - 3*x0**2)/(2*x0)
w0_0x = ((x0**2-1)/(1+x0**2)); w0_0y = (2*x0/(1+x0**2))
b0 = cancel(b0expr.subs({x: x0, g: g0, h: h0}))
P0 = {x: x0, g: g0, h: h0, a: h0, b: b0}

names = ['w1','w4','w5','w6','w7','w8']
Ws = [w1, w4, w5, w6, w7, w8]
print("=== P1: 오각형에서 벽값 (정확, 전부 0이어야) ===")
for nm, wf in zip(names, Ws):
    v = simplify(nsimplify(wf.subs(P0), [sqrt(5)]))
    print(f"  {nm}(P*) = {v}")

print("\n=== P2: grad ell (정확) ===")
gl = []
for v in (x, g, h, a, b):
    e = simplify(nsimplify(diff(ell, v).subs(P0), [sqrt(5)]))
    gl.append(e)
    print(f"  d(ell)/d{v} = {e}")

print("\n=== P3: Farkas 증서 탐색 ===")
Gm = []
for wf in Ws:
    row = [simplify(nsimplify(diff(wf, v).subs(P0), [sqrt(5)])) for v in (x, g, h, a, b)]
    Gm.append(row)

from itertools import combinations
found = None
for idxs in combinations(range(6), 5):
    A = Matrix([[Gm[j][i] for j in idxs] for i in range(5)])   # 5x5: 열 = 선택된 벽 grad — 5x5: columns = gradients of the selected walls
    try:
        lam = A.solve(Matrix(gl))
    except Exception:
        continue
    lam_s = [simplify(v) for v in lam]
    vals = [complex(v.evalf(40)).real for v in lam_s]
    if all(v >= -1e-30 for v in vals):
        found = (idxs, lam_s, vals)
        print(f"  부분집합 {[names[j] for j in idxs]}:")
        for j, lv, vn in zip(idxs, lam_s, vals):
            print(f"    lambda[{names[j]}] = {lv}   (~ {vn:.6f})")
        # 잔차 정확 검증 — exact verification of the residual
        res = [simplify(sum(lam_s[k]*Gm[idxs[k]][i] for k in range(5)) - gl[i]) for i in range(5)]
        print(f"    잔차(정확) = {res}")
        break
    else:
        print(f"  부분집합 {[names[j] for j in idxs]}: 음수 승수 {[f'{v:.3f}' for v in vals]}")
if found:
    print("  ==> Farkas 증서 존재: grad ell ∈ cone(grad w) — 오각형 1계 차단 [정확 대수]")
else:
    print("  ==> 5-부분집합에서 비음수 해 없음 — 6벡터 결합 재검토 필요")

print("\n=== P4: 슬라이스 위 d(w_j)/dx 부호 (x의 함수, 정확) ===")
gS = (1-x**2)/(2*x); hS = (1-3*x**2)/(2*x)
w0S = refl_lin((cancel(gS-hS), S(1)))((S(1),S(0)))
bS = cancel(dot(sub2((-x,S(1)), (hS,S(0))), w0S))
SL = {g: gS, h: hS, a: hS, b: bS}
for nm, wf in zip(names, Ws):
    dwdx = cancel(together(diff(wf, x).subs(SL)))
    dwdx = cancel(simplify(dwdx))
    print(f"  d({nm})/dx |slice = {dwdx}")
    print(f"      값: x=0.15: {float(dwdx.subs(x,0.15)):+.4f}, x=0.25: {float(dwdx.subs(x,0.25)):+.4f}, x=0.32: {float(dwdx.subs(x,0.32)):+.4f}")
