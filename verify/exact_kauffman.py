"""Kauffman 절단 길이 L_trunc의 오각형 1계 차단 — 정확 대수 판정.
First-order blocking of the Kauffman truncated length L_trunc at the pentagon — exact algebraic decision.

L_trunc = s_in + ell + t_out  (전부 theta의 유리함수):
L_trunc = s_in + ell + t_out  (all of them rational functions of theta):
  s_in : M1 + s(-d2) 와 직선 (G, Q1-G) [l2]의 교점 파라미터
  s_in : the intersection parameter of M1 + s(-d2) with the line (G, Q1-G) [l2]
  t_out: M3 + t v3   와 직선 y=0 의 교점 파라미터
  t_out: the intersection parameter of M3 + t v3 with the line y=0
판정: 오각형 P0에서 grad L_trunc ∈ cone{grad w_j} (8벽, taco 두 가지 각각)
decision: at the pentagon P0, grad L_trunc ∈ cone{grad w_j} (8 walls, for each of the two taco branches)
      — Caratheodory: 7벽 중 5-부분집합 전부 검사.
      — Caratheodory: every 5-subset of the 7 walls is checked.
"""
from sympy import symbols, sqrt, Rational, simplify, cancel, diff, Matrix, S, nsimplify
from itertools import combinations
import math

x, g, h, a, b = symbols('x g h a b', real=True)
R = Rational

def refl_lin(d):
    dx, dy = d; n = dx**2 + dy**2
    return lambda v: (cancel(((dx**2-dy**2)*v[0] + 2*dx*dy*v[1])/n),
                      cancel((2*dx*dy*v[0] + (dy**2-dx**2)*v[1])/n))
def cross(u,v): return u[0]*v[1]-u[1]*v[0]
def dot(u,v): return u[0]*v[0]+u[1]*v[1]
def sub2(u,v): return (u[0]-v[0], u[1]-v[1])
def add2(u,v): return (u[0]+v[0], u[1]+v[1])

E=(S(0),S(0)); F=(-x,S(1)); G=(g,S(1)); H=(h,S(0))
d2 = refl_lin(sub2(F,E))((S(1),S(0)))
w0 = refl_lin(sub2(G,H))((S(1),S(0)))
Q1 = (g + a*w0[0], 1 + a*w0[1]); Q0 = (h + b*w0[0], b*w0[1])
rl3 = refl_lin(sub2(Q0,Q1))
rho3 = lambda P: add2(Q1, rl3(sub2(P,Q1)))
dD = rl3(w0); e1p = rho3(G)
v3 = dD

M1 = (-x/2, R(1,2)); M2 = ((g+h)/2, R(1,2))
s23 = cancel(cross(sub2(Q1,M2), sub2(Q0,Q1)) / cross(w0, sub2(Q0,Q1)))
M3 = (M2[0] + s23*w0[0], M2[1] + s23*w0[1])
ell = (g+h+x)/2 + s23

# s_in: M1 + s*(-d2) = G + u*(Q1-G)
dQ = sub2(Q1, G)
s_in = cancel(cross(sub2(G, M1), dQ) / cross((-d2[0], -d2[1]), dQ))
# t_out: (M3 + t*v3).y = 0  ->  t = -M3y / v3y
t_out = cancel(-M3[1] / v3[1])
Ltr = s_in + ell + t_out

# 벽 (exact_farkas2와 동일) — walls (same as in exact_farkas2)
import sympy
def _num(e, X): return float(e.subs(X).evalf(30))
x0n = math.sqrt(1-2/math.sqrt(5)); g0n=(1-x0n**2)/(2*x0n); h0n=(1-3*x0n**2)/(2*x0n)
b0expr = dot(sub2(F,H), w0)
X0n = {x:x0n, g:g0n, h:h0n, a:h0n, b: float(b0expr.subs({x:x0n,g:g0n,h:h0n}).evalf(30))}
sA = 1 if _num(cross(sub2(F,E), d2), X0n) > 0 else -1
pD = (M3[0]+v3[0], M3[1]+v3[1])
sD = 1 if _num(cross(sub2(pD,e1p), dD), X0n) > 0 else -1
w1 = -sA*cross(sub2(G,E), d2)
w7 =  sA*cross(sub2(Q0,E), d2) - 1
w4 = -sD*cross(sub2(H,e1p), dD)
w8 =  sD*cross(sub2(F,e1p), dD) - 1
sI = cancel(cross(sub2(F,G), d2)/cross(w0, d2))
w5 = a - sI
Xl = {x:0.25, g:(1-0.0625)/0.5, h:(1-3*0.0625)/0.5, a:(1-3*0.0625)/0.5,
      b: float(b0expr.subs({x:0.25,g:(1-0.0625)/0.5,h:(1-3*0.0625)/0.5}).evalf(30))}
w6raw = cross(sub2(E,rho3(H)), dD)
s6 = 1 if _num(w6raw, Xl) > 0 else -1
w6 = s6*w6raw
den = cross(sub2(F,E), sub2(Q0,Q1))
w9  = cancel(cross(sub2(Q1,E), sub2(Q0,Q1))/den) - 1
w10 = cancel(cross(sub2(Q1,E), sub2(F,E))/den) - 1

s5 = sqrt(5); x0 = sqrt(1 - 2/s5)
g0 = (1-x0**2)/(2*x0); h0 = (1-3*x0**2)/(2*x0)
b0 = cancel(b0expr.subs({x:x0,g:g0,h:h0}))
P0 = {x:x0, g:g0, h:h0, a:h0, b:b0}

print("grad L_trunc (정확):")
gl = []
for v in (x,g,h,a,b):
    e = simplify(nsimplify(diff(Ltr, v).subs(P0), [sqrt(5)]))
    gl.append(e); print(f"  d/d{v} = {e}")

names = ['w1','w4','w5','w6','w7','w8','w9','w10']
Ws = [w1,w4,w5,w6,w7,w8,w9,w10]
print("\n벽 기울기 계산 중 (정확)...")
Gm = {nm: [simplify(nsimplify(diff(wf, vv).subs(P0), [sqrt(5)])) for vv in (x,g,h,a,b)]
      for nm, wf in zip(names, Ws)}

def farkas(wall_names):
    for idxs in combinations(wall_names, 5):
        A = Matrix([[Gm[nm][i] for nm in idxs] for i in range(5)])
        try: lam = A.solve(Matrix(gl))
        except Exception: continue
        lam_s = [simplify(v) for v in lam]
        vals = [complex(v.evalf(40)).real for v in lam_s]
        if all(v >= -1e-30 for v in vals):
            res = [simplify(sum(lam_s[k]*Gm[idxs[k]][i] for k in range(5)) - gl[i]) for i in range(5)]
            if all(r == 0 for r in res):
                return idxs, lam_s, vals
    return None

base6 = ['w1','w4','w5','w6','w7','w8']
allok = True
for extra in ('w9','w10'):
    out = farkas(base6+[extra])
    print(f"\nbranch-{extra}: " + ("증서 발견!" if out else "비음수 해 없음 — 차단 실패"))
    if out:
        idxs, lam_s, vals = out
        for nm, lv, vn in zip(idxs, lam_s, vals):
            print(f"  lambda[{nm}] = {lv}  (~{vn:.6f})")
    else:
        allok = False
print("\n====", "정리: 오각형은 L_trunc의 1계 국소 최소 [정확대수]" if allok else "차단 실패 — 추가 분석 필요", "====")
