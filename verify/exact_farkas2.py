"""
07 검증 2: taco 벽(선언적/비볼록)을 포함한 오각형 차단의 정확 판정.
07 verification 2: exact decision of the pentagon's blocking, including the taco wall (disjunctive/non-convex).

추가 벽 (크리스1 x 크리스3 선분 교차 금지 — taco-taco):
additional wall (crease1 x crease3 segment crossing is forbidden — taco-taco):
  교점 파라미터: E + t(F-E) = Q1 + u(Q0-Q1).
  intersection parameters: E + t(F-E) = Q1 + u(Q0-Q1).
  금지 = {t<1 그리고 u<1} (국소적으로, 오각형 근방).  실현 = {t>=1} ∪ {u>=1}  — 합집합!
  forbidden = {t<1 and u<1} (locally, near the pentagon).  realizable = {t>=1} ∪ {u>=1}  — a union!
  => 1계 차단 판정은 두 가지로 분해:
  => the first-order blocking test therefore splits into two branches:
     branch-9  : 벽 {w1,w4,w5,w6,w7,w8, w9=t-1}   에서 Farkas
     branch-9  : Farkas on the walls {w1,w4,w5,w6,w7,w8, w9=t-1}
     branch-10 : 벽 {w1,w4,w5,w6,w7,w8, w10=u-1}  에서 Farkas
     branch-10 : Farkas on the walls {w1,w4,w5,w6,w7,w8, w10=u-1}
  둘 다 증서가 있으면 오각형은 1계 차단 [정확 대수].
  if both admit a certificate, the pentagon is first-order blocked [exact algebra].

또한: 6벽 계의 탈출 방향 delta*를 명시적으로 구해 w9,w10을 위반하는지 확인
also: compute the escape direction delta* of the 6-wall system explicitly and check whether it violates w9, w10
(= 지난 턴 F2/F3가 '차단'으로 오판한 탈출로가 정확히 taco 벽으로 막힘을 보임).
(= shows that the escape route which F2/F3 misjudged as 'blocked' in the previous turn is in fact blocked precisely by the taco wall).
"""
from sympy import (symbols, sqrt, Rational, simplify, cancel, together, diff,
                   Matrix, S, nsimplify)
from itertools import combinations
import math

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
d2 = refl_lin(sub2(F, E))((S(1), S(0)))
w0 = refl_lin(sub2(G, H))((S(1), S(0)))
Q1 = (g + a*w0[0], 1 + a*w0[1])
Q0 = (h + b*w0[0], b*w0[1])
rl3 = refl_lin(sub2(Q0, Q1))
rho3 = lambda P: add2(Q1, rl3(sub2(P, Q1)))
dD = rl3(w0)
e1p = rho3(G)

x0n = math.sqrt(1 - 2/math.sqrt(5))
g0n = (1-x0n**2)/(2*x0n); h0n = (1-3*x0n**2)/(2*x0n)
b0expr = dot(sub2(F, H), w0)

def _num(e, X): return float(e.subs(X).evalf(30))
X0n = {x: x0n, g: g0n, h: h0n, a: h0n,
       b: float(b0expr.subs({x: x0n, g: g0n, h: h0n}).evalf(30))}
sA = 1 if _num(cross(sub2(F, E), d2), X0n) > 0 else -1
pD = add2(Q1, dD)
sD = 1 if _num(cross(sub2(pD, e1p), dD), X0n) > 0 else -1

w1 = -sA * cross(sub2(G, E), d2)
w7 =  sA * cross(sub2(Q0, E), d2) - 1
w4 = -sD * cross(sub2(H, e1p), dD)
w8 =  sD * cross(sub2(F, e1p), dD) - 1
sI = cancel(cross(sub2(F, G), d2) / cross(w0, d2))
w5 = a - sI
Xl = {x:0.25, g:(1-0.0625)/0.5, h:(1-3*0.0625)/0.5, a:(1-3*0.0625)/0.5,
      b: float(b0expr.subs({x:0.25, g:(1-0.0625)/0.5, h:(1-3*0.0625)/0.5}).evalf(30))}
w6raw = cross(sub2(E, rho3(H)), dD)
s6 = 1 if _num(w6raw, Xl) > 0 else -1
w6 = s6 * w6raw

# taco 벽: E + t(F-E) = Q1 + u(Q0-Q1) — taco wall: E + t(F-E) = Q1 + u(Q0-Q1)
den = cross(sub2(F, E), sub2(Q0, Q1))
t9  = cancel(cross(sub2(Q1, E), sub2(Q0, Q1)) / den)
u10 = cancel(cross(sub2(Q1, E), sub2(F, E)) / den)
w9  = t9 - 1
w10 = u10 - 1

M2 = ((g+h)/2, R(1,2))
s23 = cancel(cross(sub2(Q1, M2), sub2(Q0, Q1)) / cross(w0, sub2(Q0, Q1)))
ell = (g + h + x)/2 + s23

s5 = sqrt(5)
x0 = sqrt(1 - 2/s5)
g0 = (1 - x0**2)/(2*x0); h0 = (1 - 3*x0**2)/(2*x0)
b0 = cancel(b0expr.subs({x: x0, g: g0, h: h0}))
P0 = {x: x0, g: g0, h: h0, a: h0, b: b0}

names = ['w1','w4','w5','w6','w7','w8','w9','w10']
Ws = [w1, w4, w5, w6, w7, w8, w9, w10]
print("=== 오각형에서 벽값 (w9, w10 포함; 전부 0이어야) ===")
for nm, wf in zip(names, Ws):
    print(f"  {nm}(P*) = {simplify(nsimplify(wf.subs(P0), [sqrt(5)]))}")

print("\n기울기 계산 중 (정확)...")
gl = [simplify(nsimplify(diff(ell, v).subs(P0), [sqrt(5)])) for v in (x,g,h,a,b)]
Gm = {nm: [simplify(nsimplify(diff(wf, v).subs(P0), [sqrt(5)])) for v in (x,g,h,a,b)]
      for nm, wf in zip(names, Ws)}

def farkas(wall_names):
    """벽 집합에 대해 grad ell = sum lambda grad w, lambda>=0 정확해 탐색. — For the given wall set, search for an exact solution of grad ell = sum lambda grad w with lambda>=0."""
    for idxs in combinations(wall_names, 5):
        A = Matrix([[Gm[nm][i] for nm in idxs] for i in range(5)])
        try:
            lam = A.solve(Matrix(gl))
        except Exception:
            continue
        lam_s = [simplify(v) for v in lam]
        vals = [complex(v.evalf(40)).real for v in lam_s]
        if all(v >= -1e-30 for v in vals):
            res = [simplify(sum(lam_s[k]*Gm[idxs[k]][i] for k in range(5)) - gl[i])
                   for i in range(5)]
            return idxs, lam_s, vals, res
    return None

base6 = ['w1','w4','w5','w6','w7','w8']
for extra in ('w9', 'w10'):
    print(f"\n=== branch-{extra}: 벽 {base6 + [extra]} 에서 Farkas ===")
    out = farkas(base6 + [extra])
    if out:
        idxs, lam_s, vals, res = out
        print(f"  증서 발견! 부분집합 {list(idxs)}")
        for nm, lv, vn in zip(idxs, lam_s, vals):
            print(f"    lambda[{nm}] = {lv}   (~ {vn:.6f})")
        print(f"  잔차(정확) = {res}")
    else:
        print("  비음수 해 없음 — 이 가지에서 1계 탈출 가능성 잔존")

# ==== 6벽 계의 탈출 방향과 w9/w10 위반 확인 ==== — escape direction of the 6-wall system and check for w9/w10 violation ====
print("\n=== 6벽 탈출 방향 delta* 와 taco 벽의 역할 ===")
found_dir = None
for eqs in combinations(base6, 4):
    A = Matrix([[Gm[nm][i] for i in range(5)] for nm in eqs])
    ker = A.nullspace()
    if len(ker) != 1: continue
    for sgn in (1, -1):
        d = [sgn*simplify(v) for v in ker[0]]
        dn = [complex(v.evalf(40)).real for v in d]
        gln = [complex(v.evalf(40)).real for v in gl]
        if sum(gln[i]*dn[i] for i in range(5)) >= -1e-20: continue
        ok = True
        for nm in base6:
            gv = sum(complex(Gm[nm][i].evalf(40)).real*dn[i] for i in range(5))
            if gv < -1e-18: ok = False; break
        if ok:
            found_dir = (eqs, d, dn); break
    if found_dir: break
if found_dir:
    eqs, d, dn = found_dir
    nd = math.sqrt(sum(v*v for v in dn)); dn = [v/nd for v in dn]
    print(f"  delta* (등식 벽 {list(eqs)}): {[f'{v:+.4f}' for v in dn]}")
    for nm in ('w9','w10'):
        gv = sum(complex(Gm[nm][i].evalf(40)).real*dn[i] for i in range(5))
        print(f"  grad {nm} · delta* = {gv:+.6f}   {'(위반 방향!)' if gv < 0 else ''}")
    print("  => 6벽 탈출로가 taco 벽(w9, w10)을 침범하는지 확인")
else:
    print("  6벽 탈출 방향을 찾지 못함(수치 임계) — P3 결과와 대조 필요")
