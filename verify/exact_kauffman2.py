"""L_trunc 오각형 차단의 정확 판정 (2단: 수치 식별 -> 성분별 정확 검증 -> Farkas).
Exact decision of the blocking of L_trunc at the pentagon (2-stage: numerical identification -> component-wise exact verification -> Farkas).
전략: 전체 simplify 대신, 50자리 수치로 Q(sqrt5) 후보를 찾고
strategy: instead of a global simplify, locate the Q(sqrt5) candidates from 50-digit numerics and
      각 후보만 등식 검증(차이의 유리화)한다.
      verify the equality for each candidate only (by rationalizing the difference)."""
from sympy import (symbols, sqrt, Rational, simplify, cancel, diff, Matrix, S,
                   nsimplify, together, fraction, expand)
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
dD = rl3(w0); e1p = rho3(G); v3 = dD
M1 = (-x/2, R(1,2)); M2 = ((g+h)/2, R(1,2))
s23 = cancel(cross(sub2(Q1,M2), sub2(Q0,Q1)) / cross(w0, sub2(Q0,Q1)))
M3 = (M2[0] + s23*w0[0], M2[1] + s23*w0[1])
ell = (g+h+x)/2 + s23
dQ = sub2(Q1, G)
s_in = cancel(cross(sub2(G, M1), dQ) / cross((-d2[0], -d2[1]), dQ))
t_out = cancel(-M3[1] / v3[1])
Ltr = s_in + ell + t_out

s5 = sqrt(5); x0 = sqrt(1 - 2/s5)
g0 = (1-x0**2)/(2*x0); h0 = (1-3*x0**2)/(2*x0)
b0expr = dot(sub2(F,H), w0)
b0 = cancel(b0expr.subs({x:x0, g:g0, h:h0}))
P0 = {x:x0, g:g0, h:h0, a:h0, b:b0}

print("=== 1단: grad L_trunc 수치(50자리) -> Q(sqrt5) 후보 식별 ===")
import mpmath
mpmath.mp.dps = 60
S5 = mpmath.sqrt(5)
X0m = mpmath.sqrt(1 - 2/S5)
BASIS_M = [mpmath.mpf(1), S5, X0m, S5*X0m]
from sympy import Integer
x0sym = sqrt(1 - 2/sqrt(5))
BASIS_S = [Integer(1), sqrt(5), x0sym, sqrt(5)*x0sym]
def identify(val):
    """val ∈ Q(x0) = span_Q{1, sqrt5, x0, sqrt5*x0} 를 PSLQ로 정확 식별. — Identify val ∈ Q(x0) = span_Q{1, sqrt5, x0, sqrt5*x0} exactly, via PSLQ."""
    v = mpmath.mpf(val)
    if abs(v) < mpmath.mpf('1e-45'): return S(0)
    rel = mpmath.pslq(BASIS_M + [v], maxcoeff=10**10, maxsteps=10**6)
    if rel is None or rel[4] == 0: raise ValueError(f"식별 실패: {val}")
    c = rel
    return sum(Rational(-c[i], c[4])*BASIS_S[i] for i in range(4))
cands = []
for v in (x, g, h, a, b):
    d_expr = diff(Ltr, v)
    num = mpmath.mpf(str(d_expr.subs(P0).evalf(52)))
    cand = identify(num)
    cands.append(cand)
    print(f"  d/d{v}: {float(num):+.12f}  ->  후보 {cand}")

print("\n=== 2단: 성분별 정확 검증 (차이 == 0) ===")
gl = []
okall = True
for v, cand in zip((x,g,h,a,b), cands):
    dif = diff(Ltr, v).subs(P0) - cand
    z = simplify(nsimplify(dif, [sqrt(5)]))
    ok = (z == 0)
    okall &= ok
    gl.append(cand)
    print(f"  d/d{v} = {cand}: {'[정확 확인]' if ok else f'[불일치! 잔차={z}]'}")
assert okall

print("\n=== 3단: 8벽(선언적) Farkas — 정확 ===")
def _num(e, X): return float(e.subs(X).evalf(30))
x0n = math.sqrt(1-2/math.sqrt(5)); g0n=(1-x0n**2)/(2*x0n); h0n=(1-3*x0n**2)/(2*x0n)
X0n = {x:x0n, g:g0n, h:h0n, a:h0n, b: float(b0expr.subs({x:x0n,g:g0n,h:h0n}).evalf(30))}
sA = 1 if _num(cross(sub2(F,E), d2), X0n) > 0 else -1
pDp = (M3[0]+v3[0], M3[1]+v3[1])
sD = 1 if _num(cross(sub2(pDp,e1p), dD), X0n) > 0 else -1
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

names = ['w1','w4','w5','w6','w7','w8','w9','w10']
Ws = [w1,w4,w5,w6,w7,w8,w9,w10]
print("  벽 기울기: 수치 식별 -> 정확 검증 방식")
Gm = {}
for nm, wf in zip(names, Ws):
    row = []
    for v in (x,g,h,a,b):
        d_expr = diff(wf, v)
        num = mpmath.mpf(str(d_expr.subs(P0).evalf(52)))
        cand = identify(num)
        z = simplify(nsimplify(d_expr.subs(P0) - cand, [sqrt(5)]))
        assert z == 0, f"{nm} d/d{v} 식별 실패: {z}"
        row.append(cand)
    Gm[nm] = row
    print(f"    {nm}: {[str(c) for c in row]}")

def farkas(wall_names):
    for idxs in combinations(wall_names, 5):
        A = Matrix([[Gm[nm][i] for nm in idxs] for i in range(5)])
        try: lam = A.solve(Matrix(gl))
        except Exception: continue
        lam_s = [simplify(v2) for v2 in lam]
        vals = [complex(v2.evalf(40)).real for v2 in lam_s]
        if all(v2 >= -1e-30 for v2 in vals):
            res = [simplify(sum(lam_s[k]*Gm[idxs[k]][i] for k in range(5)) - gl[i]) for i in range(5)]
            if all(r == 0 for r in res):
                return idxs, lam_s, vals
    return None

base6 = ['w1','w4','w5','w6','w7','w8']
allok = True
for extra in ('w9','w10'):
    out = farkas(base6+[extra])
    print(f"\n  branch-{extra}: " + ("증서 발견!" if out else "비음수 해 없음"))
    if out:
        idxs, lam_s, vals = out
        for nm, lv, vn in zip(idxs, lam_s, vals):
            if abs(vn) > 1e-25:
                print(f"    lambda[{nm}] = {lv}  (~{vn:.6f})")
    else:
        allok = False
print("\n====", "정리 확립: 오각형은 L_trunc의 1계 국소 최소 [정확대수]" if allok else "차단 실패", "====")
