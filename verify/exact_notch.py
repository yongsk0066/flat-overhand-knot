# exact_notch.py — 노치 삼량(입구·다리·깊이)과 반사각의 Q(x†) 정확화 (body.md TODO 해소) — Q(x†)-exactification of the notch triple (mouth·leg·depth) and of the reflection angle (resolves the body.md TODO)
#
# 명제 9.14의 수치 삼량을 몫체 Q[x]/(p_x)의 닫힌 원소로 승격: — promote the numerical triple of Prop. 9.14 to closed elements of the quotient field Q[x]/(p_x):
#   X = Q0,  Y = e_B ∩ (D-밴드 먼 모서리 e0' = ρ3(H)-변),  F = (−x, 1) — X = Q0,  Y = e_B ∩ (far edge e0' of the D-band = the ρ3(H)-edge),  F = (−x, 1)
#   mouth² = |XF|²,  leg² = |XY|² (= |YF|² — 미러 항등 검증),  depth² = cross²/mouth², — mouth² = |XF|²,  leg² = |XY|² (= |YF|² — verifies the mirror identity),  depth² = cross²/mouth²,
#   cos∠Y = (X−Y)·(F−Y)/leg²  (|XY|=|YF| 덕분에 유리!) — cos∠Y = (X−Y)·(F−Y)/leg²  (rational thanks to |XY|=|YF|!)
# 전부 유리 연산 ⟹ 각 양이 차수 ≤ 11 다항식으로 표현된다. 구간 수치로 기존 값 대조. — every operation is rational ⟹ each quantity is expressed by a polynomial of degree ≤ 11. The previously known values are cross-checked by interval numerics.

import sympy as sp
from mpmath import mp, mpf, findroot, nstr

x = sp.symbols('x')
px = sp.expand(x**12 + 4*x**10 - 7*x**4 + 20*x**2 - 2)

def red(e):
    return sp.rem(sp.expand(e), px, x)

def inv(e):
    return sp.invert(red(e), px, x)

def fmul(a, b): return red(a*b)
def fdiv(a, b): return red(a*inv(b))

# ---- θ† 의 체 표현 ---- — field representation of θ† ----
h_f = sp.Rational(1, 12)*(2*x**11 + 8*x**9 - x**7 - 4*x**5 - 13*x**3 + 40*x)
g_f = fdiv(1 - x**2, 2*x)
a_f = h_f
b_f = fdiv(x**2 + 1, 2*x)
q_f = red(g_f - h_f)
nW0 = red(q_f**2 + 1)

def vec(px_, py_): return (red(px_), red(py_))
def vsub(u, v): return (red(u[0] - v[0]), red(u[1] - v[1]))
def vdot(u, v): return red(u[0]*v[0] + u[1]*v[1])
def vcross(u, v): return red(u[0]*v[1] - u[1]*v[0])

E = vec(sp.Integer(0), sp.Integer(0))
F = vec(-x, sp.Integer(1))
G = vec(g_f, sp.Integer(1))
H = vec(h_f, sp.Integer(0))
W0 = vec(red(q_f**2 - 1), red(2*q_f))          # |W0| = nW0
D2 = vec(x**2 - 1, -2*x)                        # |D2| = x²+1
W0h = (fdiv(W0[0], nW0), fdiv(W0[1], nW0))      # 단위 — unit
Q1 = (red(G[0] + a_f*W0h[0]), red(G[1] + a_f*W0h[1]))
Q0 = (red(H[0] + b_f*W0h[0]), red(H[1] + b_f*W0h[1]))
U = vsub(Q0, Q1)
UU = vdot(U, U)

def refl_pt(P, d, Qp):
    """직선(P, 방향 d)에 대한 점 반사 — 유리. — Point reflection about the line (P, direction d) — rational."""
    v = vsub(Qp, P)
    t = fdiv(red(2*vdot(v, d)), vdot(d, d))
    return (red(P[0] + t*d[0] - v[0]), red(P[1] + t*d[1] - v[1]))

def refl_lin(d, v):
    t = fdiv(red(2*vdot(v, d)), vdot(d, d))
    return (red(t*d[0] - v[0]), red(t*d[1] - v[1]))

e0p = refl_pt(Q1, U, H)          # ρ3(H): D-밴드 먼 모서리의 기준점 — ρ3(H): the reference point of the far edge of the D-band
dD = refl_lin(U, W0)             # D 방향 (비정규 — 직선용) — the D direction (unnormalized — used only for the line)

# Y = line(F, D2) ∩ line(e0p, dD):  F + t·D2 = e0p + s·dD
den = vcross(D2, dD)
t_Y = fdiv(vcross(vsub(e0p, F), dD), den)
Y = (red(F[0] + t_Y*D2[0]), red(F[1] + t_Y*D2[1]))
X = Q0

# ---- 삼량 ---- — the three quantities (mouth, leg, depth) ----
XF = vsub(X, F)
mouth2 = vdot(XF, XF)
XY = vsub(X, Y); YF = vsub(F, Y)
leg2_1 = vdot(XY, XY)
leg2_2 = vdot(YF, YF)
mirror_id = red(leg2_1 - leg2_2)
crossA = vcross(XF, vsub(Y, F))          # 부호 있는 2·면적 — signed 2·area
depth2 = fdiv(red(crossA**2), mouth2)
cosY = fdiv(vdot(XY, YF), leg2_1)        # |XY|=|YF| ⟹ 유리! — |XY|=|YF| ⟹ rational!

print("[정확 항등] |XY|² − |YF|² ≡", mirror_id, " (0 = 미러 대칭의 대수적 확인)")

print("\n[닫힌 형태 (기저 1..x^11, 유리 계수)]")
for name, e in [("mouth²", mouth2), ("leg²", leg2_1), ("depth²", depth2), ("cos∠Y", cosY)]:
    print(f"  {name} = {sp.nsimplify(e)}")

# ---- 수치 대조 (40자리) ---- — numerical cross-check (40 digits) ----
mp.dps = 50
xr = findroot(lambda t: t**12 + 4*t**10 - 7*t**4 + 20*t**2 - 2,
              mpf('0.3221274973112271792397142057801235136978'))
def ev(e):
    cs = sp.Poly(e, x).all_coeffs()
    v = mpf(0)
    for c in cs:
        v = v*xr + mpf(c.p)/mpf(c.q)
    return v

import math
m2, l2, d2v, cy = ev(mouth2), ev(leg2_1), ev(depth2), ev(cosY)
from mpmath import sqrt as msqrt, acos, degrees
print("\n[수치 대조 (기존 명제 9.14 값)]")
print(f"  mouth = {nstr(msqrt(m2), 12)}   (기존 0.084488)")
print(f"  leg   = {nstr(msqrt(l2), 12)}   (기존 0.051419)")
print(f"  depth = {nstr(msqrt(d2v), 12)}   (기존 0.029349…)")
ang = degrees(acos(cy))
print(f"  ∠XYF  = {nstr(ang, 10)}°  ⇒ 반사 내각 = {nstr(360 - ang, 10)}°   (기존 249.517°)")
assert abs(msqrt(m2) - mpf('0.084488')) < 1e-5
assert abs(msqrt(l2) - mpf('0.051419')) < 1e-5
assert abs(360 - ang - mpf('249.517')) < 1e-2
print("\n결론: 노치 삼량·반사각 전부 Q(x†) 닫힌 원소로 승격 [정확대수] — 명제 9.14의 양들.")
