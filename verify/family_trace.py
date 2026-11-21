"""Trace the second symmetric family at high precision and identify candidate algebraic coordinates.

구조 (모든 양이 theta의 유리함수 — 반사가 노름 보존이라 단위화 불필요):
structure (every quantity is a rational function of theta — reflections preserve the norm, so no normalization is needed):
  g = (1-x^2)/(2x)                    [w1=0 이 강제 — 표준 슬라이스와 동일 식]
  g = (1-x^2)/(2x)                    [forced by w1=0 — the same expression as on the standard slice]
  a = sI(x,h)                         [w5=0: 크리스3 l2-끝점 = 코너 I]
  a = sI(x,h)                         [w5=0: the l2-endpoint of crease3 = the corner I]
  (h,b) = Newton 해 of (w4, w6) = 0   [나머지 두 벽]
  (h,b) = Newton solution of (w4, w6) = 0   [the remaining two walls]
  => x 하나로 매개되는 족. ell(x) 최소화 -> x†.
  => a family parameterized by x alone. Minimizing ell(x) -> x†.

그 후 mpmath PSLQ(findpoly)로 x†, ell† 의 최소다항식 후보 탐색.
Afterwards, mpmath PSLQ (findpoly) searches for minimal-polynomial candidates for x† and ell†.
"""
from mpmath import mp, mpf, findpoly, nstr, sqrt as msqrt

mp.dps = 60

def RL(d):
    dx, dy = d
    n = dx*dx + dy*dy
    return lambda v: (((dx*dx-dy*dy)*v[0] + 2*dx*dy*v[1])/n,
                      (2*dx*dy*v[0] + (dy*dy-dx*dx)*v[1])/n)

def cross(u,v): return u[0]*v[1]-u[1]*v[0]
def dot(u,v):  return u[0]*v[0]+u[1]*v[1]
def sub(u,v):  return (u[0]-v[0], u[1]-v[1])
def add(u,v):  return (u[0]+v[0], u[1]+v[1])
def mul(u,s):  return (u[0]*s, u[1]*s)

def state(x, h, b):
    g = (1 - x*x)/(2*x)
    E = (mpf(0), mpf(0)); F = (-x, mpf(1)); G = (g, mpf(1)); H = (h, mpf(0))
    d2 = RL(sub(F,E))((mpf(1), mpf(0)))
    w0 = RL(sub(G,H))((mpf(1), mpf(0)))
    a  = cross(sub(F,G), d2)/cross(w0, d2)          # w5=0
    Q1 = add(G, mul(w0, a)); Q0 = add(H, mul(w0, b))
    rl3 = RL(sub(Q0,Q1))
    rho3 = lambda P: add(Q1, rl3(sub(P,Q1)))
    dD = rl3(w0)
    e1p = rho3(G)
    w4 = cross(sub(H, e1p), dD)
    w6 = cross(sub(E, rho3(H)), dD)
    M2 = ((g+h)/2, mpf('0.5'))
    s23 = cross(sub(Q1,M2), sub(Q0,Q1))/cross(w0, sub(Q0,Q1))
    ell = (g + h + x)/2 + s23
    return w4, w6, ell, g, a

def solve_hb(x, h0, b0, tol=mpf('1e-45')):
    h, b = mpf(h0), mpf(b0)
    for _ in range(60):
        w4, w6, _, _, _ = state(x, h, b)
        if abs(w4) < tol and abs(w6) < tol: break
        eps = mpf('1e-25')
        w4h, w6h, _,_,_ = state(x, h+eps, b)
        w4b, w6b, _,_,_ = state(x, h, b+eps)
        J = [[(w4h-w4)/eps, (w4b-w4)/eps],
             [(w6h-w6)/eps, (w6b-w6)/eps]]
        det = J[0][0]*J[1][1] - J[0][1]*J[1][0]
        dh = ( J[1][1]*w4 - J[0][1]*w6)/det
        db = (-J[1][0]*w4 + J[0][0]*w6)/det
        h, b = h - dh, b - db
    return h, b

def ell_of_x(x, hb_guess):
    h, b = solve_hb(x, *hb_guess)
    _, _, e, g, a = state(x, h, b)
    return e, (h, b), g, a

# ---- 족 추적: 오각형(x0)에서 아래로 ---- — family tracking: downward from the pentagon (x0) ----
x0 = msqrt(1 - 2/msqrt(5))
print(f"x0 (오각형) = {nstr(x0, 30)}")
guess = (mpf('1.0367'), mpf('1.7129'))
print("\n x         ell(x)                      h")
xs = [x0 - mpf(k)/1000 for k in (0, 1, 2, 3, 4, 5)]
vals = []
for xv in xs:
    e, hb, g, a = ell_of_x(xv, guess)
    guess = hb
    vals.append((xv, e))
    print(f" {nstr(xv,8)}  {nstr(e,20)}   {nstr(hb[0],12)}")

# ---- ell'(x) = 0 뉴턴 (시컨트) ---- — Newton (secant) for ell'(x) = 0 ----
def dell(x, guess):
    eps = mpf('1e-20')
    e1,_,_,_ = ell_of_x(x+eps, guess)
    e2,_,_,_ = ell_of_x(x-eps, guess)
    return (e1-e2)/(2*eps)

xa, xb_ = x0 - mpf('0.001'), x0 - mpf('0.004')
fa, fb = dell(xa, guess), dell(xb_, guess)
for _ in range(80):
    xm = xb_ - fb*(xb_-xa)/(fb-fa)
    xa, fa = xb_, fb
    xb_ = xm
    fb = dell(xm, guess)
    if abs(fb) < mpf('1e-40'): break
xd = xb_
ed, hbd, gd, ad = ell_of_x(xd, guess)
hd, bd = hbd
print(f"\n=== theta† (고정밀) ===")
print(f" x† = {nstr(xd, 40)}")
print(f" h† = {nstr(hd, 40)}")
print(f" b† = {nstr(bd, 40)}")
print(f" g† = {nstr(gd, 40)}")
print(f" a† = {nstr(ad, 40)}")
print(f" ell† = {nstr(ed, 40)}")
print(f" 오각형 ell = {nstr(1/x0 - x0, 40)}")
print(f" 차이 = {nstr(ed - (1/x0-x0), 20)}")

# ---- PSLQ: 최소다항식 후보 ---- — PSLQ: minimal-polynomial candidates ----
print("\n=== PSLQ 최소다항식 탐색 (계수 한도 10^8) ===")
mp.dps = 50
for name, val in (("x†", xd), ("ell†", ed), ("h†", hd), ("x†^2", xd*xd)):
    found = None
    for deg in range(2, 9):
        p = findpoly(val, deg, maxcoeff=10**8, tol=mpf(10)**(-40))
        if p:
            found = (deg, p); break
    print(f"  {name}: " + (f"차수 {found[0]}, 계수 {found[1]}" if found else "차수<=8, 계수<=1e8 에서 없음"))
