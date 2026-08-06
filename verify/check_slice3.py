"""
3차 검증: 확정된 기하 사전(dictionary)의 최종 확인.
Third verification: final confirmation of the fixed geometric dictionary.

확정 대상:
Items being pinned down:
  P1. 크리스2 = 직선 {y=1} (pass-1 위 모서리에 밀착): rho_{y=1}(edgeA) = (G, r).
  P1. crease2 = the line {y=1} (flush against the upper edge of pass-1): rho_{y=1}(edgeA) = (G, r).
  P2. edgeB ∩ {y=1} = F  →  크리스2의 현 = 선분 FG.
  P2. edgeB ∩ {y=1} = F  →  the chord of crease2 = segment FG.
  P3. pass-3 = 두 평행선 (G,r), (F,r) 사이의 폭 1 띠이고, H = (F,r) ∩ {y=0}.
  P3. pass-3 = the width-1 band between the two parallel lines (G,r) and (F,r), and H = (F,r) ∩ {y=0}.
  P4. 타이트(x0)에서 크리스3 방향이 d2와 평행 (pass-2 위 모서리 edgeB에 밀착),
      느슨한 x에서는 비평행.  [Delta=0 <=> 밀착의 수치 확인]
  P4. at the tight state (x0) the direction of crease3 is parallel to d2 (flush against edgeB, the
      upper edge of pass-2), while for loose x it is non-parallel.  [numerical check of Delta=0 <=> flush]
  P5. 타이트에서 오각형 다섯 변의 방향이 모두 36도의 배수.
  P5. at the tight state, the directions of the pentagon's five edges are all multiples of 36 degrees.
"""
import math

TOL = 1e-9

def sub(a, b): return (a[0]-b[0], a[1]-b[1])
def add(a, b): return (a[0]+b[0], a[1]+b[1])
def mul(a, s): return (a[0]*s, a[1]*s)
def dot(a, b): return a[0]*b[0] + a[1]*b[1]
def cross(a, b): return a[0]*b[1] - a[1]*b[0]
def norm(a): return math.hypot(a[0], a[1])
def unit(a):
    n = norm(a); return (a[0]/n, a[1]/n)

def refl(P, d):
    dh = unit(d)
    def rho(Q):
        v = sub(Q, P)
        w = sub(mul(dh, 2*dot(v, dh)), v)
        return add(P, w)
    return rho

def sdist(Q, P, d):
    dh = unit(d); v = sub(Q, P)
    return cross(v, dh) * -1.0

def check(label, val, tol=TOL):
    ok = abs(val) < tol
    print(f"  [{'OK' if ok else 'FAIL'}] {label}: {val:.3e}")
    return ok

def run(x, tight=False):
    print(f"\n=== x = {x:.7f} {'(tight)' if tight else ''} ===")
    ok = True
    E = (0.0, 0.0); F = (-x, 1.0)
    G = ((1 - x*x) / (2*x), 1.0)
    H = ((1 - 3*x*x) / (2*x), 0.0)
    r = ((x*x - 1) / (1 + x*x), 2*x / (1 + x*x))
    k = (1 - 3*x*x) / (2*x)
    K = add(G, mul(r, k))

    rho1 = refl(E, sub(F, E))
    o1 = rho1((0.0, 0.0))
    d2 = sub(rho1((1.0, 0.0)), o1)

    # P1: rho_{y=1} 가 edgeA=(E,d2)를 (G,r)로 보내는가 — does rho_{y=1} send edgeA=(E,d2) to (G,r)?
    rho2 = lambda Q: (Q[0], 2.0 - Q[1])
    for t in (-1.3, 0.7, 2.9):
        pt = rho2(add(E, mul(d2, t)))
        ok &= check(f"P1: rho_y=1(edgeA점 t={t})가 (G,r) 위", cross(sub(pt, G), r))

    # P2: edgeB=(F,d2) ∩ {y=1} = F  (d2가 수평이 아니므로 유일) — unique since d2 is not horizontal
    ok &= check("P2: F가 y=1 위 (자명)", F[1] - 1.0)
    ok &= check("P2: d2가 수평이 아님 (역수 유한)", 0.0 if abs(d2[1]) > 1e-12 else 1.0)

    # P3: pass-3 폭 = |sdist(F, G, r)| = 1, H = (F,r) ∩ {y=0} — 폭 = width
    ok &= check("P3a: 평행선 (G,r),(F,r) 사이 거리 = 1", abs(sdist(F, G, r)) - 1.0)
    tH = -F[1] / r[1]
    Hcalc = add(F, mul(r, tH))
    ok &= check("P3b: (F,r) ∩ {y=0} = H", norm(sub(Hcalc, H)))

    # P4: 크리스3 = F-K 방향 vs d2 (edgeB 방향) — crease3 = direction F-K vs d2 (edgeB direction)
    c3 = unit(sub(K, F))
    par = cross(c3, unit(d2))
    lam = (1 + x*x)**2 / (4*x*(1 - x*x))
    Delta = k - lam
    print(f"  P4: |크리스3 x edgeB| = {abs(par):.6f},  Delta = {Delta:.6f}")
    if tight:
        ok &= check("P4-tight: 크리스3이 edgeB에 밀착(평행)", par)
    else:
        ok &= check("P4-loose: 비평행이어야 함 (0이면 FAIL)", 0.0 if abs(par) > 1e-6 else 1.0)

    # P5: 타이트에서 변 방향이 36도의 배수 — at the tight state, edge directions are multiples of 36 degrees
    if tight:
        P = [E, H, G, K, F]
        for i in range(5):
            d = sub(P[(i+1) % 5], P[i])
            ang = math.degrees(math.atan2(d[1], d[0])) % 180.0
            rem = min(ang % 36.0, 36.0 - ang % 36.0)
            ok &= check(f"P5: 변 {i} 방향 {ang:.4f}deg ∈ 36Z", rem, tol=1e-6)
    return ok

allok = True
for xv in (0.15, 0.25, 0.30):
    allok &= run(xv)
allok &= run(math.sqrt(1 - 2/math.sqrt(5)), tight=True)
print("\n==== 종합:", "ALL OK" if allok else "실패 있음", "====")
