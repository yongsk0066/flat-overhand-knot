"""
2차 검증: 1차(check_slice.py)에서 확정된 사실
  - G = rho1({y=0}) ∩ {y=1}  (pass-2 아래모서리와 pass-1 위모서리의 교점)
  - H는 pass-2 모서리 위가 아님
를 바탕으로, 크리스2와 H의 정체에 대한 수정 가설을 검정한다.
Second verification: building on the facts established in the first pass (check_slice.py)
  - G = rho1({y=0}) ∩ {y=1}  (intersection of the lower edge of pass-2 with the upper edge of pass-1)
  - H does not lie on an edge of pass-2
this tests the revised hypotheses about the identity of crease2 and H.

수정 가설:
Revised hypotheses:
  W1. 크리스2 선 = G와 H를 잇는 직선 (G는 그 위 끝점, H는 선의 y=0 교점).
  W1. the crease2 line = the straight line joining G and H (G is its upper endpoint, H is that line's y=0 intersection).
  W2. l2 := rho2(edgeA) 는 G를 지나고 방향이 r과 평행. (edgeA = rho1({y=0}), G ∈ edgeA)
  W2. l2 := rho2(edgeA) passes through G and its direction is parallel to r. (edgeA = rho1({y=0}), G ∈ edgeA)
  W3. I = edgeB ∩ l2 또는 edgeA ∩ l2 = G + lambda r  중 어느 쪽인지 판별.
      (edgeB = rho1({y=1}), F 통과)
  W3. decide which of edgeB ∩ l2 and edgeA ∩ l2 is the one equal to I = G + lambda r.
      (edgeB = rho1({y=1}), passes through F)
  W4. 크리스2의 다른 끝점 P2 = (GH선) ∩ edgeB 의 위치 — y=0 위인가? (후보 접촉)
  W4. the location of crease2's other endpoint P2 = (the GH line) ∩ edgeB — is it on y=0? (candidate contact)
  W5. rho2(edgeB) (pass-3의 다른 모서리)가 H를 지나는가? E를 지나는가?
  W5. does rho2(edgeB) (the other edge of pass-3) pass through H? Through E?
  W6. H = rho2(edgeB) ∩ {y=0} 인가?
  W6. is H = rho2(edgeB) ∩ {y=0}?
  W7. rho3(H)의 정체: H가 pass-3 영역의 꼭짓점인가 — H가 l2와 rho2(edgeB) 중
      어디에 있는지, 그리고 크리스3 현의 끝점 K가 l2 위(자동)인지 확인.
  W7. the identity of rho3(H): is H a vertex of the pass-3 region — check which of l2 and
      rho2(edgeB) H lies on, and whether K, the endpoint of the crease3 chord, lies on l2 (automatic).
"""
import math

TOL = 1e-9

def sub(a, b): return (a[0]-b[0], a[1]-b[1])
def add(a, b): return (a[0]+b[0], a[1]+b[1])
def mul(a, s): return (a[0]*s, a[1]*s)
def dot(a, b): return a[0]*b[0] + a[1]*b[1]
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
    return v[0]*dh[1] - v[1]*dh[0]

def intersect(P1, d1, P2, d2):
    det = d1[0]*(-d2[1]) - d1[1]*(-d2[0])
    rhs = sub(P2, P1)
    t = (rhs[0]*(-d2[1]) - rhs[1]*(-d2[0])) / det
    return add(P1, mul(d1, t))

def report(label, val, tol=TOL):
    print(f"  [{'OK  ' if abs(val) < tol else 'no  '}] {label}: {val:.3e}")

def run(x):
    print(f"\n=== x = {x:.6f} ===")
    E = (0.0, 0.0); F = (-x, 1.0)
    G = ((1 - x*x) / (2*x), 1.0)
    H = ((1 - 3*x*x) / (2*x), 0.0)
    r = ((x*x - 1) / (1 + x*x), 2*x / (1 + x*x))
    lam = (1 + x*x)**2 / (4*x*(1 - x*x))
    k = (1 - 3*x*x) / (2*x)
    K = add(G, mul(r, k))

    rho1 = refl(E, sub(F, E))
    o1 = rho1((0.0, 0.0))
    d2 = sub(rho1((1.0, 0.0)), o1)
    edgeA = (E, d2)   # rho1({y=0}),  E 통과 — 1차에서 G∈edgeA 확인됨 — passes through E; the first pass confirmed G∈edgeA
    edgeB = (F, d2)   # rho1({y=1}),  F 통과 — passes through F

    # W1/W2: 크리스2 = GH선, l2 = rho2(edgeA) — crease2 = the GH line, l2 = rho2(edgeA)
    rho2 = refl(G, sub(H, G))
    o2 = rho2((0.0, 0.0))
    lin2 = lambda v: sub(rho2(v), o2)
    l2_dir = lin2(d2)
    report("W2a: rho2(edgeA) 방향 x r (평행성)", l2_dir[0]*r[1] - l2_dir[1]*r[0])
    report("W2b: rho2(E)가 (G, r)-직선 위", sdist(rho2(E), G, r))

    # W3: I가 어느 모서리와 l2의 교점인지 — which edge's intersection with l2 gives I
    IA = intersect(E, d2, G, r)
    IB = intersect(F, d2, G, r)
    tIA = dot(sub(IA, G), r); tIB = dot(sub(IB, G), r)
    print(f"  W3: edgeA∩l2 좌표 t={tIA:.6f}, edgeB∩l2 좌표 t={tIB:.6f}, lambda={lam:.6f}")
    report("W3a: edgeA∩l2 = G+lambda r", tIA - lam)
    report("W3b: edgeB∩l2 = G+lambda r", tIB - lam)

    # W4: 크리스2의 다른 끝점 = (GH선) ∩ edgeB — the other endpoint of crease2 = (the GH line) ∩ edgeB
    P2 = intersect(G, sub(H, G), F, d2)
    print(f"  W4: 크리스2 ∩ edgeB = ({P2[0]:.6f}, {P2[1]:.6f})")
    report("W4a: 그 점이 y=0 위인가", P2[1])
    report("W4b: 그 점이 H인가", norm(sub(P2, H)))

    # W5/W6: rho2(edgeB) 가 H 또는 E를 지나는가; {y=0}와의 교점이 H인가 — does rho2(edgeB) pass through H or E; is its intersection with {y=0} equal to H
    edgeB3 = (rho2(F), lin2(d2))   # pass-3의 다른 모서리 — the other edge of pass-3
    report("W5a: rho2(edgeB)가 H를 지남", sdist(H, *edgeB3))
    report("W5b: rho2(edgeB)가 E를 지남", sdist(E, *edgeB3))
    Q = intersect(edgeB3[0], edgeB3[1], (0.0, 0.0), (1.0, 0.0))
    print(f"  W6: rho2(edgeB) ∩ {{y=0}} = ({Q[0]:.6f}, {Q[1]:.6f}),  H = ({H[0]:.6f}, 0)")
    report("W6: 그 교점 = H", norm(sub(Q, H)))

    # W7: H의 소속 — H가 l2 위인가, edgeB3 위인가; K는 l2 위(자동인지) — where H belongs: is H on l2 or on edgeB3; is K on l2 (automatically)
    report("W7a: H가 l2 위", sdist(H, G, r))
    report("W7b: K가 l2 위 (정의상 자동)", sdist(K, G, r))
    # 크리스3 = F-K 선이 pass-3 띠의 완전 현이 되려면 다른 끝점이 edgeB3 위여야 함 — for crease3 = the F-K line to be a full chord of the pass-3 band, its other endpoint must lie on edgeB3
    rho3 = refl(F, sub(K, F))
    P3 = intersect(F, sub(K, F), edgeB3[0], edgeB3[1])
    print(f"  W7c: 크리스3 ∩ rho2(edgeB) = ({P3[0]:.6f}, {P3[1]:.6f}), F=({F[0]:.4f},1)")
    report("W7d: F가 rho2(edgeB) 위 (크리스3 끝점=F?)", sdist(F, *edgeB3))

for xv in (0.15, 0.25, 0.324920):
    run(xv)
