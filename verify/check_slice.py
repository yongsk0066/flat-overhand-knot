"""
03-moduli.md의 구조적 주장에 대한 수치 검증.
Numerical verification of the structural claims in 03-moduli.md.

검증 대상 (asset-1.md의 표준 옭매듭 가지, 폭 1 정규화):
Verification targets (the standard overhand-knot branch of asset-1.md, width-1 normalization):
  V1. 크리스2 = 선분 GH 가설: H는 pass-2의 한 모서리 위, G는 다른 모서리 위.
  V1. hypothesis crease2 = segment GH: H lies on one edge of pass-2, G on the other edge.
  V2. l2 (pass-3 모서리) = G + R*r 이고 방향이 r과 평행.
  V2. l2 (an edge of pass-3) = G + R*r, and its direction is parallel to r.
  V3. I = (pass-2 모서리) ∩ l2 = G + lambda(x) r,  lambda = (1+x^2)^2 / (4x(1-x^2)).
  V3. I = (an edge of pass-2) ∩ l2 = G + lambda(x) r,  lambda = (1+x^2)^2 / (4x(1-x^2)).
  V4. 크리스3 = F와 K=G+k(x)r를 잇는 선, k = (1-3x^2)/(2x)이고,
      rho3(H) = G + t r, t = (1-x^2)/x  (즉 rho3(H)가 l2 위).
  V4. crease3 = the line joining F and K=G+k(x)r, with k = (1-3x^2)/(2x), and
      rho3(H) = G + t r, t = (1-x^2)/x  (i.e. rho3(H) lies on l2).
  V5. Delta(x) = k - lambda = (5x^4-10x^2+1)/(4x(1-x^2)).
  V6. x0 = sqrt(1 - 2/sqrt(5)) 에서 Delta=0이고 오각형 E,H,G,K,F가 정오각형
      (다섯 변 길이 동일, 내각 108도).
  V6. at x0 = sqrt(1 - 2/sqrt(5)), Delta=0 and the pentagon E,H,G,K,F is regular
      (five equal side lengths, interior angles of 108 degrees).
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
    """P를 지나고 방향 d인 직선에 대한 반사.
    Reflection in the line through P with direction d."""
    dh = unit(d)
    def rho(Q):
        v = sub(Q, P)
        w = sub(mul(dh, 2*dot(v, dh)), v)
        return add(P, w)
    return rho

def sdist(Q, P, d):
    """점 Q에서 직선(P, 방향 d)까지의 부호 있는 거리.
    Signed distance from the point Q to the line (P, direction d)."""
    dh = unit(d); v = sub(Q, P)
    return v[0]*dh[1] - v[1]*dh[0]

def intersect(P1, d1, P2, d2):
    """P1+t d1 = P2+s d2 풀기. 반환: 교점.
    Solve P1+t d1 = P2+s d2. Returns: the intersection point."""
    det = d1[0]*(-d2[1]) - d1[1]*(-d2[0])
    rhs = sub(P2, P1)
    t = (rhs[0]*(-d2[1]) - rhs[1]*(-d2[0])) / det
    return add(P1, mul(d1, t))

def check(label, val, expect=0.0, tol=TOL):
    ok = abs(val - expect) < tol
    print(f"  [{'OK' if ok else 'FAIL'}] {label}: {val:.3e}" + ("" if expect == 0 else f" (expect {expect:.6f})"))
    return ok

def run(x, tight=False):
    print(f"\n=== x = {x:.6f} {'(tight x0)' if tight else ''} ===")
    allok = True
    E = (0.0, 0.0)
    F = (-x, 1.0)
    G = ((1 - x*x) / (2*x), 1.0)
    H = ((1 - 3*x*x) / (2*x), 0.0)
    r = ((x*x - 1) / (1 + x*x), 2*x / (1 + x*x))

    # pass 1 = 띠 {0<=y<=1}. 크리스1 = 직선 EF. rho1 = 그 반사. — pass 1 = the band {0<=y<=1}. crease1 = the line EF. rho1 = reflection in it.
    rho1 = refl(E, sub(F, E))
    o = rho1((0.0, 0.0))
    d2 = sub(rho1((1.0, 0.0)), o)          # rho1의 선형부가 보낸 (1,0) = pass-2 방향 — (1,0) carried by the linear part of rho1 = the pass-2 direction
    # pass-2 모서리: edgeA = rho1({y=0}) (E 통과), edgeB = rho1({y=1}) (F 통과) — edges of pass-2: edgeA = rho1({y=0}) (through E), edgeB = rho1({y=1}) (through F)
    edgeA = (E, d2)
    edgeB = (F, d2)

    # V1: H가 edgeA 또는 edgeB 위, G가 나머지 위 — H lies on edgeA or edgeB, and G on the remaining one
    dHA, dHB = sdist(H, *edgeA), sdist(H, *edgeB)
    dGA, dGB = sdist(G, *edgeA), sdist(G, *edgeB)
    print(f"  H->edgeA {dHA:.3e}, H->edgeB {dHB:.3e}; G->edgeA {dGA:.3e}, G->edgeB {dGB:.3e}")
    if abs(dHA) < TOL and abs(dGB) < TOL:
        print("  [OK] V1: H는 edgeA(=rho1(y=0)) 위, G는 edgeB(=rho1(y=1)) 위 -> 크리스2 = GH는 pass-2의 완전 현")
        H_edge, G_edge = edgeA, edgeB
    elif abs(dHB) < TOL and abs(dGA) < TOL:
        print("  [OK] V1: H는 edgeB 위, G는 edgeA 위 -> 크리스2 = GH는 pass-2의 완전 현")
        H_edge, G_edge = edgeB, edgeA
    else:
        print("  [FAIL] V1: 크리스2=GH 가설 위배"); allok = False
        H_edge, G_edge = edgeA, edgeB

    # V2: rho2 = GH 반사. pass-3 모서리 = rho2(pass-2 모서리). G를 지나는 쪽이 l2 = G + R r ? — rho2 = reflection in GH. edges of pass-3 = rho2(edges of pass-2). Is the one through G equal to l2 = G + R r?
    rho2 = refl(H, sub(G, H))
    o2 = rho2((0.0, 0.0))
    lin2 = lambda v: sub(rho2(v), o2)
    edgeG_img = (rho2(G_edge[0]), lin2(G_edge[1]))   # G가 놓인 pass-2 모서리의 상 — the image of the pass-2 edge on which G lies
    cross = edgeG_img[1][0]*r[1] - edgeG_img[1][1]*r[0]
    allok &= check("V2a: rho2(G쪽 모서리) 방향 x r (평행성)", cross)
    allok &= check("V2b: G가 그 상 위에 있는가", sdist(G, *edgeG_img))
    l2 = (G, r)

    # V3: I = (H쪽 pass-2 모서리) ∩ l2 가 G + lambda r 인가 — is I = (the H-side pass-2 edge) ∩ l2 equal to G + lambda r
    lam = (1 + x*x)**2 / (4*x*(1 - x*x))
    I = intersect(H_edge[0], H_edge[1], G, r)
    tI = dot(sub(I, G), unit(r)) / 1.0
    allok &= check("V3: I의 l2-좌표 = lambda", tI - lam * (1.0 if dot(r, r) else 1.0), tol=1e-8)

    # V4: 크리스3 = F--K, K = G + k r. rho3(H)가 l2 위이고 t = (1-x^2)/x 인가 — crease3 = F--K, K = G + k r. Is rho3(H) on l2 with t = (1-x^2)/x
    k = (1 - 3*x*x) / (2*x)
    K = add(G, mul(r, k / norm(r)))  # 주의: r는 단위벡터(검증) — note: r is a unit vector (verified)
    allok &= check("r는 단위벡터", norm(r) - 1.0)
    K = add(G, mul(r, k))
    rho3 = refl(F, sub(K, F))
    Hp = rho3(H)
    allok &= check("V4a: rho3(H)가 l2 위", sdist(Hp, G, r))
    tH = dot(sub(Hp, G), r)
    allok &= check("V4b: rho3(H)의 l2-좌표 = (1-x^2)/x", tH - (1 - x*x)/x, tol=1e-8)
    allok &= check("V4c: |F rho3(H)| = |F H| (반사 정합)", norm(sub(Hp, F)) - norm(sub(H, F)))

    # V5: Delta = k - lambda 공식 — the formula Delta = k - lambda
    Delta = k - lam
    Delta_formula = (5*x**4 - 10*x*x + 1) / (4*x*(1 - x*x))
    allok &= check("V5: Delta = (5x^4-10x^2+1)/(4x(1-x^2))", Delta - Delta_formula)
    print(f"  Delta(x) = {Delta:.6f}")

    # V6: 타이트면 정오각형 — if tight, a regular pentagon
    if tight:
        allok &= check("V6a: Delta(x0) = 0", Delta, tol=1e-9)
        P = [E, H, G, K, F]  # 변: EH, HG, GK, KF, FE — sides: EH, HG, GK, KF, FE
        sides = [norm(sub(P[(i+1) % 5], P[i])) for i in range(5)]
        print(f"  변 길이: {[f'{s:.6f}' for s in sides]}")
        allok &= check("V6b: 변 길이 편차", max(sides) - min(sides))
        for i in range(5):
            a, b, c = P[(i-1) % 5], P[i], P[(i+1) % 5]
            v1, v2 = unit(sub(a, b)), unit(sub(c, b))
            ang = math.degrees(math.acos(max(-1, min(1, dot(v1, v2)))))
            allok &= check(f"V6c: 내각@{i} = 108", ang, expect=108.0, tol=1e-6)
    return allok

ok = True
for xv in (0.15, 0.25, 0.30):
    ok &= run(xv)
x0 = math.sqrt(1 - 2/math.sqrt(5))
ok &= run(x0, tight=True)
print("\n==== 종합:", "ALL OK" if ok else "실패 항목 있음 — 해석 수정 필요", "====")
