"""
4차 검증: 재료 위상(material topology)의 완전 재구성 — 순방향 걷기.
Fourth verification: complete reconstruction of the material topology — the forward walk.

발견된 문제: 03의 크리스 목록 (EF, FG⊂{y=1}, FK)로는 재료 순서가 모순됨
  (크리스 FK의 현이 패스3의 재료 영역 밖에 놓임 — 순방향 도달 불가).
Problem found: with the crease list of 03 (EF, FG⊂{y=1}, FK) the material order is contradictory
  (the chord of crease FK lands outside the material region of pass-3 — unreachable going forward).

수정 가설 (이 스크립트가 검정):
Revised hypothesis (what this script tests):
  실제 크리스 = { EF, HG, FK }  (FG는 크리스가 아님; y=1 밀착 해석 철회)
  the actual creases = { EF, HG, FK }  (FG is not a crease; the "flush against y=1" reading is withdrawn)
  재료 순서:
  material order:
    레이인 [밴드(e_A,e_B), 위-오른쪽에서 진행방향 d2로 진입]
    ray-in [the band (e_A,e_B), entering from the upper right along the travel direction d2]
    -> 크리스 EF -> 중간 패널 [밴드1 {0<=y<=1}, 현 EF와 HG 사이; 오각형 영역을 덮음]
    -> crease EF -> middle panel [band1 {0<=y<=1}, between the chords EF and HG; covers the pentagon region]
    -> 크리스 HG -> 패널 pi2 [밴드((G,r),(F,r)), 진행방향 r]
    -> crease HG -> panel pi2 [the band ((G,r),(F,r)), travel direction r]
    -> 크리스 FK -> 레이아웃.
    -> crease FK -> ray-out.
  미러 대칭: E<->H, F<->G (양 끝 교환).
  mirror symmetry: E<->H, F<->G (the two ends swap).

검정 항목:
Items tested:
  T1. C1: G ∈ e_A (레이인 모서리가 중간 크리스 끝점 G를 지남) — 슬라이스 등식.
  T1. C1: G ∈ e_A (the ray-in edge passes through G, an endpoint of the middle crease) — the slice identity.
  T2. C2: F ∈ (H,r)-직선 (= rho_HG({y=0}) = pi2/레이아웃 쪽 모서리가 F를 지남) — 미러 등식.
  T2. C2: F ∈ the line (H,r) (= rho_HG({y=0}) = the pi2/ray-out-side edge passes through F) — the mirror identity.
  T3. 걷기: rho_EF(d2) = (1,0);  M1=(-x/2,1/2) -> M2=중점(H,G), |M1M2|=(1-x^2)/(2x);
      rho_HG(1,0) = r;  M3 = (M2 + s r) ∩ 직선FK, s>0, s=(1-x^2)/(2x), M3가 현 FK 내부(0<tau<1).
  T3. the walk: rho_EF(d2) = (1,0);  M1=(-x/2,1/2) -> M2=midpoint(H,G), |M1M2|=(1-x^2)/(2x);
      rho_HG(1,0) = r;  M3 = (M2 + s r) ∩ the line FK, s>0, s=(1-x^2)/(2x), M3 interior to the chord FK (0<tau<1).
  T4. ell_core(x) = |M1M2|+|M2M3| = 1/x - x  (순감소: 도함수 -1/x^2-1 < 0).
  T4. ell_core(x) = |M1M2|+|M2M3| = 1/x - x  (strictly decreasing: derivative -1/x^2-1 < 0).
  T5. 출구 방향 rho_FK(r): 단위벡터; x0에서 꼬리 사이 각도 보고.
  T5. exit direction rho_FK(r): a unit vector; reports the angle between the tails at x0.
  T6. x0에서 오각형 변 = {크리스EF, 크리스HG, 크리스FK, 모서리 y=0, 모서리 l2} 위 (구조 서술 확인).
  T6. at x0 the pentagon's sides lie on {crease EF, crease HG, crease FK, edge y=0, edge l2} (confirms the structural description).
"""
import math

TOL = 1e-9

def sub(a,b): return (a[0]-b[0], a[1]-b[1])
def add(a,b): return (a[0]+b[0], a[1]+b[1])
def mul(a,s): return (a[0]*s, a[1]*s)
def dot(a,b): return a[0]*b[0]+a[1]*b[1]
def cross(a,b): return a[0]*b[1]-a[1]*b[0]
def norm(a): return math.hypot(a[0],a[1])
def unit(a):
    n = norm(a); return (a[0]/n, a[1]/n)

def refl_lin(d):
    dh = unit(d)
    return lambda v: sub(mul(dh, 2*dot(v,dh)), v)

def intersect(P1,d1,P2,d2):
    det = d1[0]*(-d2[1]) - d1[1]*(-d2[0])
    rhs = sub(P2,P1)
    t = (rhs[0]*(-d2[1]) - rhs[1]*(-d2[0])) / det
    return add(P1, mul(d1,t)), t

def check(label, val, tol=TOL):
    ok = abs(val) < tol
    print(f"  [{'OK' if ok else 'FAIL'}] {label}: {val:.3e}")
    return ok

def run(x, tight=False):
    print(f"\n=== x = {x:.7f} {'(tight)' if tight else ''} ===")
    ok = True
    E=(0.0,0.0); F=(-x,1.0)
    G=((1-x*x)/(2*x), 1.0); H=((1-3*x*x)/(2*x), 0.0)
    r=((x*x-1)/(1+x*x), 2*x/(1+x*x))
    k=(1-3*x*x)/(2*x); K=add(G, mul(r,k))
    rl1 = refl_lin(sub(F,E))
    d2 = rl1((1.0,0.0))

    # T1, T2
    ok &= check("T1/C1: G ∈ e_A = (E,d2)", cross(sub(G,E), unit(d2)))
    ok &= check("T2/C2: F ∈ (H,r)  [= rho_HG(y=0) 모서리]", cross(sub(F,H), unit(r)))

    # T3: 걷기 — the walk
    cont1 = rl1(d2)
    ok &= check("T3a: rho_EF(d2) = (1,0)", norm(sub(unit(cont1),(1.0,0.0))))
    M1 = (-x/2, 0.5)
    M2 = ((H[0]+G[0])/2, 0.5)
    s12 = norm(sub(M2,M1))
    ok &= check("T3b: |M1M2| = (1-x^2)/(2x)", s12 - (1-x*x)/(2*x))
    rlHG = refl_lin(sub(G,H))
    contHG = rlHG((1.0,0.0))
    ok &= check("T3c: rho_HG(1,0) = r", norm(sub(unit(contHG), unit(r))))
    M3, s23 = intersect(M2, unit(r), F, sub(K,F))
    ok &= check("T3d: M3 전방 도달 (s>0)", 0.0 if s23 > 0 else 1.0)
    ok &= check("T3e: |M2M3| = (1-x^2)/(2x)", s23 - (1-x*x)/(2*x))
    tau = dot(sub(M3,F), sub(K,F))/dot(sub(K,F),sub(K,F))
    ok &= check("T3f: M3가 현 FK 내부 (0<tau<1)", 0.0 if 0 < tau < 1 else 1.0)
    print(f"        tau(M3 on FK) = {tau:.6f}  (대칭이면 0.5)")

    # T4: ell_core
    ell = s12 + s23
    ok &= check("T4: ell_core = 1/x - x", ell - (1/x - x))
    print(f"        ell_core({x:.4f}) = {ell:.6f},  d(ell)/dx = {-(1/x**2+1):.4f} < 0")

    # T5: 출구 방향 — exit direction
    rlFK = refl_lin(sub(K,F))
    out = rlFK(unit(r))
    ok &= check("T5: |rho_FK(r)| = 1", norm(out)-1.0)
    if tight:
        ang_in  = math.degrees(math.atan2(d2[1], d2[0])) % 360
        ang_out = math.degrees(math.atan2(out[1], out[0])) % 360
        print(f"        레이인 진행각 {ang_in:.3f}°, 레이아웃 진행각 {ang_out:.3f}°, 사이각 {abs(ang_in-ang_out)%360:.3f}°")
        # T6: 오각형 변의 지지선 — the supporting lines of the pentagon's sides
        P = {'E':E,'H':H,'G':G,'K':K,'F':F}
        sides = [('E','H','y=0 (밴드1 모서리)'), ('H','G','크리스 HG'),
                 ('G','K','l2=(G,r) (pi2 모서리)'), ('K','F','크리스 FK'), ('F','E','크리스 EF')]
        lines = {'y=0 (밴드1 모서리)': (E,(1.0,0.0)), '크리스 HG': (H, sub(G,H)),
                 'l2=(G,r) (pi2 모서리)': (G,r), '크리스 FK': (F, sub(K,F)), '크리스 EF': (E, sub(F,E))}
        for a,b,ln in sides:
            P0,d0 = lines[ln]
            v1 = cross(sub(P[a],P0), unit(d0)); v2 = cross(sub(P[b],P0), unit(d0))
            ok &= check(f"T6: 변 {a}{b} ⊂ {ln}", abs(v1)+abs(v2))
    return ok

allok = True
for xv in (0.15, 0.25, 0.30):
    allok &= run(xv)
allok &= run(math.sqrt(1 - 2/math.sqrt(5)), tight=True)
print("\n==== 종합:", "ALL OK — 수정 재료 위상 확정" if allok else "실패 있음", "====")
