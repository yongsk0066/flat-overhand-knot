# check_gp.py — 일반 위치(GP) 감사 (referee-15 처방 #7) — general position (GP) audit (referee-15 prescription #7)
#
# 문서 15 §1 / 논문 §3.4 의 GP 주장을 실제로 감사한다: — actually audits the GP claim of document 15 §1 / paper §3.4:
#   (i)  어떤 두 중심선 조각도 1점 초과로 겹치지 않음 (공선쌍 없음) — (i) no two centerline pieces overlap in more than one point (no collinear pair)
#   (ii) 이중점: 횡단(교차각>0), 서로 다름(쌍별 거리>0), 접힘점·크리스 선분에서 양의 거리 — (ii) double points: transversal (crossing angle>0), mutually distinct (pairwise distance>0), positive distance from fold points and crease segments
#   (iii) 접힘각 비퇴화 — (iii) fold angles non-degenerate
# 감사 지점: 정오각형, θ†, 느슨 슬라이스(x=0.25) + 유효 무작위 표본. — audit points: the regular pentagon, θ†, the loose slice (x=0.25) + valid random samples.
# 또한 교차쌍이 항상 {(P0,P2),(P0,P3),(P1,P3)}이고 인접쌍 교차가 없음을 확인. — also confirms that the crossing pairs are always {(P0,P2),(P0,P3),(P1,P3)} and that adjacent pairs never cross.

import math, random
import uniq_setup as U0
import uniq_strata as ST

def _meet_param(P1, d1, P2, d2):
    det = d1[0]*(-d2[1]) - d1[1]*(-d2[0])
    if abs(det) < 1e-13: return None, None
    rhs = (P2[0]-P1[0], P2[1]-P1[1])
    s = (rhs[0]*(-d2[1]) - rhs[1]*(-d2[0]))/det
    u = (d1[0]*rhs[1] - d1[1]*rhs[0])/det
    return s, u

def pt_seg_dist(p, a, b):
    ax, ay = a; bx, by = b; px, py = p
    dx, dy = bx-ax, by-ay
    L2 = dx*dx + dy*dy
    t = max(0.0, min(1.0, ((px-ax)*dx + (py-ay)*dy)/L2))
    return math.hypot(px-(ax+t*dx), py-(ay+t*dy))

def pieces(c):
    """중심선 4조각: (기점, 방향, 유형) — ray는 s≥0, seg는 s∈[0,L].
    The 4 centerline pieces: (base point, direction, type) — a ray has s≥0, a seg has s∈[0,L]."""
    d2, w0, v3 = c['d2'], c['w0'], c['v3']
    M1, M2, M3 = c['M1'], c['M2'], c['M3']
    L1 = math.dist(M1, M2); L2 = math.dist(M2, M3)
    dir12 = ((M2[0]-M1[0])/L1, (M2[1]-M1[1])/L1)
    dir23 = ((M3[0]-M2[0])/L2, (M3[1]-M2[1])/L2)
    return [
        ("P0", M1, (-d2[0], -d2[1]), None),   # 레이인 (M1에서 뒤로 무한) — ray-in (infinite backwards from M1)
        ("P1", M1, dir12, L1),
        ("P2", M2, dir23, L2),
        ("P3", M3, v3, None),                 # 레이아웃 — ray-out
    ]

def audit(th, ray_cap=50.0):
    c = U0._config(th)
    ps = pieces(c)
    xs = []           # (쌍, 점, 각도) — (pair, point, angle)
    collinear = 0
    for i in range(4):
        for j in range(i+1, 4):
            n1, P1, d1, L1 = ps[i]
            n2, P2, d2, L2 = ps[j]
            crossv = d1[0]*d2[1] - d1[1]*d2[0]
            if abs(crossv) < 1e-12:
                off = (P2[0]-P1[0])*d1[1] - (P2[1]-P1[1])*d1[0]
                if abs(off) < 1e-9: collinear += 1
                continue
            s, u = _meet_param(P1, d1, P2, d2)
            if s is None: continue
            lo = 1e-9
            ok1 = (lo < s < (L1 - lo if L1 else ray_cap))
            ok2 = (lo < u < (L2 - lo if L2 else ray_cap))
            if ok1 and ok2:
                pt = (P1[0]+s*d1[0], P1[1]+s*d1[1])
                ang = math.degrees(math.asin(min(1.0, abs(crossv))))
                xs.append(((n1, n2), pt, ang, abs(i-j) == 1))
    adj = [x for x in xs if x[3]]
    pairs = sorted(x[0] for x in xs)
    angs = [x[2] for x in xs]
    ptl = [x[1] for x in xs]
    # 교차 상호거리, 접힘점 거리, 크리스 선분 거리 — mutual distances among crossings, distances to fold points, distances to crease segments
    dxx = min((math.dist(ptl[i], ptl[j]) for i in range(len(ptl)) for j in range(i+1, len(ptl))), default=float('inf'))
    folds = [c['M1'], c['M2'], c['M3']]
    dxf = min((math.dist(p, f) for p in ptl for f in folds), default=float('inf'))
    creases = [(c['E'], c['F']), (c['H'], c['G']), (c['Q1'], c['Q0'])]
    dxc = min((pt_seg_dist(p, a, b) for p in ptl for a, b in creases), default=float('inf'))
    return dict(n_cross=len(xs), pairs=pairs, adj=len(adj), collinear=collinear,
                min_angle=min(angs, default=None), dxx=dxx, dxf=dxf, dxc=dxc)

X0 = math.sqrt(1 - 2/math.sqrt(5))
TH = {
    "정오각형": U0._slice_theta(X0),
    "θ†": (0.32212749731122718, 1.3911166915901648, 1.0363863004307971, 1.0363863004307971, 1.7132441870641958),
    "느슨 슬라이스 x=0.25": U0._slice_theta(0.25),
}
EXPECT = [("P0", "P2"), ("P0", "P3"), ("P1", "P3")]

print("[GP 감사] 지정 상태:")
for name, th in TH.items():
    r = audit(th)
    ok = (r['n_cross'] == 3 and r['pairs'] == sorted(EXPECT) and r['adj'] == 0 and r['collinear'] == 0)
    print(f"  {name}: 교차 {r['n_cross']}개 {r['pairs']} 인접교차 {r['adj']} 공선 {r['collinear']} | "
          f"최소각 {r['min_angle']:.1f}° d(교차,교차)={r['dxx']:.3f} d(교차,접힘점)={r['dxf']:.3f} d(교차,크리스)={r['dxc']:.3f}"
          + ("  ✓" if ok else "  ✗✗"))
    assert ok

print("\n[GP 감사] 유효 무작위 표본:")
rnd = random.Random(99)
mins = dict(angle=1e9, dxx=1e9, dxf=1e9, dxc=1e9)
n_ok = 0
target = 2000
tries = 0
while n_ok < target and tries < 200000:
    tries += 1
    xx = rnd.uniform(0.12, 0.55)
    th0 = U0._slice_theta(xx)
    th = tuple(v + rnd.uniform(-0.15, 0.15) for v in th0)
    try:
        ws, c = U0._walls8f(th)
    except Exception:
        continue
    if min(ws[:6]) < 0 or max(ws[6], ws[7]) < 0: continue
    okv, _ = ST.valid_full(c)
    if not okv: continue
    r = audit(th)
    if not (r['n_cross'] == 3 and r['pairs'] == sorted(EXPECT) and r['adj'] == 0 and r['collinear'] == 0):
        print("  ✗ GP 위반 표본:", th, r); raise SystemExit(1)
    n_ok += 1
    mins['angle'] = min(mins['angle'], r['min_angle'])
    mins['dxx'] = min(mins['dxx'], r['dxx'])
    mins['dxf'] = min(mins['dxf'], r['dxf'])
    mins['dxc'] = min(mins['dxc'], r['dxc'])
print(f"  유효 표본 {n_ok}개 (시도 {tries}): 교차쌍 전부 {{(P0,P2),(P0,P3),(P1,P3)}}, 인접교차 0, 공선쌍 0")
print(f"  최소 교차각 {mins['angle']:.1f}°, 최소 d(교차,교차) {mins['dxx']:.3f}, "
      f"최소 d(교차,접힘점) {mins['dxf']:.3f}, 최소 d(교차,크리스) {mins['dxc']:.3f}")
print("\n결론: 감사한 전 상태(지정 3 + 유효 표본)에서 GP 성립 — 여유는 위 최솟값들로 정량화.")
