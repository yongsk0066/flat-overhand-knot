"""Explore the three-crease branch with multistart descent and search for other blocked states.

설계 (2단계):
Design (2 stages):
  0. 시작점: 표준 슬라이스 x in {0.10,...,0.32} 에서 좌표별 +-m (m in {0.03,0.06,0.10})
     균등 섭동 후 유효한 점만 채택 (6벽 >= -1e-9, taco max(w9,w10) >= -1e-9, valid_full).
  0. start points: from the standard slice x in {0.10,...,0.32}, perturb each coordinate uniformly
     by +-m (m in {0.03,0.06,0.10}) and keep only the valid points (6 walls >= -1e-9,
     taco max(w9,w10) >= -1e-9, valid_full).
  1. [탐사] 무작위 실현가능 하강: 무작위 단위방향, 한 스텝 크기에서 연속 120회 실패 시
     절반 (2e-3 -> 1e-5, 성공 시 1.3배 재성장, 상한 1e-2). 매 수락 스텝마다 valid_full.
     -> 어느 활성면(벽 집합)에 떨어지는지 무편향 탐사.
  1. [exploration] random feasible descent: random unit directions; after 120 consecutive failures
     at a given step size, halve it (2e-3 -> 1e-5, regrown by 1.3x on success, capped at 1e-2).
     valid_full at every accepted step.
     -> an unbiased exploration of which active face (wall set) it falls onto.
  2. [정련] 활성면 사영 하강: 정체점의 활성 집합 S(|w|<1e-4)로 뉴턴 사영 후,
     {w_j=0, j in S} 접공간에 사영한 -grad ell 방향으로 백트래킹 하강. 벽 충돌 시
     활성 집합 확장. |사영 grad| < 1e-7 이면 면-KKT 점으로 정지.
     (파일럿 관찰: 1단계 무작위 하강은 4-등식 면 위에서 접선 원뿔이 얇아 조기 정체 —
      전부 제2족 위의 stragglers 였다. 2단계가 이를 면-최소점까지 끌고 간다.)
  2. [refinement] projected descent on the active face: Newton-project onto the active set
     S (|w|<1e-4) of the stagnation point, then backtracking descent along -grad ell projected onto
     the tangent space of {w_j=0, j in S}. On a wall collision, enlarge the active set.
     Stop at a face-KKT point once |projected grad| < 1e-7.
     (pilot observation: the stage-1 random descent stagnates early on the 4-equation face, where the
      tangent cone is thin — they were all stragglers on the second family. Stage 2 drags them down
      to the face minimum.)
  3. 종점 클러스터링 (theta 최대노름 2e-3), 활성 집합(|w|<1e-6) 기록,
     알려진 어트랙터 theta†/오각형과 비교.
  3. cluster the endpoints (max-norm 2e-3 in theta), record the active set (|w|<1e-6), and
     compare against the known attractors theta†/the pentagon.
  4. theta†/오각형 외 수렴 종점 발견 시 check_branch 방식 합법성 감사
     (조각 교차 조사 + 크리스쌍 + 199점 공리).
  4. if a converged endpoint other than theta†/the pentagon turns up, audit its legitimacy the way
     check_branch does (survey of piece crossings + crease pairs + the 199-point axiom).

check_edge.py 의 config/walls8f/valid_full/slice_theta/gradv/covers_* 를 재사용
(모듈 하단의 추적 루프는 잘라내고 함수 정의부만 exec).
Reuses config/walls8f/valid_full/slice_theta/gradv/covers_* from check_edge.py
(the tracing loop at the bottom of the module is cut away; only the function definitions are exec'd).
"""
import math, random, time, os, sys, json

HERE = os.path.dirname(os.path.abspath(__file__))
_src = open(os.path.join(HERE, 'check_edge.py')).read()
_head = _src.split('\nx0 = math.sqrt')[0]          # 함수 정의부만 (추적 루프 제외) — the function definitions only (the tracing loop is excluded)
NS = {}
exec(compile(_head, 'check_edge.py[head]', 'exec'), NS)
config      = NS['config']
walls8f     = NS['walls8f']
valid_full  = NS['valid_full']
slice_theta = NS['slice_theta']
gradv       = NS['gradv']
covers_A    = NS['covers_A']
covers_D    = NS['covers_D']
sub, add, mul, dot, cross, norm, unit = (NS[k] for k in
    ('sub','add','mul','dot','cross','norm','unit'))

# ---------- 파라미터 (조절 내역은 stdout 로그에 남김) ---------- — parameters (the tuning record is left in the stdout log)
SLICE_XS   = [0.10, 0.15, 0.20, 0.25, 0.30, 0.32]
MAGS       = ([float(z) for z in os.environ['MAGS'].split(',')]
              if 'MAGS' in os.environ else [0.03, 0.06, 0.10])
KEEP_CELL  = int(os.environ.get('KEEP_CELL', 5))      # 셀당 채택 시작점 (6*3*5 = 90) — start points kept per cell (6*3*5 = 90)
MAX_TRY    = 2000    # 셀당 섭동 시도 상한 ([조절] 오각형 인접 셀 x=0.32,m=0.10의
                     # 채택률이 <1/400 로 낮아 상향; 비용은 셀당 ~1s)
                     # — cap on perturbation attempts per cell ([tuning] raised because the cell
                     # adjacent to the pentagon, x=0.32, m=0.10, has an acceptance rate below 1/400;
                     # the cost is ~1s per cell)
STEP0      = 2e-3
STEP_MIN   = 1e-5
STEP_MAX   = 1e-2    # [조절] 원거리 시작점(x=0.10, ell~9.9)의 수렴을 위해 성공 시 1.3배
                     # 재성장을 1e-2까지 허용. 축소 스케줄(120연속실패 -> 절반, 하한 1e-5)은 명세대로.
                     # — [tuning] to let far-off start points (x=0.10, ell~9.9) converge in time,
                     # regrowth by 1.3x on success is allowed up to 1e-2. The shrink schedule
                     # (120 consecutive failures -> halve, floor 1e-5) is as specified.
FAILS_LVL  = 120
MAX_STEPS  = 12000
BUDGET1    = float(os.environ.get('BUDGET1', 300.0))   # 1단계 총 예산(초) — total budget for stage 1 (seconds)
BUDGET2    = float(os.environ.get('BUDGET2', 400.0))   # 2단계 총 예산(초) — total budget for stage 2 (seconds)
SEED       = int(os.environ.get('SEED', 20260803))
CLUST_TOL  = 2e-3
ACT_TOL    = 1e-6    # 최종 활성 판정 — final activity decision
PG_TOL     = 1e-7    # 사영 기울기 정지 (유한차분 잡음 바닥 ~1e-8 고려) — projected-gradient stopping (accounting for the ~1e-8 finite-difference noise floor)

WNAMES = ['w1','w4','w5','w6','w7','w8']

# 알려진 어트랙터 — the known attractors
XD = 0.3221274973112271792397142057801235136978
TH_DAG = (XD, (1-XD*XD)/(2*XD),
          1.036386300430797104445252761398636908419,
          1.036386300430797104445252761398636908419,
          1.713244187064195830219559769216923187227)
ELL_DAG = 2.749630487494992934664812530615560095646
X_PENT = math.sqrt(1-2/math.sqrt(5))
TH_PENT = slice_theta(X_PENT)

def dmax(a, b): return max(abs(a[k]-b[k]) for k in range(5))

def ell_c(c):
    return norm(sub(c['M2'],c['M1'])) + norm(sub(c['M3'],c['M2']))

def ell_fn(th):
    try:
        c = config(th)
        if c['M3'] is None: return float('inf')
        return ell_c(c)
    except Exception:
        return float('inf')

def feas(th, tol=-1e-9):
    """벽 실현가능성 + ell. valid_full은 호출자가 별도로 (비싸므로 나중에).
    Wall feasibility + ell. valid_full is left to the caller separately (expensive, so deferred)."""
    try:
        c = config(th)
    except Exception:
        return False, None, None, None
    if c['M3'] is None or c['s23'] is None or c['s23'] <= 0:
        return False, None, None, None
    try:
        ws = walls8f(th, c)
    except Exception:
        return False, None, None, None
    if min(ws[:6]) < tol or max(ws[6], ws[7]) < tol:
        return False, c, ws, None
    return True, c, ws, ell_c(c)

def fully_valid(th):
    ok, c, ws, e = feas(th)
    if not ok: return False, c, ws, e
    okv, _why = valid_full(c)
    if not okv: return False, c, ws, e
    return True, c, ws, e

# ---------- 선형대수 소도구 ---------- — small linear-algebra utilities
def solve_k(M, rhs):
    k = len(rhs)
    A = [M[i][:] + [rhs[i]] for i in range(k)]
    for col in range(k):
        p = max(range(col, k), key=lambda r: abs(A[r][col]))
        if abs(A[p][col]) < 1e-14: return None
        A[col], A[p] = A[p], A[col]
        for r in range(k):
            if r != col:
                f = A[r][col]/A[col][col]
                A[r] = [A[r][m] - f*A[col][m] for m in range(k+1)]
    return [A[i][k]/A[i][i] for i in range(k)]

def jacS(th, S):
    return [gradv(lambda t, j=j: walls8f(t)[j], th) for j in S]

def project_null(v, J):
    """v 를 rows(J) 의 영공간에 사영.
    Project v onto the null space of rows(J)."""
    k = len(J)
    if k == 0: return list(v)
    JJT = [[sum(J[i][m]*J[j][m] for m in range(5)) for j in range(k)] for i in range(k)]
    rhs = [sum(J[i][m]*v[m] for m in range(5)) for i in range(k)]
    y = solve_k(JJT, rhs)
    if y is None: return None
    return [v[m] - sum(J[i][m]*y[i] for i in range(k)) for m in range(5)]

def newton_polish(th, S, iters=14):
    """{w_j = 0 : j in S} 로 최소노름 뉴턴 사영.
    Minimum-norm Newton projection onto {w_j = 0 : j in S}."""
    if not S: return th, True
    for _ in range(iters):
        try: ws = walls8f(th)
        except Exception: return th, False
        v = [ws[j] for j in S]
        if max(abs(z) for z in v) < 1e-13: return th, True
        J = jacS(th, S)
        k = len(S)
        JJT = [[sum(J[i][m]*J[j][m] for m in range(5)) for j in range(k)] for i in range(k)]
        y = solve_k(JJT, v)
        if y is None: return th, False
        dth = [sum(J[i][m]*y[i] for i in range(k)) for m in range(5)]
        th = tuple(th[m] - dth[m] for m in range(5))
    try: ws = walls8f(th)
    except Exception: return th, False
    return th, max(abs(ws[j]) for j in S) < 1e-10

def active_set_names(ws, tol):
    S = [WNAMES[j] for j in range(6) if abs(ws[j]) < tol]
    if abs(max(ws[6], ws[7])) < tol: S.append('taco')
    return S

# ---------- 0. 시작점 생성 ---------- — generating the start points
def gen_starts(rng):
    starts, log = [], []
    for x in SLICE_XS:
        base = slice_theta(x)
        okb, _, _, eb = fully_valid(base)
        log.append(f"  슬라이스 x={x:.2f}: 기준점 유효={okb}" +
                   (f", ell={eb:.6f}" if okb else ""))
        for m in MAGS:
            kept, tried = 0, 0
            while kept < KEEP_CELL and tried < MAX_TRY:
                tried += 1
                th0 = tuple(base[k] + rng.uniform(-m, m) for k in range(5))
                ok, _, _, e0 = fully_valid(th0)
                if ok:
                    starts.append(dict(x=x, m=m, theta0=th0, ell0=e0))
                    kept += 1
            log.append(f"    m={m:.2f}: {kept}/{KEEP_CELL} 채택 ({tried}회 시도)")
    return starts, log

# ---------- 1. 무작위 실현가능 하강 (탐사) ---------- — random feasible descent (exploration)
def descend_random(th0, e0, budget_s, rng):
    th, prev = th0, e0
    step, fails, accepted = STEP0, 0, 0
    t0 = time.time()
    while accepted < MAX_STEPS:
        if time.time() - t0 > budget_s:
            return th, prev, accepted, 'budget'
        d = [rng.gauss(0.0, 1.0) for _ in range(5)]
        n = math.sqrt(sum(z*z for z in d))
        d = [z/n for z in d]
        th_try = tuple(th[k] + step*d[k] for k in range(5))
        ok, c, ws, e = feas(th_try)
        if ok and e < prev - 1e-15:
            okv, _why = valid_full(c)          # 매 수락 스텝마다 유효성 검사 — a validity check at every accepted step
            if okv:
                th, prev = th_try, e
                accepted += 1; fails = 0
                step = min(step*1.3, STEP_MAX)
                continue
        fails += 1
        if fails >= FAILS_LVL:
            step *= 0.5; fails = 0
            if step < STEP_MIN:
                return th, prev, accepted, 'converged'
    return th, prev, accepted, 'maxsteps'

# ---------- 2. 활성면 사영 하강 (정련) ---------- — projected descent on the active face (refinement)
def descend_face(th, budget_s, max_iter=3000):
    t0 = time.time()
    ok, c, ws, e = feas(th)
    if not ok: return th, ell_fn(th), [], 'infeas0'
    S = [j for j in range(6) if abs(ws[j]) < 1e-4]
    thp, okp = newton_polish(th, S)
    okv, cp, wsp, ep_ = fully_valid(thp)
    if okp and okv:
        th, ws, e = thp, wsp, ep_
    else:
        S = [j for j in range(6) if abs(ws[j]) < 1e-6]   # 보수적 재시도 — a conservative retry
        thp, okp = newton_polish(th, S)
        okv, cp, wsp, ep_ = fully_valid(thp)
        if okp and okv: th, ws, e = thp, wsp, ep_
        else: return th, e, S, 'polish-fail'
    alpha0 = 2e-3
    status = 'maxiter'
    for it in range(max_iter):
        if time.time() - t0 > budget_s: status = 'budget'; break
        if len(S) >= 5: status = 'vertex'; break
        J = jacS(th, S)
        g = gradv(ell_fn, th)
        gp = project_null(g, J)
        if gp is None: status = 'singular'; break
        npg = math.sqrt(sum(z*z for z in gp))
        if npg < PG_TOL: status = 'kkt'; break
        dirv = [-z/npg for z in gp]
        alpha, moved = alpha0, False
        while alpha >= 1e-8:
            cand0 = tuple(th[m] + alpha*dirv[m] for m in range(5))
            cand, okp = newton_polish(cand0, S)
            if okp:
                ok2, c2, ws2, e2 = feas(cand)
                if ok2 and e2 < e - 1e-15 and valid_full(c2)[0]:
                    th, ws, e = cand, ws2, e2
                    alpha0 = min(alpha*2.0, 1e-2)
                    moved = True
                    break
            alpha *= 0.5
        if moved:
            # 새로 활성화된 벽 편입 — absorb the newly activated walls
            newS = sorted(set(S) | {j for j in range(6) if abs(ws[j]) < 1e-9})
            if newS != S:
                thp, okp = newton_polish(th, newS)
                okv2, _, wsp2, ep2 = fully_valid(thp)
                if okp and okv2:
                    S, th, ws, e = newS, thp, wsp2, ep2
            continue
        # 백트래킹 전멸: 눌린 벽을 활성 집합에 편입 시도 — backtracking exhausted: try absorbing the pressed walls into the active set
        press = sorted(set(S) | {j for j in range(6) if abs(ws[j]) < 1e-5})
        if press != S and len(press) <= 5:
            thp, okp = newton_polish(th, press)
            okv2, _, wsp2, ep2 = fully_valid(thp)
            if okp and okv2 and ep2 <= e + 1e-12:
                S, th, ws, e = press, thp, wsp2, ep2
                alpha0 = 2e-3
                continue
        if max(ws[6], ws[7]) < 1e-5:
            status = 'taco-floor'; break
        status = 'floor'; break
    return th, e, S, status

# ---------- 4. 합법성 감사 (check_branch 방식, 일반 theta) ---------- — legitimacy audit (the check_branch approach, for a general theta)
def meet_params(P1, d1, P2, d2):
    det = d1[0]*(-d2[1]) - d1[1]*(-d2[0])
    if abs(det) < 1e-15: return None, None
    rhs = sub(P2, P1)
    s = (rhs[0]*(-d2[1]) - rhs[1]*(-d2[0]))/det
    u = (d1[0]*rhs[1] - d1[1]*rhs[0])/det
    return s, u

def audit(th, N=199):
    """B1 조각 전 쌍별 교차 + B2 순서 + B3 크리스쌍 + B4 199점 공리.
    B1 crossings over all pairs of pieces + B2 order + B3 crease pairs + B4 the 199-point axiom."""
    res = dict(ok=True, detail=[])
    def chk(label, cond, extra=""):
        res['detail'].append(f"[{'OK' if cond else 'FAIL'}] {label} {extra}".strip())
        if not cond: res['ok'] = False
    c = config(th)
    if c['M3'] is None or c['s23'] is None or c['s23'] <= 0:
        res['ok'] = False; res['detail'].append("[FAIL] M3/s23 퇴화"); return res
    INF = 1e6
    pieces = {'A': (c['M1'], mul(c['d2'], -1.0), 0.0, INF),
              'B': (c['M1'], unit(sub(c['M2'], c['M1'])), 0.0, norm(sub(c['M2'], c['M1']))),
              'C': (c['M2'], c['w0'], 0.0, c['s23']),
              'D': (c['M3'], c['v3'], 0.0, INF)}
    expect1 = {('A','C'), ('A','D'), ('B','D')}
    names = ['A','B','C','D']; crossings = {}
    for i in range(4):
        for j in range(i+1, 4):
            n1, n2 = names[i], names[j]
            P1, dd1, lo1, hi1 = pieces[n1]; P2, dd2, lo2, hi2 = pieces[n2]
            s, u = meet_params(P1, dd1, P2, dd2)
            cnt = 1 if (s is not None and lo1+1e-9 < s < hi1-1e-9
                        and lo2+1e-9 < u < hi2-1e-9) else 0
            crossings[(n1, n2)] = (cnt, s, u)
            want = 1 if (n1, n2) in expect1 else 0
            chk(f"B1 {n1}x{n2} 교차 {cnt} (기대 {want})", cnt == want)
    if res['ok']:
        sAC = crossings[('A','C')][1]; sAD = crossings[('A','D')][1]
        tDA = crossings[('A','D')][2]; tDB = crossings[('B','D')][2]
        chk("B2 A에서 xC 먼저", sAC > sAD, f"sAC={sAC:.4f} sAD={sAD:.4f}")
        chk("B2 D에서 xA 먼저", tDA < tDB, f"tDA={tDA:.4f} tDB={tDB:.4f}")
    creases = {'c1': (c['E'], c['F']), 'c2': (c['H'], c['G']), 'c3': (c['Q1'], c['Q0'])}
    for (na, nb) in (('c1','c2'), ('c2','c3'), ('c1','c3')):
        A1, A2 = creases[na]; B1, B2 = creases[nb]
        s, u = meet_params(A1, sub(A2, A1), B1, sub(B2, B1))
        inside = (s is not None and 1e-9 < s < 1-1e-9 and 1e-9 < u < 1-1e-9)
        chk(f"B3 {na}x{nb} 내부교차 없음", not inside)
    viol = []
    for i in range(1, N+1):
        t = i/(N+1)
        p1 = add(c['E'], mul(sub(c['F'], c['E']), t))
        p2 = add(c['H'], mul(sub(c['G'], c['H']), t))
        p3 = add(c['Q1'], mul(sub(c['Q0'], c['Q1']), t))
        if covers_D(c, p1): viol.append(f"D@c1({t:.3f})")
        if covers_A(c, p2): viol.append(f"A@c2({t:.3f})")
        if covers_D(c, p2): viol.append(f"D@c2({t:.3f})")
        if covers_A(c, p3): viol.append(f"A@c3({t:.3f})")
    chk(f"B4 크리스 이웃 공리 {N}점", not viol, str(viol[:4]))
    return res

# ---------- 메인 ---------- — main
def main():
    rng = random.Random(SEED)
    t_all = time.time()
    print(f"=== 전역 유일성 스윕 (seed={SEED}, 예산: 1단계 {BUDGET1:.0f}s + 2단계 {BUDGET2:.0f}s) ===")
    print(f"theta† = {tuple(round(v,6) for v in TH_DAG)}, ell† = {ELL_DAG:.9f}")
    print(f"오각형 = {tuple(round(v,6) for v in TH_PENT)}, ell = {ell_fn(TH_PENT):.9f}")

    print("\n--- 0. 시작점 생성 ---")
    starts, slog = gen_starts(rng)
    for L in slog: print(L)
    print(f"총 시작점: {len(starts)}")

    print("\n--- 1. 무작위 실현가능 하강 (탐사) ---")
    mids = []
    t1 = time.time()
    for i, st in enumerate(starts):
        remain = BUDGET1 - (time.time() - t1)
        per = max(3.0, remain/max(1, len(starts)-i))
        th, e, acc, status = descend_random(st['theta0'], st['ell0'], per, rng)
        _, _, ws, _ = feas(th)
        mids.append(dict(idx=i, x=st['x'], m=st['m'], theta1=th, ell1=e,
                         steps=acc, status1=status,
                         act1=active_set_names(ws, 1e-4) if ws else []))
    print(f"  1단계 완료: {len(mids)}개, 경과 {time.time()-t1:.0f}s, "
          f"상태 분포: { {s: sum(1 for m in mids if m['status1']==s) for s in set(m['status1'] for m in mids)} }")
    from collections import Counter
    print(f"  1단계 정체점 활성집합(1e-4) 분포: {dict(Counter(tuple(m['act1']) for m in mids))}")

    print("\n--- 2. 활성면 사영 하강 (정련) ---")
    endpoints = []
    t2 = time.time()
    for i, md in enumerate(mids):
        remain = BUDGET2 - (time.time() - t2)
        per = max(3.0, remain/max(1, len(mids)-i))
        th, e, S, status = descend_face(md['theta1'], per)
        ok, c, ws, _ = feas(th)
        ep = dict(**md, theta=th, ell=e, S=[WNAMES[j] for j in S], status2=status,
                  walls=ws, act=active_set_names(ws, ACT_TOL) if ws else [])
        endpoints.append(ep)
        if (i+1) % 12 == 0 or i == len(mids)-1:
            print(f"  [{i+1:3d}/{len(mids)}] 최근: ell={e:.9f} ({status}) "
                  f"S={ep['S']} d(theta†)={dmax(th,TH_DAG):.2e} 경과={time.time()-t2:.0f}s")
    print(f"  2단계 상태 분포: { {s: sum(1 for m in endpoints if m['status2']==s) for s in set(m['status2'] for m in endpoints)} }")

    print("\n--- 3. 클러스터링 ---")
    clusters = []
    for ep in sorted(endpoints, key=lambda z: z['ell']):
        for cl in clusters:
            if dmax(ep['theta'], cl['rep']['theta']) < CLUST_TOL and \
               abs(ep['ell'] - cl['rep']['ell']) < 1e-3:
                cl['members'].append(ep); break
        else:
            clusters.append(dict(rep=ep, members=[ep]))
    print(f"클러스터 수: {len(clusters)}")

    table = []
    for ci, cl in enumerate(clusters):
        rep = cl['rep']; th = rep['theta']
        row = dict(cluster=ci, count=len(cl['members']),
                   theta=list(th), ell=rep['ell'],
                   act=rep['act'], S=rep['S'],
                   d_dag=dmax(th, TH_DAG), d_pent=dmax(th, TH_PENT),
                   st2={s: sum(1 for m in cl['members'] if m['status2']==s)
                        for s in set(m['status2'] for m in cl['members'])},
                   spread=max(dmax(m['theta'], th) for m in cl['members']))
        table.append(row)

    print("\n--- 4. 보고: 클러스터 표 ---")
    hdr = (f"{'#':>2} {'도달':>4} {'ell':>14} {'d(theta†)':>10} {'d(오각형)':>10} "
           f"{'활성(1e-6)':<16} {'2단계상태':<22} {'퍼짐':>8}  대표 theta")
    print(hdr); print('-'*len(hdr))
    for r in table:
        print(f"{r['cluster']:>2} {r['count']:>4} {r['ell']:>14.9f} "
              f"{r['d_dag']:>10.2e} {r['d_pent']:>10.2e} "
              f"{','.join(r['act']) or '-':<16} {str(r['st2']):<22} {r['spread']:>8.1e}  "
              f"({', '.join(f'{v:.6f}' for v in r['theta'])})")

    print("\n--- 5. 판정 ---")
    new_minima = []
    for r in table:
        kkt = r['st2'].get('kkt', 0)
        if r['d_dag'] < 1e-3:
            print(f"클러스터 {r['cluster']}: theta† 도달 ({r['count']}회, 그중 면-KKT {kkt})")
        elif r['d_pent'] < 2e-3:
            print(f"클러스터 {r['cluster']}: 오각형 코너 근방 ({r['count']}회) — "
                  f"알려진 안장(정리 7.2), 알고리즘 정체로 해석")
        elif kkt == 0:
            print(f"클러스터 {r['cluster']}: 비-KKT 정체 ({r['count']}회, {r['st2']}) — "
                  f"최소 후보 아님 (예산/스텝/특이 정지)")
        else:
            print(f"클러스터 {r['cluster']}: *** theta† 외 면-KKT 종점 *** ell={r['ell']:.9f} — 감사 수행")
            new_minima.append(r)
    if not new_minima:
        print("=> theta† 외의 합법적 국소 최소 후보 없음 (이번 스윕 범위 내)")

    for r in new_minima:
        print(f"\n--- 6. 신규 후보 감사 (클러스터 {r['cluster']}) ---")
        res = audit(tuple(r['theta']))
        for L in res['detail']: print("  " + L)
        print(f"  감사 종합: {'전부 통과 — 합법적 신규 타이트 후보!' if res['ok'] else '실패 있음 — 비합법/가지 이탈'}")
        r['audit_ok'] = res['ok']; r['audit'] = res['detail']

    out = dict(seed=SEED, budget1=BUDGET1, budget2=BUDGET2, n_starts=len(starts),
               step0=STEP0, step_min=STEP_MIN, step_max=STEP_MAX,
               fails_per_level=FAILS_LVL, pg_tol=PG_TOL,
               th_dag=list(TH_DAG), ell_dag=ELL_DAG,
               starts=[dict(x=s['x'], m=s['m'], theta0=list(s['theta0']),
                            ell0=s['ell0']) for s in starts],
               endpoints=[dict(idx=e['idx'], x=e['x'], m=e['m'],
                               theta1=list(e['theta1']), ell1=e['ell1'],
                               steps=e['steps'], status1=e['status1'], act1=e['act1'],
                               theta=list(e['theta']), ell=e['ell'],
                               S=e['S'], status2=e['status2'], act=e['act'])
                          for e in endpoints],
               clusters=table)
    outpath = os.environ.get('OUT', os.path.join(HERE, 'global_sweep_results.json'))
    with open(outpath, 'w') as f: json.dump(out, f, indent=1)
    print(f"\n결과 JSON: {outpath}")
    print(f"총 소요: {time.time()-t_all:.0f}s")
    print(f"\n[조절 내역] 시작점 셀당 {KEEP_CELL}개(6x3셀). 1단계 무작위 하강은 명세 스텝 스케줄"
          f"(2e-3, 120연속실패시 절반, 하한 1e-5)에 성공 시 1.3배 재성장(상한 1e-2)을 추가 — "
          f"원거리 시작점(x=0.10, ell~9.9)의 시간 내 수렴용. 파일럿에서 무작위 하강이 "
          f"4-등식 활성면 위에서 조기 정체함을 관찰, 활성면 사영 하강(2단계)을 추가해 "
          f"면-KKT 점까지 정련. 예산: 1단계 {BUDGET1:.0f}s, 2단계 {BUDGET2:.0f}s "
          f"(시작점별 잔여시간/잔여시작점 균등 분배, 최소 3s).")

if __name__ == '__main__':
    main()
