# uniq_strata.py — 전역 유일성(16) 2단계: 활성층 전수 KKT 해집(harvest) — Step 2 of global uniqueness (16): exhaustive KKT harvest over the active strata
#
# 대상: 층 S ⊆ {1,4,5,6,7,8} × taco ∈ {∅,{9},{10},{9,10}}, |S|≤5 (LICQ 상한). — Scope: strata S ⊆ {1,4,5,6,7,8} × taco ∈ {∅,{9},{10},{9,10}}, |S|≤5 (LICQ upper bound).
# 각 층의 정사각 다항계  { W_j = 0 (j∈S),  ∇ell − Σ λ_j ∇W_j = 0 (5식) } — The square polynomial system of each stratum { W_j = 0 (j∈S), ∇ell − Σ λ_j ∇W_j = 0 (5 equations) }
# 를 다중시작 뉴턴으로 전수 수확한 뒤 필터: — is harvested exhaustively by multistart Newton, then filtered:
#   (F1) 나머지 벽 실현측:  w_i ≥ −tol (i∉S, i≤8),  taco는 선언적 max(w9,w10) ≥ −tol — (F1) realizable side of the remaining walls: w_i ≥ −tol (i∉S, i≤8); taco is declarative, max(w9,w10) ≥ −tol
#   (F2) KKT 부호: λ_j ≥ −tol — (F2) KKT sign: λ_j ≥ −tol
#   (F3) 유효성(valid_full: 전방성, 재료 정합, 크리스 이웃 공리 29점) — (F3) validity (valid_full: forwardness, material consistency, crease-neighborhood axiom at 29 points)
# 통과점 = 도달가능 타이트 후보. 기대: {1,4,5,6}의 θ† (±거울)뿐. — Passing points = accessible tight candidates. Expected: only θ† of {1,4,5,6} (± mirror).
#
# 부수 스캔: ≥6벽 동시 영점(과결정 코너) 탐색 — 오각형(8벽) 외 존재 여부. — Side scan: search for ≥6 walls vanishing simultaneously (overdetermined corners) — whether any exists besides the pentagon (8 walls).

import math, random, itertools, sys, time
import sympy as sp
import uniq_setup as U

x, g, h, a, b = U.TH
NAMES = ["W1", "W4", "W5", "W6", "W7", "W8", "W9", "W10"]
WALLS = [U.WALLS[n] for n in NAMES]      # 부호 고정판 (실현측 ≥ 0) — sign-fixed version (realizable side ≥ 0)
ELL = U.ELL

# ---- 심볼릭 기울기 → 수치 함수 (simplify 금지: 원시 트리 그대로) ---- — symbolic gradients → numeric functions (no simplify: raw expression tree as is)
dELL = [sp.diff(ELL, v) for v in U.TH]
dW = [[sp.diff(w_, v) for v in U.TH] for w_ in WALLS]
f_W = [sp.lambdify(U.TH, w_, 'math', cse=True) for w_ in WALLS]
f_dE = [sp.lambdify(U.TH, e_, 'math', cse=True) for e_ in dELL]
f_dW = [[sp.lambdify(U.TH, e_, 'math', cse=True) for e_ in row] for row in dW]
f_ELL = sp.lambdify(U.TH, ELL, 'math', cse=True)

# ---- 유효성 (check_edge.valid_full 이식; uniq_setup의 수치 헬퍼 재사용) ---- — validity (ported from check_edge.valid_full; reuses the numeric helpers of uniq_setup)
_un, _cr, _sub, _add, _mul, _dot = U._unit, U._cross, U._sub, U._add, U._mul, U._dot

def covers_A(c, p):
    du = c['d2']
    if _cr(_sub(p, c['E']), du) * _cr(_sub(p, c['F']), du) > 1e-12: return False
    pa = _sub(c['M1'], _mul(du, 1.0))
    return _cr(_sub(c['F'], c['E']), _sub(p, c['E'])) * _cr(_sub(c['F'], c['E']), _sub(pa, c['E'])) > 0

def covers_D(c, p):
    e1p = c['rho3'](c['G']); e0p = c['rho3'](c['H']); dD = c['rl3'](c['w0'])
    if _cr(_sub(p, e1p), _un(dD)) * _cr(_sub(p, e0p), _un(dD)) > 1e-12: return False
    pd = _add(c['M3'], _mul(c['v3'], 1.0))
    s = _cr(_sub(c['Q0'], c['Q1']), _sub(p, c['Q1']))
    sD = _cr(_sub(c['Q0'], c['Q1']), _sub(pd, c['Q1']))
    return s * sD > 0

def _meet(P1, d1, P2, d2v):
    det = d1[0]*(-d2v[1]) - d1[1]*(-d2v[0])
    if abs(det) < 1e-14: return None, None
    rhs = _sub(P2, P1)
    s = (rhs[0]*(-d2v[1]) - rhs[1]*(-d2v[0]))/det
    return _add(P1, _mul(d1, s)), s

def valid_full(c, nsamp=29):
    if c['M3'] is None or c['s23'] is None or c['s23'] <= 0: return False, "M3"
    if c['M2'][0] <= c['M1'][0]: return False, "M2<M1"
    _, sAC = _meet(c['M1'], _mul(c['d2'], -1.0), c['M2'], c['w0'])
    _, sAD = _meet(c['M1'], _mul(c['d2'], -1.0), c['M3'], c['v3'])
    _, tDA = _meet(c['M3'], c['v3'], c['M1'], _mul(c['d2'], -1.0))
    _, tDB = _meet(c['M3'], c['v3'], (0.0, 0.5), (1.0, 0.0))
    if None in (sAC, sAD, tDA, tDB): return False, "extent-none"
    if not (sAC > 0 and sAD > 0 and tDA > 0 and tDB > 0): return False, "extent-sign"
    xB = c['M3'][0] + tDB*c['v3'][0]
    if not (c['M1'][0] < xB < c['M2'][0]): return False, "BxD"
    _, uCA = _meet(c['M2'], c['w0'], c['M1'], _mul(c['d2'], -1.0))
    if uCA is None or not (0 < uCA < c['s23']): return False, "CxA"
    for i in range(1, nsamp + 1):
        t = i/(nsamp + 1)
        p1 = _add(c['E'], _mul(_sub(c['F'], c['E']), t))
        p2 = _add(c['H'], _mul(_sub(c['G'], c['H']), t))
        p3 = _add(c['Q1'], _mul(_sub(c['Q0'], c['Q1']), t))
        if covers_D(c, p1): return False, "ax-c1"
        if covers_A(c, p2) or covers_D(c, p2): return False, "ax-c2"
        if covers_A(c, p3): return False, "ax-c3"
    return True, ""

# ---- 뉴턴 (유한차분 야코비, float) ---- — Newton (finite-difference Jacobian, float)
def newton(Ffun, z0, iters=30, tol=1e-11):
    z = list(z0); n = len(z)
    for _ in range(iters):
        try:
            Fv = Ffun(z)
        except Exception:
            return None
        if any(math.isnan(v) or math.isinf(v) for v in Fv): return None
        res = max(abs(v) for v in Fv)
        if res < tol: return z
        J = []
        for i in range(n):
            e = 1e-7 * max(1.0, abs(z[i]))
            zp = z[:]; zp[i] += e
            try:
                Fp = Ffun(zp)
            except Exception:
                return None
            J.append([(Fp[k] - Fv[k])/e for k in range(n)])
        # solve J^T dz = -F  (행 우선 구성이라 전치 주의: J[i][k] = dF_k/dz_i) — built row-major, so mind the transpose: J[i][k] = dF_k/dz_i
        M = [[J[i][k] for i in range(n)] + [-Fv[k]] for k in range(n)]
        for col in range(n):
            p = max(range(col, n), key=lambda r0: abs(M[r0][col]))
            if abs(M[p][col]) < 1e-14: return None
            M[col], M[p] = M[p], M[col]
            piv = M[col][col]
            for r0 in range(n):
                if r0 != col:
                    f0 = M[r0][col]/piv
                    M[r0] = [M[r0][kk] - f0*M[col][kk] for kk in range(n + 1)]
        dz = [M[i][n]/M[i][i] for i in range(n)]
        step = 1.0
        nrm = max(abs(v) for v in dz)
        if nrm > 1.0: step = 1.0/nrm          # 감쇠 — damping
        z = [z[i] + step*dz[i] for i in range(n)]
        if max(abs(v) for v in z) > 50: return None
    return None

def make_system(S):
    k = len(S)
    def Ffun(z):
        th = tuple(z[:5]); lam = z[5:]
        out = [f_W[j](*th) for j in S]
        for comp in range(5):
            v = f_dE[comp](*th)
            for i, j in enumerate(S):
                v -= lam[i]*f_dW[j][comp](*th)
            out.append(v)
        return out
    return Ffun

def lam_seed(S, th):
    # 최소자승 λ 시드: ∇W_S^T λ ≈ ∇ell — least-squares λ seed: ∇W_S^T λ ≈ ∇ell
    k = len(S)
    A = [[f_dW[j][comp](*th) for j in S] for comp in range(5)]
    bv = [f_dE[comp](*th) for comp in range(5)]
    # normal equations k×k
    N = [[sum(A[c][i]*A[c][jj] for c in range(5)) for jj in range(k)] for i in range(k)]
    y = [sum(A[c][i]*bv[c] for c in range(5)) for i in range(k)]
    M = [N[i][:] + [y[i]] for i in range(k)]
    for col in range(k):
        p = max(range(col, k), key=lambda r0: abs(M[r0][col]))
        if abs(M[p][col]) < 1e-12: return [0.5]*k
        M[col], M[p] = M[p], M[col]
        for r0 in range(k):
            if r0 != col:
                f0 = M[r0][col]/M[col][col]
                M[r0] = [M[r0][kk] - f0*M[col][kk] for kk in range(k + 1)]
    return [M[i][k]/M[i][i] for i in range(k)]

# ---- 시드 풀 ---- — seed pool
rnd = random.Random(20260803)
def seeds(nn=40):
    out = []
    for _ in range(nn):
        xx = rnd.uniform(0.10, 0.56)
        th0 = U._slice_theta(xx)
        out.append(tuple(v + rnd.uniform(-0.25, 0.25) for v in th0))
    return out

TOL = 1e-7
IDX = {n_: i for i, n_ in enumerate(NAMES)}

def check_point(th, lam, S):
    """필터 F1–F3. 반환: (통과여부, 사유)
    Filters F1–F3. Returns: (passed or not, reason)"""
    try:
        ws, c = U._walls8f(th)
    except Exception:
        return False, "config"
    for i in range(6):                     # w1..w8 (수치 배열 순서: w1,w4,w5,w6,w7,w8) — (order of the numeric array: w1,w4,w5,w6,w7,w8)
        if i not in [ [0,1,2,3,4,5][k] for k in range(6)]: pass
    numeric_names = ["W1","W4","W5","W6","W7","W8","W9","W10"]
    for i, n_ in enumerate(numeric_names[:6]):
        if IDX[n_] not in S and ws[i] < -TOL:
            return False, f"{n_}<0"
    # taco 선언적 — taco, declarative
    t9a, t10a = (IDX["W9"] in S), (IDX["W10"] in S)
    if not (t9a or t10a):
        if max(ws[6], ws[7]) < -TOL: return False, "taco"
    elif t9a and not t10a:
        if abs(ws[6]) > 1e-5: return False, "W9!=0"
    elif t10a and not t9a:
        if abs(ws[7]) > 1e-5: return False, "W10!=0"
    if any(l < -TOL for l in lam): return False, "λ<0"
    ok, why = valid_full(c)
    if not ok: return False, f"valid:{why}"
    return True, ""

def strata_list():
    base = list(range(6))                 # 인덱스 0..5 = W1,W4,W5,W6,W7,W8 — indices 0..5 = W1,W4,W5,W6,W7,W8
    taco_opts = [(), (6,), (7,), (6, 7)]  # 6=W9, 7=W10
    out = []
    for kbase in range(0, 6 + 1):
        for Sb in itertools.combinations(base, kbase):
            for tc in taco_opts:
                S = tuple(sorted(Sb + tc))
                if len(S) == 0 or len(S) > 5:      # |S|=0: 명제 6.1 (내부 임계점 부재) — Proposition 6.1 (no interior critical point)
                    continue
                if set(S) <= {1, 3}:               # 보조정리 A: S ⊆ {W4,W6} 은 정확 배제 — Lemma A: S ⊆ {W4,W6} is excluded exactly
                    continue
                out.append(S)
    return out

def harvest_one(S):
    Ffun = make_system(S)
    found = {}
    for th0 in SEEDPOOL:
        lam0 = lam_seed(S, th0)
        z = newton(Ffun, list(th0) + lam0)
        if z is None: continue
        key = tuple(round(v, 6) for v in z[:5])
        if key not in found:
            found[key] = (tuple(z[:5]), list(z[5:]))
    return S, list(found.values())

SEEDPOOL = seeds(30)

def run():
    import multiprocessing as mp_
    t0 = time.time()
    all_pass, all_sols = [], 0
    strata_with_sols = 0
    report = []
    SL = strata_list()
    print(f"층 {len(SL)}개 (보조정리 A 배제 후), 시드 {len(SEEDPOOL)}개/층, 8-way")
    with mp_.get_context("fork").Pool(8) as pool:
        for S, sols in pool.imap_unordered(harvest_one, SL):
            if not sols: continue
            strata_with_sols += 1
            for th, lam in sols:
                all_sols += 1
                ok, why = check_point(th, lam, S)
                tag = "PASS" if ok else why
                if ok or why.startswith("λ") or why.startswith("valid"):
                    report.append((S, th, lam, f_ELL(*th), tag))
                if ok:
                    all_pass.append((S, th, lam))
    print(f"층 수확 완료: 방정식 해 {all_sols}개 (해 보유 층 {strata_with_sols}개), 경과 {time.time()-t0:.0f}s")
    print("\n관심 해 (PASS 또는 λ/유효성 필터 근처 탈락):")
    for S, th, lam, e_, tag in report:
        Sn = ",".join(NAMES[j] for j in S)
        print(f"  S={{{Sn}}}: ell={e_:.7f}, x={th[0]:.6f}, tag={tag}")
        print(f"      θ=({', '.join(f'{v:.6f}' for v in th)}), λ=({', '.join(f'{v:.4f}' for v in lam)})")
    print(f"\n전 필터 통과(도달가능 타이트 후보): {len(all_pass)}개")
    for S, th, lam in all_pass:
        Sn = ",".join(NAMES[j] for j in S)
        print(f"  ★ S={{{Sn}}}: θ=({', '.join(f'{v:.8f}' for v in th)}), ell={f_ELL(*th):.9f}")

    # ---- 과결정 코너 스캔: ≥6벽 동시 영점 ---- — overdetermined-corner scan: ≥6 walls vanishing simultaneously
    print("\n[코너 스캔] Σ w_j² 최소화 (8벽 전부):")
    def sumsq(th):
        try:
            ws, _ = U._walls8f(th)
        except Exception:
            return None
        return sum(v*v for v in ws[:6]) + min(ws[6], ws[7])**2 * 0  # taco는 코너 판단 별도 — taco is judged as a corner separately
    corners = {}
    for th0 in seeds(60):
        th = list(th0)
        for _ in range(60):     # 좌표 순환 국소 최소화 (거친 스캔) — cyclic-coordinate local minimization (coarse scan)
            improved = False
            s0 = sumsq(tuple(th))
            if s0 is None: break
            for i in range(5):
                for d in (1e-3, -1e-3, 1e-4, -1e-4):
                    th2 = th[:]; th2[i] += d
                    s2 = sumsq(tuple(th2))
                    if s2 is not None and s2 < s0:
                        th, s0, improved = th2, s2, True
            if not improved: break
        s_f = sumsq(tuple(th))
        if s_f is not None and s_f < 1e-8:
            key = tuple(round(v, 4) for v in th)
            corners[key] = s_f
    print(f"  6벽(슬리버) 동시 영점 후보: {len(corners)}개")
    for key in corners:
        print(f"    θ≈{key}  (오각형 x0=0.324920 와 비교)")

if __name__ == "__main__":
    run()
