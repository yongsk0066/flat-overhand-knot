# uniq_census.py — 전역 유일성(16) 3단계: 다항계 기반 **전수 센서스** (독립 2차 구현) — Step 3 of global uniqueness (16): **exhaustive census** based on the polynomial system (independent second implementation)
#
# uniq_strata.py(megatree 유리계)와 독립으로, uniq_poly.py 의 구조형 벽 + 유리 ell 로 — Independently of uniq_strata.py (the megatree rational system), using the structural walls of uniq_poly.py + the rational ell,
# 같은 215개 층의 KKT 정사각계를 100 시드 다중시작 수확하고, 모든 해를 분류해 — this harvests the KKT square systems of the same 215 strata by 100-seed multistart, classifies every solution,
# uniq_census.json 으로 덤프한다. 두 구현의 일치 = 수확의 교차 검증. — and dumps them to uniq_census.json. Agreement of the two implementations = cross-validation of the harvest.
#
# 추가: 벽-전용 5-부분집합 C(8,5)=56 정사각계로 과결정 코너(≥6벽 동시 영점) 전수. — Additionally: an exhaustive sweep for overdetermined corners (≥6 walls vanishing simultaneously) via the wall-only 5-subset C(8,5)=56 square systems.

import math, random, itertools, json, time
import sympy as sp
import uniq_poly as P
import uniq_setup as U0

TH = U0.TH
NAMES = P.NAMES
STRUCT = [P.STRUCT[n] for n in NAMES]
Np, Dp = P.Np, P.Dp

# ---- 수치 함수 ---- — numeric functions
f_W = [sp.lambdify(TH, e, 'math', cse=True) for e in STRUCT]
dW = [[sp.diff(e, v) for v in TH] for e in STRUCT]
f_dW = [[sp.lambdify(TH, e, 'math', cse=True) for e in row] for row in dW]
dNp = [sp.diff(Np, v) for v in TH]
dDp = [sp.diff(Dp, v) for v in TH]
f_Np = sp.lambdify(TH, Np, 'math', cse=True)
f_Dp = sp.lambdify(TH, Dp, 'math', cse=True)
f_dNp = [sp.lambdify(TH, e, 'math', cse=True) for e in dNp]
f_dDp = [sp.lambdify(TH, e, 'math', cse=True) for e in dDp]
LIN = [0.5, 0.5, 0.5, 0.0, 0.0]   # ∇(g+h+x)/2

def grad_ell(th):
    Npv, Dpv = f_Np(*th), f_Dp(*th)
    return [LIN[k] + (f_dNp[k](*th)*Dpv - Npv*f_dDp[k](*th))/Dpv**2 for k in range(5)]

def ell(th):
    return (th[0] + th[1] + th[2])/2 + f_Np(*th)/f_Dp(*th)

# ---- 뉴턴 (uniq_strata와 동일 구조) ---- — Newton (same structure as in uniq_strata)
def newton(Ffun, z0, iters=30, tol=1e-11):
    z = list(z0); n = len(z)
    for _ in range(iters):
        try:
            Fv = Ffun(z)
        except Exception:
            return None
        if any(math.isnan(v) or math.isinf(v) for v in Fv): return None
        if max(abs(v) for v in Fv) < tol: return z
        J = []
        for i in range(n):
            e = 1e-7*max(1.0, abs(z[i]))
            zp = z[:]; zp[i] += e
            try:
                Fp = Ffun(zp)
            except Exception:
                return None
            J.append([(Fp[k] - Fv[k])/e for k in range(n)])
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
        st = 1.0
        nrm = max(abs(v) for v in dz)
        if nrm > 1.0: st = 1.0/nrm
        z = [z[i] + st*dz[i] for i in range(n)]
        if max(abs(v) for v in z) > 50: return None
    return None

def make_kkt(S):
    def Ffun(z):
        th = tuple(z[:5]); lam = z[5:]
        out = [f_W[j](*th) for j in S]
        gE = grad_ell(th)
        for k in range(5):
            v = gE[k]
            for i, j in enumerate(S):
                v -= lam[i]*f_dW[j][k](*th)
            out.append(v)
        return out
    return Ffun

def lam_seed(S, th):
    k = len(S)
    A = [[f_dW[j][c](*th) for j in S] for c in range(5)]
    bv = grad_ell(th)
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

# ---- 층 목록 (uniq_strata와 동일 규약: 인덱스 0..5=W1,W4,W5,W6,W7,W8; 6=W9,7=W10) ---- — stratum list (same convention as uniq_strata: indices 0..5=W1,W4,W5,W6,W7,W8; 6=W9, 7=W10)
def strata_list():
    out = []
    for kb in range(0, 7):
        for Sb in itertools.combinations(range(6), kb):
            for tc in [(), (6,), (7,), (6, 7)]:
                S = tuple(sorted(Sb + tc))
                if len(S) == 0 or len(S) > 5: continue
                if set(S) <= {1, 3}: continue        # 보조정리 A (W4,W6 x-독립) — Lemma A (W4, W6 are x-independent)
                out.append(S)
    return out

TH_DAG = (0.32212749731122718, 1.3911166915901648, 1.0363863004307971, 1.0363863004307971, 1.7132441870641958)
X0 = math.sqrt(1 - 2/math.sqrt(5))
TH_PEN = None  # 슬라이스에서 계산 — computed on the slice

def classify(th, lam, S):
    ws, c = U0._walls8f(th)
    d_dag = max(abs(th[i] - TH_DAG[i]) for i in range(5))
    d_pen = max(abs(th[i] - TH_PEN[i]) for i in range(5))
    kind = "theta_dagger" if d_dag < 1e-4 else ("pentagon" if d_pen < 1e-4 else "other")
    fails = []
    for i, n_ in enumerate(NAMES[:6]):
        if NAMES.index(n_) not in S and ws[i] < -1e-7: fails.append(f"{n_}<0({ws[i]:.2e})")
    t9a, t10a = (6 in S), (7 in S)
    if not (t9a or t10a) and max(ws[6], ws[7]) < -1e-7: fails.append(f"taco({ws[6]:.2e},{ws[7]:.2e})")
    negl = [i for i, l in enumerate(lam) if l < -1e-7]
    if negl: fails.append("λ<0@" + ",".join(NAMES[S[i]] for i in negl))
    import uniq_strata as ST
    ok, why = ST.valid_full(c)
    if not ok: fails.append(f"valid:{why}")
    return kind, fails

SEEDS = None
def _harvest_one(S):
    Ffun = make_kkt(S)
    found = {}
    for th0 in SEEDS:
        z = newton(Ffun, list(th0) + lam_seed(S, th0))
        if z is None: continue
        key = tuple(round(v, 6) for v in z[:5])
        found.setdefault(key, (tuple(z[:5]), list(z[5:])))
    return S, list(found.values())

def run():
    global TH_PEN
    TH_PEN = U0._slice_theta(X0)
    rnd = random.Random(41)
    seeds = []
    for _ in range(100):
        xx = rnd.uniform(0.10, 0.56)
        th0 = U0._slice_theta(xx)
        seeds.append(tuple(v + rnd.uniform(-0.3, 0.3) for v in th0))
    t0 = time.time()
    census = []
    SL = strata_list()
    global SEEDS
    SEEDS = seeds
    import multiprocessing as mp_
    with mp_.get_context("fork").Pool(8) as pool:
        for S, sols in pool.imap_unordered(_harvest_one, SL):
            for th, lam in sols:
                try:
                    kind, fails = classify(th, lam, S)
                except Exception as ex:
                    kind, fails = "config-fail", [str(ex)[:40]]
                census.append(dict(S=[NAMES[j] for j in S], theta=list(th), lam=lam,
                                   ell=ell(th), kind=kind, fails=fails))
    print(f"[센서스] 층 {len(SL)}, 시드 100/층: 해 {len(census)}개, {time.time()-t0:.0f}s")
    kinds = {}
    for r in census:
        kk = r["kind"] + ("|PASS" if not r["fails"] else "")
        kinds[kk] = kinds.get(kk, 0) + 1
    print("[센서스] 분류:", kinds)
    passes = [r for r in census if not r["fails"]]
    others = [r for r in census if r["kind"] == "other"]
    print(f"[센서스] 전 필터 통과: {len(passes)}개")
    for r in passes:
        print(f"   ★ S={{{','.join(r['S'])}}} θ=({', '.join(f'{v:.8f}' for v in r['theta'])}) ell={r['ell']:.9f}")
    print(f"[센서스] 비(오각형/θ†) 해: {len(others)}개 — 전부 필터 탈락이어야 함")
    for r in others[:20]:
        print(f"   S={{{','.join(r['S'])}}} x={r['theta'][0]:.4f} ell={r['ell']:.5f} fails={r['fails'][:3]}")

    # ---- 벽-전용 코너 전수: C(8,5)=56 정사각계 ---- — exhaustive wall-only corner sweep: C(8,5)=56 square systems
    corners = {}
    for S in itertools.combinations(range(8), 5):
        def Fw(z, S=S):
            return [f_W[j](*tuple(z)) for j in S]
        for th0 in seeds[:40]:
            z = newton(Fw, list(th0))
            if z is None: continue
            th = tuple(z)
            try:
                ws, _ = U0._walls8f(th)
            except Exception:
                continue          # 퇴화 배치 (평행선 등) — 코너 아님 — degenerate configuration (parallel lines etc.) — not a corner
            nact = sum(1 for v in ws if abs(v) < 1e-7)
            if nact >= 6:
                key = tuple(round(v, 5) for v in th)
                corners.setdefault(key, (nact, th))
    print(f"[코너] ≥6벽 동시 영점: {len(corners)}개")
    for key, (nact, th) in corners.items():
        d_pen = max(abs(th[i] - TH_PEN[i]) for i in range(5))
        print(f"   {nact}벽 활성 θ≈{key}  {'= 오각형' if d_pen < 1e-4 else '★ 신규!'}")

    with open('verify/uniq_census.json', 'w') as f:
        json.dump(dict(census=census,
                       corners=[dict(nact=n_, theta=list(t_)) for n_, t_ in corners.values()],
                       theta_dagger=list(TH_DAG), pentagon=list(TH_PEN)), f, indent=1)
    print("[센서스] uniq_census.json 저장")

if __name__ == "__main__":
    run()
