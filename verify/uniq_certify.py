# uniq_certify.py — 전역 유일성(16) 3단계: 수확된 KKT 후보의 구간 증서 — Step 3 of global uniqueness (16): interval certificates for the harvested KKT candidates
#
# 용도: uniq_strata.py 가 찾은 층 S의 방정식 해 z* = (theta*, lambda*) 에 대해 — Purpose: for an equation solution z* = (theta*, lambda*) of a stratum S found by uniq_strata.py,
#   (1) Krawczyk 연산자로 z* 근방 상자 X 에 정확히 한 개의 실근 존재를 인증하고 — (1) certify with the Krawczyk operator that a box X around z* contains exactly one real root, and
#   (2) 그 근방 위에서 "죽이는 양"(어떤 λ_i, 어떤 벽 W_i, 또는 s23 등)의 부호를 — (2) on that neighbourhood, decide the sign of the "killing quantity" (some λ_i, some wall W_i, or s23, etc.)
#       구간 산술로 엄밀 판정한다. — rigorously by interval arithmetic.
# 모든 식이 theta의 유리식이므로 lambdify('mpmath') 결과가 iv 산술에서 그대로 닫힌다. — Since every expression is rational in theta, the lambdify('mpmath') result closes directly under iv arithmetic.
#
# 사용: from uniq_certify import certify;  certify(S, zstar, kill=('lam', i) or ('wall', j) or ('expr', sympy_expr)) — 사용 = Usage

import sympy as sp
from mpmath import iv, mp, mpf, matrix, lu_solve
import uniq_setup as U

mp.dps = 40
iv.dps = 40

x, g, h, a, b = U.TH
NAMES = ["W1", "W4", "W5", "W6", "W7", "W8", "W9", "W10"]
WALLS = [U.WALLS[n] for n in NAMES]
LAM = sp.symbols('l0 l1 l2 l3 l4')

_dELL = [sp.diff(U.ELL, v) for v in U.TH]
_dW = [[sp.diff(w_, v) for v in U.TH] for w_ in WALLS]

_cache = {}

def _system(S):
    """층 S의 정사각계 EQS(z), 야코비 JAC(z), 변수 목록. z = (x,g,h,a,b, l0..l_{k-1}).
    Square system EQS(z) of stratum S, Jacobian JAC(z), and variable list. z = (x,g,h,a,b, l0..l_{k-1})."""
    S = tuple(S)
    if S in _cache: return _cache[S]
    k = len(S)
    lam = LAM[:k]
    eqs = [WALLS[j] for j in S]
    for comp in range(5):
        e = _dELL[comp]
        for i, j in enumerate(S):
            e = e - lam[i]*_dW[j][comp]
        eqs.append(e)
    vars_ = list(U.TH) + list(lam)
    jac = [[sp.diff(e, v) for v in vars_] for e in eqs]
    f_eqs = [sp.lambdify(vars_, e, 'mpmath', cse=True) for e in eqs]
    f_jac = [[sp.lambdify(vars_, e, 'mpmath', cse=True) for e in row] for row in jac]
    _cache[S] = (f_eqs, f_jac, vars_)
    return _cache[S]

def _iv_box(z, rad):
    return [iv.mpf([mpf(v) - rad, mpf(v) + rad]) for v in z]

def _endpoints(I):
    lo, hi = I.a, I.b
    lo = lo.a if hasattr(lo, 'a') else lo
    hi = hi.b if hasattr(hi, 'b') else hi
    return mpf(lo), mpf(hi)

def _contains(X, K):
    for Xi, Ki in zip(X, K):
        xa, xb = _endpoints(Xi); ka, kb = _endpoints(Ki)
        if not (xa < ka and kb < xb): return False
    return True

def krawczyk(S, zstar, rad=mpf('1e-10'), refine=6):
    """Krawczyk 인증: 상자 X ∋ 유일 실근. 성공 시 (True, 근 상자) 반환.
    Krawczyk certification: box X ∋ a unique real root. On success returns (True, root box)."""
    f_eqs, f_jac, _ = _system(S)
    n = len(zstar)
    # 고정밀 뉴턴 정련 (mpf) — high-precision Newton refinement (mpf)
    z = [mpf(v) for v in zstar]
    for _ in range(refine):
        Fv = matrix([f(*z) for f in f_eqs])
        Jv = matrix([[f_jac[i][j2](*z) for j2 in range(n)] for i in range(n)])
        try:
            dz = lu_solve(Jv, -Fv)
        except Exception:
            return False, None
        z = [z[i] + dz[i] for i in range(n)]
    X = _iv_box(z, rad)
    m = [iv.mpf(v) for v in z]
    Fm = [f(*m) for f in f_eqs]
    Jmid = matrix([[mpf(f_jac[i][j2](*z)) for j2 in range(n)] for i in range(n)])
    try:
        Y = Jmid**-1
    except Exception:
        return False, None
    JX = [[f_jac[i][j2](*X) for j2 in range(n)] for i in range(n)]
    # K = m − Y·F(m) + (I − Y·J(X))·(X − m)
    K = []
    for i in range(n):
        acc = m[i]
        for j2 in range(n):
            acc = acc - iv.mpf(Y[i, j2])*Fm[j2]
        for j2 in range(n):
            coef = (iv.mpf(1) if i == j2 else iv.mpf(0))
            s_ = coef
            for k2 in range(n):
                s_ = s_ - iv.mpf(Y[i, k2])*JX[k2][j2]
            acc = acc + s_*(X[j2] - m[j2])
        K.append(acc)
    if _contains(X, K):
        return True, X
    return False, None

def sign_on(expr_or_kill, S, box):
    """근 상자 위에서 양의 배율 무관 부호 판정. 반환 '+', '-', '?'.
    Sign decision on the root box, independent of any positive scale factor. Returns '+', '-', '?'."""
    _, _, vars_ = _system(S)
    if isinstance(expr_or_kill, tuple):
        kind, idx = expr_or_kill
        if kind == 'lam':
            I = box[5 + idx]
            lo, hi = _endpoints(I)
            return '+' if lo > 0 else ('-' if hi < 0 else '?')
        elif kind == 'wall':
            e = WALLS[idx]
        else:
            raise ValueError(kind)
    else:
        e = expr_or_kill
    f = sp.lambdify(U.TH, e, 'mpmath', cse=True)
    I = f(*box[:5])
    lo, hi = _endpoints(I)
    return '+' if lo > 0 else ('-' if hi < 0 else '?')

def certify(S, zstar, kills):
    """kills: [('lam', i) | ('wall', j) | ('expr', e)] — 하나라도 '-' 인증되면 해당 후보는 배제.
    반환: (인증 성공 여부, 상세 dict)
    kills: [('lam', i) | ('wall', j) | ('expr', e)] — if any one of them is certified '-', that candidate is excluded.
    Returns: (whether certification succeeded, detail dict)"""
    ok, box = krawczyk(tuple(S), zstar)
    if not ok:
        return False, {"stage": "krawczyk", "msg": "근 인증 실패 (상자 축약 안 됨)"}
    out = {"stage": "signs", "box_mid": [str(_endpoints(I)[0]) for I in box], "kills": []}
    for kill in kills:
        s_ = sign_on(kill if not (isinstance(kill, tuple) and kill[0] == 'expr') else kill[1], tuple(S), box)
        out["kills"].append((kill if isinstance(kill, tuple) else str(kill), s_))
        if s_ == '-':
            out["verdict"] = f"배제 인증: {kill} < 0 (구간 엄밀)"
            return True, out
    out["verdict"] = "배제 실패 — 모든 kill 양수/불확정"
    return False, out

if __name__ == "__main__":
    # 자기 시험: 층 {W1,W4,W5,W6} 에서 θ† 인증 (배제가 아니라 존재+λ>0 확인) — self-test: certify θ† on the stratum {W1,W4,W5,W6} (confirming existence + λ>0, not exclusion)
    S = (0, 1, 2, 3)   # 인덱스: W1,W4,W5,W6 — indices: W1,W4,W5,W6
    th_dag = (0.3221274973112272, 1.3906913792102, 1.0363863004308, 1.0363863004308, 1.7132441870642)
    # λ 시드는 대략값 — 뉴턴 정련이 잡아줌 — the λ seeds are rough values — the Newton refinement pins them down
    z0 = list(th_dag) + [1.0, 1.0, 1.0, 1.0]
    ok, box = krawczyk(S, z0)
    print("[자기 시험] 층 {W1,W4,W5,W6}: Krawczyk", "성공" if ok else "실패")
    if ok:
        for i in range(4):
            lo, hi = _endpoints(box[5 + i])
            print(f"   λ_{i} ∈ [{lo}, {hi}]  ⇒ {'양수 인증' if lo > 0 else '???'}")
        lo, hi = _endpoints(box[0])
        print(f"   x ∈ [{lo}, {hi}]  (x† = 0.32212749731122717…)")
