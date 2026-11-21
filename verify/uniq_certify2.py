# uniq_certify2.py — 전역 유일성(16) 4단계: **다항계** Krawczyk 후보 인증 — Step 4 of global uniqueness (16): Krawczyk candidate certification on the **polynomial system**
#
# uniq_census.json 의 각 해 (층 S, θ*, λ*) 에 대해: — For each solution (stratum S, θ*, λ*) in uniq_census.json:
#   시스템 (전부 ℤ-계수 다항): — the system (all polynomial with ℤ coefficients):
#     W_j^poly(θ) = 0                          (j ∈ S)
#     2·Dp²·LIN_k + 2(dNp_k·Dp − Np·dDp_k) − 2·Σ λ_j Dp²·dW_jk^poly = 0   (k=1..5)
#   ← ∇ell = Σ λ ∇W 에 2Dp² 를 곱한 것 (Dp≠0: 가지에서 부호 + 고정 — uniq_poly 검증). — ← this is ∇ell = Σ λ ∇W multiplied by 2Dp² (Dp≠0: the sign is fixed to + on the branch — verified in uniq_poly).
#   Krawczyk(dps 40)로 근 상자 인증 후, 죽이는 양(λ_i<0 / 벽 위반)의 부호를 구간 엄밀 판정. — After certifying the root box with Krawczyk (dps 40), the sign of the killing quantity (λ_i<0 / wall violation) is decided rigorously by intervals.
#
# 출력: 후보별 인증 결과 표 + uniq_certified.json — Output: a table of per-candidate certification results + uniq_certified.json

import json, itertools
import sympy as sp
from mpmath import iv, mp, mpf, matrix
import uniq_poly as P
import uniq_setup as U0

mp.dps = 40
iv.dps = 40

TH = U0.TH
NAMES = P.NAMES
POLY = [sp.expand(P.POLY[n]) for n in NAMES]
Np, Dp = sp.expand(P.Np), sp.expand(P.Dp)
LAM = sp.symbols('l0 l1 l2 l3 l4')

dW = [[sp.expand(sp.diff(e, v)) for v in TH] for e in POLY]
# 정수 계수화: Np 의 절반 계수를 2Np 치환으로 제거 (iv-lambdify가 mpf 상수를 만들지 않도록) — Integerizing the coefficients: the half coefficients of Np are removed by substituting 2Np (so that iv-lambdify does not create mpf constants)
Np2 = sp.expand(2*Np)
dNp2 = [sp.expand(sp.diff(Np2, v)) for v in TH]
dDp = [sp.expand(sp.diff(Dp, v)) for v in TH]
LIN2 = [1, 1, 1, 0, 0]     # 2·∇(g+h+x)/2

_sys_cache = {}

def system(S):
    S = tuple(S)
    if S in _sys_cache: return _sys_cache[S]
    k = len(S)
    lam = LAM[:k]
    Dp2 = Dp**2
    eqs = [POLY[j] for j in S]
    for c in range(5):
        # 2Dp²·(KKT_c):  LIN2_c·Dp² + dNp2_c·Dp − Np2·dDp_c − Σ 2λ_j Dp² dW_jc  — 전 계수 정수 — all coefficients are integers
        e = LIN2[c]*Dp2 + dNp2[c]*Dp - Np2*dDp[c]
        for i, j in enumerate(S):
            e = e - 2*lam[i]*Dp2*dW[j][c]
        eqs.append(e)                      # 전개 없이 구조 유지 (iv 안정 + 미분만 전개됨) — structure kept without expanding (iv-stable + only the derivatives are expanded)
    vars_ = list(TH) + list(lam)
    jac = [[sp.diff(e, v) for v in vars_] for e in eqs]
    f_eqs = [sp.lambdify(vars_, e, 'mpmath', cse=True) for e in eqs]
    f_jac = [[sp.lambdify(vars_, e, 'mpmath', cse=True) for e in row] for row in jac]
    _sys_cache[S] = (f_eqs, f_jac)
    return _sys_cache[S]

def _ep(I):
    lo, hi = I.a, I.b
    lo = lo.a if hasattr(lo, 'a') else lo
    hi = hi.b if hasattr(hi, 'b') else hi
    return mpf(lo), mpf(hi)

def krawczyk(S, z0, rad=mpf('1e-12'), refine=8):
    f_eqs, f_jac = system(S)
    n = len(z0)
    z = [mpf(repr(v)) for v in z0]
    for _ in range(refine):
        Fv = matrix([f(*z) for f in f_eqs])
        Jv = matrix([[f_jac[i][j](*z) for j in range(n)] for i in range(n)])
        try:
            from mpmath import lu_solve
            dz = lu_solve(Jv, -Fv)
        except Exception:
            return None
        z = [z[i] + dz[i] for i in range(n)]
    X = [iv.mpf([v - rad, v + rad]) for v in z]
    m = [iv.mpf(v) for v in z]
    Fm = [f(*m) for f in f_eqs]
    Jm = matrix([[mpf(f_jac[i][j](*z)) for j in range(n)] for i in range(n)])
    try:
        Y = Jm**-1
    except Exception:
        return None
    JX = [[f_jac[i][j](*X) for j in range(n)] for i in range(n)]
    K = []
    for i in range(n):
        acc = m[i]
        for j in range(n):
            acc = acc - iv.mpf(Y[i, j])*Fm[j]
        for j in range(n):
            s_ = (iv.mpf(1) if i == j else iv.mpf(0))
            for k2 in range(n):
                s_ = s_ - iv.mpf(Y[i, k2])*JX[k2][j]
            acc = acc + s_*(X[j] - m[j])
        K.append(acc)
    for Xi, Ki in zip(X, K):
        xa, xb = _ep(Xi); ka, kb = _ep(Ki)
        if not (xa < ka and kb < xb): return None
    return X

_fivP = {n_: sp.lambdify(TH, e_, 'mpmath', cse=True) for n_, e_ in zip(NAMES, POLY)}

# W_j = c_j w_j. A scaled wall certifies the true sign only when c_j is positive.
# c1 = x^2+1 and c7 = (q^2+1)(x^2+1) are positive throughout the chart.
# The remaining unconditional factors use positive powers of nW0 and the squared crease norm.
#   c5 = (q+x)(1-qx);  c9 = c10 ∝ cross(F, Ū)  — 가지 밖에서 부호 반전 가능 — the sign can flip outside the branch
SAFE_WALLS = ("W1", "W4", "W6", "W7", "W8")
x_, g_, h_, a_, b_ = TH
_q = g_ - h_
CFACT = {
    "W5": (_q + x_)*(1 - _q*x_),
    "W9": (-x_)*P.Ub[1] - P.Ub[0],       # cross(F, Ū), F=(-x,1)
    "W10": (-x_)*P.Ub[1] - P.Ub[0],
}
_fivC = {n_: sp.lambdify(TH, e_, 'mpmath', cse=True) for n_, e_ in CFACT.items()}

def wall_sign_on_box(name, box):
    I = _fivP[name](*box[:5])
    lo, hi = _ep(I)
    return '+' if lo > 0 else ('-' if hi < 0 else '?')

def cfact_positive_on_box(name, box):
    """조건부 배율 벽/λ 증서에 필수인 c_j > 0 공동 인증. 무조건 벽은 True.
    Joint certification of c_j > 0, which is required for the certificates of conditional-scale-factor walls / λ. Unconditional walls give True."""
    if name in SAFE_WALLS or name == "W1": return True
    I = _fivC[name](*box[:5])
    lo, _hi = _ep(I)
    return lo > 0

# 최종 티어: 참-벽(유리, uniq_setup.WALLS) 구간 평가 — 배율이 x²+1 / q²+1 / 1 로 전부 — Final tier: interval evaluation of the true walls (rational, uniq_setup.WALLS) — the scale factors are x²+1 / q²+1 / 1, all
# 무조건 양수이므로 부호가 참 벽과 항상 일치. 분모 극점(예: w5의 q+x=0)은 iv 나눗셈이 — unconditionally positive, so the sign always agrees with the true wall. At a denominator pole (e.g. q+x=0 for w5), the iv division
# 0-포함 구간을 만나면 폭발 → '?' 로 자동 스킵된다. — blows up when it meets a zero-containing interval → it is automatically skipped as '?'.
_fivR = {n_: sp.lambdify(TH, U0.WALLS[n_], 'mpmath', cse=True) for n_ in NAMES}

def true_wall_sign_on_box(name, box):
    try:
        I = _fivR[name](*box[:5])
        lo, hi = _ep(I)
    except Exception:
        return '?'
    return '+' if lo > 0 else ('-' if hi < 0 else '?')

def main():
    data = json.load(open('verify/uniq_census.json'))
    rows = [r for r in data['census'] if r['kind'] == 'other']
    print(f"[인증] 'other' 후보 {len(rows)}개")
    out = []
    for r in rows:
        S = tuple(NAMES.index(n) for n in r['S'])
        z0 = list(r['theta']) + list(r['lam'])
        box = krawczyk(S, z0)
        if box is None:
            out.append(dict(r=r, verdict='KRAWCZYK-FAIL'))
            print(f"  S={{{','.join(r['S'])}}} x={r['theta'][0]:.4f}: Krawczyk 실패 (퇴화/스퓨리어스)")
            continue
        verdict = None
        # Prefer walls with unconditionally positive scale factors.
        for n_ in SAFE_WALLS:
            j = NAMES.index(n_)
            if j in S: continue
            if wall_sign_on_box(n_, box) == '-':
                verdict = f"{n_}<0 구간 인증 (무조건 배율)"
                break
        # (2) 조건부 벽 위반 — c_j>0 공동 인증 필수 — (2) conditional wall violation — joint certification of c_j>0 is required
        if verdict is None:
            for n_ in ("W5",):
                j = NAMES.index(n_)
                if j in S: continue
                if wall_sign_on_box(n_, box) == '-' and cfact_positive_on_box(n_, box):
                    verdict = f"{n_}<0 구간 인증 (+c>0 공동 인증)"
                    break
        # (3) λ 음수 인증 — 해당 활성 벽의 c_j>0 공동 인증 필수 — (3) certification of a negative λ — joint certification of c_j>0 for the corresponding active wall is required
        if verdict is None:
            for i in range(len(S)):
                lo, hi = _ep(box[5 + i])
                n_ = NAMES[S[i]]
                if hi < 0 and cfact_positive_on_box(n_, box):
                    verdict = f"λ({n_})<0 구간 인증 (+c>0 공동 인증)"
                    break
        # (4) taco 쌍 동시 위반 — c9>0 공동 인증 필수 — (4) simultaneous violation of the taco pair — joint certification of c9>0 is required
        if verdict is None and not (6 in S or 7 in S):
            s9 = wall_sign_on_box("W9", box); s10 = wall_sign_on_box("W10", box)
            if s9 == '-' and s10 == '-' and cfact_positive_on_box("W9", box):
                verdict = "taco 쌍 동시 위반 구간 인증 (+c>0)"
        # Evaluate the unscaled rational wall directly.
        if verdict is None:
            for j, n_ in enumerate(NAMES):
                if j in S or n_ in ("W9", "W10"): continue
                if true_wall_sign_on_box(n_, box) == '-':
                    verdict = f"{n_}<0 참-벽 구간 인증"
                    break
        if verdict is None and not (6 in S or 7 in S):
            if true_wall_sign_on_box("W9", box) == '-' and true_wall_sign_on_box("W10", box) == '-':
                verdict = "taco 쌍 동시 위반 참-벽 구간 인증"
        out.append(dict(r=r, verdict=verdict or '미결(추가 검사 필요)'))
        print(f"  S={{{','.join(r['S'])}}} x={r['theta'][0]:.4f} ell={r['ell']:.5f}: {verdict or '미결!'}")
    n_ok = sum(1 for o in out if o['verdict'] not in (None, '미결(추가 검사 필요)') and 'FAIL' not in o['verdict'])
    print(f"\n[인증] 배제 인증 완료 {n_ok}/{len(rows)}  (KRAWCZYK-FAIL = 판정 불능 — 부재 증명 아님; 전원 참벽 위반 보유는 별도 감사)")
    json.dump([dict(S=o['r']['S'], theta=o['r']['theta'], verdict=o['verdict']) for o in out],
              open('verify/uniq_certified.json', 'w'), indent=1)
    print("[인증] uniq_certified.json 저장")

if __name__ == "__main__":
    main()
