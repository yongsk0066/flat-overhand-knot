# uniq_poly.py — 전역 유일성(16): 벽 시스템의 **다항화** (megatree 제거) — global uniqueness (16): **polynomialization** of the wall system (megatree removed)
#
# 구성 원리: 모든 분모가 가지 위 부호 고정 인자들의 곱임을 이용해 양수 배율로 소거. — Construction principle: every denominator is a product of factors of fixed sign on the branch, so it is cleared by a positive scale factor.
#   nW0 = q^2+1 > 0,   Ū·Ū > 0 (크리스3 비퇴화),   |R̄| = (Ū·Ū)·nW0 (반사 등거리), — nW0 = q^2+1 > 0, Ū·Ū > 0 (crease 3 nondegenerate), |R̄| = (Ū·Ū)·nW0 (reflection is an isometry),
#   cross(D2,W0) = 2(q+x)(qx−1) — 가지에서 (q+x)>0, (1−qx)>0. — on the branch, (q+x)>0 and (1−qx)>0.
# 결과: W1..W10 의 다항 대표 (실현측 ≥ 0, 원 벽의 양수 배율) + ell 의 소분모 유리식. — Result: polynomial representatives of W1..W10 (realizable side ≥ 0, a positive scale factor times the original walls) + a small-denominator rational expression for ell.
#
# 정의 (전부 다항): — Definitions (all polynomial):
#   Ū  = nW0·(H−G) + (b−a)·W0          ∝ 크리스3 방향 U — ∝ crease-3 direction U
#   R̄  = 2(W0·Ū)Ū − (Ū·Ū)W0            ∝ 패스 D 방향 dD  (|R̄| = (Ū·Ū)nW0) — ∝ pass D direction dD
#   V̄4 = nW0(Ū·Ū)(H−G) − a(Ū·Ū)W0 + a R̄            ∝ H − ρ3(G)
#   V̄8 = (Ū·Ū)(nW0(F−G) − aW0) + a R̄               ∝ F − ρ3(G)
#   V̄6 = −(Ū·Ū)(nW0·G + aW0) − [2((Ū−bW0)·Ū)Ū − (Ū·Ū)(Ū−bW0)]   ∝ E − ρ3(H)
#
#   W1p  = 2gx + x² − 1
#   W5p  = a(q+x)(1−qx) − x(g+x)(q²+1)
#   W7p  = 2b(q+x)(1−qx) − (2hx+x²+1)(q²+1)
#   W4p  = σ4·cross(V̄4, R̄)
#   W6p  = σ6·cross(V̄6, R̄)
#   W8p  = σ8·cross(V̄8, R̄) − nW0²(Ū·Ū)²
#   W9p  = σ9·[cross(nW0·G + aW0, Ū) − nW0·cross(F, Ū)]
#   W10p = σ9·[cross(nW0·G + aW0, F) − cross(F, Ū)]
#   ell  = (g+h+x)/2 + Np/Dp,  Np = cross(nW0(G−M2)+aW0, Ū), Dp = cross(W0, Ū)
# σ들은 표본 대조로 고정한다. — The σ's are fixed by sample cross-checking.

import math, random
import sympy as sp
import uniq_setup as U0   # 수치 대조기(_walls8f 등) 재사용 — 임포트시 자체 시험 수행 — reuses the numerical cross-checker (_walls8f etc.) — runs its self-test on import

x, g, h, a, b = U0.TH
q = g - h
nW0 = q**2 + 1
W0 = sp.Matrix([q**2 - 1, 2*q])
E = sp.Matrix([0, 0]); F = sp.Matrix([-x, 1]); G = sp.Matrix([g, 1]); H = sp.Matrix([h, 0])
M2 = sp.Matrix([(g + h)/2, sp.Rational(1, 2)])

def cross2(u, v): return sp.expand(u[0]*v[1] - u[1]*v[0])
def cross2_nx(u, v): return u[0]*v[1] - u[1]*v[0]      # 비전개 (float 안정) — unexpanded (float-stable)
def dot2(u, v): return sp.expand(u[0]*v[0] + u[1]*v[1])

Ub = sp.Matrix([sp.expand(nW0*(H[i] - G[i]) + (b - a)*W0[i]) for i in range(2)])
UU = dot2(Ub, Ub)
Rb = sp.Matrix([sp.expand(2*dot2(W0, Ub)*Ub[i] - UU*W0[i]) for i in range(2)])

V4 = sp.Matrix([sp.expand(nW0*UU*(H[i] - G[i]) - a*UU*W0[i] + a*Rb[i]) for i in range(2)])
V8 = sp.Matrix([sp.expand(UU*(nW0*(F[i] - G[i]) - a*W0[i]) + a*Rb[i]) for i in range(2)])
_t6 = sp.Matrix([sp.expand(Ub[i] - b*W0[i]) for i in range(2)])
V6 = sp.Matrix([sp.expand(-UU*(nW0*G[i] + a*W0[i])
                          - (2*dot2(_t6, Ub)*Ub[i] - UU*_t6[i])) for i in range(2)])

W1p = sp.expand(2*g*x + x**2 - 1)
W5p = sp.expand(a*(q + x)*(1 - q*x) - x*(g + x)*(q**2 + 1))
W7p = sp.expand(2*b*(q + x)*(1 - q*x) - (2*h*x + x**2 + 1)*(q**2 + 1))
W4p_raw = cross2(V4, Rb)
W6p_raw = cross2(V6, Rb)
W8p_raw = sp.expand(cross2(V8, Rb))            # σ8 적용 전; −nW0²UU² 는 σ8 적용 후에 — before applying σ8; −nW0²UU² is added after σ8 is applied
QA = sp.Matrix([sp.expand(nW0*G[i] + a*W0[i]) for i in range(2)])   # = nW0·Q1
W9p_raw = sp.expand(cross2(QA, Ub) - nW0*cross2(F, Ub))
W10p_raw = sp.expand(cross2(QA, F) - cross2(F, Ub))

Np = cross2(sp.Matrix([sp.expand(nW0*(G[i] - M2[i]) + a*W0[i]) for i in range(2)]), Ub)
Dp = cross2(W0, Ub)
ELL_LIN = (g + h + x)/2

# ---------- 부호 고정 (수치 대조) ---------- — sign fixing (numerical cross-check)
NAMES = ["W1", "W4", "W5", "W6", "W7", "W8", "W9", "W10"]
_raws = {"W1": W1p, "W4": W4p_raw, "W5": W5p, "W6": W6p_raw,
         "W7": W7p, "W8": None, "W9": W9p_raw, "W10": W10p_raw}

_fW4 = sp.lambdify(U0.TH, W4p_raw, 'math', cse=True)
_fW6 = sp.lambdify(U0.TH, W6p_raw, 'math', cse=True)
_fW8r = sp.lambdify(U0.TH, W8p_raw, 'math', cse=True)
_fUU = sp.lambdify(U0.TH, UU, 'math', cse=True)
_fnW0 = sp.lambdify(U0.TH, nW0, 'math', cse=True)
_f1 = sp.lambdify(U0.TH, W1p, 'math', cse=True)
_f5 = sp.lambdify(U0.TH, W5p, 'math', cse=True)
_f7 = sp.lambdify(U0.TH, W7p, 'math', cse=True)
_f9 = sp.lambdify(U0.TH, W9p_raw, 'math', cse=True)
_f10 = sp.lambdify(U0.TH, W10p_raw, 'math', cse=True)
_fNp = sp.lambdify(U0.TH, Np, 'math', cse=True)
_fDp = sp.lambdify(U0.TH, Dp, 'math', cse=True)

def _detect_signs(pts):
    sigs = {}
    # W8: σ8·cross(V8,R̄) − nW0²UU² 가 수치 w8 의 양수 배율이 되도록 σ8 결정 — choose σ8 so that σ8·cross(V8,R̄) − nW0²UU² is a positive scale factor times the numeric w8
    for nm in ["W4", "W6", "W9", "W10"]:
        idx = {"W4": 1, "W6": 3, "W9": 6, "W10": 7}[nm]
        f = {"W4": _fW4, "W6": _fW6, "W9": _f9, "W10": _f10}[nm]
        rs = []
        for th, ws, c in pts:
            num = ws[idx]
            if abs(num) > 1e-8:
                v = f(*th)
                if abs(v) > 1e-12: rs.append(v/num)
        pos = sum(1 for z in rs if z > 0)
        assert pos == len(rs) or pos == 0, (nm, pos, len(rs))
        sigs[nm] = 1 if pos == len(rs) else -1
    for s8 in (1, -1):
        okc = 0; tot = 0
        for th, ws, c in pts:
            num = ws[5]
            v = s8*_fW8r(*th) - (_fnW0(*th)**2)*(_fUU(*th)**2)
            if abs(num) > 1e-8 and abs(v) > 1e-12:
                tot += 1
                if v/num > 0: okc += 1
        if tot and okc == tot:
            sigs["W8"] = s8
            break
    assert "W8" in sigs, "W8 부호 결정 실패"
    return sigs

pts = U0.sample_points(200, seed=13)
SIG = _detect_signs(pts)
print(f"[다항화] 부호 상수: {SIG}")

# 구조형(비전개: float 안정) — 수치 대조·수확용 — structural form (unexpanded: float-stable) — for numerical cross-checking and harvesting
W4s = SIG["W4"]*cross2_nx(V4, Rb)
W6s = SIG["W6"]*cross2_nx(V6, Rb)
W8s = SIG["W8"]*cross2_nx(V8, Rb) - nW0**2*UU**2
W9s = SIG["W9"]*W9p_raw
W10s = SIG["W10"]*W10p_raw
STRUCT = {"W1": W1p, "W4": W4s, "W5": W5p, "W6": W6s,
          "W7": W7p, "W8": W8s, "W9": W9s, "W10": W10s}

# 전개형(정확대수 전용: float 평가 금지 — 대규모 소거로 불안정) — expanded form (exact algebra only: no float evaluation — unstable because of massive cancellation)
W4p = sp.expand(SIG["W4"]*W4p_raw)
W6p = sp.expand(SIG["W6"]*W6p_raw)
W8p = sp.expand(SIG["W8"]*W8p_raw - nW0**2*UU**2)
W9p = sp.expand(SIG["W9"]*W9p_raw)
W10p = sp.expand(SIG["W10"]*W10p_raw)
POLY = {"W1": sp.expand(W1p), "W4": W4p, "W5": sp.expand(W5p), "W6": W6p,
        "W7": sp.expand(W7p), "W8": W8p, "W9": W9p, "W10": W10p}

# ---------- 전면 대조: 200 표본에서 부호·영점 동시 검증 ---------- — full cross-check: signs and zeros verified simultaneously on 200 samples
_fP = {n_: sp.lambdify(U0.TH, e_, 'math', cse=True) for n_, e_ in STRUCT.items()}
# 벽별 배율 c(θ) = v/num 의 범위를 강한 점(|num|>1e-3)에서 추정 — estimate the range of the per-wall scale factor c(θ) = v/num at strong points (|num|>1e-3)
scale = {}
for i, n_ in enumerate(NAMES):
    rs = [ _fP[n_](*th)/ws[i] for th, ws, c in pts if abs(ws[i]) > 1e-3 ]
    scale[n_] = (min(rs), max(rs))
    assert scale[n_][0] > 0, (n_, scale[n_])   # 배율 전부 양수 = 부호 일치 (강한 점) — all scale factors positive = signs agree (at strong points)
print("[다항화] 배율 범위 (강한 점):", {k: (f"{v[0]:.3g}", f"{v[1]:.3g}") for k, v in scale.items()})
bad = []
for th, ws, c in pts:
    for i, n_ in enumerate(NAMES):
        num = ws[i]; v = _fP[n_](*th)
        cmax = scale[n_][1]
        if abs(num) > 1e-6 and v*num <= 0:
            bad.append((n_, num, v, th))
        if abs(num) <= 1e-6 and abs(v) > 10*cmax*1e-6:
            bad.append((n_, num, v, th))
for n_, num, v, th in bad[:10]:
    print(f"   실패: {n_} num={num:.3e} poly={v:.3e} θ={tuple(round(t,4) for t in th)}")
print(f"[다항화] 부호·영점 대조 실패 수: {len(bad)} / {200*8}")
assert len(bad) == 0

# ell 대조 (Dp 부호: 가지에서 고정 — 확인) — ell cross-check (sign of Dp: fixed on the branch — confirmed)
sDp = [1 if _fDp(*th) > 0 else -1 for th, ws, c in pts]
assert len(set(sDp)) == 1, "Dp 부호 비고정!"
print(f"[다항화] ell 분모 Dp = cross(W0,Ū) 부호: {'+' if sDp[0]>0 else '−'} 고정")
errE = 0.0
for th, ws, c in pts:
    ell_num = math.dist(c['M2'], c['M1']) + math.dist(c['M3'], c['M2'])
    ell_sym = (th[1] + th[2] + th[0])/2 + _fNp(*th)/_fDp(*th)
    errE = max(errE, abs(ell_num - ell_sym))
print(f"[다항화] ell 대조 최대 오차: {errE:.2e}")
assert errE < 5e-9

# 전개형 ≡ 구조형: 고정밀(mpmath dps50) 20점 대조 — 전개 자체의 무결성 확인 — expanded form ≡ structural form: 20-point high-precision (mpmath dps50) cross-check — confirms the integrity of the expansion itself
from mpmath import mp as _mp
_mp.dps = 50
_fS_hp = {n_: sp.lambdify(U0.TH, e_, 'mpmath', cse=True) for n_, e_ in STRUCT.items()}
_fP_hp = {n_: sp.lambdify(U0.TH, e_, 'mpmath', cse=True) for n_, e_ in POLY.items()}
_worst = 0.0
for th, ws, c in pts[:20]:
    thm = tuple(_mp.mpf(repr(v)) for v in th)
    for n_ in NAMES:
        d_ = abs(_fS_hp[n_](*thm) - _fP_hp[n_](*thm))
        _worst = max(_worst, float(d_))
assert _worst < 1e-25, _worst
print(f"[다항화] 전개형 ≡ 구조형 (dps50, 20점): 최대차 {_worst:.1e}")

deg = {n_: sp.total_degree(e_) for n_, e_ in POLY.items()}
size = {n_: len(e_.as_ordered_terms()) for n_, e_ in POLY.items()}
print(f"[다항화] 차수: {deg}")
print(f"[다항화] 항 수: {size}")
print("다항 벽 시스템 준비 완료 — uniq_harvest2 / 정확 소거가 이것을 임포트한다.")
