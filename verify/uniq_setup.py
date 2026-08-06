# uniq_setup.py — 전역 유일성(16)의 1단계: 유리 벽 시스템의 심볼릭 구축과 대조 검증 — Step 1 of global uniqueness (16): symbolic construction of the rational wall system and cross-check verification
#
# 산출물: — Outputs:
#   (1) theta=(x,g,h,a,b)의 유리(분모 양수) 벽 함수 W1..W10 — check_edge.walls8f와 — (1) Rational (positive-denominator) wall functions W1..W10 of theta=(x,g,h,a,b) — verified against check_edge.walls8f
#       영점·부호가 일치함을 200+ 표본에서 검증 (정규화 양수 배율만 다름). — that zeros and signs agree on 200+ samples (they differ only by a positive normalization scale factor).
#   (2) ell = (g+h+x)/2 + r,  r^2 = R2(theta) (유리),  s23 (유리). — (2) ell = (g+h+x)/2 + r, r^2 = R2(theta) (rational), s23 (rational).
#   (3) 각 벽의 x-편미분의 부호 표본 조사 — 층 가지치기 보조정리의 데이터. — (3) Sample survey of the sign of each wall's x-partial derivative — data for the stratum-pruning lemma.
#
# 핵심 구조 발견(계산 확인 대상): — Key structural discovery (to be confirmed by computation):
#   D2 = (x^2-1, -2x),        |D2| = x^2+1      (다항!) — (polynomial!)
#   W0 = (q^2-1, 2q), q=g-h,  |W0| = q^2+1      (다항!) — (polynomial!)
#   ⇒ 단위벡터가 유리 ⇒ 모든 벽·중심선 자료가 theta의 유리식. — ⇒ the unit vectors are rational ⇒ all wall and centerline data are rational expressions in theta.

import math, random
import sympy as sp

x, g, h, a, b, r = sp.symbols('x g h a b r', real=True)
TH = (x, g, h, a, b)

q = g - h
D2 = sp.Matrix([x**2 - 1, -2*x])          # 패스 A 방향 (비정규); |D2| = x^2+1 — pass A direction (unnormalized); |D2| = x^2+1
W0 = sp.Matrix([q**2 - 1, 2*q])           # 패스 C 방향 (비정규); |W0| = q^2+1 — pass C direction (unnormalized); |W0| = q^2+1
nD2 = x**2 + 1
nW0 = q**2 + 1

E = sp.Matrix([0, 0]); F = sp.Matrix([-x, 1]); G = sp.Matrix([g, 1]); H = sp.Matrix([h, 0])
Q1 = G + (a/nW0) * W0                     # 크리스3 끝점 (a,b는 단위 스케일 파라미터) — crease-3 endpoint (a, b are unit-scale parameters)
Q0 = H + (b/nW0) * W0
U = Q0 - Q1                               # 크리스3 방향 — crease-3 direction

def cross2(u, v):
    return sp.expand(u[0]*v[1] - u[1]*v[0])

def refl_lin(dvec, v):
    """직선 방향 dvec에 대한 선형 반사 (유리): 2 (v·d) d / |d|^2 − v.
    Linear reflection about the line direction dvec (rational): 2 (v·d) d / |d|^2 − v."""
    return (2*dvec.dot(v)/dvec.dot(dvec))*dvec - v

def refl_pt(P, dvec, Qp):
    return P + refl_lin(dvec, Qp - P)

e1p = refl_pt(Q1, U, G)      # rho3(G): D 밴드 모서리 기준점 — rho3(G): reference point for the D-band edge
e0p = refl_pt(Q1, U, H)      # rho3(H)
dD = refl_lin(U, W0)                      # 패스 D 방향 (비정규; 배율 양수 무관) — pass D direction (unnormalized; independent of the positive scale factor)

M1 = sp.Matrix([-x/2, sp.Rational(1, 2)])
M2 = sp.Matrix([(h + g)/2, sp.Rational(1, 2)])
# M3 = line(M2, W0) ∩ line(Q1, U):  M2 + t*W0 = Q1 + s*U
_den23 = cross2(U, W0)
_s23 = cross2(M2 - Q1, W0) / _den23       # M2 + t·W0 = Q1 + s·U 의 s (크리스3 파라미터) — the s in M2 + t·W0 = Q1 + s·U (crease-3 parameter)
M3 = Q1 + _s23 * U
R2 = sp.expand((M3 - M2).dot(M3 - M2))    # |M3-M2|^2 (유리) — (rational)
# 수치 코드의 s23 = M2에서 단위 w0 방향 파라미터 t. 유효 영역(t>0)에서 |M3-M2| = t: — s23 in the numeric code = the parameter t along the unit w0 direction from M2. On the valid region (t>0), |M3-M2| = t:
#   ⇒ ell = (g+h+x)/2 + S23 은 theta의 **유리함수** (sqrt 불필요) — KKT 계가 다항계로 환원. — ⇒ ell = (g+h+x)/2 + S23 is a **rational function** of theta (no sqrt needed) — the KKT system reduces to a polynomial system.
S23 = nW0 * cross2(Q1 - M2, U) / cross2(W0, U)
ELL = (g + h + x)/2 + S23

# ---- 벽 (실현측 >= 0; 부호 상수는 느슨 슬라이스에서 고정) ---- — walls (realizable side >= 0; the sign constants are fixed on the loose slice)
# 원본 코드의 sA, sD, s6에 해당하는 배향은 아래에서 수치 대조로 결정·고정한다. — The orientations corresponding to sA, sD, s6 in the original code are determined and fixed below by numerical cross-checking.
W1_raw = cross2(G - E, D2)                # ∝ ±w1
W7_raw = cross2(Q0 - E, D2) - nD2         # (부호 sA 적용 전) 거리−1 을 nD2 배율로 — (before applying the sign sA) distance−1, at the nD2 scale factor
W4_raw = cross2(H - e1p, dD)
W6_raw = cross2(E - e0p, dD)
# w5 = a − λ_I,  I = line(F, D2) ∩ line(G, W0):  G + λ_I·(W0/nW0) = F + t·D2
_denI = cross2(D2, W0)
lamI = nW0 * cross2(G - F, D2) / _denI    # I = G + (lamI/nW0)·W0 ∈ e_B:  cross(I−F, D2)=0 풀이 — solving cross(I−F, D2)=0
W5_raw = a - lamI
# taco: crease1 = E + t (F−E), crease3 = Q1 + u U
_denT = cross2(F - E, U)
t9 = cross2(Q1 - E, U) / _denT
u10 = cross2(Q1 - E, F - E) / _denT
W9_raw = t9 - 1
W10_raw = u10 - 1

# |dD| = |W0| (반사는 등거리): 일반 기호로 추상 검증 — 큰 식 simplify 회피 — |dD| = |W0| (reflection is an isometry): verified abstractly with generic symbols — avoids simplify on large expressions
_u1,_u2,_w1,_w2 = sp.symbols('_u1 _u2 _w1 _w2')
_uv = sp.Matrix([_u1,_u2]); _wv = sp.Matrix([_w1,_w2])
_dgen = (2*_uv.dot(_wv)/_uv.dot(_uv))*_uv - _wv
_id = sp.simplify(_dgen.dot(_dgen) - _wv.dot(_wv))
print("[구조] |refl_u(w)|^2 − |w|^2 =", _id, " (일반 기호 항등식 ⇒ |dD| = |W0| = q^2+1)")
assert _id == 0
W8_raw = cross2(F - e1p, dD) - (q**2 + 1)

# ---------- 수치 대조: check_edge.walls8f 재구현 (독립 복사) ---------- — numerical cross-check: re-implementation of check_edge.walls8f (independent copy)
def _unit(v):
    n = math.hypot(v[0], v[1]); return (v[0]/n, v[1]/n)
def _cross(u, v): return u[0]*v[1] - u[1]*v[0]
def _sub(u, v): return (u[0]-v[0], u[1]-v[1])
def _add(u, v): return (u[0]+v[0], u[1]+v[1])
def _mul(u, s): return (u[0]*s, u[1]*s)
def _dot(u, v): return u[0]*v[0] + u[1]*v[1]
def _refl_lin(d):
    dh = _unit(d)
    return lambda v: _sub(_mul(dh, 2*_dot(v, dh)), v)
def _refl_pt(P, d):
    L = _refl_lin(d)
    return lambda Qp: _add(P, L(_sub(Qp, P)))
def _meet(P1, d1, P2, d2v):
    det = d1[0]*(-d2v[1]) - d1[1]*(-d2v[0])
    if abs(det) < 1e-14: return None, None
    rhs = _sub(P2, P1)
    s = (rhs[0]*(-d2v[1]) - rhs[1]*(-d2v[0]))/det
    return _add(P1, _mul(d1, s)), s

def _config(th):
    xx, gg, hh, aa, bb = th
    Ee=(0.0,0.0); Ff=(-xx,1.0); Gg=(gg,1.0); Hh=(hh,0.0)
    d2u = _unit(_refl_lin(_sub(Ff,Ee))((1.0,0.0)))
    w0u = _unit(_refl_lin(_sub(Gg,Hh))((1.0,0.0)))
    Q1n = _add(Gg,_mul(w0u,aa)); Q0n = _add(Hh,_mul(w0u,bb))
    rho3 = _refl_pt(Q1n,_sub(Q0n,Q1n)); rl3 = _refl_lin(_sub(Q0n,Q1n))
    M1n=(-xx/2,0.5); M2n=((hh+gg)/2,0.5)
    M3n,s23n = _meet(M2n,w0u,Q1n,_sub(Q0n,Q1n))
    v3 = _unit(rl3(w0u)) if M3n else None
    return dict(E=Ee,F=Ff,G=Gg,H=Hh,d2=d2u,w0=w0u,Q1=Q1n,Q0=Q0n,rho3=rho3,rl3=rl3,
                M1=M1n,M2=M2n,M3=M3n,s23=s23n,v3=v3)

def _walls8f(th):
    c = _config(th)
    Ee,Ff,Gg,Hh = c['E'],c['F'],c['G'],c['H']
    d2u,w0u = c['d2'],c['w0']
    sA = 1.0 if _cross(_sub(Ff,Ee), d2u) > 0 else -1.0
    e1pn = c['rho3'](Gg); dDn = c['rl3'](w0u)
    pD = _add(c['M3'], c['v3'])
    sD = 1.0 if _cross(_sub(pD,e1pn), dDn) > 0 else -1.0
    w1 = -sA*_cross(_sub(Gg,Ee), d2u)
    w7 =  sA*_cross(_sub(c['Q0'],Ee), d2u) - 1.0
    w4 = -sD*_cross(_sub(Hh,e1pn), dDn)
    w8 =  sD*_cross(_sub(Ff,e1pn), dDn) - 1.0
    I,_ = _meet(Ff, d2u, Gg, w0u)
    w5 = th[3] - _dot(_sub(I,Gg), w0u)
    den = _cross(_sub(Ff,Ee), _sub(c['Q0'],c['Q1']))
    t9n  = _cross(_sub(c['Q1'],Ee), _sub(c['Q0'],c['Q1']))/den
    u10n = _cross(_sub(c['Q1'],Ee), _sub(Ff,Ee))/den
    w6 = _S6*_cross(_sub(Ee,c['rho3'](Hh)), dDn)
    return [w1,w4,w5,w6,w7,w8,t9n-1.0,u10n-1.0], c

def _slice_theta(xx):
    gg=(1-xx*xx)/(2*xx); hh=(1-3*xx*xx)/(2*xx)
    Gg=(gg,1.0); Hh=(hh,0.0); Ff=(-xx,1.0)
    w0u=_unit(_refl_lin(_sub(Gg,Hh))((1.0,0.0)))
    return (xx,gg,hh,hh,_dot(_sub(Ff,Hh),w0u))

_cS = _config(_slice_theta(0.25))
_S6 = 1.0 if _cross(_sub(_cS['E'], _cS['rho3'](_cS['H'])), _cS['rl3'](_cS['w0'])) > 0 else -1.0

# ---------- 심볼릭 벽의 부호 상수 결정 + 대조 ---------- — determine the sign constants of the symbolic walls + cross-check
RAWS = [W1_raw, W4_raw, W5_raw, W6_raw, W7_raw, W8_raw, W9_raw, W10_raw]
NAMES = ["W1","W4","W5","W6","W7","W8","W9","W10"]
F_raw = [sp.lambdify(TH, e, 'math', cse=True) for e in RAWS]
F_R2  = sp.lambdify(TH, R2, 'math', cse=True)
F_s23 = sp.lambdify(TH, _s23, 'math', cse=True)

def sample_points(n=200, seed=7):
    rnd = random.Random(seed)
    pts = []
    while len(pts) < n:
        xx = rnd.uniform(0.12, 0.55)
        th0 = _slice_theta(xx)
        th = tuple(v + rnd.uniform(-0.12, 0.12) for v in th0)
        try:
            ws, c = _walls8f(th)
        except Exception:
            continue
        if c['M3'] is None or c['s23'] is None or c['s23'] <= 0: continue
        if c['M2'][0] <= c['M1'][0]: continue
        pts.append((th, ws, c))
    return pts

pts = sample_points()
print(f"[대조] 표본 {len(pts)}개 (슬라이스 근방 유효 배치)")

# 각 심볼릭 벽에 대해: 수치 벽과의 비율이 항상 양수인 부호 s를 결정, 영점 일치 확인 — For each symbolic wall: determine the sign s that makes the ratio to the numeric wall always positive, and confirm the zeros agree
import statistics
sign_fix = {}
for k, (name, fraw) in enumerate(zip(NAMES, F_raw)):
    ratios = []
    for th, ws, c in pts:
        num = ws[k]
        try:
            sym = fraw(*th)
        except Exception:
            sym = float('nan')
        if abs(num) > 1e-8 and abs(sym) > 1e-12 and not math.isnan(sym):
            ratios.append(sym / num)
    pos = sum(1 for z in ratios if z > 0); neg = len(ratios) - pos
    if pos and neg:
        print(f"  {name}: 부호 비일관!! (+{pos} / -{neg}) — 배율이 부호를 바꿈: 재검토 필요")
        sign_fix[name] = None
    else:
        s = 1 if pos else -1
        sign_fix[name] = s
        mag = [abs(z) for z in ratios]
        print(f"  {name}: 부호 {'+' if s>0 else '-'} 고정, |배율| 범위 [{min(mag):.3g}, {max(mag):.3g}] (양수 유리 배율)")
assert all(v is not None for v in sign_fix.values())

# r, s23, ell 대조 — cross-check of r, s23, ell
F_S23 = sp.lambdify(TH, S23, 'math', cse=True)
F_ELL = sp.lambdify(TH, ELL, 'math', cse=True)
err_s = max(abs(F_S23(*th) - c['s23']) for th, ws, c in pts)
err_e = max(abs(F_ELL(*th) - (math.dist(c['M2'],c['M1']) + math.dist(c['M3'],c['M2']))) for th, ws, c in pts)
print(f"[대조] s23(단위 w0 파라미터) 최대 오차 = {err_s:.2e},  ell 최대 오차 = {err_e:.2e}")
# (참고: 전개형 R2 는 대규모 소거로 float 평가가 불안정 — S23 항등 |M3M2|=s23 이 그 역할을 대체) — (note: the expanded R2 is unstable under float evaluation because of massive cancellation — the S23 identity |M3M2|=s23 takes over that role)
assert err_s < 5e-9 and err_e < 5e-9
print("[구조] 유효 영역에서 |M3-M2| = s23 (단위 파라미터) ⇒ ell 이 유리함수임을 확인")

# ---------- x-편미분 부호 조사 (가지치기 데이터) ---------- — survey of the x-partial-derivative signs (pruning data)
print("\n[가지치기 데이터] 각 벽의 ∂/∂x 부호 (부호 고정 후, 표본 전체):")
for name, e, s in [(n_, e_, sign_fix[n_]) for n_, e_ in zip(NAMES, RAWS)]:
    dex = sp.diff(s*e, x)
    f = sp.lambdify(TH, dex, 'math', cse=True)
    vals = [f(*th) for th, ws, c in pts]
    pos = sum(1 for v in vals if v > 1e-12); neg = sum(1 for v in vals if v < -1e-12)
    zero = len(vals) - pos - neg
    print(f"  ∂x({name}): + {pos:3d} / - {neg:3d} / ≈0 {zero:3d}"
          + ("   ← x-독립" if sp.simplify(dex) == 0 else ""))

print("\n[구조 확인] W1(부호 고정 후) =", sp.factor(sign_fix['W1']*W1_raw))
print("\n설정 완료 — 이 부호·식을 uniq_strata.py 가 임포트한다.")

# 임포트용 내보내기 (원시 트리 유지 — simplify 금지 규율) — exports for importing (raw expression tree preserved — the no-simplify discipline)
WALLS = {n_: sign_fix[n_]*e_ for n_, e_ in zip(NAMES, RAWS)}
ELL_LIN = (g + h + x)/2
