"""
결정 실험: 오각형에서 나가는 4-벽 모서리 {w1=w4=w5=w6=0}를 비선형 추적.
Decisive experiment: nonlinear tracing of the 4-wall edge {w1=w4=w5=w6=0} leaving the pentagon.

방법: 접선 = [grad w1; w4; w5; w6]의 영공간(1차원), ell 감소 방향으로 배향.
Method: the tangent = the (1-dimensional) null space of [grad w1; w4; w5; w6], oriented so that ell decreases.
  스텝 후 뉴턴 사영으로 모서리에 복귀. 매 스텝:
  After each step, return to the edge by Newton projection. At every step:
    - 8벽 전부 >= -tol (w9,w10은 선언적: max(w9,w10) >= -tol)
    - all 8 walls >= -tol (w9,w10 are disjunctive: max(w9,w10) >= -tol)
    - 크리스 이웃 공리 (각 크리스 29개 내부점 + 끝점 슬리버는 벽이 담당)
    - the crease-neighborhood axiom (29 interior points per crease + endpoint slivers are handled by the walls)
    - 재료 정합성 (전방 도달, 교차 범위)
    - material consistency (forward reach, crossing range)
  기록: ell, theta, 벽값. 정지: 유효성 실패 / 벽 충돌 / ell 정체.
  Recorded: ell, theta, wall values. Stops on: validity failure / wall collision / ell stagnation.

결과 해석:
Interpreting the result:
  - 오래 진행하며 ell이 실제로 감소 => 오각형은 이 모델의 국소 최소가 아님 (중대 발견).
  - if it runs a long way and ell genuinely decreases => the pentagon is not a local minimum of this model (a major finding).
  - 즉시 유효성 실패 => 실패 원인 = 누락 벽의 정체.
  - if validity fails immediately => the cause of the failure identifies the missing wall.
"""
import math, random

def sub(a,b): return (a[0]-b[0], a[1]-b[1])
def add(a,b): return (a[0]+b[0], a[1]+b[1])
def mul(a,s): return (a[0]*s, a[1]*s)
def dot(a,b): return a[0]*b[0]+a[1]*b[1]
def cross(a,b): return a[0]*b[1]-a[1]*b[0]
def norm(a): return math.hypot(a[0],a[1])
def unit(a):
    n=norm(a); return (a[0]/n, a[1]/n)

def refl_lin(d):
    dh = unit(d)
    return lambda v: sub(mul(dh, 2*dot(v,dh)), v)

def refl_pt(P, d):
    L = refl_lin(d)
    return lambda Q: add(P, L(sub(Q,P)))

def meet(P1,d1,P2,d2):
    det = d1[0]*(-d2[1]) - d1[1]*(-d2[0])
    if abs(det) < 1e-14: return None, None
    rhs = sub(P2,P1)
    s = (rhs[0]*(-d2[1]) - rhs[1]*(-d2[0]))/det
    return add(P1, mul(d1,s)), s

def slice_theta(x):
    g = (1-x*x)/(2*x); h = (1-3*x*x)/(2*x)
    G=(g,1.0); H=(h,0.0); F=(-x,1.0)
    w0 = unit(refl_lin(sub(G,H))((1.0,0.0)))
    return (x, g, h, h, dot(sub(F,H), w0))

def config(th):
    x,g,h,a,b = th
    E=(0.0,0.0); F=(-x,1.0); G=(g,1.0); H=(h,0.0)
    d2 = unit(refl_lin(sub(F,E))((1.0,0.0)))
    w0 = unit(refl_lin(sub(G,H))((1.0,0.0)))
    Q1 = add(G, mul(w0,a)); Q0 = add(H, mul(w0,b))
    rho3 = refl_pt(Q1, sub(Q0,Q1))
    rl3 = refl_lin(sub(Q0,Q1))
    M1 = (-x/2, 0.5); M2 = ((h+g)/2, 0.5)
    M3, s23 = meet(M2, w0, Q1, sub(Q0,Q1))
    v3 = unit(rl3(w0)) if M3 else None
    return dict(E=E,F=F,G=G,H=H,d2=d2,w0=w0,Q1=Q1,Q0=Q0,rho3=rho3,rl3=rl3,
                M1=M1,M2=M2,M3=M3,s23=s23,v3=v3)

def walls8(th, c=None):
    if c is None: c = config(th)
    E,F,G,H = c['E'],c['F'],c['G'],c['H']
    d2,w0 = c['d2'],c['w0']
    sA = 1.0 if cross(sub(F,E), d2) > 0 else -1.0
    e1p = c['rho3'](G); dD = c['rl3'](w0)
    pD = add(c['M3'], c['v3'])          # D 밴드 내부점 (퇴화 없는 프로브) — an interior point of band D (a non-degenerate probe)
    sD = 1.0 if cross(sub(pD,e1p), dD) > 0 else -1.0
    w1 = -sA*cross(sub(G,E), d2)
    w7 =  sA*cross(sub(c['Q0'],E), d2) - 1.0
    w4 = -sD*cross(sub(H,e1p), dD)
    w8 =  sD*cross(sub(F,e1p), dD) - 1.0
    I,_ = meet(F, d2, G, w0)
    w5 = th[3] - dot(sub(I,G), w0)
    w6 = cross(sub(E,c['rho3'](H)), dD)
    if w6 < 0: w6 = -w6  # 방향 상수는 연속 추적에서 고정되므로 절대값 사용 회피 위해 아래서 s6 고정 — the orientation constant stays fixed along the continuous trace, so s6 is pinned below to avoid taking an absolute value
    den = cross(sub(F,E), sub(c['Q0'],c['Q1']))
    t9  = cross(sub(c['Q1'],E), sub(c['Q0'],c['Q1']))/den
    u10 = cross(sub(c['Q1'],E), sub(F,E))/den
    return [w1,w4,w5,w6,w7,w8,t9-1.0,u10-1.0]

# w6 부호 상수 고정 (느슨 슬라이스에서 +) — pin the sign constant of w6 (+ on the loose slice)
_c = config(slice_theta(0.25))
_s6 = 1.0 if cross(sub(_c['E'],_c['rho3'](_c['H'])), _c['rl3'](_c['w0'])) > 0 else -1.0
def walls8f(th, c=None):
    if c is None: c = config(th)
    ws = walls8(th, c)
    w6 = _s6*cross(sub(c['E'],c['rho3'](c['H'])), c['rl3'](c['w0']))
    ws[3] = w6
    return ws

def covers_A(c, p):
    du = c['d2']
    if cross(sub(p,c['E']),du)*cross(sub(p,c['F']),du) > 1e-12: return False
    pa = sub(c['M1'], mul(du,1.0))
    return cross(sub(c['F'],c['E']), sub(p,c['E']))*cross(sub(c['F'],c['E']), sub(pa,c['E'])) > 0

def covers_D(c, p):
    e1p = c['rho3'](c['G']); e0p = c['rho3'](c['H']); dD = c['rl3'](c['w0'])
    if cross(sub(p,e1p),unit(dD))*cross(sub(p,e0p),unit(dD)) > 1e-12: return False
    pd = add(c['M3'], mul(c['v3'],1.0))
    s = cross(sub(c['Q0'],c['Q1']), sub(p,c['Q1']))
    sD = cross(sub(c['Q0'],c['Q1']), sub(pd,c['Q1']))
    return s*sD > 0

def valid_full(c, nsamp=29):
    if c['M3'] is None or c['s23'] is None or c['s23'] <= 0: return False, "M3"
    if c['M2'][0] <= c['M1'][0]: return False, "M2<M1"
    _, sAC = meet(c['M1'], mul(c['d2'],-1.0), c['M2'], c['w0'])
    _, sAD = meet(c['M1'], mul(c['d2'],-1.0), c['M3'], c['v3'])
    _, tDA = meet(c['M3'], c['v3'], c['M1'], mul(c['d2'],-1.0))
    _, tDB = meet(c['M3'], c['v3'], (0.0,0.5), (1.0,0.0))
    if None in (sAC,sAD,tDA,tDB): return False, "extent-none"
    if not (sAC>0 and sAD>0 and tDA>0 and tDB>0): return False, "extent-sign"
    xB = c['M3'][0] + tDB*c['v3'][0]
    if not (c['M1'][0] < xB < c['M2'][0]): return False, "BxD-extent"
    _, uCA = meet(c['M2'], c['w0'], c['M1'], mul(c['d2'],-1.0))
    if uCA is None or not (0 < uCA < c['s23']): return False, "CxA-extent"
    for i in range(1, nsamp+1):
        t = i/(nsamp+1)
        p1 = add(c['E'], mul(sub(c['F'],c['E']), t))
        p2 = add(c['H'], mul(sub(c['G'],c['H']), t))
        p3 = add(c['Q1'], mul(sub(c['Q0'],c['Q1']), t))
        if covers_D(c, p1): return False, f"D@crease1(t={t:.2f})"
        if covers_A(c, p2): return False, f"A@crease2(t={t:.2f})"
        if covers_D(c, p2): return False, f"D@crease2(t={t:.2f})"
        if covers_A(c, p3): return False, f"A@crease3(t={t:.2f})"
    return True, ""

def ell(th):
    c = config(th)
    if c['M3'] is None: return float('inf')
    return norm(sub(c['M2'],c['M1'])) + norm(sub(c['M3'],c['M2']))

def gradv(f, th, eps=1e-7):
    out=[]
    for i in range(5):
        tp=list(th); tm=list(th); tp[i]+=eps; tm[i]-=eps
        out.append((f(tuple(tp))-f(tuple(tm)))/(2*eps))
    return out

EQ = [0,1,2,3]   # w1,w4,w5,w6 를 등식으로 유지 — keep w1,w4,w5,w6 as equalities

def eq_vals(th):
    ws = walls8f(th)
    return [ws[j] for j in EQ]

def jac(th, eps=1e-7):
    return [gradv(lambda t,j=j: walls8f(t)[j], th, eps) for j in EQ]

def newton_project(th, iters=8):
    for _ in range(iters):
        v = eq_vals(th)
        if max(abs(z) for z in v) < 1e-12: break
        J = jac(th)
        # 최소자승 보정: J^T (J J^T)^{-1} v — 4x4 역행렬 수작업 — least-squares correction: J^T (J J^T)^{-1} v — the 4x4 inverse done by hand
        JJT = [[sum(J[i][k]*J[j][k] for k in range(5)) for j in range(4)] for i in range(4)]
        # 가우스 소거 — Gaussian elimination
        M = [row[:] + [v[i]] for i,row in enumerate(JJT)]
        for col in range(4):
            p = max(range(col,4), key=lambda r0: abs(M[r0][col]))
            M[col],M[p] = M[p],M[col]
            for r0 in range(4):
                if r0!=col:
                    f0 = M[r0][col]/M[col][col]
                    M[r0] = [M[r0][k]-f0*M[col][k] for k in range(5)]
        y = [M[i][4]/M[i][i] for i in range(4)]
        dth = [sum(J[i][k]*y[i] for i in range(4)) for k in range(5)]
        th = tuple(th[k]-dth[k] for k in range(5))
    return th

def tangent(th):
    J = jac(th)
    # 영공간: 5차원에서 4개 행에 수직인 벡터 — 그램슈미트 — null space: the vector in 5 dimensions orthogonal to the 4 rows — Gram-Schmidt
    basis = []
    for row in J:
        v = row[:]
        for u0 in basis:
            pr = sum(v[k]*u0[k] for k in range(5))
            v = [v[k]-pr*u0[k] for k in range(5)]
        n = math.sqrt(sum(z*z for z in v))
        if n > 1e-12: basis.append([z/n for z in v])
    # 임의 벡터를 직교화해 영공간 벡터 구성 — build the null-space vector by orthogonalizing an arbitrary vector
    for seed in ([1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]):
        v = list(map(float,seed))
        for u0 in basis:
            pr = sum(v[k]*u0[k] for k in range(5))
            v = [v[k]-pr*u0[k] for k in range(5)]
        n = math.sqrt(sum(z*z for z in v))
        if n > 1e-6:
            return [z/n for z in v]
    return None

x0 = math.sqrt(1-2/math.sqrt(5))
th = slice_theta(x0)
ell0 = ell(th)
print(f"오각형: ell = {ell0:.9f}")
print(f"초기 벽값: {[f'{v:+.2e}' for v in walls8f(th)]}")

step = 2e-3
log_every = 25
s_total = 0.0
prev_ell = ell0
for it in range(1, 2001):
    tg = tangent(th)
    if tg is None: print("접선 소실"); break
    gl = gradv(ell, th)
    if sum(gl[k]*tg[k] for k in range(5)) > 0: tg = [-z for z in tg]
    if abs(sum(gl[k]*tg[k] for k in range(5))) < 1e-10:
        print(f"it {it}: 접선 방향 하강률 0 — 모서리 위 정체점"); break
    th_try = newton_project(tuple(th[k]+step*tg[k] for k in range(5)))
    c_try = config(th_try)
    ws = walls8f(th_try, c_try)
    okwalls = (min(ws[:6]) > -1e-9) and (max(ws[6], ws[7]) > -1e-9)
    okval, why = valid_full(c_try)
    e_try = ell(th_try)
    if not okwalls or not okval or e_try >= prev_ell:
        print(f"it {it}: 정지 — 벽OK={okwalls}, 유효={okval}({why}), ell_try={e_try:.9f}")
        print(f"   벽값: {[f'{v:+.3e}' for v in ws]}")
        break
    th = th_try; prev_ell = e_try; s_total += step
    if it % log_every == 0:
        ws = walls8f(th)
        print(f"it {it:4d}: s={s_total:.3f}, ell={prev_ell:.9f} (d={prev_ell-ell0:+.2e}), "
              f"x={th[0]:.4f}, w=[{','.join(f'{v:+.1e}' for v in ws)}]")

print(f"\n최종: ell = {prev_ell:.9f}, 오각형 대비 {prev_ell-ell0:+.3e}")
print(f"theta = {[f'{v:.6f}' for v in th]}")
c = config(th)
print(f"꼭짓점: E={c['E']}, H=({c['H'][0]:.4f},0), G=({c['G'][0]:.4f},1), "
      f"Q1=({c['Q1'][0]:.4f},{c['Q1'][1]:.4f}), F=({c['F'][0]:.4f},1)")
