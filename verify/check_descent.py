"""
갈림길 판정: 오각형에서 발견된 하강 방향을 실제로 따라간다.
Fork decision: actually follow the descent direction found at the pentagon.

강화된 유효성 검사:
Strengthened validity checks:
  - 크리스 이웃 공리를 각 크리스의 내부 9개 표본점에서 검사 (중점만 X)
  - check the crease-neighborhood axiom at 9 interior sample points of each crease (not the midpoint alone)
  - 재료 정합성(전방 도달, 교차 범위) 유지
  - material consistency (forward reach, crossing range) is maintained
  (한계: taco-taco형 조건은 미구현 — 결과 해석 시 명시)
  (limitation: taco-taco type conditions are not implemented — state this when interpreting results)

절차: 오각형 theta*에서 투영 경사 하강.
Procedure: projected gradient descent from the pentagon theta*.
  각 스텝: 무작위 후보 방향 중 [ell 감소 & c5,c6 >= 0 유지(1계) & 강화 유효성 통과]
  를 고르고 이동. 정지 조건: 후보 실패 연속 다수 or 스텝 상한.
  Each step: among random candidate directions, pick one that [decreases ell & keeps c5,c6 >= 0
  (first order) & passes the strengthened validity] and move. Stopping condition: many consecutive
  candidate failures, or the step cap.
  기록: ell, x, (c1..c6), 최소 크리스 여유(clearance).
  Recorded: ell, x, (c1..c6), minimum crease margin (clearance).
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

def config(theta):
    x,g,h,a,b = theta
    E=(0.0,0.0); F=(-x,1.0); G=(g,1.0); H=(h,0.0)
    d2 = unit(refl_lin(sub(F,E))((1.0,0.0)))
    w0 = unit(refl_lin(sub(G,H))((1.0,0.0)))
    Q1 = add(G, mul(w0,a)); Q0 = add(H, mul(w0,b))
    rho3 = refl_pt(Q1, sub(Q0,Q1))
    M1 = (-x/2, 0.5); M2 = ((h+g)/2, 0.5)
    M3, s23 = meet(M2, w0, Q1, sub(Q0,Q1))
    rl3 = refl_lin(sub(Q0,Q1))
    v3 = unit(rl3(w0)) if M3 else None
    return dict(E=E,F=F,G=G,H=H,d2=d2,w0=w0,Q1=Q1,Q0=Q0,rho3=rho3,rl3=rl3,
                M1=M1,M2=M2,M3=M3,s23=s23,v3=v3)

def in_band(p, P1, P2, d):
    du = unit(d)
    return cross(sub(p,P1), du) * cross(sub(p,P2), du) <= 1e-12

def covers_A(c, p):
    if not in_band(p, c['E'], c['F'], c['d2']): return False
    pa = sub(c['M1'], mul(c['d2'], 1.0))
    sEF = cross(sub(c['F'],c['E']), sub(p,c['E']))
    sA  = cross(sub(c['F'],c['E']), sub(pa,c['E']))
    return sEF*sA > 0

def covers_D(c, p):
    e1p = c['rho3'](c['G']); e0p = c['rho3'](c['H'])
    dD = c['rl3'](c['w0'])
    if not in_band(p, e1p, e0p, dD): return False
    pd = add(c['M3'], mul(c['v3'], 1.0))
    s = cross(sub(c['Q0'],c['Q1']), sub(p,c['Q1']))
    sD = cross(sub(c['Q0'],c['Q1']), sub(pd,c['Q1']))
    return s*sD > 0

def valid_strong(c, nsamp=9):
    if c['M3'] is None or c['s23'] is None or c['s23'] <= 0: return False
    if c['M2'][0] <= c['M1'][0]: return False
    _, sAC = meet(c['M1'], mul(c['d2'],-1.0), c['M2'], c['w0'])
    _, sAD = meet(c['M1'], mul(c['d2'],-1.0), c['M3'], c['v3'])
    _, tDA = meet(c['M3'], c['v3'], c['M1'], mul(c['d2'],-1.0))
    _, tDB = meet(c['M3'], c['v3'], (0.0,0.5), (1.0,0.0))
    if None in (sAC,sAD,tDA,tDB): return False
    if not (sAC > 0 and sAD > 0 and 0 < tDA and 0 < tDB): return False
    xB = c['M3'][0] + tDB*c['v3'][0]
    if not (c['M1'][0] < xB < c['M2'][0]): return False
    _, uCA = meet(c['M2'], c['w0'], c['M1'], mul(c['d2'],-1.0))
    if uCA is None or not (0 < uCA < c['s23']): return False
    for i in range(1, nsamp+1):
        t = i/(nsamp+1)
        p1 = add(c['E'], mul(sub(c['F'],c['E']), t))
        p2 = add(c['H'], mul(sub(c['G'],c['H']), t))
        p3 = add(c['Q1'], mul(sub(c['Q0'],c['Q1']), t))
        if covers_D(c, p1): return False
        if covers_A(c, p2) or covers_D(c, p2): return False
        if covers_A(c, p3): return False
    return True

def contacts6(theta, c=None):
    if c is None: c = config(theta)
    x,g,h,a,b = theta
    c1 = cross(sub(c['G'],c['E']), c['d2'])
    c2 = cross(sub(c['F'],c['H']), c['w0'])
    c3 = cross(sub(c['F'],c['Q1']), unit(sub(c['Q0'],c['Q1'])))
    c4 = cross(sub(c['rho3'](c['H']),c['G']), c['w0'])
    I,_ = meet(c['F'], c['d2'], c['G'], c['w0'])
    c5 = a - dot(sub(I,c['G']), c['w0'])
    eD0p = c['rho3'](c['H']); dD = c['rl3'](c['w0'])
    c6 = cross(sub(c['E'],eD0p), unit(dD))
    return [c1,c2,c3,c4,c5,c6]

def ell(theta):
    c = config(theta)
    if c['M3'] is None: return float('inf')
    return norm(sub(c['M2'],c['M1'])) + norm(sub(c['M3'],c['M2']))

# ===== 하강 ===== — descent
random.seed(2)
x0 = math.sqrt(1-2/math.sqrt(5))
th = slice_theta(x0)
print(f"시작(오각형): ell = {ell(th):.6f}, x = {th[0]:.6f}")
eta = 2e-3
fails = 0
for step in range(1, 20001):
    improved = False
    for _ in range(60):
        d = [random.gauss(0,1) for _ in range(5)]
        nd = math.sqrt(sum(t*t for t in d)); d = [t/nd for t in d]
        t2 = tuple(th[i] + eta*d[i] for i in range(5))
        c2c = config(t2)
        if c2c['M3'] is None: continue
        if ell(t2) < ell(th) - 1e-12:
            cs = contacts6(t2, c2c)
            if cs[4] >= -1e-12 and cs[5] >= -1e-12 and valid_strong(c2c):
                th = t2; improved = True; break
    if not improved:
        fails += 1
        if fails >= 8:
            print(f"정지: step {step}, 개선 실패 누적")
            break
    else:
        fails = 0
    if step % 500 == 0 or step == 1:
        cs = contacts6(th)
        print(f"  step {step:5d}: ell = {ell(th):.6f}, x = {th[0]:.4f}, "
              f"c=[{','.join(f'{v:+.3f}' for v in cs)}]")

cs = contacts6(th)
cfin = config(th)
print(f"\n종점: ell = {ell(th):.6f}  (오각형 = {1/x0 - x0:.6f})")
print(f"theta = {[f'{v:.4f}' for v in th]}")
print(f"c1..c6 = {[f'{v:+.4f}' for v in cs]}")
print(f"strong valid: {valid_strong(cfin)}")
print("\n해석 가이드: ell이 오각형보다 유의미하게 작으면 '모델상 오각형 너머 조임 가능'.")
print("먼저 c5(Delta), c6가 0으로 닫히며 멈추면 '벽에서 정지 — 새 타이트 후보'.")
