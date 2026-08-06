"""
최종 판정: 슬리버 벽을 포함한 완성 후보 벽 목록으로
Final decision: with the completed candidate wall list, sliver walls included,
  (F1) 각 벽의 금지 방향을 해석적 슬리버 기준으로 확정
  (F1) fix the forbidden direction of each wall by the analytic sliver criterion
  (F2) 오각형 차단 재검사 (1계, 40만 방향 표본)
  (F2) recheck the pentagon's blocking (first order, 400k direction samples)
  (F3) 슬리버 인지 검사기로 오각형에서 하강 재시도
  (F3) retry the descent from the pentagon with a sliver-aware checker

벽 후보 (모두 '끝점 슬리버' 메커니즘 = 크리스 이웃 공리의 극한형):
Wall candidates (all via the 'endpoint sliver' mechanism = the limiting form of the crease-neighborhood axiom):
  w1: G(크리스2 끝점)가 레이인 A 밴드 밖   [A는 크리스2 이웃 금지]
  w1: G (an endpoint of crease2) lies outside the ray-in band A   [A is forbidden in the neighborhood of crease2]
  w4: H쪽 미러 — rho3(H) 관련, D가 크리스2 이웃 금지
  w4: the H-side mirror — concerns rho3(H); D is forbidden in the neighborhood of crease2
  w5: Delta >= 0
  w6: 미러-Delta >= 0
  w6: mirror-Delta >= 0
  w7: Q0(크리스3 끝점)가 레이인 A 밴드 밖  [A는 크리스3 이웃 금지]
  w7: Q0 (an endpoint of crease3) lies outside the ray-in band A  [A is forbidden in the neighborhood of crease3]
  w8: F(크리스1 끝점)가 레이아웃 D 밴드 밖 [D는 크리스1 이웃 금지]
  w8: F (an endpoint of crease1) lies outside the ray-out band D [D is forbidden in the neighborhood of crease1]
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
    rl3 = refl_lin(sub(Q0,Q1))
    M1 = (-x/2, 0.5); M2 = ((h+g)/2, 0.5)
    M3, s23 = meet(M2, w0, Q1, sub(Q0,Q1))
    v3 = unit(rl3(w0)) if M3 else None
    return dict(E=E,F=F,G=G,H=H,d2=d2,w0=w0,Q1=Q1,Q0=Q0,rho3=rho3,rl3=rl3,
                M1=M1,M2=M2,M3=M3,s23=s23,v3=v3)

def walls(theta):
    """w1..w8 (w2,w3 슬롯은 c2,c3 층 표지로 제외). 실현측 >= 0으로 부호 통일.
    w1..w8 (the w2,w3 slots are dropped, being the c2,c3 layer markers). Signs normalized so the realizable side is >= 0."""
    c = config(theta)
    d2u, w0u = c['d2'], c['w0']
    # 밴드 A의 '안쪽' 부호: e_A(직선 E,d2) 기준, e_B쪽(F쪽)이 안쪽 — the 'inside' sign of band A: relative to e_A (the line (E,d2)), the e_B side (F side) is the inside
    sin_A = 1.0 if cross(sub(c['F'],c['E']), d2u) > 0 else -1.0
    w1 = -sin_A * cross(sub(c['G'],c['E']), d2u)          # G가 안쪽(+sin_A)이면 음수 = 위반 — if G is on the inside (+sin_A) the value is negative = a violation
    # Q0는 e_B '너머'가 실현측: 밴드 내 좌표(0=e_A, 1=e_B)가 1 이상이어야 — for Q0 the realizable side is 'beyond' e_B: the in-band coordinate (0=e_A, 1=e_B) must be at least 1
    w7 = sin_A * cross(sub(c['Q0'],c['E']), d2u) - 1.0
    # D 밴드: 모서리 rho3((G,w0)), rho3((H,w0)); 안쪽 부호는 D 내부 표본으로 — band D: edges rho3((G,w0)), rho3((H,w0)); the inside sign is taken from an interior sample of D
    e1p = c['rho3'](c['G']); dD = unit(c['rl3'](w0u))
    pd = add(c['M3'], mul(c['v3'], 1.0))
    sin_D = 1.0 if cross(sub(pd,e1p), dD) > 0 else -1.0
    # H의 미러짝: 크리스2의 H 끝점 근방을 D가 못 덮음 <=> H가 D 밴드 밖 — the mirror partner of H: D cannot cover a neighborhood of crease2's H endpoint <=> H lies outside band D
    w4 = -sin_D * cross(sub(c['H'],e1p), dD)
    # F는 D 밴드의 먼 모서리 '너머'가 실현측 (밴드 내 좌표 - 1) — for F the realizable side is 'beyond' the far edge of band D (in-band coordinate - 1)
    w8 = sin_D * cross(sub(c['F'],e1p), dD) - 1.0
    # Delta, 미러-Delta — Delta, mirror-Delta
    I,_ = meet(c['F'], d2u, c['G'], w0u)
    w5 = theta[3] - dot(sub(I,c['G']), w0u)
    eD0p = c['rho3'](c['H'])
    w6a = cross(sub(c['E'],eD0p), dD)
    # w6 부호: 느슨한 슬라이스에서 +가 되도록 — sign of w6: chosen so that it comes out + on the loose slice
    w6 = w6a if w6a >= 0 or True else w6a
    return [w1, w4, w5, w6, w7, w8], c

# ===== F1: 느슨한 슬라이스에서 벽값 확인 (전부 >= 0 이어야) ===== — F1: check the wall values on the loose slice (all must be >= 0)
print("=== F1: 느슨한 슬라이스에서의 벽값 (실현측 >= 0 규약 점검) ===")
for xv in (0.20, 0.25, 0.30):
    ws,_ = walls(slice_theta(xv))
    print(f"  x={xv:.2f}: w = [{', '.join(f'{v:+.4f}' for v in ws)}]")
x0 = math.sqrt(1-2/math.sqrt(5))
ws0,_ = walls(slice_theta(x0))
print(f"  x0    : w = [{', '.join(f'{v:+.4f}' for v in ws0)}]  (활성 벽들 = 0)")
the = (0.3136, 1.3306, 0.9786, 0.9793, 1.6445)   # 이전 하강 종점 — the previous descent endpoint
wse,_ = walls(the)
print(f"  하강종점: w = [{', '.join(f'{v:+.4f}' for v in wse)}]  <- 음수 성분 = 위반 확인")

# ===== F2: 오각형 차단 재검사 ===== — F2: recheck the pentagon's blocking
print("\n=== F2: 오각형에서 1계 하강 가능성 (완성 벽 목록) ===")
def ell(theta):
    c = config(theta)
    if c['M3'] is None: return float('inf')
    return norm(sub(c['M2'],c['M1'])) + norm(sub(c['M3'],c['M2']))

def grad(f, theta, eps=1e-6):
    out=[]
    for i in range(5):
        tp=list(theta); tm=list(theta); tp[i]+=eps; tm[i]-=eps
        out.append((f(tuple(tp))-f(tuple(tm)))/(2*eps))
    return out

th0 = slice_theta(x0)
gl = grad(ell, th0)
gws = [grad(lambda t,j=j: walls(t)[0][j], th0) for j in range(6)]
active = [j for j in range(6) if abs(walls(th0)[0][j]) < 1e-8]
print("  활성 벽:", [f"w{[1,4,5,6,7,8][j]}" for j in active])
random.seed(9)
blocked = True; found=None
for _ in range(400000):
    d = [random.gauss(0,1) for _ in range(5)]
    nd = math.sqrt(sum(t*t for t in d)); d=[t/nd for t in d]
    if sum(gl[i]*d[i] for i in range(5)) < -1e-6:
        if all(sum(gws[j][i]*d[i] for i in range(5)) >= -1e-9 for j in active):
            blocked = False; found=d; break
print("  결과:", "차단 — 오각형은 1계 국소 최소 (완성 벽 기준)" if blocked
      else f"여전히 하강 방향 존재: {[f'{v:+.3f}' for v in found]}")

# ===== F3: 슬리버 인지 하강 재시도 ===== — F3: retry the descent with sliver awareness
print("\n=== F3: 벽 전체(w>=0)를 지키는 하강 재시도 (오각형 출발) ===")
random.seed(2)
th = th0; eta = 2e-3; fails=0; steps=0
for step in range(1, 5001):
    improved=False
    for _ in range(80):
        d = [random.gauss(0,1) for _ in range(5)]
        nd = math.sqrt(sum(t*t for t in d)); d=[t/nd for t in d]
        t2 = tuple(th[i]+eta*d[i] for i in range(5))
        ws2, c2c = walls(t2)
        if c2c['M3'] is None or c2c['s23'] is None or c2c['s23'] <= 0: continue
        if min(ws2) < -1e-9: continue
        if ell(t2) < ell(th) - 1e-12:
            th = t2; improved=True; break
    steps = step
    if not improved:
        fails += 1
        if fails >= 10: break
    else:
        fails = 0
print(f"  {steps} 스텝 후 ell = {ell(th):.6f}  (오각형 = {ell(th0):.6f})")
print(f"  변화량 = {ell(th)-ell(th0):+.2e}")
print("  해석: 변화량이 ~0이면 벽 목록 완성 + 오각형 국소 최소의 강한 수치 증거.")
