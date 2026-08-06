"""
벽 분류 결정판: 유효 상태들의 부호 패턴으로 벽의 방향을 경험적으로 확정.
Definitive wall classification: fix the walls' directions empirically from the sign patterns of valid states.

방법 (W1): 느슨한 슬라이스 점 theta(0.25) 주변 무작위 섭동 수천 개를 생성,
  크리스 이웃 공리 + 재료 정합성(전방 도달, 교차 범위)을 만족하는 것만 채택,
  채택된 상태들의 (c1,...,c6) 부호 패턴을 집계한다.
Method (W1): generate thousands of random perturbations around the loose slice point theta(0.25),
  keep only those satisfying the crease-neighborhood axiom + material consistency (forward reach,
  crossing range), and tabulate the (c1,...,c6) sign patterns of the accepted states.
  - 어떤 c_j가 유효 표본에서 한 부호만 가지면 => c_j는 벽 (금지 방향 = 반대 부호).
  - if some c_j takes only one sign over the valid samples => c_j is a wall (forbidden direction = the opposite sign).
  - 양쪽 부호 모두 나타나면 => c_j는 벽이 아니라 층 표지.
  - if both signs occur => c_j is not a wall but a layer marker.

(W2): 확정된 부호로 오각형 점에서 1계 하강 가능성 재검사.
(W2): with the fixed signs, recheck first-order descent feasibility at the pentagon point.
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
    v3 = unit(rl3(w0))
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

def valid(c):
    """크리스 이웃 공리 + 재료 정합성.
    Crease-neighborhood axiom + material consistency."""
    if c['M3'] is None or c['s23'] is None or c['s23'] <= 0: return False
    if c['M2'][0] <= c['M1'][0]: return False
    # 교차 범위: A×C, A×D, B×D가 실제 재료 범위 안에서 일어나는가 — crossing range: do A×C, A×D, B×D occur within the actual material range
    _, sAC = meet(c['M1'], mul(c['d2'],-1.0), c['M2'], c['w0'])
    _, sAD = meet(c['M1'], mul(c['d2'],-1.0), c['M3'], c['v3'])
    _, tDA = meet(c['M3'], c['v3'], c['M1'], mul(c['d2'],-1.0))
    _, tDB = meet(c['M3'], c['v3'], (0.0,0.5), (1.0,0.0))
    if None in (sAC,sAD,tDA,tDB): return False
    if not (sAC > 0 and sAD > 0 and 0 < tDA and 0 < tDB): return False
    xB = c['M3'][0] + tDB*c['v3'][0]
    if not (c['M1'][0] < xB < c['M2'][0]): return False          # B×D가 B 범위 안 — B×D lies within B's range
    _, uCA = meet(c['M2'], c['w0'], c['M1'], mul(c['d2'],-1.0))
    if uCA is None or not (0 < uCA < c['s23']): return False     # C×A가 C 범위 안 — C×A lies within C's range
    # 크리스 이웃 공리 — crease-neighborhood axiom
    mc1 = mul(add(c['E'],c['F']),0.5)
    mc2 = mul(add(c['H'],c['G']),0.5)
    mc3 = mul(add(c['Q0'],c['Q1']),0.5)
    if covers_D(c, mc1): return False
    if covers_A(c, mc2) or covers_D(c, mc2): return False
    if covers_A(c, mc3): return False
    return True

def contacts6(theta):
    c = config(theta)
    x,g,h,a,b = theta
    c1 = cross(sub(c['G'],c['E']), c['d2'])
    c2 = cross(sub(c['F'],c['H']), c['w0'])
    c3 = cross(sub(c['F'],c['Q1']), unit(sub(c['Q0'],c['Q1'])))
    c4 = cross(sub(c['rho3'](c['H']),c['G']), c['w0'])
    I,_ = meet(c['F'], c['d2'], c['G'], c['w0'])
    c5 = a - dot(sub(I,c['G']), c['w0'])
    eD0p = c['rho3'](c['H']); dD = c['rl3'](c['w0'])
    c6 = cross(sub(c['E'],eD0p), unit(dD))
    return [c1,c2,c3,c4,c5,c6], c

# ======== W1: 유효 표본의 부호 집계 ======== — W1: tabulate the signs over the valid samples
print("=== W1: 무작위 섭동 중 유효 상태의 (c1..c6) 부호 패턴 ===")
random.seed(42)
base = slice_theta(0.25)
count = {j: {+1:0, -1:0} for j in range(6)}
nvalid = 0; N = 20000
for _ in range(N):
    t2 = tuple(base[i] + random.uniform(-0.06, 0.06) for i in range(5))
    cs, c = contacts6(t2)
    if valid(c):
        nvalid += 1
        for j in range(6):
            if abs(cs[j]) > 1e-9:
                count[j][+1 if cs[j] > 0 else -1] += 1
print(f"  유효 표본: {nvalid}/{N}")
feas_sign = {}
for j in range(6):
    p, m = count[j][+1], count[j][-1]
    if p and not m: verdict, feas_sign[j] = "벽: 실현측 = (+), 금지 = (-)", +1
    elif m and not p: verdict, feas_sign[j] = "벽: 실현측 = (-), 금지 = (+)", -1
    elif p and m: verdict, feas_sign[j] = f"양쪽 실현(+{p}/-{m}) => 벽 아님(층 표지)", 0
    else: verdict, feas_sign[j] = "표본 없음", 0
    print(f"  c{j+1}: +{p} / -{m}   => {verdict}")

# ======== W2: 오각형 차단 재검사 (확정 부호 사용) ======== — W2: recheck the pentagon's blocking (using the fixed signs)
print("\n=== W2: 오각형에서 1계 하강 가능성 (경험 부호 사용) ===")
x0 = math.sqrt(1-2/math.sqrt(5))
th0 = slice_theta(x0)

def ell(theta):
    c = config(theta)
    return norm(sub(c['M2'],c['M1'])) + norm(sub(c['M3'],c['M2']))

def grad(f, theta, eps=1e-6):
    out=[]
    for i in range(5):
        tp=list(theta); tm=list(theta); tp[i]+=eps; tm[i]-=eps
        out.append((f(tuple(tp))-f(tuple(tm)))/(2*eps))
    return out

gl = grad(ell, th0)
gcs = [grad(lambda t,j=j: contacts6(t)[0][j], th0) for j in range(6)]
walls = [j for j in range(6) if feas_sign[j] != 0]
print("  벽으로 확정된 접촉:", [f"c{j+1}(실현측 {'+' if feas_sign[j]>0 else '-'})" for j in walls])
random.seed(5)
blocked = True; found = None
for _ in range(400000):
    d = [random.gauss(0,1) for _ in range(5)]
    nd = math.sqrt(sum(t*t for t in d)); d = [t/nd for t in d]
    if sum(gl[i]*d[i] for i in range(5)) < -1e-6:
        if all(feas_sign[j]*sum(gcs[j][i]*d[i] for i in range(5)) >= -1e-9 for j in walls):
            blocked = False; found = d; break
print("  결과:", "차단됨 — 오각형은 (확정 벽 기준) 1계 국소 최소" if blocked
      else f"하강 방향 존재: {[f'{v:+.3f}' for v in found]}")

# ======== W3: 하강 방향이 있다면 실제 유효 경로인지 2차 확인 ======== — W3: if a descent direction exists, a second check that it is an actually valid path
if not blocked and found:
    print("\n=== W3: 그 방향으로 작게 이동한 점의 유효성/ell 변화 ===")
    for step in (1e-3, 5e-3, 1e-2):
        t2 = tuple(th0[i] + step*found[i] for i in range(5))
        cs, c = contacts6(t2)
        print(f"  step={step:g}: valid={valid(c)}, ell 변화={ell(t2)-ell(th0):+.2e}, "
              f"c부호=[{','.join('+' if v>1e-12 else ('-' if v<-1e-12 else '0') for v in cs)}]")
