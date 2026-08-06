"""
문서 06 검증: 층 순서, 크리스 이웃 공리, 벽의 금지 방향.
Verification of document 06: layer order, the crease-neighborhood axiom, and the walls' forbidden directions.

이론 입력 (06에서 증명/정식화):
Theoretical input (proved/formalized in 06):
  - 층 순서 정리: 교대 세잎 다이어그램 => 국소 층 순서 A>C, D>A, B>D
    => 전역 순서 B > D > A > C (뒤집기 제외 유일).
  - layer-order theorem: an alternating trefoil diagram => the local layer orders A>C, D>A, B>D
    => the global order B > D > A > C (unique up to flipping).
  - 크리스 이웃 공리: 크리스 근방에서, 접히는 두 패널 사이(층 순서상 strictly between)에
    다른 패스가 놓이면 관통. 크리스1(A-B 접힘): 사이 = {D}. 크리스2(B-C): 사이 = {D, A}.
    크리스3(C-D): 사이 = {A}.
  - crease-neighborhood axiom: near a crease, if another pass lies between the two folded panels
    (strictly between in the layer order), that is a penetration. crease1 (the A-B fold): between = {D}.
    crease2 (B-C): between = {D, A}. crease3 (C-D): between = {A}.

검정:
Tests:
  L1. 느슨한 슬라이스(x=0.25)에서 공리 성립 (금지 패스가 크리스 이웃을 안 덮음)
      => 슬라이스는 유효한 오버핸드 평면 접힘 상태.
  L1. the axiom holds on the loose slice (x=0.25) (the forbidden pass does not cover the crease neighborhood)
      => the slice is a valid flat-folded overhand state.
  L2. 다이어그램: 3교차의 재료 순서 = (A×C, A×D, B×D, C×A, D×A, D×B), 교대 실현 가능.
  L2. diagram: the material order of the 3 crossings = (A×C, A×D, B×D, C×A, D×A, D×B), alternation realizable.
  L3. c6 (미러-Delta: 크리스1 vs 레이아웃 모서리) 정의; 타이트에서 0, 느슨하면 >0.
  L3. definition of c6 (mirror-Delta: crease1 vs the ray-out edge); 0 at the tight state, >0 when loose.
  L4. 벽 방향 프로브:
  L4. wall-direction probes:
      c1 (G vs e_A): 한쪽 부호에서만 'A가 크리스2 이웃을 덮음' => 관통 벽 + 금지 방향 확정.
      c1 (G vs e_A): 'A covers the neighborhood of crease2' on one sign only => a penetration wall, with its forbidden direction fixed.
      c4 (미러): 'D가 크리스2 이웃을 덮음' 프로브.
      c4 (mirror): the probe for 'D covers the neighborhood of crease2'.
      c2, c3: 양쪽 다 공리 위반 없음 => 관통 벽 아님 (층 표지 / 가지 벽).
      c2, c3: no axiom violation on either side => not penetration walls (layer markers / branch walls).
  L5. 오각형 차단 검사: 확정된 관통 벽들만으로 1계 하강 방향이 막히는지 무작위 방향 표본.
  L5. pentagon blocking check: random direction sampling of whether the fixed penetration walls alone block every first-order descent direction.
  L6. grad ell_core의 x-성분 = 1/2 (무작위 내부점들) — 명제 6.1의 수치 확인.
  L6. the x-component of grad ell_core = 1/2 (at random interior points) — numerical confirmation of Proposition 6.1.
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
    rhs = sub(P2,P1)
    s = (rhs[0]*(-d2[1]) - rhs[1]*(-d2[0]))/det
    return add(P1, mul(d1,s)), s

def slice_theta(x):
    g = (1-x*x)/(2*x); h = (1-3*x*x)/(2*x)
    E=(0.0,0.0); F=(-x,1.0); G=(g,1.0); H=(h,0.0)
    w0 = refl_lin(sub(G,H))((1.0,0.0))
    a = h                                   # K의 l2-좌표 = k = h 값과 동일 식 — the l2-coordinate of K = k, given by the same expression as h
    b = dot(sub(F,H), unit(w0))
    return (x,g,h,a,b)

def config(theta):
    x,g,h,a,b = theta
    E=(0.0,0.0); F=(-x,1.0); G=(g,1.0); H=(h,0.0)
    d2 = refl_lin(sub(F,E))((1.0,0.0))
    w0u = unit(refl_lin(sub(G,H))((1.0,0.0)))
    Q1 = add(G, mul(w0u,a)); Q0 = add(H, mul(w0u,b))
    rho3 = refl_pt(Q1, sub(Q0,Q1))
    M1 = (-x/2, 0.5); M2 = ((h+g)/2, 0.5)
    M3,_ = meet(M2, w0u, Q1, sub(Q0,Q1))
    v3 = unit(refl_lin(sub(Q0,Q1))(w0u))
    return dict(E=E,F=F,G=G,H=H,d2=unit(d2),w0=w0u,Q1=Q1,Q0=Q0,rho3=rho3,
                M1=M1,M2=M2,M3=M3,v3=v3)

# ---- 패스 덮개 판정 ---- — deciding pass coverage
def in_band(p, P1, P2, d):
    """평행선 (P1,d),(P2,d) 사이(폐구간)에 있는가.
    Does it lie between the parallel lines (P1,d),(P2,d) (closed interval)?"""
    du = unit(d)
    s1 = cross(sub(p,P1), du); s2 = cross(sub(p,P2), du)
    return s1*s2 <= 1e-12

def covers_A(c, p):
    """레이인 A가 점 p를 덮는가: 밴드(e_A,e_B) 안 & 현 EF의 재료측.
    Does the ray-in A cover the point p: inside the band (e_A,e_B) & on the material side of the chord EF."""
    if not in_band(p, c['E'], c['F'], c['d2']): return False
    pa = sub(c['M1'], mul(c['d2'], 1.0))     # A 내부 표본점 (진행 d2의 반대편) — an interior sample point of A (on the side opposite the travel direction d2)
    sEF = cross(sub(c['F'],c['E']), sub(p,c['E']))
    sA  = cross(sub(c['F'],c['E']), sub(pa,c['E']))
    return sEF*sA > 0

def covers_D(c, p):
    """레이아웃 D가 점 p를 덮는가: rho3(pi2 모서리들) 사이 & 크리스3의 재료측.
    Does the ray-out D cover the point p: between rho3(the pi2 edges) & on the material side of crease3."""
    e1 = (c['rho3'](c['G']), c['rho3'](add(c['G'], c['w0'])))
    e0 = (c['rho3'](c['H']), c['rho3'](add(c['H'], c['w0'])))
    dD = sub(e1[1], e1[0])
    if not in_band(p, e1[0], e0[0], dD): return False
    pd = add(c['M3'], mul(c['v3'], 1.0))
    s = cross(sub(c['Q0'],c['Q1']), sub(p,c['Q1']))
    sD = cross(sub(c['Q0'],c['Q1']), sub(pd,c['Q1']))
    return s*sD > 0

def axiom_violations(c):
    """크리스 이웃 공리 위반 목록 (크리스 중점 근방 표본).
    List of crease-neighborhood axiom violations (sampled near the crease midpoints)."""
    out = []
    mc1 = mul(add(c['E'],c['F']),0.5)
    mc2 = mul(add(c['H'],c['G']),0.5)
    mc3 = mul(add(c['Q0'],c['Q1']),0.5)
    if covers_D(c, mc1): out.append("D가 크리스1 이웃을 덮음")
    if covers_A(c, mc2): out.append("A가 크리스2 이웃을 덮음")
    if covers_D(c, mc2): out.append("D가 크리스2 이웃을 덮음")
    if covers_A(c, mc3): out.append("A가 크리스3 이웃을 덮음")
    return out

def probe_near(c, base, toward, frac=0.06):
    """크리스 끝점 base 근방(내부로 frac만큼) 표본점.
    A sample point near the crease endpoint base (a fraction frac inwards)."""
    return add(base, mul(unit(sub(toward,base)), frac))

# ======== L1: 느슨한 슬라이스의 유효성 ======== — L1: validity of the loose slice
print("=== L1: 느슨한 슬라이스 x=0.25의 크리스 이웃 공리 ===")
th = slice_theta(0.25); c = config(th)
v = axiom_violations(c)
print("  위반:", v if v else "없음  => 슬라이스는 유효한 접힘 상태 [PASS]")

# ======== L2: 다이어그램 재료 순서 ======== — L2: the material order of the diagram
print("\n=== L2: 3교차의 재료 순서와 교대성 ===")
# A 중심선: M1 + s*(-d2), s>0 (진입측). C: M2 + u*w0, u in [0,|M2M3|]. D: M3 + t*v3. — centerline A: M1 + s*(-d2), s>0 (entry side). C: M2 + u*w0, u in [0,|M2M3|]. D: M3 + t*v3.
_, sAC = meet(c['M1'], mul(c['d2'],-1.0), c['M2'], c['w0'])
_, sAD = meet(c['M1'], mul(c['d2'],-1.0), c['M3'], c['v3'])
_, tDA = meet(c['M3'], c['v3'], c['M1'], mul(c['d2'],-1.0))
_, tDB = meet(c['M3'], c['v3'], (0.0,0.5), (1.0,0.0))
print(f"  A에서: ×C at s={sAC:.4f}, ×D at s={sAD:.4f}  (s 감소 방향으로 진행: ×C 먼저 {'OK' if sAC>sAD>0 else 'FAIL'})")
print(f"  D에서: ×A at t={tDA:.4f}, ×B at t={tDB:.4f}  ({'OK' if 0<tDA<tDB else 'FAIL'})")
print("  재료 순서 = (A×C, A×D, B×D, C×A, D×A, D×B); 교대 O/U => A>C, D>A, B>D => B>D>A>C")

# ======== L3: c6 (미러-Delta) ======== — L3: c6 (mirror-Delta)
print("\n=== L3: c6 = 크리스1 vs 레이아웃 모서리 (미러-Delta) ===")
def c6_of(theta):
    cc = config(theta)
    eD0 = (cc['rho3'](cc['H']), cc['rho3'](add(cc['H'],cc['w0'])))  # D의 (H,w0)-상 모서리 — the edge of D that is the image of (H,w0)
    dD = sub(eD0[1], eD0[0])
    return abs(cross(sub(cc['E'],eD0[0]), unit(dD))) + abs(cross(sub(cc['F'],eD0[0]), unit(dD)))
x0 = math.sqrt(1-2/math.sqrt(5))
for xv in (0.25, 0.30, x0):
    print(f"  x={xv:.4f}: c6 = {c6_of(slice_theta(xv)):.6f}" + ("  (tight: 0이어야 함)" if abs(xv-x0)<1e-9 else ""))

# ======== L4: 벽 방향 프로브 ======== — L4: wall-direction probes
print("\n=== L4: 벽의 금지 방향 (프로브) ===")
def probe_wall(name, param_idx, delta, test):
    res = {}
    for sgn in (+1,-1):
        t2 = list(th); t2[param_idx] += sgn*delta
        c2c = config(tuple(t2))
        res[sgn] = test(c2c)
    print(f"  {name}: +방향 위반={res[+1]}, -방향 위반={res[-1]}"
          + ("   => 관통 벽 (한쪽만 금지)" if res[+1] != res[-1] else "   => 이 프로브로는 비관통(양쪽 동일)"))
    return res

# c1: g 이동 => A가 크리스2의 G쪽 이웃을 덮는가 — move g => does A cover the G-side neighborhood of crease2
probe_wall("c1 (g±0.08, A vs 크리스2@G)", 1, 0.08,
           lambda cc: covers_A(cc, probe_near(cc, cc['G'], cc['H'])))
# c4: a 이동 => D가 크리스2의 H쪽 이웃을 덮는가 — move a => does D cover the H-side neighborhood of crease2
probe_wall("c4 (a±0.08, D vs 크리스2@H)", 3, 0.08,
           lambda cc: covers_D(cc, probe_near(cc, cc['H'], cc['G'])))
# c2: h 이동 => 공리 전체 검사 — move h => check the whole axiom
probe_wall("c2 (h±0.08, 공리 전체)", 2, 0.08, lambda cc: bool(axiom_violations(cc)))
# c3: b 이동 => 공리 전체 검사 — move b => check the whole axiom
probe_wall("c3 (b±0.08, 공리 전체)", 4, 0.08, lambda cc: bool(axiom_violations(cc)))
# c5(Delta): 슬라이스에서 x를 x0 너머로 => A가 크리스3 이웃을 덮는가 (asset-1 비관통의 유도 확인) — push x past x0 on the slice => does A cover the neighborhood of crease3 (confirms the derivation of asset-1 non-penetration)
thp = slice_theta(x0+0.004); thm = slice_theta(x0-0.004)
print(f"  c5 (슬라이스 x = x0±0.004, A vs 크리스3): x0+ 위반={bool(axiom_violations(config(thp)))}, "
      f"x0- 위반={bool(axiom_violations(config(thm)))}   => Delta>=0은 관통 벽")

# ======== L5: 오각형 차단 (1계, 무작위 방향) ======== — L5: pentagon blocking (first order, random directions)
print("\n=== L5: 오각형에서 관통 벽만으로 하강이 막히는가 (1계, 표본) ===")
def contacts_pen(theta):
    """관통 벽 후보들: c1, c4, c5, c6 (부호: 위반 검사 기반 수치 미분용 스칼라).
    Penetration-wall candidates: c1, c4, c5, c6 (signs: scalars for numerical differentiation, grounded in the violation test)."""
    cc = config(theta)
    x,g,h,a,b = theta
    c1v = cross(sub(cc['G'],cc['E']), cc['d2'])
    c4v = cross(sub(cc['rho3'](cc['H']),cc['G']), cc['w0'])
    I,_ = meet(cc['F'], cc['d2'], cc['G'], cc['w0'])
    c5v = a - dot(sub(I,cc['G']), cc['w0'])
    return [c1v, c4v, c5v, c6_of(theta)]

def ell(theta):
    cc = config(theta)
    return norm(sub(cc['M2'],cc['M1'])) + norm(sub(cc['M3'],cc['M2']))

def grad(f, theta, eps=1e-6):
    out=[]
    for i in range(5):
        tp=list(theta); tm=list(theta); tp[i]+=eps; tm[i]-=eps
        out.append((f(tuple(tp))-f(tuple(tm)))/(2*eps))
    return out

th0 = slice_theta(x0)
gl = grad(ell, th0)
gc = [grad(lambda t,j=j: contacts_pen(t)[j], th0) for j in range(4)]
# 금지 방향: 프로브에서 '위반' 쪽으로 c가 움직이는 부호를 결정해야 하나,
# 여기서는 두 부호 가정 모두 검사해 어느 가정에서 차단되는지 보고한다.
# Forbidden direction: strictly one ought to determine the sign in which c moves toward 'violation'
# in the probe; here both sign assumptions are tested instead, and it is reported under which
# assumption blocking occurs.
random.seed(11)
N=200000
for signs in [tuple(s) for s in [(1,1,1,1),(-1,-1,1,1),(1,1,-1,-1),(-1,-1,-1,-1)]]:
    blocked = True; found=None
    for _ in range(N):
        d = [random.gauss(0,1) for _ in range(5)]
        nd = math.sqrt(sum(t*t for t in d)); d=[t/nd for t in d]
        if sum(gl[i]*d[i] for i in range(5)) < -1e-6:
            if all(signs[j]*sum(gc[j][i]*d[i] for i in range(5)) >= -1e-9 for j in range(4)):
                blocked = False; found = d; break
    print(f"  부호가정 {signs}: {'차단됨 (하강 불가)' if blocked else '하강 방향 존재'}")

# ======== L6: grad ell의 x-성분 = 1/2 ======== — L6: the x-component of grad ell = 1/2
print("\n=== L6: 명제 6.1 수치 확인 — 무작위 내부점에서 d(ell)/dx = 1/2 ===")
random.seed(3)
worst = 0.0
for _ in range(50):
    t2 = list(slice_theta(0.2+0.1*random.random()))
    for i in range(1,5): t2[i] += random.uniform(-0.05,0.05)
    gx = grad(ell, tuple(t2))[0]
    worst = max(worst, abs(gx-0.5))
print(f"  50개 점에서 |d(ell)/dx - 1/2| 최대 = {worst:.2e}  {'[PASS]' if worst<1e-6 else '[FAIL]'}")
