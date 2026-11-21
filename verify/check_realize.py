"""Sample separation conditions and the height construction for Theorem F.

  R1. 여유(클리어런스) 계산기:
  R1. margin (clearance) calculator:
      c_nb = min over {(k1,R_D),(k2,R_A),(k2,R_D),(k3,R_A)} 의 세그먼트-볼록영역 거리
      c_nb = the min over {(k1,R_D),(k2,R_A),(k2,R_D),(k3,R_A)} of the segment-to-convex-region distance
      c_cr = d(k1, k3)  (비인접 크리스쌍; 인접쌍은 강체성으로 자동 분리)
      c_cr = d(k1, k3)  (the non-adjacent crease pair; adjacent pairs are separated automatically by rigidity)
      m_mat = min(d(k1,k2), d(k2,k3))  (강체성에 의해 재료 현 간격과 동일)
      m_mat = min(d(k1,k2), d(k2,k3))  (equal to the material chord spacing, by rigidity)
      Lambda = max 크리스 현 길이,  eps0 = (2/(pi*Lambda)) * min(m_mat, 2*c_nb, c_cr)
      Lambda = the maximum crease chord length,  eps0 = (2/(pi*Lambda)) * min(m_mat, 2*c_nb, c_cr)
  R2. theta†: 활성 벽 {w1,w4,w5,w6} 기울기 랭크 = 4, 개방 방향 delta 존재
      (J delta = 1 최소노름 해), theta† + t*delta 가 열린 상태(8벽 전부 양수,
      유효성 통과, c_nb > 0)임을 t = 1e-3, 1e-2 에서 확인.
  R2. theta†: the gradients of the active walls {w1,w4,w5,w6} have rank 4, and an opening direction
      delta exists (the minimum-norm solution of J delta = 1); confirms at t = 1e-3, 1e-2 that
      theta† + t*delta is an open state (all 8 walls positive, validity passes, c_nb > 0).
  R3. 오각형 P*: 8벽 전부 활성. 최소자승 방향으로 전 벽 동시 개방 가능한지 확인
      (성립하면 P*도 열린 층의 극한 => 극한 실현 가능).
  R3. the pentagon P*: all 8 walls are active. Checks whether every wall can be opened simultaneously
      along the least-squares direction (if so, P* too is a limit of open strata => realizable in the limit).
  R4. 여유 표: theta†, P*, 개방점들에서 c_nb 4쌍 각각 + c_cr + m_mat + eps0.
  R4. margin table: each of the 4 c_nb pairs + c_cr + m_mat + eps0, at theta†, P*, and the opened points.

  자기 검사: theta†에서 활성 벽 4개가 정확히 여유 0 접촉과 대응하는지
      (d(G,R_A)=0, d(H,R_D)=0, d(Q1,R_A)=0, d(E,R_D)=0).
  Self-check: whether the 4 active walls at theta† correspond exactly to the zero-margin contacts
      (d(G,R_A)=0, d(H,R_D)=0, d(Q1,R_A)=0, d(E,R_D)=0).
"""
import math

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
    pD = add(c['M3'], c['v3'])
    sD = 1.0 if cross(sub(pD,e1p), dD) > 0 else -1.0
    w1 = -sA*cross(sub(G,E), d2)
    w7 =  sA*cross(sub(c['Q0'],E), d2) - 1.0
    w4 = -sD*cross(sub(H,e1p), dD)
    w8 =  sD*cross(sub(F,e1p), dD) - 1.0
    I,_ = meet(F, d2, G, w0)
    w5 = th[3] - dot(sub(I,G), w0)
    w6 = cross(sub(E,c['rho3'](H)), dD)
    den = cross(sub(F,E), sub(c['Q0'],c['Q1']))
    t9  = cross(sub(c['Q1'],E), sub(c['Q0'],c['Q1']))/den
    u10 = cross(sub(c['Q1'],E), sub(F,E))/den
    return [w1,w4,w5,w6,w7,w8,t9-1.0,u10-1.0]

_c0 = config(slice_theta(0.25))
_s6 = 1.0 if cross(sub(_c0['E'],_c0['rho3'](_c0['H'])), _c0['rl3'](_c0['w0'])) > 0 else -1.0
def walls8f(th, c=None):
    if c is None: c = config(th)
    ws = walls8(th, c)
    ws[3] = _s6*cross(sub(c['E'],c['rho3'](c['H'])), c['rl3'](c['w0']))
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
        if covers_D(c, p1): return False, f"D@c1({t:.2f})"
        if covers_A(c, p2): return False, f"A@c2({t:.2f})"
        if covers_D(c, p2): return False, f"D@c2({t:.2f})"
        if covers_A(c, p3): return False, f"A@c3({t:.2f})"
    return True, ""

# ---------- 볼록 영역 거리 ---------- — distance to convex regions
def region_A(c):
    # 면: (E,d2) 모서리 e_A, (F,d2) 모서리 e_B, (E,F-E) 현 EF의 꼬리쪽 — faces: the edge e_A (E,d2), the edge e_B (F,d2), and (E,F-E) the tail side of the chord EF
    pref = sub(c['M1'], mul(c['d2'], 2.0))
    faces = [(c['E'], c['d2']), (c['F'], c['d2']), (c['E'], sub(c['F'],c['E']))]
    verts = [c['E'], c['F']]
    return faces, verts, pref

def region_D(c):
    e1p = c['rho3'](c['G']); e0p = c['rho3'](c['H']); dD = unit(c['rl3'](c['w0']))
    pref = add(c['M3'], mul(c['v3'], 2.0))
    faces = [(e1p, dD), (e0p, dD), (c['Q1'], sub(c['Q0'],c['Q1']))]
    verts = [c['Q1'], c['Q0']]
    return faces, verts, pref

def dist_region(p, reg):
    faces, verts, pref = reg
    sgn = []
    for (P,v) in faces:
        vh = unit(v)
        s = cross(vh, sub(pref,P))
        sgn.append(1.0 if s > 0 else -1.0)
    def g(q):
        return [sgn[i]*cross(unit(faces[i][1]), sub(q,faces[i][0])) for i in range(len(faces))]
    if min(g(p)) >= -1e-12: return 0.0
    best = min(norm(sub(p,V)) for V in verts)
    for i,(P,v) in enumerate(faces):
        vh = unit(v)
        q = add(P, mul(vh, dot(sub(p,P), vh)))
        if min(g(q)) >= -1e-9:
            best = min(best, norm(sub(p,q)))
    return best

def seg_region_dist(A, B, reg, n=2001):
    return min(dist_region(add(A, mul(sub(B,A), i/(n-1))), reg) for i in range(n))

def seg_seg_dist(A,B,C,D):
    def pt_seg(p, S, T):
        d = sub(T,S); L2 = dot(d,d)
        t = 0.0 if L2 == 0 else max(0.0, min(1.0, dot(sub(p,S), d)/L2))
        return norm(sub(p, add(S, mul(d,t))))
    # 내부 교차? — interior crossing?
    def orient(P,Q,R): return cross(sub(Q,P), sub(R,P))
    o1,o2 = orient(A,B,C), orient(A,B,D)
    o3,o4 = orient(C,D,A), orient(C,D,B)
    if o1*o2 < 0 and o3*o4 < 0: return 0.0
    return min(pt_seg(A,C,D), pt_seg(B,C,D), pt_seg(C,A,B), pt_seg(D,A,B))

def clearances(th):
    c = config(th)
    k1 = (c['E'], c['F']); k2 = (c['H'], c['G']); k3 = (c['Q1'], c['Q0'])
    RA = region_A(c); RD = region_D(c)
    d_k1_RD = seg_region_dist(*k1, RD)
    d_k2_RA = seg_region_dist(*k2, RA)
    d_k2_RD = seg_region_dist(*k2, RD)
    d_k3_RA = seg_region_dist(*k3, RA)
    c_nb = min(d_k1_RD, d_k2_RA, d_k2_RD, d_k3_RA)
    c_cr = seg_seg_dist(*k1, *k3)
    m_mat = min(seg_seg_dist(*k1, *k2), seg_seg_dist(*k2, *k3))
    Lam = max(norm(sub(k[1],k[0])) for k in (k1,k2,k3))
    eps0 = (2.0/(math.pi*Lam))*min(m_mat, 2*c_nb, c_cr)
    return dict(pairs=(d_k1_RD, d_k2_RA, d_k2_RD, d_k3_RA),
                c_nb=c_nb, c_cr=c_cr, m_mat=m_mat, Lam=Lam, eps0=eps0)

# ---------- 선형대수 소도구 ---------- — small linear-algebra utilities
def gradv(f, th, eps=1e-7):
    out=[]
    for i in range(5):
        tp=list(th); tm=list(th); tp[i]+=eps; tm[i]-=eps
        out.append((f(tuple(tp))-f(tuple(tm)))/(2*eps))
    return out

def jac_rows(th, idxs):
    return [gradv(lambda t,j=j: walls8f(t)[j], th) for j in idxs]

def solve_sym(M, b):
    n = len(b)
    A = [row[:]+[b[i]] for i,row in enumerate(M)]
    for col in range(n):
        p = max(range(col,n), key=lambda r: abs(A[r][col]))
        if abs(A[p][col]) < 1e-14: return None
        A[col],A[p] = A[p],A[col]
        for r in range(n):
            if r != col:
                f0 = A[r][col]/A[col][col]
                A[r] = [A[r][k]-f0*A[col][k] for k in range(n+1)]
    return [A[i][n]/A[i][i] for i in range(n)]

def min_norm_sol(J, rhs):
    m = len(J)
    JJT = [[sum(J[i][k]*J[j][k] for k in range(5)) for j in range(m)] for i in range(m)]
    y = solve_sym(JJT, rhs)
    if y is None: return None
    return [sum(J[i][k]*y[i] for i in range(m)) for k in range(5)]

def lsq_sol(J, rhs):
    JTJ = [[sum(J[i][a]*J[i][b] for i in range(len(J))) for b in range(5)] for a in range(5)]
    JTr = [sum(J[i][a]*rhs[i] for i in range(len(J))) for a in range(5)]
    return solve_sym(JTJ, JTr)

def num_rank(J, tol=1e-8):
    A = [row[:] for row in J]; m=len(A); r=0; piv=[]
    for col in range(5):
        p = max(range(r,m), key=lambda i: abs(A[i][col])) if r < m else None
        if p is None or abs(A[p][col]) < tol: continue
        A[r],A[p] = A[p],A[r]; piv.append(abs(A[r][col]))
        for i in range(m):
            if i != r:
                f0 = A[i][col]/A[r][col]
                A[i] = [A[i][k]-f0*A[r][k] for k in range(5)]
        r += 1
    return r, piv

def ell(th):
    c = config(th)
    return norm(sub(c['M2'],c['M1'])) + norm(sub(c['M3'],c['M2']))

ok = True
def check(label, cond, detail=""):
    global ok
    print(f"  [{'OK' if cond else 'FAIL'}] {label}" + (f"  {detail}" if detail else ""))
    ok = ok and cond

# ============ theta† ============
xd = 0.3221274973112271792
hd = 1.0363863004307971044
bd = 1.7132441870641958302
gd = (1-xd*xd)/(2*xd)
cB = config((xd,gd,hd,hd,bd))  # a는 아래서 정확 재계산 — a is recomputed exactly below
ad = cross(sub(cB['F'],cB['G']), cB['d2'])/cross(cB['w0'], cB['d2'])
thD = (xd, gd, hd, ad, bd)

print("=== R0: 자기 검사 — theta†의 활성 벽 = 여유 0 접촉 ===")
c = config(thD)
RA = region_A(c); RD = region_D(c)
check("w1=0  <-> d(G, R_A)=0", dist_region(c['G'], RA) < 1e-10, f"{dist_region(c['G'],RA):.2e}")
check("w4=0  <-> d(H, R_D)=0", dist_region(c['H'], RD) < 1e-10, f"{dist_region(c['H'],RD):.2e}")
check("w5=0  <-> d(Q1, R_A)=0", dist_region(c['Q1'], RA) < 1e-9, f"{dist_region(c['Q1'],RA):.2e}")
check("w6=0  <-> d(E, R_D)=0", dist_region(c['E'], RD) < 1e-9, f"{dist_region(c['E'],RD):.2e}")

def report(name, th):
    cl = clearances(th)
    ws = walls8f(th)
    cfg = config(th)
    v, why = valid_full(cfg)
    print(f"  {name}: valid={v}{('/'+why) if why else ''}")
    print(f"    벽 w=[{', '.join(f'{w:+.3e}' for w in ws)}]")
    p = cl['pairs']
    print(f"    d(k1,R_D)={p[0]:.4e}  d(k2,R_A)={p[1]:.4e}  d(k2,R_D)={p[2]:.4e}  d(k3,R_A)={p[3]:.4e}")
    print(f"    c_nb={cl['c_nb']:.4e}  c_cr={cl['c_cr']:.4e}  m_mat={cl['m_mat']:.4e}  "
          f"Lambda={cl['Lam']:.4f}  eps0={cl['eps0']:.4e}")
    return cl

print("\n=== R1: 여유 표 — 닫힌(타이트) 상태 ===")
clD = report("theta†", thD)
check("theta†: c_nb = 0 (4쌍 전부 접촉)", max(clD['pairs']) < 1e-8)
check("theta†: c_cr > 0 (taco 여유 잔존)", clD['c_cr'] > 1e-3, f"{clD['c_cr']:.4f}")

x0 = math.sqrt(1-2/math.sqrt(5))
thP = slice_theta(x0)
clP = report("오각형 P*", thP)
check("P*: c_nb = 0 그리고 c_cr = 0 (8벽 코너)", max(clP['pairs']) < 1e-8 and clP['c_cr'] < 1e-8)

print("\n=== R2: theta†의 개방 — 활성 4벽 동시 완화 방향 ===")
J4 = jac_rows(thD, [0,1,2,3])
r4, piv4 = num_rank(J4)
check("rank(J[w1,w4,w5,w6]) = 4", r4 == 4, f"피벗 {['%.3f'%p for p in piv4]}")
delta = min_norm_sol(J4, [1.0,1.0,1.0,1.0])
check("개방 방향 delta 존재 (J delta = 1)", delta is not None,
      f"delta={['%+.4f'%d for d in delta]}" if delta else "")
for t in (1e-3, 1e-2):
    tho = tuple(thD[k] + t*delta[k] for k in range(5))
    ws = walls8f(tho)
    cfg = config(tho); v,_ = valid_full(cfg)
    cl = clearances(tho)
    check(f"t={t:g}: 8벽 전부 > 0, 유효, c_nb > 0",
          min(ws) > 0 and v and cl['c_nb'] > 0,
          f"min w={min(ws):.2e}, c_nb={cl['c_nb']:.3e}, eps0={cl['eps0']:.3e}, "
          f"ell={ell(tho):.6f}")
    if t == 1e-2: report(f"theta†+{t:g}·delta", tho)

print("\n=== R3: 오각형 P*의 개방 — 8벽 동시 완화 ===")
J8 = jac_rows(thP, list(range(8)))
r8, piv8 = num_rank(J8)
print(f"  rank(J[w1..w10]) = {r8} (5차원이므로 최대 5)")
dl = lsq_sol(J8, [1.0]*8)
prods = [sum(J8[i][k]*dl[k] for k in range(5)) for i in range(8)]
print(f"  최소자승 방향의 벽 방향도함수: {['%+.4f'%p for p in prods]}")
check("전 벽 방향도함수 > 0 (P* 개방 가능)", min(prods) > 0)
for t in (1e-3, 1e-2):
    tho = tuple(thP[k] + t*dl[k] for k in range(5))
    ws = walls8f(tho)
    cfg = config(tho); v,_ = valid_full(cfg)
    cl = clearances(tho)
    check(f"t={t:g}: 8벽 전부 > 0, 유효, c_nb > 0, c_cr > 0",
          min(ws) > 0 and v and cl['c_nb'] > 0 and cl['c_cr'] > 0,
          f"min w={min(ws):.2e}, c_nb={cl['c_nb']:.3e}, c_cr={cl['c_cr']:.3e}, "
          f"eps0={cl['eps0']:.3e}")
    if t == 1e-2: report(f"P*+{t:g}·delta", tho)

print("\n=== R4: 개방 경로의 극한 수렴 (정리 9.8용) ===")
# t -> 0 에서 개방점의 평면 상태가 타이트 상태로 수렴하는 속도 (구성의 연속성) — the rate at which the planar state of the opened point converges to the tight state as t -> 0 (continuity of the construction)
for t in (1e-2, 1e-3, 1e-4):
    tho = tuple(thD[k] + t*delta[k] for k in range(5))
    co = config(tho)
    dev = max(norm(sub(co[key], config(thD)[key])) for key in ('E','F','G','H','Q1','Q0','M3'))
    print(f"  t={t:g}: 꼭짓점 최대 이동 {dev:.3e}, eps0={clearances(tho)['eps0']:.3e}")

print("\n====", "ALL OK" if ok else "실패 있음", "====")
