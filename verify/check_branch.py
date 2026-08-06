"""
08 검증 2: theta†의 합법성 최종 사냥.
08 verification 2: final hunt for the legitimacy of theta†.

  B1. 중심선 조각 A,B,C,D의 '전' 쌍별 교차 조사 (범위 포함):
  B1. survey of crossings over 'all' pairs of the centerline pieces A,B,C,D (ranges included):
      비인접 (A,C),(A,D),(B,D) = 각 1회, 인접 = 0회 이어야 가지 B_tref 유지.
      non-adjacent (A,C),(A,D),(B,D) = once each and adjacent = zero, for the branch B_tref to hold.
  B2. 재료 통과 순서와 교대성: (A×C, A×D, B×D, C×A, D×A, D×B) 유지 확인.
  B2. material passage order and alternation: confirm (A×C, A×D, B×D, C×A, D×A, D×B) is preserved.
  B3. 크리스 선분쌍 교차: EF×HG, HG×크리스3, EF×크리스3 — 전부 내부 교차 없음이어야.
  B3. crossings of crease segment pairs: EF×HG, HG×crease3, EF×crease3 — all must have no interior crossing.
  B4. 크리스 이웃 공리: 각 크리스 199개 내부점.
  B4. crease-neighborhood axiom: 199 interior points on each crease.
  B5. 8벽 값; Q0가 pi2 재료측인지; b† = dot(F-H, w0) 안사츠 검사.
  B5. the 8 wall values; whether Q0 lies on the pi2 material side; check of the ansatz b† = dot(F-H, w0).
"""
import math

def sub(a,b): return (a[0]-b[0], a[1]-b[1])
def add(a,b): return (a[0]+b[0], a[1]+b[1])
def mul(a,s): return (a[0]*s, a[1]*s)
def dot(a,b): return a[0]*b[0]+a[1]*b[1]
def cross(a,b): return a[0]*b[1]-a[1]*b[0]
def norm(a): return math.hypot(*a)
def unit(a):
    n=norm(a); return (a[0]/n,a[1]/n)

def RL(d):
    dx,dy=d; n=dx*dx+dy*dy
    return lambda v: (((dx*dx-dy*dy)*v[0]+2*dx*dy*v[1])/n,
                      (2*dx*dy*v[0]+(dy*dy-dx*dx)*v[1])/n)

def meet_params(P1,d1,P2,d2):
    det = d1[0]*(-d2[1]) - d1[1]*(-d2[0])
    if abs(det) < 1e-15: return None, None
    rhs = sub(P2,P1)
    s = (rhs[0]*(-d2[1]) - rhs[1]*(-d2[0]))/det
    u = (d1[0]*rhs[1] - d1[1]*rhs[0])/det
    return s, u

# theta† (family_trace 고정밀 값) — theta† (high-precision value from family_trace)
xd = 0.3221274973112271792
hd = 1.0363863004307971044
bd = 1.7132441870641958302
gd = (1-xd*xd)/(2*xd)

E=(0.0,0.0); F=(-xd,1.0); G=(gd,1.0); H=(hd,0.0)
d2 = RL(sub(F,E))((1.0,0.0))
w0 = RL(sub(G,H))((1.0,0.0))
ad = cross(sub(F,G), d2)/cross(w0, d2)
Q1 = add(G, mul(w0, ad)); Q0 = add(H, mul(w0, bd))
rl3 = RL(sub(Q0,Q1))
rho3 = lambda P: add(Q1, rl3(sub(P,Q1)))
dD = rl3(w0); v3 = dD
M1 = (-xd/2, 0.5); M2 = ((hd+gd)/2, 0.5)
sC, _ = meet_params(M2, w0, Q1, sub(Q0,Q1))
M3 = add(M2, mul(w0, sC))

ok = True
def check(label, cond, detail=""):
    global ok
    print(f"  [{'OK' if cond else 'FAIL'}] {label}" + (f"  {detail}" if detail else ""))
    ok = ok and cond

print("=== B1: 중심선 조각 전 쌍별 교차 조사 ===")
# 조각: A: M1 + s*(-d2), s in (0, INF); B: M1->M2 (선분); C: M2 + u*w0, u in (0, sC); D: M3 + t*v3, t in (0, INF) — 조각 = pieces, 선분 = segment
INF = 1e6
pieces = {
    'A': (M1, mul(d2,-1.0), 0.0, INF),
    'B': (M1, unit(sub(M2,M1)), 0.0, norm(sub(M2,M1))),
    'C': (M2, w0, 0.0, sC),
    'D': (M3, v3, 0.0, INF),
}
adj = {('A','B'),('B','C'),('C','D')}
expect1 = {('A','C'),('A','D'),('B','D')}
names = ['A','B','C','D']
crossings = {}
for i in range(4):
    for j in range(i+1,4):
        n1,n2 = names[i],names[j]
        P1,dd1,lo1,hi1 = pieces[n1]; P2,dd2,lo2,hi2 = pieces[n2]
        s,u = meet_params(P1,dd1,P2,dd2)
        cnt = 0
        if s is not None and lo1+1e-9 < s < hi1-1e-9 and lo2+1e-9 < u < hi2-1e-9:
            cnt = 1
        crossings[(n1,n2)] = (cnt, s, u)
        want = 1 if (n1,n2) in expect1 else 0
        check(f"{n1}×{n2}: 교차 {cnt}회 (기대 {want})", cnt == want,
              f"s={s if s is None else round(s,4)}, u={u if u is None else round(u,4)}")

print("\n=== B2: 재료 통과 순서와 교대성 ===")
sAC = crossings[('A','C')][1]; sAD = crossings[('A','D')][1]
tDA = crossings[('A','D')][2]; tDB = crossings[('B','D')][2]
check("A에서 ×C 먼저 (sAC > sAD)", sAC > sAD, f"sAC={sAC:.4f}, sAD={sAD:.4f}")
check("D에서 ×A 먼저 (tDA < tDB)", tDA < tDB, f"tDA={tDA:.4f}, tDB={tDB:.4f}")
print("  => 통과 순서 (A×C, A×D, B×D, C×A, D×A, D×B) 유지 — 교대 배정 B>D>A>C 그대로")

print("\n=== B3: 크리스 선분쌍 내부 교차 없음 ===")
creases = {'c1': (E, F), 'c2': (H, G), 'c3': (Q1, Q0)}
for (na, nb) in (('c1','c2'), ('c2','c3'), ('c1','c3')):
    P1,P2 = creases[na]; Q1_,Q2_ = creases[nb]
    s,u = meet_params(P1, sub(P2,P1), Q1_, sub(Q2_,Q1_))
    inside = (s is not None and 1e-9 < s < 1-1e-9 and 1e-9 < u < 1-1e-9)
    check(f"{na}×{nb} 내부 교차 없음", not inside,
          f"s={s if s is None else round(s,4)}, u={u if u is None else round(u,4)}")

print("\n=== B4: 크리스 이웃 공리 (199점) ===")
def covers_A(p):
    if cross(sub(p,E),d2)*cross(sub(p,F),d2) > 1e-13: return False
    pa = sub(M1, d2)
    return cross(sub(F,E), sub(p,E))*cross(sub(F,E), sub(pa,E)) > 0
def covers_D(p):
    e1p = rho3(G); e0p = rho3(H)
    if cross(sub(p,e1p),dD)*cross(sub(p,e0p),dD) > 1e-13: return False
    pd = add(M3, v3)
    return cross(sub(Q0,Q1), sub(p,Q1))*cross(sub(Q0,Q1), sub(pd,Q1)) > 0
viol = []
N=199
for i in range(1,N+1):
    t=i/(N+1)
    p1 = add(E, mul(sub(F,E),t)); p2 = add(H, mul(sub(G,H),t)); p3 = add(Q1, mul(sub(Q0,Q1),t))
    if covers_D(p1): viol.append(f"D@c1({t:.3f})")
    if covers_A(p2): viol.append(f"A@c2({t:.3f})")
    if covers_D(p2): viol.append(f"D@c2({t:.3f})")
    if covers_A(p3): viol.append(f"A@c3({t:.3f})")
check("위반 없음", not viol, str(viol[:5]))

print("\n=== B5: 벽·재료측·안사츠 ===")
sA_ = 1.0 if cross(sub(F,E), d2) > 0 else -1.0
w1v = -sA_*cross(sub(G,E), d2)
w7v = sA_*cross(sub(Q0,E), d2) - 1.0
e1p = rho3(G)
pdD = add(M3, v3)
sD_ = 1.0 if cross(sub(pdD,e1p), dD) > 0 else -1.0
w4v = -sD_*cross(sub(H,e1p), dD)
w8v = sD_*cross(sub(F,e1p), dD) - 1.0
sI_ = cross(sub(F,G), d2)/cross(w0, d2)
w5v = ad - sI_
w6v = cross(sub(E,rho3(H)), dD)
den = cross(sub(F,E), sub(Q0,Q1))
w9v = cross(sub(Q1,E), sub(Q0,Q1))/den - 1.0
w10v = cross(sub(Q1,E), sub(F,E))/den - 1.0
print(f"  w = [w1 {w1v:+.2e}, w4 {w4v:+.2e}, w5 {w5v:+.2e}, w6 {abs(w6v):+.2e}, "
      f"w7 {w7v:+.4f}, w8 {w8v:+.4f}, w9 {w9v:+.4f}, w10 {w10v:+.4f}]")
check("활성 {w1,w4,w5,w6}=0, {w7..w10}>0",
      abs(w1v)<1e-10 and abs(w4v)<1e-10 and abs(w5v)<1e-10 and abs(w6v)<1e-10
      and min(w7v,w8v,w9v,w10v) > 0)
fwd = cross(sub(G,H), sub(add(M2,w0), H)) * cross(sub(G,H), sub(Q0, H))
check("Q0가 크리스2(HG) 기준 pi2 재료측", fwd > 0)
b_ansatz = dot(sub(F,H), w0)
print(f"  b† = {bd:.12f} vs dot(F-H,w0) = {b_ansatz:.12f}  (차 {bd-b_ansatz:+.2e})")

print("\n====", "ALL OK — theta†는 현행 전 검사를 통과 (가지 유지 포함)" if ok else "실패 있음 — 원인 = 새 제약/가지 이탈", "====")
