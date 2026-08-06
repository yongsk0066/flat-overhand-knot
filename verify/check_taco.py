"""
최종 판정 프로브: 하강 종점의 합법성.
Final decision probe: legitimacy of the descent endpoint.

T1. 하강 종점(check_descent.py의 theta)에서 크리스1(EF)과 크리스3(Q1Q0) '선분'이
    교차하는가? 교차하면: 층 순서 B>D>A>C에서 접힘쌍 (A,B)와 (C,D)가 인터리브
    => taco-taco 관통 (검사기 미구현 조건) => 하강은 불법 이동이었다.
    교차점 근방에서 네 패스가 실제로 겹치는지도 확인.
T1. at the descent endpoint (the theta of check_descent.py), do the 'segments' crease1 (EF) and
    crease3 (Q1Q0) cross? If they do: under the layer order B>D>A>C the fold pairs (A,B) and (C,D)
    interleave => taco-taco penetration (a condition the checker does not implement) => the descent
    was an illegal move.
    Also checks whether the four passes really do overlap near the crossing point.

T2. c1의 끝점 슬리버: 크리스2의 G쪽 극단(t=0.98) 표본으로 g±프로브 재실행
    => c1이 사실은 한쪽 금지(벽)인지 판정.
T2. endpoint sliver of c1: rerun the g± probe with a sample at the G-side extreme (t=0.98) of crease2
    => decide whether c1 is in fact one-sided forbidden (a wall).

T3. 슬라이스에서 c3 프로브: b±로 크리스1×크리스3 선분 교차 발생 여부
    => phi-벽(taco-taco)의 금지 방향 확정.
T3. c3 probe on the slice: whether a crease1×crease3 segment crossing occurs under b±
    => fix the forbidden direction of the phi-wall (taco-taco).
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

def seg_cross(P1,P2,Q1,Q2):
    """선분 P1P2와 Q1Q2의 내부 교차 여부와 매개변수.
    Whether the segments P1P2 and Q1Q2 cross in their interiors, together with the parameters."""
    d1 = sub(P2,P1); d2 = sub(Q2,Q1)
    det = d1[0]*(-d2[1]) - d1[1]*(-d2[0])
    if abs(det) < 1e-14: return False, None, None
    rhs = sub(Q1,P1)
    t = (rhs[0]*(-d2[1]) - rhs[1]*(-d2[0]))/det
    u = (d1[0]*rhs[1] - d1[1]*rhs[0])/det
    return (0 < t < 1 and 0 < u < 1), t, u

def slice_theta(x):
    g = (1-x*x)/(2*x); h = (1-3*x*x)/(2*x)
    G=(g,1.0); H=(h,0.0); F=(-x,1.0)
    w0 = unit(refl_lin(sub(G,H))((1.0,0.0)))
    return (x, g, h, h, dot(sub(F,H), w0))

def geo(theta):
    x,g,h,a,b = theta
    E=(0.0,0.0); F=(-x,1.0); G=(g,1.0); H=(h,0.0)
    d2 = unit(refl_lin(sub(F,E))((1.0,0.0)))
    w0 = unit(refl_lin(sub(G,H))((1.0,0.0)))
    Q1 = add(G, mul(w0,a)); Q0 = add(H, mul(w0,b))
    return E,F,G,H,d2,w0,Q1,Q0

# ===== T1: 하강 종점 ===== — descent endpoint
theta_end = (0.3136, 1.3306, 0.9786, 0.9793, 1.6445)
E,F,G,H,d2,w0,Q1,Q0 = geo(theta_end)
hit, t, u = seg_cross(E,F,Q1,Q0)
print("=== T1: 하강 종점에서 크리스1 × 크리스3 선분 교차? ===")
print(f"  교차 = {hit}" + (f",  t(크리스1) = {t:.4f}, u(크리스3) = {u:.4f}" if t is not None else ""))
if hit:
    p = add(E, mul(sub(F,E), t))
    print(f"  교차점 = ({p[0]:.4f},{p[1]:.4f})")
    print("  => 접힘쌍 (A,B)@크리스1 과 (C,D)@크리스3 이 같은 점에서 힌지.")
    print("     층 순서 B>D>A>C: [A,B] 구간이 D를 포함, [C,D] 구간이 A를 포함 => 인터리브")
    print("     => taco-taco 관통. 하강 종점은 불법 상태 (검사기의 구멍으로 통과).")

# 오각형/슬라이스에서는 t=1 (F에서 접촉)이어야 함 — for the pentagon/slice, t must be 1 (contact at F)
x0 = math.sqrt(1-2/math.sqrt(5))
E,F,G,H,d2,w0,Q1,Q0 = geo(slice_theta(x0))
_, t0, u0 = seg_cross(E,F,Q1,Q0)
print(f"  (참고) 오각형: t = {t0:.6f} (1.0 = 끝점 F 접촉이 경계)")

# ===== T2: c1 끝점 슬리버 ===== — c1 endpoint sliver
print("\n=== T2: c1 재판정 — 크리스2의 G쪽 극단 표본 (t=0.98) ===")
def covers_A_pt(theta, p):
    E,F,G,H,d2,w0,Q1,Q0 = geo(theta)
    du = d2
    s1 = cross(sub(p,E), du); s2 = cross(sub(p,F), du)
    if s1*s2 > 1e-12: return False
    M1 = (-theta[0]/2, 0.5)
    pa = sub(M1, mul(d2,1.0))
    if cross(sub(F,E), sub(p,E)) * cross(sub(F,E), sub(pa,E)) <= 0: return False
    return True

base = slice_theta(0.25)
for dg in (+0.05, -0.05):
    th2 = (base[0], base[1]+dg, base[2], base[3], base[4])
    _,_,G2,H2,_,_,_,_ = geo(th2)
    p = add(H2, mul(sub(G2,H2), 0.98))
    c1v = cross(sub(G2,(0,0)), geo(th2)[4])
    print(f"  dg={dg:+.2f}: c1 = {c1v:+.4f},  A가 크리스2@t=0.98 덮음 = {covers_A_pt(th2, p)}")

# ===== T3: 슬라이스에서 b± 로 크리스1×크리스3 교차 발생 방향 ===== — on the slice, the b± direction in which the crease1×crease3 crossing arises
print("\n=== T3: phi-벽(크리스1×크리스3)의 금지 방향 — b± 프로브 ===")
for db in (+0.05, -0.05):
    th2 = (base[0], base[1], base[2], base[3], base[4]+db)
    E,F,G,H,d2,w0,Q1,Q0 = geo(th2)
    hit, t, u = seg_cross(E,F,Q1,Q0)
    # c3 부호도 같이 — the sign of c3 as well
    c3v = cross(sub(F,Q1), unit(sub(Q0,Q1)))
    print(f"  db={db:+.2f}: 선분 교차 = {hit} (t={t if t is None else round(t,4)}), c3 = {c3v:+.4f}")
print("\n결론 후보: 교차가 발생하는 부호 = 금지측. phi = t-1 >= 0 (교점이 F 너머)가 실현측.")
