"""
잼 문제(04 추측 4.7)를 향한 첫 계산: 오각형 점에서의 KKT 검사.
First computation toward the jam problem (04 Conjecture 4.7): a KKT check at the pentagon point.

5차원 전좌표 (03 명제 3.1의 자유도):
The 5-dimensional full coordinates (the degrees of freedom of 03 Proposition 3.1):
  theta = (x, g, h, a, b)
    x    : 크리스1 = E(0,0)-F(-x,1) 의 기울기 매개변수 (게이지: E 고정)
    x    : the slope parameter of crease1 = E(0,0)-F(-x,1) (gauge: E fixed)
    g, h : 크리스2 = H(h,0)-G(g,1)  (밴드1의 현)
    g, h : crease2 = H(h,0)-G(g,1)  (the chord of band1)
    a, b : 크리스3 = Q1-Q0, Q1 = G + a*w0 ∈ rho_HG({y=1}), Q0 = H + b*w0 ∈ rho_HG({y=0})
    a, b : crease3 = Q1-Q0, Q1 = G + a*w0 ∈ rho_HG({y=1}), Q0 = H + b*w0 ∈ rho_HG({y=0})
           (w0 = rho_HG의 선형부가 보낸 (1,0) — pi2의 진행방향)
           (w0 = (1,0) carried by the linear part of rho_HG — the travel direction of pi2)

접촉 함수 (03 정리 3.2의 C1-C4를 부등식으로 완화 + C5 = Delta):
Contact functions (C1-C4 of 03 Theorem 3.2 relaxed to inequalities + C5 = Delta):
  c1 = cross(G - E, unit(d2))          [G가 레이인 모서리 e_A 위 <=> 0]
  c1 = cross(G - E, unit(d2))          [G lies on the ray-in edge e_A <=> 0]
  c2 = cross(F - H, unit(w0))          [F가 pi2 모서리 rho_HG({y=0}) 위 <=> 0]
  c2 = cross(F - H, unit(w0))          [F lies on the pi2 edge rho_HG({y=0}) <=> 0]
  c3 = cross(F - Q1, unit(Q0 - Q1))    [크리스3이 F를 지남 <=> 0]
  c3 = cross(F - Q1, unit(Q0 - Q1))    [crease3 passes through F <=> 0]
  c4 = cross(rho3(H) - G, unit(w0))    [접힘3이 H를 l2 위로 <=> 0]
  c4 = cross(rho3(H) - G, unit(w0))    [fold 3 carries H onto l2 <=> 0]
  c5 = a - lambda~                     [Delta: 크리스3 l2-끝점 vs 모서리 교점 I]
  c5 = a - lambda~                     [Delta: the l2-endpoint of crease3 vs the edge intersection I]

목적 함수: ell_core(theta) = |M1M2| + |M2M3| (중심선 크리스점 사이 거리).
Objective function: ell_core(theta) = |M1M2| + |M2M3| (distance between the centerline crease points).

검사: 오각형 점 theta*에서
Checks: at the pentagon point theta*
  (i)   c1..c5 = 0 (5접촉 활성)
  (i)   c1..c5 = 0 (all 5 contacts active)
  (ii)  {grad c_i} 의 랭크 = 5 (횡단성 — 03 열린 질문 2의 수치 판정)
  (ii)  rank of {grad c_i} = 5 (transversality — a numerical decision on open question 2 of 03)
  (iii) grad ell = sum mu_i grad c_i 의 승수 mu_i (부호는 벽 방향 확정 후 해석)
  (iii) the multipliers mu_i in grad ell = sum mu_i grad c_i (signs to be interpreted once the wall directions are fixed)
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
    rhs = sub(P2,P1)
    s = (rhs[0]*(-d2[1]) - rhs[1]*(-d2[0]))/det
    return add(P1, mul(d1,s)), s

def state(theta):
    x,g,h,a,b = theta
    E=(0.0,0.0); F=(-x,1.0); G=(g,1.0); H=(h,0.0)
    d2 = refl_lin(sub(F,E))((1.0,0.0))
    w0 = refl_lin(sub(G,H))((1.0,0.0))
    Q1 = add(G, mul(w0,a)); Q0 = add(H, mul(w0,b))
    return E,F,G,H,d2,w0,Q1,Q0

def contacts(theta):
    E,F,G,H,d2,w0,Q1,Q0 = state(theta)
    rho3 = refl_pt(Q1, sub(Q0,Q1))
    c1 = cross(sub(G,E), unit(d2))
    c2 = cross(sub(F,H), unit(w0))
    c3 = cross(sub(F,Q1), unit(sub(Q0,Q1)))
    c4 = cross(sub(rho3(H),G), unit(w0))
    I,_ = meet(F, d2, G, w0)          # e_B ∩ l2
    lam = dot(sub(I,G), unit(w0))
    c5 = theta[3] - lam               # a - lambda
    return [c1,c2,c3,c4,c5]

def ell_core(theta):
    x,g,h,a,b = theta
    E,F,G,H,d2,w0,Q1,Q0 = state(theta)
    M1 = (-x/2, 0.5)
    M2 = ((h+g)/2, 0.5)
    M3,_ = meet(M2, unit(w0), Q1, sub(Q0,Q1))
    return norm(sub(M2,M1)) + norm(sub(M3,M2))

def grad(f, theta, eps=1e-7):
    out=[]
    for i in range(5):
        tp = list(theta); tm = list(theta)
        tp[i]+=eps; tm[i]-=eps
        out.append((f(tuple(tp))-f(tuple(tm)))/(2*eps))
    return out

# ---- 오각형 점 theta* ---- — the pentagon point theta*
x0 = math.sqrt(1-2/math.sqrt(5))
g0 = (1-x0*x0)/(2*x0)
h0 = (1-3*x0*x0)/(2*x0)
k0 = (1-3*x0*x0)/(2*x0)
# b0: F = H + b0*w0  =>  b0 = dot(F-H, w0)
E,F,G,H,d2,w0,_,_ = state((x0,g0,h0,0.0,0.0))
b0 = dot(sub(F,H), unit(w0))
theta0 = (x0, g0, h0, k0, b0)

print("=== 오각형 점 theta* ===")
print("theta* =", [f"{v:.6f}" for v in theta0])
cs = contacts(theta0)
print("(i) 접촉값 c1..c5:", [f"{c:.2e}" for c in cs], "-> 전부 0이어야 함")

print("\n(ii) 접촉 기울기 행렬과 랭크")
Gm = [grad(lambda t, j=j: contacts(t)[j], theta0) for j in range(5)]
for j,row in enumerate(Gm):
    print(f"  grad c{j+1} =", [f"{v:+.4f}" for v in row])
# 랭크: 가우스 소거 — rank: Gaussian elimination
import copy
A = copy.deepcopy(Gm); n=5; rank=0; piv_tol=1e-8
for col in range(n):
    p = max(range(rank,n), key=lambda r0: abs(A[r0][col]))
    if abs(A[p][col]) < piv_tol: continue
    A[rank],A[p] = A[p],A[rank]
    for r0 in range(n):
        if r0!=rank and abs(A[r0][col])>0:
            f0 = A[r0][col]/A[rank][col]
            A[r0] = [A[r0][i]-f0*A[rank][i] for i in range(n)]
    rank+=1
print("  rank =", rank, "(5면 횡단성 성립: 5접촉이 독립 — 타이트 점은 0차원)")

print("\n(iii) KKT 승수: grad ell = sum mu_j grad c_j")
gl = grad(ell_core, theta0)
print("  grad ell =", [f"{v:+.4f}" for v in gl])
# 5x5 선형계 풀기 (Gm^T mu = gl) — solve the 5x5 linear system (Gm^T mu = gl)
M = [[Gm[j][i] for j in range(5)] for i in range(5)]
v = list(gl)
mu = [0.0]*5
# 가우스 소거 with partial pivoting — Gaussian elimination
Ab = [M[i]+[v[i]] for i in range(5)]
for col in range(5):
    p = max(range(col,5), key=lambda r0: abs(Ab[r0][col]))
    Ab[col],Ab[p] = Ab[p],Ab[col]
    for r0 in range(5):
        if r0!=col:
            f0 = Ab[r0][col]/Ab[col][col]
            Ab[r0] = [Ab[r0][i]-f0*Ab[col][i] for i in range(6)]
for i in range(5): mu[i] = Ab[i][5]/Ab[i][i]
print("  mu =", [f"{m:+.6f}" for m in mu])
res = [sum(Gm[j][i]*mu[j] for j in range(5)) - gl[i] for i in range(5)]
print("  잔차 |res| =", f"{math.hypot(*res[:2])+math.hypot(*res[2:4])+abs(res[4]):.2e}")
print("\n해석: rank 5 + 유일 승수 => 오각형에서 KKT는 자동 성립(0차원 꼭짓점).")
print("      승수의 부호가 '벽의 금지 방향'과 일치하는지가 잼 없음의 관건 —")
print("      부호 해석은 벽 분류(층 순서 공리) 확정 후. mu 값은 그 대조를 위한 기록.")

# ---- 보너스: 슬라이스 밖 무작위 인근 점에서 ell_core 하강 방향 존재 확인 ---- — bonus: confirm that a descent direction for ell_core exists at random nearby points off the slice
print("\n=== 인근 내부점에서 하강 가능성 표본 검사 ===")
import random
random.seed(7)
cnt_desc = 0; N = 200
for _ in range(N):
    d = [random.gauss(0,1) for _ in range(5)]
    nd = math.sqrt(sum(t*t for t in d)); d=[t/nd for t in d]
    tpert = tuple(theta0[i] - 0.02*abs(d[i]) if i==0 else theta0[i] + 0.01*d[i] for i in range(5))
    gl2 = grad(ell_core, tpert)
    ngl2 = math.sqrt(sum(t*t for t in gl2))
    if ngl2 > 1e-6: cnt_desc += 1
print(f"무작위 인근점 {N}개 중 grad ell != 0 인 점: {cnt_desc}/{N}")
print("(내부 임계점 부재의 표본 증거 — grad ell = 0인 내부점이 없으면 잼은 경계에서만 가능)")
