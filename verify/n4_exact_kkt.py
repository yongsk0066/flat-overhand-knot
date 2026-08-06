# n4_exact_kkt.py — θ*_{4_1}의 KKT 승수 정확 승격 + L_trunc 정확 계산 (전부 Q(√15)) — exact promotion of the KKT multipliers of θ*_{4_1} + exact computation of L_trunc (all in Q(√15))
# 핵심: 현 벡터 r_i u_i = (d_i + d_{i+1}) / (2 sin τ_{i+1})  — 반각 소거 — key: chord vector r_i u_i = (d_i + d_{i+1}) / (2 sin τ_{i+1}) — half-angle elimination
import sympy as sp

t1,t2,t3,t4,l1,l2,l3=sp.symbols('t1 t2 t3 t4 l1 l2 l3', real=True)
TH=[t1,t2,t3,t4,l1,l2,l3]
phi=[0,t1,t1+t2,t1+t2+t3,t1+t2+t3+t4]
d=[sp.Matrix([sp.cos(p),sp.sin(p)]) for p in phi]
M=[sp.Matrix([0,0])]
for i in range(3): M.append(M[-1]+TH[4+i]*d[i+1])
CV=[(d[i]+d[i+1])/(2*sp.sin(TH[i])) if False else (d[i]+d[i+1])/(2*sp.sin([t1,t2,t3,t4][i])) for i in range(4)]
E=[(M[i]-CV[i], M[i]+CV[i]) for i in range(4)]   # (E-, E+)
lat=lambda pt,P,dd: -(pt[0]-P[0])*dd[1]+(pt[1]-P[1])*dd[0]

r15=sp.sqrt(15)
SUB={sp.cos(t1):-sp.Rational(1,4), sp.sin(t1):r15/4,
     sp.cos(t2):-sp.Rational(7,8), sp.sin(t2):r15/8,
     sp.cos(t3):-sp.Rational(7,8), sp.sin(t3):-r15/8,
     sp.cos(t4):-sp.Rational(1,4), sp.sin(t4):-r15/4,
     l1:2*r15/5, l2:8*r15/15, l3:2*r15/5}
def X(e): return sp.simplify(sp.expand_trig(e).subs(SUB))

# 접촉 함수 8개 (n4_kkt와 동일 부호; E±는 수치 대조로 라벨 확정) — the 8 contact functions (same signs as n4_kkt; the E± labels are fixed by numerical cross-check)
E1m,E1p=E[0]; E2m,E2p=E[1]; E3m,E3p=E[2]; E4m,E4p=E[3]
G=[('g1',  -lat(E1p,M[3],d[4])-sp.Rational(1,2)),
   ('g2',   lat(E1m,M[1],d[2])-sp.Rational(1,2)),
   ('g3a',  lat(E2p,M[3],d[4])-sp.Rational(1,2)),
   ('g3b',  lat(E2m,M[3],d[4])-sp.Rational(1,2)),
   ("g1p",  lat(E4m,M[0],d[0])-sp.Rational(1,2)),
   ("g2p", -lat(E4p,M[1],d[2])-sp.Rational(1,2)),
   ("g3ap",-lat(E3p,M[0],d[0])-sp.Rational(1,2)),
   ("g3bp",-lat(E3m,M[0],d[0])-sp.Rational(1,2))]
print("[활성 검산 — 전부 0이어야]")
vals={}
for n_,g in G:
    v=X(g); vals[n_]=v
    print(f"  {n_:5}: {v}")
# 라벨 스왑 필요 여부: 크리스3,4의 sin(τ/2)<0 ⟹ CV 부호 반전 → E± 교환됨. — whether a label swap is needed: for creases 3,4 sin(τ/2)<0 ⟹ the CV sign flips → E± get exchanged.
# 위 출력에서 0이 아닌 항목이 있으면 해당 크리스의 ± 교환으로 재시도 — if any entry above is nonzero, retry with the ± exchange for that crease
swap={}
for idx,(cr,(nm1,nm2)) in enumerate([(2,("g3ap","g3bp")),(3,("g1p","g2p"))]):
    if vals[nm1]!=0 or vals[nm2]!=0: swap[cr]=True
if swap:
    print("  → 라벨 스왑:", swap)
    if 2 in swap: E3m,E3p=E3p,E3m
    if 3 in swap: E4m,E4p=E4p,E4m
    G=[('g1',  -lat(E1p,M[3],d[4])-sp.Rational(1,2)),
       ('g2',   lat(E1m,M[1],d[2])-sp.Rational(1,2)),
       ('g3a',  lat(E2p,M[3],d[4])-sp.Rational(1,2)),
       ('g3b',  lat(E2m,M[3],d[4])-sp.Rational(1,2)),
       ("g1p",  lat(E4m,M[0],d[0])-sp.Rational(1,2)),
       ("g2p", -lat(E4p,M[1],d[2])-sp.Rational(1,2)),
       ("g3ap",-lat(E3p,M[0],d[0])-sp.Rational(1,2)),
       ("g3bp",-lat(E3m,M[0],d[0])-sp.Rational(1,2))]
    for n_,g in G: print(f"  {n_:5}: {X(g)}")
# 야코비안 (정확) — Jacobian (exact)
names=[n_ for n_,_ in G]
J=sp.Matrix([[X(sp.diff(g,v)) for v in TH] for _,g in G])  # 8x7
grad_ell=sp.Matrix([0,0,0,0,1,1,1])
# n4_kkt가 찾은 7-생성원 부분집합: g2 제외 — the 7-generator subset found by n4_kkt: g2 excluded
sel=[0,2,3,4,5,6,7]
A=J[sel,:].T  # 7x7
lam=A.solve(grad_ell)
lam=sp.simplify(lam)
print("\n[정확 KKT 승수] (g1,g3a,g3b,g1',g2',g3a',g3b')")
# 라벨 규약 [referee-H F5]: 크리스3은 sin(τ/2)<0라 CV-규약 E±가 config-규약과 반대. — label convention [referee-H F5]: for crease 3, sin(τ/2)<0, so the CV-convention E± is opposite to the config convention.
# 본 스크립트의 (g3ap,g3bp) = config의 (E3-,E3+). 후보를 그 순서로 기재: (...,17/3,5/3). — this script's (g3ap,g3bp) = config's (E3-,E3+). The candidates are listed in that order: (...,17/3,5/3).
cand=[7,1,7,sp.Rational(13,3),sp.Rational(4,3),sp.Rational(17,3),sp.Rational(5,3)]
for k,(nm,lv) in enumerate(zip([names[i] for i in sel],lam)):
    match = sp.simplify(lv-cand[k]/r15)==0
    print(f"  λ_{nm:5} = {lv}  = {sp.nsimplify(lv, [r15])}   후보일치: {match}  양수: {lv>0}")
res=sp.simplify(A*lam-grad_ell)
print("  잔차:", list(res), " (전부 0이어야)")

# L_trunc 정확 계산 — 규약: 문서 11의 교정된 규칙(중심선이 코어 외곽을 "이탈"하는 점까지의 스텁) — exact computation of L_trunc — convention: the corrected rule of doc 11 (the stub up to the point where the centerline "exits" the core outline)
# [정정 — referee-H F1] 이전판은 최외곽 중심선 교차 규칙(오류). n=3에서 두 규칙이 다르고 — [correction — referee-H F1] the previous version used the outermost centerline-intersection rule (wrong). At n=3 the two rules differ and
# 외곽-이탈 규칙만 Kauffman 값 4/√(5−2√5)를 재현함. — only the outline-exit rule reproduces Kauffman's value 4/√(5−2√5).
print("\n[L_trunc — 외곽-이탈 규약 (referee-H F1 정정)]")
def meet(P1,D1,P2,D2):
    s=sp.Matrix([[D1[0],-D2[0]],[D1[1],-D2[1]]]).solve(P2-P1)
    return sp.simplify(sp.expand_trig(s).subs(SUB))
# 육각형 꼭짓점 (인증 완료 좌표) — hexagon vertices (certified coordinates)
V=[sp.Matrix([-13*r15/30,sp.Rational(1,2)]), sp.Matrix([-7*r15/30,sp.Rational(3,2)]),
   sp.Matrix([r15/30,sp.Rational(3,2)]),     sp.Matrix([3*r15/30,sp.Rational(1,2)]),
   sp.Matrix([-3*r15/30,-sp.Rational(1,2)]), sp.Matrix([-11*r15/30,-sp.Rational(1,2)])]
def ray_exit(P0,D0):
    """반직선 P0+tD0 (t>0)가 볼록 육각형을 이탈하는 t (변별 교차 중 내부 통과 후 최초 이탈)
    the t at which the ray P0+tD0 (t>0) exits the convex hexagon (among the per-edge intersections, the first exit after passing through the interior)"""
    ts=[]
    for k in range(6):
        A_,B_=V[k],V[(k+1)%6]
        e=B_-A_
        den=sp.simplify(D0[0]*(-e[1])+D0[1]*e[0])
        if den==0: continue
        # cross(e, P0+tD0-A_)=0
        t=sp.simplify((e[0]*(P0[1]-A_[1])-e[1]*(P0[0]-A_[0]))/(D0[0]*e[1]-D0[1]*e[0]))
        u_pt=P0+t*D0
        # 변 위인지 (파라미터 [0,1]) — whether it lies on the edge (parameter in [0,1])
        uden=sp.simplify(e.dot(e))
        uu=sp.simplify(e.dot(u_pt-A_)/uden)
        if t.is_positive and (uu>=0)==True and (uu<=1)==True:
            ts.append(sp.nsimplify(t,[r15]))
    return sp.Max(*ts) if ts else None
Min=sp.Matrix([0,0]); Din=sp.Matrix([-1,0])           # 들어오는 꼬리: M0에서 -d0 방향 — incoming tail: from M0 in the -d0 direction
Mout=sp.Matrix([-r15/3,1]); Dout=sp.Matrix([1,0])     # 나가는 꼬리: M3에서 d4 방향 — outgoing tail: from M3 in the d4 direction
s_in=ray_exit(Min,Din); s_out=ray_exit(Mout,Dout)
Lc=2*r15/5*2+8*r15/15
Lt=sp.simplify(s_in+Lc+s_out)
print("  스텁(입구,출구) =", s_in, ",", s_out, "  (= L1 = L3 =", sp.nsimplify(2*r15/5,[r15]), ")")
print("  검증 스텁-L1 =", sp.simplify(s_in-2*r15/5), sp.simplify(s_out-2*r15/5))
print("  코어 =", sp.nsimplify(Lc,[r15]))
print("  L_trunc =", sp.nsimplify(Lt,[r15]), "=", sp.N(Lt,30))
print("  비율 L_trunc/ell_core =", sp.nsimplify(sp.simplify(Lt/Lc),[r15]))

# ---- 접근성 방향의 정확화: A7 v = 1 을 Q(√15)에서 풀고 8번째 방향도함수 부호 확인 ---- — making the accessibility direction exact: solve A7 v = 1 in Q(√15) and check the sign of the 8th directional derivative
print("\n[접근성 방향 — 정확]")
A7=J[sel,:]           # 7x7 (g2 제외) — 7x7 (g2 excluded)
ones=sp.Matrix([1]*7)
v_int=sp.simplify(A7.solve(ones))
d8=sp.simplify((J[1,:]*v_int)[0,0])   # ∇g2 · v
print("  v (Q(√15)) =", [sp.nsimplify(vi,[r15]) for vi in v_int])
print("  ∇g2·v =", sp.nsimplify(d8,[r15]), " >0:", sp.simplify(d8)>0)
print("  ⟹ 8개 방향도함수 = (1,1,1,1,1,1,1,", sp.nsimplify(d8,[r15]), ") — 전부 정확 양수" if sp.simplify(d8)>0 else "  조정 필요")
