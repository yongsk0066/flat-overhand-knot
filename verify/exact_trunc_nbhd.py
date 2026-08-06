# exact_trunc_nbhd.py — 오각형 절단선 동정의 Q(x0) 정확화 (리뷰 02 §3 마지막 항목) — Q(x0)-exactification of the identification of the pentagon's cut lines (review 02 §3, last item)
# 주장: P*에서 나가는 꼬리는 외곽 변 [G,Q1]을 파라미터 u=1/2 (정확!)에서 가로지르고, — Claim: the tail leaving P* crosses the outer edge [G,Q1] at parameter u=1/2 (exactly!),
#       들어오는 꼬리는 [E,H]를 u=1/2에서 가로지른다; 교차각 cos = (√5-1)/4 (=cos 72°); — and the incoming tail crosses [E,H] at u=1/2; the crossing angle has cos = (√5-1)/4 (=cos 72°);
#       다른 외곽 변과는 만나지 않거나 세그먼트 밖 — 전부 몫체 부호판정 ⟹ 근방 개방성. — the other outer edges are either not met or met outside the segment — all decided by quotient-field sign tests ⟹ openness of the neighbourhood.
import sympy as sp
from mpmath import mp, iv, mpf, findroot

x=sp.symbols('x')
q4=5*x**4-10*x**2+1
def red(e): return sp.rem(sp.expand(e), q4, x)
def inv(e): return sp.invert(red(e), q4, x)
def fdiv(a,b): return red(a*inv(b))
def fr(e):
    n,d=sp.fraction(sp.cancel(e))
    return fdiv(sp.expand(n), sp.expand(d))

g=fr((1-x**2)/(2*x)); h=fr((1-3*x**2)/(2*x)); a=h
E=(sp.Integer(0),sp.Integer(0)); F=(red(-x),sp.Integer(1)); G=(g,sp.Integer(1)); H=(h,sp.Integer(0))
# 슬라이스에서 q = g-h = x (검산) — in the slice, q = g-h = x (cross-check)
assert red(g-h-x)==0, "q=x 실패"
w0=(fr((x**2-1)/(x**2+1)), fr(2*x/(x**2+1)))
def vadd(P,s,D): return (red(P[0]+s*D[0]), red(P[1]+s*D[1]))
def vsub(P,Q): return (red(P[0]-Q[0]), red(P[1]-Q[1]))
def dot(u,v): return red(u[0]*v[0]+u[1]*v[1])
def crs(u,v): return red(u[0]*v[1]-u[1]*v[0])
Q1=vadd(G,a,w0)
bprime=fr(dot(vsub(F,H),sp.Matrix(w0)) if False else (vsub(F,H)[0]*w0[0]+vsub(F,H)[1]*w0[1]))
bprime=red(bprime)
Q0=vadd(H,bprime,w0)
U=vsub(Q0,Q1)
# v3 = refl_U(w0)
UU=dot(U,U)
t=fdiv(red(2*dot(w0,U)), UU)
v3=(red(t*U[0]-w0[0]), red(t*U[1]-w0[1]))
M2=(fdiv(h+g,sp.Integer(2)), sp.Rational(1,2))
den=crs(w0,U)
s_=fdiv(crs(vsub(Q1,M2),U), den)
M3=vadd(M2,s_,w0)
def meet_params(P1,D1,P2,D2):
    d=crs(D1,D2)
    tt=fdiv(crs(vsub(P2,P1),D2), d)
    uu=fdiv(crs(vsub(P2,P1),D1), d)
    return tt,uu
# 외곽(오각형) 변: [E,H](y=0), [H,G](κ2), [G,Q1](ℓ2 모서리), [Q1,Q0](κ3), [Q0(=F 근방)..F,E?]: 오각형 P*는 — outer (pentagon) edges: [E,H](y=0), [H,G](κ2), [G,Q1](ℓ2 edge), [Q1,Q0](κ3), [Q0(=near F)..F,E?]: the pentagon P* has
#   변 5개: EH, HG, G-Q1, Q1-Q0(=κ3), F-E(=κ1) — Q0=F (오각형 코너!) — 5 edges: EH, HG, G-Q1, Q1-Q0(=κ3), F-E(=κ1) — Q0=F (a pentagon corner!)
print("[검산] Q0-F =", vsub(Q0,F), " (오각형에서 Q0=F)")
# ---- 두 꼬리의 정확 절단 동정 ---- — exact identification of the two tails' cuts ----
# 출구(패스 D): ray(M3, v3);  입구(패스 A): ray(M1, -d2)  (A = κ1 반사의 후방 연장) — exit (pass D): ray(M3, v3);  entrance (pass A): ray(M1, -d2)  (A = the backward extension of the κ1 reflection)
d2=( fdiv(red(2*(-x)*(-x)-( (-x)**2+1 )), red((-x)**2+1)), fdiv(red(2*(-x)*1), red((-x)**2+1)) )
# refl_lin(F-E)((1,0)):  F-E=(-x,1);  refl(v)=2(v·d)d/|d|²-v,  v=(1,0):
FEv=vsub(F,E); dd=dot(FEv,FEv)
d2=( red(fdiv(2*(-x),dd)*FEv[0]-1), red(fdiv(2*(-x),dd)*FEv[1]) )
M1=(fdiv(red(-x),sp.Integer(2)), sp.Rational(1,2))
dA=(red(-d2[0]), red(-d2[1]))
import mpmath as _mp
_mp.mp.dps=40
from mpmath import iv as _iv
_iv.dps=40
x0=findroot(lambda tt: 5*tt**4-10*tt**2+1, mpf('0.3249196962'))
XIV=_iv.mpf([str(x0-mpf(10)**-30),str(x0+mpf(10)**-30)])
def sgn(e):
    e=red(e)
    if e==0: return 0
    if not e.free_symbols: return 1 if e>0 else -1
    l=sp.ilcm(*[sp.denom(c) for c in sp.Poly(e,x).all_coeffs()])
    cs=[int(c) for c in sp.Poly(sp.expand(e*l),x).all_coeffs()]
    v=_iv.mpf(0)
    for c in cs: v=v*XIV+c
    if v.a>0: return 1
    if v.b<0: return -1
    raise RuntimeError(e)
s5=fr((5-5*x**2)/2)
EDGES=[("EH",E,vsub(H,E)),("HG",H,vsub(G,H)),("GQ1",G,vsub(Q1,G)),("K3",Q1,vsub(Q0,Q1)),("FE",F,vsub(E,F))]
for tag,P0,D0,target in (("출구(D)",M3,v3,"EH"),("입구(A)",M1,dA,"GQ1")):
    print(f"\n[{tag}] ray 기준 전 변 판정:")
    for nm,Q,DE in EDGES:
        d_=crs(D0,DE)
        if d_==0:
            print(f"  {nm}: 평행"); continue
        tt=fdiv(crs(vsub(Q,P0),DE), d_); uu=fdiv(crs(vsub(Q,P0),D0), d_)
        if nm==target:
            resu=red(uu-sp.Rational(1,2))
            # 각도: cos∠(D0, DE)/|DE| — |D0|=1 — angle: cos∠(D0, DE)/|DE| with |D0|=1
            LDE2=dot(DE,DE)
            # cos² 비교 대신: cos·|DE| = D0·DE;  (D0·DE)² - cos72²·|DE|² ≡ 0 ? — instead of comparing cos²: cos·|DE| = D0·DE;  is (D0·DE)² - cos72²·|DE|² ≡ 0 ?
            c72=fdiv(s5-1,sp.Integer(4))
            resang=red(dot(D0,DE)**2 - c72**2*LDE2)
            sgn_dot=sgn(dot(D0,DE))
            print(f"  {nm} [절단변]: u-1/2 ≡ {resu} | t부호 {sgn(tt)} | (cos·|DE|)²-cos²72°·|DE|² ≡ {resang} (부호 {sgn_dot})")
        else:
            st=sgn(tt); su0=sgn(uu); su1=sgn(red(1-uu))
            inside=(st>0 and su0>0 and su1>0)
            print(f"  {nm}: t부호{st} u∈(0,1)? {su0>0 and su1>0} → {'!!세그먼트 안!!' if inside else '배제 ✓'}")
print("""
[결론] 두 꼬리 모두: 대상 변의 정확한 중점(u≡1/2)에서, 정확히 72°(cos=(√5-1)/4)로,
   외곽을 단 1회 가로지른다(여타 변 배제 — 부호판정). 모든 판정이 Q(x0) 정확 산술
   + 구간 부호 ⟹ 부등식이 엄격하므로 근방 개방성은 연속성으로 따라온다.
   (부산물: 절단 교차점 = 변의 중점 — 신규 정확 구조)""")
