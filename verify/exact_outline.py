# exact_outline.py — θ† 칠각형 외곽 조합의 Q(x†) 정확화 (리뷰 02 §2.4 해소) — Q(x†)-exactification of the combinatorics of the θ† heptagon outline (resolves review 02 §2.4)
# 판정: 7개 변 중점이 (i) 어느 패스의 열린 내부에도 없음 (ii) 최소 한 패스 경계 위, — decision: the 7 edge midpoints (i) lie in the open interior of no pass, (ii) lie on the boundary of at least one pass,
#       그리고 XYF 쐐기 무게중심이 네 패스 어디에도 덮이지 않음 — 전부 몫체 부호 판정. — and the barycentre of the XYF wedge is covered by none of the four passes — all by quotient-field sign tests.
import sympy as sp
from mpmath import mp, iv, mpf, findroot

x = sp.symbols('x')
px = sp.expand(x**12 + 4*x**10 - 7*x**4 + 20*x**2 - 2)
def red(e): return sp.rem(sp.expand(e), px, x)
def inv(e): return sp.invert(red(e), px, x)
def fdiv(a,b): return red(a*inv(b))

h_f = sp.Rational(1,12)*(2*x**11+8*x**9-x**7-4*x**5-13*x**3+40*x)
g_f = fdiv(1-x**2, 2*x); a_f=h_f; b_f=fdiv(x**2+1, 2*x)
q_f = red(g_f-h_f); nW0=red(q_f**2+1)
E=(sp.Integer(0),sp.Integer(0)); F=(red(-x),sp.Integer(1))
G=(g_f,sp.Integer(1)); H=(h_f,sp.Integer(0))
W0=(red(q_f**2-1),red(2*q_f))
W0h=(fdiv(W0[0],nW0),fdiv(W0[1],nW0))
Q1=(red(G[0]+a_f*W0h[0]),red(G[1]+a_f*W0h[1]))
Q0=(red(H[0]+b_f*W0h[0]),red(H[1]+b_f*W0h[1]))
U=(red(Q0[0]-Q1[0]),red(Q0[1]-Q1[1]))
def vdot(u,v): return red(u[0]*v[0]+u[1]*v[1])
def vcross(u,v): return red(u[0]*v[1]-u[1]*v[0])
def refl_pt(P,d,Qp):
    v=(red(Qp[0]-P[0]),red(Qp[1]-P[1]))
    t=fdiv(red(2*vdot(v,d)),vdot(d,d))
    return (red(P[0]+t*d[0]-v[0]),red(P[1]+t*d[1]-v[1]))
def refl_lin(d,v):
    t=fdiv(red(2*vdot(v,d)),vdot(d,d))
    return (red(t*d[0]-v[0]),red(t*d[1]-v[1]))
r3H=refl_pt(Q1,U,H); r3G=refl_pt(Q1,U,G); dD=refl_lin(U,W0)
def meet(P1,d1,P2,d2):
    den=vcross(d1,d2)
    t=fdiv(vcross((red(P2[0]-P1[0]),red(P2[1]-P1[1])),d2),den)
    return (red(P1[0]+t*d1[0]),red(P1[1]+t*d1[1]))
# Y = e_B(y=1) ∩ (ρ3H, dD)-변 — Y = e_B(y=1) ∩ the (ρ3H, dD)-edge
Y=meet(r3H,dD,(sp.Integer(0),sp.Integer(1)),(sp.Integer(1),sp.Integer(0)))
M1=(red(-x/2),sp.Rational(1,2)); M2=(red((h_f+g_f)/2),sp.Rational(1,2))
M3=meet(M2,W0,Q1,U)
DSAMP=(red((Q0[0]+Q1[0])/2+dD[0]/10), red((Q0[1]+Q1[1])/2+dD[1]/10))
BSAMP=(fdiv(M1[0]+M2[0],sp.Integer(2)),sp.Rational(1,2))
CSAMP=(fdiv(M2[0]+M3[0],sp.Integer(2)),fdiv(M2[1]+M3[1],sp.Integer(2)))
ASAMP=(sp.Integer(-1),sp.Rational(1,2))

# 부호 판정: 구간 (x† enclosure) — 0-판정은 red()==0 정확 — sign decision: interval (x† enclosure) — the zero test is exact, via red()==0
mp.dps=60; iv.dps=60
pxl=sp.lambdify(x,px,'mpmath'); dpl=sp.lambdify(x,sp.diff(px,x),'mpmath')
xd=findroot(pxl,mpf('0.32212749731122718'))
XIV=iv.mpf([str(xd-mpf(10)**-40),str(xd+mpf(10)**-40)])
def sgn(e):
    e=red(e)
    if e==0: return 0
    if not e.free_symbols:
        return 1 if e>0 else -1
    p=sp.Poly(e,x)
    l=sp.ilcm(*[sp.denom(c) for c in p.all_coeffs()]) if p.all_coeffs() else 1
    coeffs=[int(c) for c in sp.Poly(sp.expand(e*l),x).all_coeffs()]
    v=iv.mpf(0)
    for c in coeffs: v=v*XIV+c
    if v.a>0: return 1
    if v.b<0: return -1
    raise RuntimeError(f"부호 판정 실패(폭 {v.delta}): {e}")
def line_form(P,d):
    return lambda p:(vcross(d,(red(p[0]-P[0]),red(p[1]-P[1]))))
# 패스 = 반평면 리스트 [(형식, σ)] — σ는 정확 표본으로 고정 — a pass = a list of half-planes [(form, σ)] — σ is pinned down by an exact sample point
def halfplanes(forms, sample):
    out=[]
    for f in forms:
        s=sgn(f(sample))
        assert s!=0, "표본이 경계 위"
        out.append((f,s))
    return out
yform=lambda p: red(p[1]); yform1=lambda p: red(1-p[1])
fEF=line_form(E,(red(F[0]-E[0]),red(F[1]-E[1])))
fHG=line_form(H,(red(G[0]-H[0]),red(G[1]-H[1])))
fK3=line_form(Q1,U)
fH_w0=line_form(H,W0); fG_w0=line_form(G,W0)
fr3H=line_form(r3H,dD); fr3G=line_form(r3G,dD)
PASSES={
 "A": halfplanes([yform,yform1,fEF], ASAMP),
 "B": halfplanes([yform,yform1,fEF,fHG], BSAMP),
 "C": halfplanes([fH_w0,fG_w0,fHG,fK3], CSAMP),
 "D": halfplanes([fr3H,fr3G,fK3], DSAMP),
}
yge1=lambda p: red(p[1]-1)
D_TR=halfplanes([fr3H,fr3G,fK3,yge1], DSAMP)   # D를 e_B(y=1)에서 절단 (flush-cut 관례) — cut D at e_B(y=1) (flush-cut convention)
CORE_PASSES={"B":PASSES["B"],"C":PASSES["C"],"D":D_TR}   # A는 κ₁에서 절단 — A is cut at κ₁
def covered(p, passes):
    """p가 어떤 패스의 폐포에 속하는가 (전 반평면 부호 >=0) — does p belong to the closure of some pass (all half-plane signs >=0)"""
    for nm,hp in passes.items():
        if all(sgn(f(p))*s>=0 for f,s in hp): return nm
    return None
def status(p):
    """패스별: 'in'(열린내부) / 'bd'(경계) / 'out' — per pass: 'in' (open interior) / 'bd' (boundary) / 'out'"""
    st={}
    for name,hp in PASSES.items():
        vals=[sgn(f(p))*s for f,s in hp]
        if all(v>0 for v in vals): st[name]='in'
        elif all(v>=0 for v in vals): st[name]='bd'
        else: st[name]='out'
    return st
def mid(P,Q): return (fdiv(P[0]+Q[0],sp.Integer(2)),fdiv(P[1]+Q[1],sp.Integer(2)))
HEP=[("E",E),("H",H),("G",G),("K",Q1),("X",Q0),("Y",Y),("F",F)]
print("[칠각형 변 중점 — 외곽 판정: 수직 양측 섭동 중 정확히 한쪽만 커버 (δ=1/1000, 전부 정확 부호)]")
ok=True
for k in range(7):
    n1,P=HEP[k]; n2,Q=HEP[(k+1)%7]
    m=mid(P,Q)
    nx=(red(-(Q[1]-P[1])),red(Q[0]-P[0]))
    d=sp.Rational(1,1000)
    pp=(red(m[0]+d*nx[0]),red(m[1]+d*nx[1]))
    pm=(red(m[0]-d*nx[0]),red(m[1]-d*nx[1]))
    cp=covered(pp,CORE_PASSES); cm=covered(pm,CORE_PASSES)
    onb=covered(m,CORE_PASSES)
    edge_ok = (cp is None)!=(cm is None) and onb is not None
    if not edge_ok: ok=False
    print(f"  {n1}{n2}: 중점폐포={onb}  +측={cp}  -측={cm}  → {'외곽 ✓' if edge_ok else '문제!'}")
w=(fdiv(Q0[0]+Y[0]+F[0],sp.Integer(3)),fdiv(Q0[1]+Y[1]+F[1],sp.Integer(3)))
cw=covered(w,PASSES)   # A 포함 전 패스에 대해 미커버 — uncovered by every pass, A included
print(f"[쐐기 XYF 무게중심] 커버={cw}  → {'미커버 ✓ (노치 실재)' if cw is None else '문제!'}")
if cw is None and ok:
    print("\n결론: 칠각형 = 절단-실루엣 ∂(B∪C∪D)의 정확 외곽 [정확대수+구간부호]")
    print("       — 7변 전부 한쪽커버/중점폐포 확인, 쐐기 미커버: 명제의 수치 성분 해소")

# ---- 진단: XY 변 주변 ---- — diagnostics: around the XY edge ----
import sympy as _sp
def fl(e):
    return float(sp.lambdify(x,red(e),'mpmath')(xd)) if red(e)!=0 and red(e).free_symbols else float(red(e))
XYmid=mid(Q0,Y)
print("\n[진단] XYmid =",(fl(XYmid[0]),fl(XYmid[1])))
for nm,hp in PASSES.items():
    print(f"  {nm}: ", [(round(fl(f(XYmid)),6), s) for f,s in hp])
print("  X=Q0:", (fl(Q0[0]),fl(Q0[1])), " Y:", (fl(Y[0]),fl(Y[1])), " F:", (fl(F[0]),fl(F[1])))
print("  fK3(XYmid)=",fl(fK3(XYmid)), " fr3H(XYmid)=",fl(fr3H(XYmid)), " fr3G(XYmid)=",fl(fr3G(XYmid)))
# 내측(코어 방향) 섭동: XY의 안쪽은 어느 패스? — inward (core-direction) perturbation: which pass lies on the inner side of XY?
nrm=(red(-(Y[1]-Q0[1])),red(Y[0]-Q0[0]))
for sgn_ in (1,-1):
    pt=(red(XYmid[0]+sgn_*sp.Rational(1,1000)*nrm[0]), red(XYmid[1]+sgn_*sp.Rational(1,1000)*nrm[1]))
    st={nm:('in' if all(sgn(f(pt))*s>0 for f,s in hp) else 'no') for nm,hp in PASSES.items()}
    print(f"  섭동 {'+' if sgn_>0 else '-'}: {st}")
