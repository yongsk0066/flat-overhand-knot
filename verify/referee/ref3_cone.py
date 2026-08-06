# ref3_cone.py — REFEREE: is theta* an ACCESSIBLE boundary point?
#  (1) exact linearized feasible cone C = {v : J v >= 0} at theta*  (8 constraints, 7 dims)
#  (2) does ell strictly increase along every extreme ray?
#  (3) high-precision *nonlinear* feasibility along cone rays: is theta* a limit of valid states?
#  (4) why does the authors' float oracle reject theta*?
import sympy as sp, itertools, math, mpmath as mp
mp.mp.dps = 50
r15 = sp.sqrt(15)
T = sp.symbols('T1 T2 T3 T4', real=True); LL = sp.symbols('L1 L2 L3', positive=True)
TH = list(T)+list(LL)
phi=[sp.Integer(0)]
for t in T: phi.append(phi[-1]+t)
d=[sp.Matrix([sp.cos(p),sp.sin(p)]) for p in phi]
M=[sp.Matrix([0,0])]
for i in range(3): M.append(M[-1]+LL[i]*d[i+1])
SIG=[1,1,-1,-1]; CH=[]
for i in range(4):
    v=d[i]+d[i+1]; n=sp.Matrix([-d[i][1],d[i][0]]); hv=v/(2*SIG[i]*v.dot(n))
    CH.append((M[i]-hv, M[i]+hv))
SUB={sp.cos(T[0]):sp.Rational(-1,4),sp.sin(T[0]):r15/4,sp.cos(T[1]):sp.Rational(-7,8),sp.sin(T[1]):r15/8,
     sp.cos(T[2]):sp.Rational(-7,8),sp.sin(T[2]):-r15/8,sp.cos(T[3]):sp.Rational(-1,4),sp.sin(T[3]):-r15/4,
     LL[0]:6/r15,LL[1]:8/r15,LL[2]:6/r15}
EV=lambda e: sp.radsimp(sp.simplify(sp.expand_trig(sp.expand(e)).subs(SUB)))
PASS={0:(M[0],-d[0],sp.Integer(60)),1:(M[0],d[1],LL[0]),2:(M[1],d[2],LL[1]),3:(M[2],d[3],LL[2]),4:(M[3],d[4],sp.Integer(60))}
lat=lambda p,P,dd: -(p[0]-P[0])*dd[1]+(p[1]-P[1])*dd[0]
ORDER=[1,4,2,0,3]; rank={p:i for i,p in enumerate(ORDER)}
G=[]
for i in range(4):
    ra,rb=rank[i],rank[i+1]; lo,hi=min(ra,rb),max(ra,rb)
    for k in [ORDER[r] for r in range(lo+1,hi)]:
        P,dd,ln=PASS[k]
        s=1 if float(sp.N(EV(lat(CH[i][0],P,dd))))>0 else -1
        for e in (0,1):
            g=s*lat(CH[i][e],P,dd)-sp.Rational(1,2)
            if EV(g)==0: G.append((f"c{i+1}/P{k}/e{e}", g))
J=sp.Matrix([[EV(sp.diff(g,v)) for v in TH] for _,g in G])
gl=sp.Matrix([0,0,0,0,1,1,1])
print("active:", [n for n,_ in G], " J:", J.shape, " rank", J.rank())

# ---- (1) exact cone -------------------------------------------------------
S=[0,2,3,4,5,6,7]                       # the 7-subset the paper uses (drops c1/P2/e0)
Js=J[S,:]; W=Js.inv()                    # columns = extreme rays of {v: J_S v >= 0}
rays=[sp.simplify(W[:,k]) for k in range(7)]
k8=[k for k in range(8) if k not in S][0]
print(f"\n[cone] extreme rays of the 7-subset cone, tested against the 8th constraint {G[k8][0]}:")
a8=[sp.simplify((J[k8,:]*rays[k])[0]) for k in range(7)]
dl=[sp.simplify((gl.T*rays[k])[0]) for k in range(7)]
for k in range(7):
    print(f"  ray {k} (increases only {G[S[k]][0]:14}):  grad_g8.ray = {sp.nsimplify(sp.radsimp(a8[k]))}"
          f"   d(ell) = {sp.nsimplify(sp.radsimp(dl[k]))}")
neg=[k for k in range(7) if a8[k]<0]; pos=[k for k in range(7) if a8[k]>0]
print(f"  -> 8th constraint is violated by rays {neg}, satisfied by {pos}")
if not pos and neg: print("  *** CONE IS {0}: theta* is a first-order ISOLATED feasible point ***")
else: print(f"  cone is FULL-DIMENSIONAL (nonempty interior): {len(pos)>0 and len(neg)<7}")
print("  d(ell) along every extreme ray of the 7-cone:",
      "ALL STRICTLY POSITIVE" if all(x>0 for x in dl) else f"NOT ALL POSITIVE: {dl}")

# an interior direction of the full 8-constraint cone (if any)
vint=None
for combo in itertools.product(range(3),repeat=7):
    if sum(combo)==0: continue
    v=sum((sp.Integer(c)*rays[k] for k,c in enumerate(combo)), sp.zeros(7,1))
    if sp.simplify((J[k8,:]*v)[0])>0:
        vint=sp.simplify(v); break
print("  interior direction of the FULL cone found:", vint is not None)
if vint is not None:
    print("   v =", [sp.N(x,8) for x in vint], "  J v =", [sp.N(sp.simplify((J[k,:]*vint)[0]),6) for k in range(8)])
    print("   d(ell) along it =", sp.N(sp.simplify((gl.T*vint)[0]),10))

# ---- (3) high-precision NONLINEAR feasibility along that ray -------------
print("\n[nonlinear] exact-arithmetic axiom test at theta* + r*v  (mpmath 50 digits)")
tst=mp.acos(mp.mpf(-1)/4); sst=2*mp.pi-2*tst
TH0=[tst,sst,-sst,-tst,6/mp.sqrt(15),8/mp.sqrt(15),6/mp.sqrt(15)]
def axiom_mp(th):
    """exact-arithmetic version of axiom_exact_ok: no float tolerance games."""
    phi=[mp.mpf(0)]
    for t in th[:4]: phi.append(phi[-1]+t)
    dd=[(mp.cos(p),mp.sin(p)) for p in phi]
    MM=[(mp.mpf(0),mp.mpf(0))]
    for i in range(3): MM.append((MM[-1][0]+th[4+i]*dd[i+1][0], MM[-1][1]+th[4+i]*dd[i+1][1]))
    CHm=[]
    for i in range(4):
        vx,vy=dd[i][0]+dd[i+1][0], dd[i][1]+dd[i+1][1]
        nx,ny=-dd[i][1],dd[i][0]; den=2*(vx*nx+vy*ny)
        hx,hy=vx/abs(den)*mp.sign(den)*mp.sign(den), vy/abs(den)  # keep magnitude convention
        hx,hy=vx/abs(den), vy/abs(den)
        CHm.append(((MM[i][0]-hx,MM[i][1]-hy),(MM[i][0]+hx,MM[i][1]+hy)))
    PS={0:(MM[0],(-dd[0][0],-dd[0][1]),mp.mpf(60)),1:(MM[0],dd[1],th[4]),2:(MM[1],dd[2],th[5]),
        3:(MM[2],dd[3],th[6]),4:(MM[3],dd[4],mp.mpf(60))}
    worst=mp.mpf(0)
    for i in range(4):
        ra,rb=rank[i],rank[i+1]; lo,hi=min(ra,rb),max(ra,rb)
        for k in [ORDER[r] for r in range(lo+1,hi)]:
            P,ddk,ln=PS[k]
            A,B=CHm[i]
            la=lambda p: -(p[0]-P[0])*ddk[1]+(p[1]-P[1])*ddk[0]
            lo_=lambda p: (p[0]-P[0])*ddk[0]+(p[1]-P[1])*ddk[1]
            lA,lB=la(A),la(B); nA,nB=lo_(A),lo_(B)
            # depth of penetration into the OPEN band, sampled on a fine grid + endpoints
            for j in range(0,201):
                f=mp.mpf(j)/200
                L=lA+(lB-lA)*f; N=nA+(nB-nA)*f
                if 0<N<ln:
                    dep=mp.mpf(0.5)-abs(L)
                    if dep>worst: worst=dep
    return worst   # >0 means penetration depth (violation), <=0 means OK
print(f"   at theta* itself: max penetration depth = {mp.nstr(axiom_mp(TH0),5)}  (<=0 means axiom OK)")
if vint is not None:
    vv=[mp.mpf(str(sp.N(x,40))) for x in vint]; nv=mp.sqrt(sum(x*x for x in vv)); vv=[x/nv for x in vv]
    for r in ('1e-2','1e-3','1e-4','1e-6','1e-8'):
        rr=mp.mpf(r); th2=[TH0[j]+rr*vv[j] for j in range(7)]
        dep=axiom_mp(th2); dell=sum(th2[4:])-sum(TH0[4:])
        print(f"   r={r}: max penetration = {mp.nstr(dep,4):>12}   d(ell) = {mp.nstr(dell,6)}"
              f"   {'FEASIBLE' if dep<=0 else 'violates'}")

# ---- (4) why the float oracle rejects theta* -----------------------------
print("\n[float-oracle forensics] authors' axiom_exact_ok at theta*")
import sys; sys.path.insert(0,'/Users/yongseok/dev/pentagon/verify')
from n4_tighten2 import config, pieces
t0=math.acos(-0.25); s0=2*math.pi-2*t0
TH0f=(t0,s0,-s0,-t0,6/math.sqrt(15),8/math.sqrt(15),6/math.sqrt(15))
c=config(TH0f); segs=pieces(c)
for i in range(4):
    ra,rb=rank[i],rank[i+1]; lo,hi=min(ra,rb),max(ra,rb)
    for k in [ORDER[r] for r in range(lo+1,hi)]:
        P,ddk,ln=segs[k]; (ax,ay),(bx,by)=c['CH'][i]
        n=(-ddk[1],ddk[0])
        lA=(ax-P[0])*n[0]+(ay-P[1])*n[1]; lB=(bx-P[0])*n[0]+(by-P[1])*n[1]
        print(f"   crease{i+1} vs P{k}: latA={lA!r} latB={lB!r}  |dlat|={abs(lB-lA):.3e}"
              f"  {'DEGENERATE-BRANCH' if abs(lB-lA)<1e-15 else ''}"
              f"  {'<0.5 -> flagged as INSIDE slab' if abs(lB-lA)<1e-15 and -0.5<lA<0.5 else ''}")
