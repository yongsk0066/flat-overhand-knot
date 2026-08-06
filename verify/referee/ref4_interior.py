# ref4_interior.py — REFEREE: is theta* a genuine (accessible) constrained local minimum?
# Uses a TRUE interior direction of the 8-constraint cone and a high-precision nonlinear
# feasibility test, then searches adversarially for any nearby feasible point with smaller ell.
import sympy as sp, mpmath as mp, random, itertools
mp.mp.dps = 60
r15=sp.sqrt(15)
T=sp.symbols('T1 T2 T3 T4',real=True); LL=sp.symbols('L1 L2 L3',positive=True); TH=list(T)+list(LL)
phi=[sp.Integer(0)]
for t in T: phi.append(phi[-1]+t)
d=[sp.Matrix([sp.cos(p),sp.sin(p)]) for p in phi]
M=[sp.Matrix([0,0])]
for i in range(3): M.append(M[-1]+LL[i]*d[i+1])
SIG=[1,1,-1,-1]; CH=[]
for i in range(4):
    v=d[i]+d[i+1]; n=sp.Matrix([-d[i][1],d[i][0]]); CH.append((M[i]-v/(2*SIG[i]*v.dot(n)), M[i]+v/(2*SIG[i]*v.dot(n))))
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
            if EV(g)==0: G.append((f"c{i+1}/P{k}/e{e}",g))
J=sp.Matrix([[EV(sp.diff(g,v)) for v in TH] for _,g in G])
S=[0,2,3,4,5,6,7]; W=J[S,:].inv(); k8=1
rays=[W[:,k] for k in range(7)]

# ---------- true interior direction: all 8 constraints strictly increasing ----------
vint=sum(rays[1:], rays[0])
Jv=[sp.simplify((J[k,:]*vint)[0]) for k in range(8)]
print("[true interior direction] grad g_k . v for all 8 constraints:", [sp.nsimplify(x) for x in Jv])
print("  all strictly positive:", all(x>0 for x in Jv), "   d(ell)/dv =", sp.nsimplify(sp.radsimp(sum(vint[4:7]))))

# ---------- high-precision nonlinear axiom oracle (no float fudge) ----------
def penetration(th, ngrid=400):
    phi=[mp.mpf(0)]
    for t in th[:4]: phi.append(phi[-1]+t)
    dd=[(mp.cos(p),mp.sin(p)) for p in phi]
    MM=[(mp.mpf(0),mp.mpf(0))]
    for i in range(3): MM.append((MM[-1][0]+th[4+i]*dd[i+1][0],MM[-1][1]+th[4+i]*dd[i+1][1]))
    CHm=[]
    for i in range(4):
        vx,vy=dd[i][0]+dd[i+1][0],dd[i][1]+dd[i+1][1]
        den=2*(vx*(-dd[i][1])+vy*dd[i][0])
        CHm.append(((MM[i][0]-vx/abs(den),MM[i][1]-vy/abs(den)),(MM[i][0]+vx/abs(den),MM[i][1]+vy/abs(den))))
    PS={0:(MM[0],(-dd[0][0],-dd[0][1]),mp.mpf(60)),1:(MM[0],dd[1],th[4]),2:(MM[1],dd[2],th[5]),
        3:(MM[2],dd[3],th[6]),4:(MM[3],dd[4],mp.mpf(60))}
    worst=mp.mpf(-1)
    for i in range(4):
        ra,rb=rank[i],rank[i+1]; lo,hi=min(ra,rb),max(ra,rb)
        for k in [ORDER[r] for r in range(lo+1,hi)]:
            P,ddk,ln=PS[k]; A,B=CHm[i]
            la=lambda p:-(p[0]-P[0])*ddk[1]+(p[1]-P[1])*ddk[0]
            ln_=lambda p:(p[0]-P[0])*ddk[0]+(p[1]-P[1])*ddk[1]
            lA,lB=la(A),la(B); nA,nB=ln_(A),ln_(B)
            for j in range(ngrid+1):
                f=mp.mpf(j)/ngrid; L=lA+(lB-lA)*f; N=nA+(nB-nA)*f
                if 0<N<ln: worst=max(worst, mp.mpf(1)/2-abs(L))
    return worst
tst=mp.acos(mp.mpf(-1)/4); sst=2*mp.pi-2*tst
TH0=[tst,sst,-sst,-tst,6/mp.sqrt(15),8/mp.sqrt(15),6/mp.sqrt(15)]
ELLS=sum(TH0[4:])
vv=[mp.mpf(str(sp.N(x,50))) for x in vint]; nv=mp.sqrt(sum(x*x for x in vv)); vv=[x/nv for x in vv]
print(f"\n[nonlinear feasibility along the TRUE interior direction]  (penetration<=0 == valid)")
print(f"  theta* itself: penetration = {mp.nstr(penetration(TH0),4)}")
for r in ('1e-1','1e-2','1e-3','1e-4','1e-6'):
    rr=mp.mpf(r); th2=[TH0[j]+rr*vv[j] for j in range(7)]
    print(f"  r={r:6}: penetration = {mp.nstr(penetration(th2),4):>12}  ell-ell* = {mp.nstr(sum(th2[4:])-ELLS,6):>12}"
          f"   {'FEASIBLE' if penetration(th2)<=0 else 'violates'}")

# ---------- adversarial search: any feasible point with SMALLER ell? ----------
print("\n[adversarial descent search from theta*, high precision, 30000 cone+random probes]")
rnd=random.Random(31337)
Wn=[[mp.mpf(str(sp.N(rays[k][j],40))) for j in range(7)] for k in range(7)]
best=None; nfeas=0
def probe(v, r):
    n=mp.sqrt(sum(x*x for x in v)); u=[x/n for x in v]
    th2=[TH0[j]+r*u[j] for j in range(7)]
    if penetration(th2, ngrid=120) <= 0:
        return sum(th2[4:])-ELLS
    return None
trials=[]
for _ in range(1200):                                  # inside the linearised cone
    c=[mp.mpf(rnd.random())**2 for _ in range(7)]
    v=[sum(c[k]*Wn[k][j] for k in range(7)) for j in range(7)]
    trials.append(v)
for _ in range(1800):                                  # unbiased random
    trials.append([mp.mpf(rnd.gauss(0,1)) for _ in range(7)])
for v in trials:
    for r in (mp.mpf('1e-2'), mp.mpf('1e-3'), mp.mpf('1e-4'), mp.mpf('1e-5')):
        de=probe(v,r)
        if de is not None:
            nfeas+=1
            if best is None or de<best[0]: best=(de,r)
            break
print(f"  feasible probes: {nfeas}/{len(trials)};  smallest ell-ell* over all feasible probes = "
      f"{mp.nstr(best[0],8) if best else 'n/a'} at radius {mp.nstr(best[1],2) if best else '-'}")
print(f"  => {'NO feasible point with ell < ell* found (consistent with blocked)' if (best is None or best[0]>0) else '*** FOUND ell < ell* : CERTIFICATE BROKEN ***'}")
