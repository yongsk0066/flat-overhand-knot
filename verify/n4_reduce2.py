# n4_reduce2.py — 접촉 4등식 완전계 뉴턴 (0차원) + PSLQ — Newton on the complete system of 4 contact equations (0-dimensional) + PSLQ
# g1: E1+ ∈ P4 하모서리 (⟺ M3.y=1),  g2: E1- ∈ P2 모서리, — g1: E1+ on the lower edge of P4 (⟺ M3.y=1); g2: E1- on the P2 edge;
# g3a,g3b: E2± ∈ P4 상모서리 (⟺ 크리스2 평행밀착: s=2π-2t & L1 sin t=3/2) — g3a,g3b: E2± on the upper edge of P4 (⟺ crease 2 flush and parallel: s=2π-2t & L1 sin t=3/2)
import mpmath as mp, json
mp.mp.dps=60

def F(v):
    t,s,L1,L2=v
    d1=[mp.cos(t),mp.sin(t)]; d2=[mp.cos(t+s),mp.sin(t+s)]
    M1=[L1*d1[0],L1*d1[1]]; M2=[M1[0]+L2*d2[0],M1[1]+L2*d2[1]]
    M3=[M2[0]+L1*d1[0],M2[1]+L1*d1[1]]
    rt=1/(2*mp.sin(t/2)); rs=1/(2*mp.sin(s/2))
    et=[mp.cos(t/2),mp.sin(t/2)]; es=[mp.cos(t+s/2),mp.sin(t+s/2)]
    E1p=[rt*et[0],rt*et[1]]; E1m=[-rt*et[0],-rt*et[1]]
    E2p=[M1[0]+rs*es[0],M1[1]+rs*es[1]]; E2m=[M1[0]-rs*es[0],M1[1]-rs*es[1]]
    lat2=lambda pt: -(pt[0]-M1[0])*d2[1]+(pt[1]-M1[1])*d2[0]
    g1=-(E1p[1]-M3[1])-mp.mpf(1)/2          # s1=-1: M3.y-E1p.y=1/2
    g2=lat2(E1m)-mp.mpf(1)/2
    g3a=(E2p[1]-M3[1])-mp.mpf(1)/2
    g3b=(E2m[1]-M3[1])-mp.mpf(1)/2
    return mp.matrix([g1,g2,g3a,g3b])

st=json.load(open("verify/n4_tight41.json"))
th=st["theta"]
v0=(mp.mpf((th[0]-th[3])/2), mp.mpf((th[1]-th[2])/2), mp.mpf((th[4]+th[6])/2), mp.mpf(th[5]))
print("초기 잔차:",[mp.nstr(x,3) for x in F(v0)])
sol=mp.findroot(lambda a,b,c,d: F((a,b,c,d)), v0, tol=mp.mpf(10)**(-50))
t,s,L1,L2=sol
ell=2*L1+L2
print("\n[정확 결정계 해] (dps60)")
for n_,v_ in (("t*",t),("s*",s),("L1*",L1),("L2*",L2),("ell*",ell)):
    print(f"  {n_:4}= {mp.nstr(v_,45)}")
print("\n[구조 항등 검증]")
print("  s*-(2π-2t*)     =",mp.nstr(s-(2*mp.pi-2*t),8))
print("  L1* sin t* -3/2 =",mp.nstr(L1*mp.sin(t)-mp.mpf(3)/2,8))
print("  M3.y-1          =",mp.nstr(2*L1*mp.sin(t)+L2*mp.sin(t+s)-1,8))
print("\n[대수 동정]")
ct=mp.cos(t)
print("  cos t*          =",mp.nstr(ct,40))
print("  cos t*+1/4      =",mp.nstr(ct+mp.mpf(1)/4,8))
for name,val in [("cos t*",ct),("L1*",L1),("L2*",L2),("ell*",ell),("sin t*",mp.sin(t))]:
    for deg in (1,2,3,4,6,8):
        rel=mp.pslq([val**k for k in range(deg+1)],maxcoeff=10**10,maxsteps=10**6)
        if rel:
            print(f"  {name}: 최소다항 후보 차수{deg}: {rel}")
            break
    else: print(f"  {name}: 관계 미발견")
json.dump({k:str(v) for k,v in [("t",t),("s",s),("L1",L1),("L2",L2),("ell",ell)]},
          open("verify/n4_exact_theta.json","w"))
