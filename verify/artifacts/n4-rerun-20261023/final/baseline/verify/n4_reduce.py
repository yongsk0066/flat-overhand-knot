# n4_reduce.py — C2 대칭 축소 (t,s,L1,L2) 고정밀 해 + PSLQ 동정 — C2-symmetric reduction: high-precision solution in (t,s,L1,L2) + PSLQ identification
# 접촉(수치 상태에서 판독): 6개 = 대칭 3쌍 — contacts (read off from the numerical state): 6 = 3 symmetric pairs
#   g1: 크리스1 f=1 끝점 vs P4 밴드 (lat-0.5) — g1: f=1 endpoint of crease 1 vs the P4 band (lat-0.5)
#   g2: 크리스1 f=0 끝점 vs P2 밴드 (lat-0.5) — g2: f=0 endpoint of crease 1 vs the P2 band (lat-0.5)
#   g3: 크리스2 f=1 끝점 vs P4 밴드 (lat-0.5) — g3: f=1 endpoint of crease 2 vs the P4 band (lat-0.5)
import mpmath as mp, json, sys
mp.mp.dps=50

def geom(t,s,L1,L2):
    d0=mp.matrix([1,0]); d1=mp.matrix([mp.cos(t),mp.sin(t)])
    d2=mp.matrix([mp.cos(t+s),mp.sin(t+s)])
    M0=mp.matrix([0,0]); M1=M0+L1*d1; M2=M1+L2*d2; M3=M2+L1*d1
    rt=1/(2*mp.sin(t/2)); rs=1/(2*mp.sin(s/2))
    et=mp.matrix([mp.cos(t/2),mp.sin(t/2)])
    es=mp.matrix([mp.cos(t+s/2),mp.sin(t+s/2)])
    return d0,d1,d2,M0,M1,M2,M3,rt,rs,et,es

def lat(pt,P,dd):
    return -(pt[0]-P[0])*dd[1]+(pt[1]-P[1])*dd[0]

def residuals(t,s,L1,L2, br):
    d0,d1,d2,M0,M1,M2,M3,rt,rs,et,es=geom(t,s,L1,L2)
    E1p=M0+br['e1p']*rt*et    # 크리스1의 P4쪽 끝점 — P4-side endpoint of crease 1
    E1m=M0+br['e1m']*rt*et    # 크리스1의 P2쪽 끝점 — P2-side endpoint of crease 1
    E2p=M1+br['e2p']*rs*es    # 크리스2의 P4쪽 끝점 — P4-side endpoint of crease 2
    d4=mp.matrix([1,0])       # P4 방향(축소 대칭이면 heading 0) — P4 direction (heading 0 under the reduced symmetry)
    g1=br['s1']*lat(E1p,M3,d4)-mp.mpf(1)/2
    g2=br['s2']*lat(E1m,M1,d2)-mp.mpf(1)/2
    g3=br['s3']*lat(E2p,M3,d4)-mp.mpf(1)/2
    return mp.matrix([g1,g2,g3])

def detect_branches(th):
    """수치 타이트 상태에서 끝점 부호·lat 부호 브랜치 판독
    Read off the endpoint-sign and lat-sign branches from the numerically tight state."""
    import math
    t,s = th[0], th[1]; L1=(th[4]+th[6])/2; L2=th[5]
    d0,d1,d2,M0,M1,M2,M3,rt,rs,et,es=geom(t,s,L1,L2)
    d4=mp.matrix([1,0])
    br={}
    # 크리스1 끝점 두 개: ±rt et. P4에 가까운 쪽(lat |.|-0.5 최소)이 e1p — the two endpoints of crease 1 are ±rt et; the one closer to P4 (minimal |lat|-0.5) is e1p
    c1a=M0+rt*et; c1b=M0-rt*et
    la=abs(abs(lat(c1a,M3,d4))-mp.mpf(1)/2); lb=abs(abs(lat(c1b,M3,d4))-mp.mpf(1)/2)
    br['e1p']=1 if la<lb else -1
    E1p=M0+br['e1p']*rt*et; br['s1']=1 if lat(E1p,M3,d4)>0 else -1
    br['e1m']=-br['e1p']
    E1m=M0+br['e1m']*rt*et; br['s2']=1 if lat(E1m,M1,d2)>0 else -1
    c2a=M1+rs*es; c2b=M1-rs*es
    la=abs(abs(lat(c2a,M3,d4))-mp.mpf(1)/2); lb=abs(abs(lat(c2b,M3,d4))-mp.mpf(1)/2)
    br['e2p']=1 if la<lb else -1
    E2p=M1+br['e2p']*rs*es; br['s3']=1 if lat(E2p,M3,d4)>0 else -1
    return br

def solve_sL(t, br, guess):
    """t 고정: g1=g2=g3=0 을 (s,L1,L2)에 대해 뉴턴
    t held fixed: Newton-solve g1=g2=g3=0 for (s,L1,L2)."""
    f=lambda s,L1,L2: residuals(t,s,L1,L2,br)
    return mp.findroot(f, guess, tol=mp.mpf(10)**(-40))

if __name__=="__main__":
    st=json.load(open("verify/n4_tight41.json"))
    th=st["theta"]
    br=detect_branches(th)
    print("브랜치:",br)
    t0=mp.mpf((th[0]-th[3])/2)   # 대칭화 초기값 — symmetrized initial value
    s0=mp.mpf((th[1]-th[2])/2)
    L10=mp.mpf((th[4]+th[6])/2); L20=mp.mpf(th[5])
    # 검산: 초기 잔차 — sanity check: initial residuals
    print("초기 잔차:",[mp.nstr(v,3) for v in residuals(t0,s0,L10,L20,br)])
    ell=lambda t: None
    def ell_of_t(t):
        s,L1,L2=solve_sL(t,br,(s0,L10,L20))
        return 2*L1+L2,(s,L1,L2)
    # 황금분할 1차원 최소화 — golden-section 1-D minimization
    gr=(mp.sqrt(5)-1)/2
    a,b=t0-mp.mpf('0.02'), t0+mp.mpf('0.02')
    c=b-gr*(b-a); d=a+gr*(b-a)
    fc,_=ell_of_t(c); fd,_=ell_of_t(d)
    for it in range(200):
        if fc<fd: b,d,fd=d,c,fc; c=b-gr*(b-a); fc,_=ell_of_t(c)
        else: a,c,fc=c,d,fd; d=a+gr*(b-a); fd,_=ell_of_t(d)
        if b-a<mp.mpf(10)**(-42): break
    tstar=(a+b)/2
    estar,(sstar,L1s,L2s)=ell_of_t(tstar)
    print(f"\n[고정밀 해] (dps50)")
    print("t*  =",mp.nstr(tstar,40))
    print("s*  =",mp.nstr(sstar,40))
    print("L1* =",mp.nstr(L1s,40))
    print("L2* =",mp.nstr(L2s,40))
    print("ell*=",mp.nstr(estar,40))
    print("\n[가설 검증]")
    print("cos t* +1/4      =",mp.nstr(mp.cos(tstar)+mp.mpf(1)/4,10))
    print("s* -(2pi-2t*)    =",mp.nstr(sstar-(2*mp.pi-2*tstar),10))
    print("M3.y - 1         =",mp.nstr(2*L1s*mp.sin(tstar)+L2s*mp.sin(tstar+sstar)-1,10))
    print("L2*-L1* - 2/√15  =",mp.nstr(L2s-L1s-2/mp.sqrt(15),10))
    print("\n[PSLQ 동정]")
    for name,val in [("cos t*",mp.cos(tstar)),("L1*",L1s),("L2*",L2s),("ell*",estar),
                     ("sin t*",mp.sin(tstar)),("M3.y",2*L1s*mp.sin(tstar)+L2s*mp.sin(tstar+sstar))]:
        # 대수 관계: 최소다항식 시도 (차수 <=6) — algebraic relation: try for a minimal polynomial (degree <=6)
        for deg in (1,2,3,4,6,8):
            rel=mp.pslq([val**k for k in range(deg+1)], maxcoeff=10**12, maxsteps=10**6)
            if rel:
                print(f"  {name}: 차수{deg} 계수 {rel}")
                break
        else:
            print(f"  {name}: 미발견 (<=8)")
    json.dump({"t":str(tstar),"s":str(sstar),"L1":str(L1s),"L2":str(L2s),"ell":str(estar)},
              open("verify/n4_exact_theta.json","w"))
