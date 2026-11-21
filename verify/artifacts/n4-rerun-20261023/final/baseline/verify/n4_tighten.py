# n4_tighten.py — 프로그램 2단계 본편: B_{4_1}(8자 매듭)의 폭-1 유효성 + 조임 [수치 1차] — main stage 2 of the program: width-1 validity + tightening for B_{4_1} (figure-8 knot) [first numerical pass]
# 차트: theta = (tau1..tau4, L1..L3), dim 7 = 2n-1.  ell_core = L1+L2+L3 (선형!). — chart: theta = (tau1..tau4, L1..L3), dim 7 = 2n-1; ell_core = L1+L2+L3 (linear!)
# 층 순서: P1 > P3 > P4 > P0 > P2  (n4_witness). 코드: canon = ABCADCBD. — layer order: P1 > P3 > P4 > P0 > P2 (n4_witness); code: canon = ABCADCBD
import math, random, sys

LAYER = [1,3,4,0,2]            # 위 -> 아래 — top -> bottom
RANK = {p:i for i,p in enumerate(LAYER)}
TARGET_CANON = "ABCADCBD"

def canon(code):
    n=len(code); best=None
    for r in range(n):
        rot=code[r:]+code[:r]; lab={}; out=[]
        for ch in rot:
            if ch not in lab: lab[ch]=chr(ord('A')+len(lab))
            out.append(lab[ch])
        cand="".join(out)
        best=cand if best is None or cand<best else best
    return best

def config(th):
    tau=th[:4]; L=th[4:]
    cs=[sum(tau[:j]) for j in range(5)]
    d=[(math.cos(c),math.sin(c)) for c in cs]
    M=[(0.0,0.0)]
    for i in range(3): M.append((M[-1][0]+L[i]*d[i+1][0], M[-1][1]+L[i]*d[i+1][1]))
    # 크리스 방향 u_i ∥ d_i + d_{i+1} (반사 항등 — 세잎 슬라이스에서 검증) — crease direction u_i ∥ d_i + d_{i+1} (reflection identity, verified on the trefoil slice)
    U=[]
    for i in range(4):
        ux,uy=d[i][0]+d[i+1][0], d[i][1]+d[i+1][1]
        nn=math.hypot(ux,uy)
        if nn<1e-9: return None
        U.append((ux/nn,uy/nn))
    # 크리스 현: line(M_i, U_i)를 폭-1 밴드(중심선±1/2)로 절단. — crease chord: truncate line(M_i, U_i) to the width-1 band (centerline±1/2)
    # 현 반길이 t: |t (U·n̂)| = 1/2, n̂ = perp(d_i) — chord half-length t: |t (U·n̂)| = 1/2, n̂ = perp(d_i)
    CH=[]
    for i in range(4):
        n_=(-d[i][1],d[i][0])
        s=abs(U[i][0]*n_[0]+U[i][1]*n_[1])
        if s<1e-6: return None
        t=0.5/s
        CH.append(((M[i][0]-t*U[i][0],M[i][1]-t*U[i][1]),(M[i][0]+t*U[i][0],M[i][1]+t*U[i][1])))
    return dict(d=d,M=M,U=U,CH=CH,L=list(L))

RAY=60.0
def pieces(c):
    d,M=c['d'],c['M']
    return {0:(M[0],(-d[0][0],-d[0][1]),RAY),1:(M[0],d[1],c['L'][0]),
            2:(M[1],d[2],c['L'][1]),3:(M[2],d[3],c['L'][2]),4:(M[3],d[4],RAY)}

def crossings(c):
    segs=pieces(c); X=[]
    for i in range(5):
        for j in range(i+2,5):
            P1,d1,l1=segs[i]; P2,d2,l2=segs[j]
            det=d1[0]*d2[1]-d1[1]*d2[0]
            if abs(det)<1e-12: continue
            rx,ry=P2[0]-P1[0],P2[1]-P1[1]
            s=(rx*d2[1]-ry*d2[0])/det; u=(rx*d1[1]-ry*d1[0])/det
            if 1e-9<s<l1-1e-9 and 1e-9<u<l2-1e-9: X.append(((i,j),s,u))
    return X

def code_of(X):
    if len(X)!=4: return None
    lab={p:chr(ord('A')+k) for k,(p,_,_) in enumerate(X)}
    ev=[]
    for p,s,u in X:
        i,j=p
        ev.append((i,-s if i==0 else s,lab[p]))
        ev.append((j,-u if j==0 else u,lab[p]))
    ev.sort()
    return "".join(t[2] for t in ev)

def covers(c,k,pt,tol=0.0):
    """패스 k(폭-1 밴드, 종방향 범위 포함)가 점 pt를 (tol 여유로) 덮는가
    Does pass k (the width-1 band, longitudinal range included) cover the point pt, within margin tol?"""
    P,dd,l=pieces(c)[k]
    rx,ry=pt[0]-P[0],pt[1]-P[1]
    lon=rx*dd[0]+ry*dd[1]; lat=abs(-rx*dd[1]+ry*dd[0])
    return (-tol<lon<l+tol) and lat<0.5+tol

def axiom_margin(c,nsamp=25):
    """크리스 이웃 공리: 크리스 i(P_{i-1}↔P_i 연결)의 내부점을, 층에서 둘 사이에 끼는
    패스가 덮으면 위반. margin = 위반 시 음수(깊이), 아니면 최소 이탈 거리의 근사(>0 표기 생략, 0.0).
    Crease-neighbour axiom: it is a violation if an interior point of crease i (which joins
    P_{i-1}↔P_i) is covered by a pass sandwiched between the two in the layer order.
    margin = negative (the depth) on violation, otherwise an approximation of the minimal
    clearance distance (the >0 value is not reported, 0.0 is used)."""
    worst=1.0
    for i in range(4):
        a,bnd=i,i+1
        ra,rb=RANK[a],RANK[bnd]
        lo,hi=min(ra,rb),max(ra,rb)
        between=[LAYER[r] for r in range(lo+1,hi)]
        if not between: continue
        (x1,y1),(x2,y2)=c['CH'][i]
        for t in range(1,nsamp+1):
            f=t/(nsamp+1)
            pt=(x1+f*(x2-x1), y1+f*(y2-y1))
            for k in between:
                if covers(c,k,pt,tol=-1e-9):
                    # 위반 깊이: 밴드 경계까지 거리(근사) — violation depth: distance to the band boundary (approximate)
                    P,dd,l=pieces(c)[k]
                    rx,ry=pt[0]-P[0],pt[1]-P[1]
                    lat=abs(-rx*dd[1]+ry*dd[0]); lon=rx*dd[0]+ry*dd[1]
                    depth=min(0.5-lat, lon if lon<l/2 else l-lon)
                    worst=min(worst,-max(depth,1e-6))
    return worst

def valid(th):
    c=config(th)
    if c is None: return None,None,"config"
    X=crossings(c)
    cd=code_of(X)
    if cd is None or canon(cd)!=TARGET_CANON: return None,None,"code"
    m=axiom_margin(c)
    if m<0: return None,None,f"axiom{m:.3f}"
    return c,X,""

def ell(th): return th[4]+th[5]+th[6]

# ---- 유효 시작점 수확 ---- — harvest valid starting points
rnd=random.Random(4141)
starts=[]
tries=0
while len(starts)<40 and tries<600000:
    tries+=1
    rho=rnd.choice([1,-1])
    th=tuple([rho*rnd.uniform(0.6,3.05) for _ in range(3)]+[rho*rnd.uniform(0.3,3.05)]+
             [rnd.uniform(0.4,3.5) for _ in range(3)])
    c,X,why=valid(th)
    if c: starts.append(th)
print(f"[유효 시작점] {len(starts)}개 (시도 {tries})  — 폭-1 공리까지 통과한 8자 평면 접힘 상태의 실재")
if starts:
    e0=[round(ell(t),3) for t in starts[:8]]
    print(f"  초기 ell 예: {e0}")

# ---- 조임: 좌표 하강 + 무작위 방향, 유효성 거부 ---- — tightening: coordinate descent + random directions, with rejection on invalidity
def tighten(th0, iters=4000):
    th=list(th0); e=ell(th); step=0.08
    for it in range(iters):
        moved=False
        for _ in range(14):
            dvec=[rnd.gauss(0,1) for _ in range(7)]
            # ell 감소 편향: L 성분에 음의 성분 강제 혼합 — bias toward decreasing ell: force a negative component into the L coordinates
            for k in (4,5,6): dvec[k]-=abs(rnd.gauss(0,1))*0.7
            nn=math.sqrt(sum(v*v for v in dvec))
            cand=[th[k]+step*dvec[k]/nn for k in range(7)]
            if any(cand[k]<0.05 for k in (4,5,6)): continue
            c,X,why=valid(tuple(cand))
            if c and ell(cand)<e-1e-12:
                th,e=cand,ell(cand); moved=True; break
        if not moved:
            step*=0.5
            if step<1e-7: break
    return tuple(th),e

results=[]
for th0 in starts[:25]:
    thf,ef=tighten(th0)
    results.append((ef,thf))
results.sort()
print(f"\n[조임 결과] {len(results)}개 하강 종점, ell 분포:")
import collections
cl=collections.Counter(round(e,3) for e,_ in results)
print("  ", dict(sorted(cl.items())[:8]))
if results:
    ebest,thbest=results[0]
    print(f"\n[최소] ell* = {ebest:.6f}")
    print(f"  tau = {[round(v,5) for v in thbest[:4]]}")
    print(f"  L   = {[round(v,5) for v in thbest[4:]]}")
    c,X,_=valid(thbest)
    print(f"  교차: {[(p,round(s,3),round(u,3)) for p,s,u in X]}")
    print(f"  공리 margin: {axiom_margin(c):.4f}")
    import json
    json.dump({"ell":ebest,"theta":list(thbest)}, open("verify/n4_tight.json","w"))

# ---- 2단계: 경계 활강 폴리시 (정체 탈출: ell 비증가 측면 이동 허용) ---- — stage 2: boundary-gliding polish (escape stalls: allow lateral moves that do not increase ell)
def polish(th0, rounds=6):
    best=list(th0); ebest=ell(best)
    th=list(th0)
    for rd in range(rounds):
        step=0.05/(2**rd)
        stall=0
        for it in range(1500):
            e=ell(th); improved=False
            for _ in range(30):
                dvec=[rnd.gauss(0,1) for _ in range(7)]
                for k in (4,5,6): dvec[k]-=abs(rnd.gauss(0,1))*0.5
                nn=math.sqrt(sum(v*v for v in dvec))
                cand=[th[k]+step*dvec[k]/nn for k in range(7)]
                if any(cand[k]<0.05 for k in (4,5,6)): continue
                c,X,why=valid(tuple(cand))
                if not c: continue
                ec=ell(cand)
                if ec<e-1e-12:
                    th=cand; improved=True; break
                if ec<e+step*0.05 and rnd.random()<0.35:   # 측면 활강 — lateral glide
                    th=cand; break
            if ell(th)<ebest-1e-12: ebest=ell(th); best=list(th); stall=0
            else:
                stall+=1
                if stall>250: break
    return tuple(best), ebest

if "--polish" in sys.argv:
    import json
    seed_pts=[th for _,th in results[:10]]+starts[:6]
    fin=[]
    for th0 in seed_pts:
        thf,ef=polish(th0)
        fin.append((ef,thf))
        print(f"  polish: {ell(th0):.4f} -> {ef:.6f}")
    fin.sort()
    eb,tb=fin[0]
    print(f"\n[폴리시 최소] ell* = {eb:.7f}")
    print(f"  tau = {[round(v,6) for v in tb[:4]]}")
    print(f"  L   = {[round(v,6) for v in tb[4:]]}")
    c,X,_=valid(tb)
    print(f"  교차 파라미터: {[(p,round(s,4),round(u,4)) for p,s,u in X]}")
    # 근접 종점 클러스터 — cluster of nearby endpoints
    close=[t for e,t in fin if e<eb+0.02]
    print(f"  ell*+0.02 이내 종점: {len(close)}개")
    json.dump({"ell":eb,"theta":list(tb)}, open("verify/n4_tight.json","w"))
