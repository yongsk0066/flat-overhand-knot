# n4_tighten2.py — B_{4_1} 조임 본편 v2: 매듭형 오라클(Fox5) + 층 순서 탐색 + 층 고정 하강 — main B_{4_1} tightening v2: knot-type oracle (Fox5) + layer-order search + layer-fixed descent
# 교훈(v1): 그림자 canon 만으로는 매듭형이 보존되지 않는다 (Fox5=5 로 판명, 언노트 누출). — lesson from v1: the shadow canon alone does not preserve the knot type (Fox5=5 was found, i.e. unknots leaked in)
import math, random, sys, json, itertools

TARGET_CANON="ABCADCBD"
# RAY: 반직선 절단 길이. 수용 상태에는 영향 없음 [referee-H F6 증명 요지]: — RAY: the ray truncation length; it has no effect on accepted states [referee-H F6 proof sketch]:
# 코어 반경 <= ΣL <= 12 ⟹ 반직선-코어 교차 파라미터 <= 24 < 60; 장거리 교차 후보는 — core radius <= ΣL <= 12 ⟹ the ray-core intersection parameter is <= 24 < 60; the only long-range crossing candidate is
# (0,4) 반직선쌍뿐인데 유효 세포에서 M3.y ≈ 1 > 0 ⟹ P4가 y=0에 닿는 x가 양수라 — the (0,4) ray pair, but in a valid cell M3.y ≈ 1 > 0 ⟹ the x at which P4 meets y=0 is positive, so
# P0 ⊂ {x<=0}과 만나지 않는다. (검증: RAY 60→1e7 변경 시 수용 상태 변화 0/총계) — it never meets P0 ⊂ {x<=0}. (check: changing RAY from 60 to 1e7 changes 0 of all accepted states)
RAY=60.0

def canon(code):
    best=None
    for r in range(len(code)):
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
    U=[]; CH=[]
    for i in range(4):
        ux,uy=d[i][0]+d[i+1][0], d[i][1]+d[i+1][1]
        nn=math.hypot(ux,uy)
        if nn<1e-9: return None
        u=(ux/nn,uy/nn); U.append(u)
        n_=(-d[i][1],d[i][0]); s=abs(u[0]*n_[0]+u[1]*n_[1])
        if s<1e-6: return None
        t=0.5/s
        CH.append(((M[i][0]-t*u[0],M[i][1]-t*u[1]),(M[i][0]+t*u[0],M[i][1]+t*u[1])))
    return dict(d=d,M=M,U=U,CH=CH,L=list(L))

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

def events(X):
    ev=[]
    for cid,((i,j),s,u) in enumerate(X):
        ev.append((i,-s if i==0 else s,cid,'i'))
        ev.append((j,-u if j==0 else u,cid,'j'))
    ev.sort()
    return ev

def code_str(X):
    if len(X)!=4: return None
    lab={cid:chr(ord('A')+cid) for cid in range(4)}
    return "".join(lab[cid] for _,_,cid,_ in events(X))

def fox5_events(ev_over):
    """ev_over: [(cid, is_over)] 재료 순서. 폐포 후 Fox 5-채색 수.
    ev_over: [(cid, is_over)] in material order. Number of Fox 5-colorings after closing up."""
    A=sum(1 for _,ov in ev_over if not ov)
    if A==0: return 5
    arc=0; visits=[]
    for cid,ov in ev_over:
        visits.append((cid,ov,arc))
        if not ov: arc+=1
    rows=[]; over_arc={}; und={}
    for k,(cid,ov,a0) in enumerate(visits):
        if ov: over_arc[cid]=a0%A
        else: und.setdefault(cid,[]).append((a0%A,(a0+1)%A))
    for cid,oa in over_arc.items():
        (ain,aout)=und[cid][0]
        r=[0]*A; r[oa]+=2; r[ain]-=1; r[aout]-=1
        rows.append([x%5 for x in r])
    m=[r[:] for r in rows]; rank=0; col=0; R=len(m)
    while rank<R and col<A:
        piv=next((rr for rr in range(rank,R) if m[rr][col]%5),None)
        if piv is None: col+=1; continue
        m[rank],m[piv]=m[piv],m[rank]
        inv=pow(m[rank][col],3,5)
        m[rank]=[(v*inv)%5 for v in m[rank]]
        for rr in range(R):
            if rr!=rank and m[rr][col]%5:
                f=m[rr][col]; m[rr]=[(m[rr][k]-f*m[rank][k])%5 for k in range(A)]
        rank+=1; col+=1
    return 5**(A-rank)

def fox3_events(ev_over):
    A=sum(1 for _,ov in ev_over if not ov)
    if A==0: return 3
    arc=0; visits=[]
    for cid,ov in ev_over:
        visits.append((cid,ov,arc))
        if not ov: arc+=1
    rows=[]; over_arc={}; und={}
    for cid,ov,a0 in visits:
        if ov: over_arc[cid]=a0%A
        else: und.setdefault(cid,[]).append((a0%A,(a0+1)%A))
    for cid,oa in over_arc.items():
        ain,aout=und[cid][0]
        r=[0]*A; r[oa]+=2; r[ain]-=1; r[aout]-=1
        rows.append([x%3 for x in r])
    m=[r[:] for r in rows]; rank=0; col=0; R=len(m)
    while rank<R and col<A:
        piv=next((rr for rr in range(rank,R) if m[rr][col]%3),None)
        if piv is None: col+=1; continue
        m[rank],m[piv]=m[piv],m[rank]
        inv=m[rank][col]%3  # 1->1, 2->2 (2*2=4=1)
        m[rank]=[(v*inv)%3 for v in m[rank]]
        for rr in range(R):
            if rr!=rank and m[rr][col]%3:
                f=m[rr][col]; m[rr]=[(m[rr][k]-f*m[rank][k])%3 for k in range(A)]
        rank+=1; col+=1
    return 3**(A-rank)

def overs_from_order(X, order):
    """order: 패스 리스트 위->아래. 각 교차 cid: over 방문은 상위 패스.
    order: list of passes, top -> bottom. For each crossing cid, the over visit is the higher pass."""
    rank={p:i for i,p in enumerate(order)}
    ov={}
    for cid,((i,j),_,_) in enumerate(X):
        ov[cid]='i' if rank[i]<rank[j] else 'j'
    return ov

def ev_over_of(X, order):
    ov=overs_from_order(X,order)
    return [(cid, side==ov[cid]) for _,_,cid,side in events(X)]

def linear_extensions(pairs, nodes=(0,1,2,3,4)):
    """pairs: (top,bottom) 제약. 모든 선형 확장.
    pairs: (top,bottom) constraints. Returns all linear extensions."""
    from itertools import permutations
    out=[]
    for perm in permutations(nodes):
        pos={p:i for i,p in enumerate(perm)}
        if all(pos[a]<pos[b] for a,b in pairs): out.append(list(perm))
    return out

def orders_making_41(X):
    """4교차 실현에서: 4_1 (Fox5=25)을 만드는 over/under 배정(비순환) 전부와 그 선형 확장들.
    For a 4-crossing realization: all (acyclic) over/under assignments yielding 4_1 (Fox5=25), with their linear extensions."""
    res=[]
    for mask in range(16):
        pairs=[]
        for cid,((i,j),_,_) in enumerate(X):
            top,bot=(i,j) if (mask>>cid)&1 else (j,i)
            pairs.append((top,bot))
        exts=linear_extensions(pairs)
        if not exts: continue
        evo=[(cid, (side=='i')==((mask>>cid)&1==1) ) for _,_,cid,side in events(X)]
        if fox5_events(evo)==25:
            res.append((mask,exts))
    return res

def covers(c,k,pt,tol=0.0):
    P,dd,l=pieces(c)[k]
    rx,ry=pt[0]-P[0],pt[1]-P[1]
    lon=rx*dd[0]+ry*dd[1]; lat=abs(-rx*dd[1]+ry*dd[0])
    return (-tol<lon<l+tol) and lat<0.5+tol

def axiom_ok(c, order, nsamp=31):
    rank={p:i for i,p in enumerate(order)}
    for i in range(4):
        ra,rb=rank[i],rank[i+1]
        lo,hi=min(ra,rb),max(ra,rb)
        between=[order[r] for r in range(lo+1,hi)]
        if not between: continue
        (x1,y1),(x2,y2)=c['CH'][i]
        for t in range(1,nsamp+1):
            f=t/(nsamp+1); pt=(x1+f*(x2-x1),y1+f*(y2-y1))
            for k in between:
                if covers(c,k,pt,tol=-1e-9): return False
    return True

def valid_any(th):
    c=config(th)
    if c is None: return None
    X=crossings(c)
    cd=code_str(X)
    if cd is None or canon(cd)!=TARGET_CANON: return None
    good=[]
    for mask,exts in orders_making_41(X):
        for order in exts:
            if axiom_ok(c,order): good.append(order)
    if not good: return None
    return dict(c=c,X=X,orders=good)

def valid_fixed(th, order):
    c=config(th)
    if c is None: return None
    X=crossings(c)
    cd=code_str(X)
    if cd is None or canon(cd)!=TARGET_CANON: return None
    if fox5_events(ev_over_of(X,order))!=25: return None
    if not axiom_ok(c,order): return None
    return dict(c=c,X=X)

def ell(th): return th[4]+th[5]+th[6]

# ---------------- 단위검증 ---------------- — unit verification
def _test():
    # 표준 교대 코드: 세잎 O1U2O3U1O2U3, 8자 O1U2O3U4O2U1O4U3, 언노트 O1U1 — standard alternating codes: trefoil O1U2O3U1O2U3, figure-8 O1U2O3U4O2U1O4U3, unknot O1U1
    tre=[(0,True),(1,False),(2,True),(0,False),(1,True),(2,False)]
    f8 =[(0,True),(1,False),(2,True),(3,False),(1,True),(0,False),(3,True),(2,False)]
    unk=[(0,True),(0,False)]
    print("세잎: fox3 =",fox3_events(tre),"(기대 9*? => 3^(1+null)=9? det3 => 27? )","fox5 =",fox5_events(tre))
    print("8자 : fox3 =",fox3_events(f8), "fox5 =",fox5_events(f8),"(기대 25)")
    print("언노트: fox3 =",fox3_events(unk),"fox5 =",fox5_events(unk),"(기대 3, 5)")

if __name__=="__main__":
    if "--test" in sys.argv: _test(); sys.exit()

# ---------------- 수확 (진단 포함) ---------------- — harvest (with diagnostics)
def _worker(seed_n):
    seed,n=seed_n
    rnd=random.Random(seed)
    diag={"canon":0}
    pairstat={}   # pairing-tuple -> [seen, fox_acyclic, axiom_valid]
    hits=[]
    for _ in range(n):
        if rnd.random()<0.5:
            rho=rnd.choice([1,-1]); sg=[rho]*4
        else:
            sg=[rnd.choice([1,-1]) for _ in range(4)]
        th=tuple([sg[k]*rnd.uniform(0.3,3.1) for k in range(4)]+
                 [rnd.uniform(0.3,4.0) for _ in range(3)])
        c=config(th)
        if c is None: continue
        X=crossings(c)
        cd=code_str(X)
        if cd is None or canon(cd)!=TARGET_CANON: continue
        diag["canon"]+=1
        key=tuple(sorted(p for p,_,_ in X))
        st=pairstat.setdefault(key,[0,0,0])
        st[0]+=1
        om=orders_making_41(X)
        if not om: continue
        st[1]+=1
        ok=[o for m,exts in om for o in exts if axiom_ok(c,o)]
        if ok:
            st[2]+=1
            hits.append((list(th),ok[0]))
    return diag,pairstat,hits

def harvest(total=3_000_000, procs=8):
    import multiprocessing as mp
    per=total//procs
    with mp.get_context("fork").Pool(procs) as pool:
        parts=pool.map(_worker,[(1000+i,per) for i in range(procs)])
    canon_n=sum(d["canon"] for d,_,_ in parts)
    agg={}
    hits=[]
    for _,ps,h in parts:
        hits+=h
        for k,v in ps.items():
            a=agg.setdefault(k,[0,0,0])
            for t in range(3): a[t]+=v[t]
    print(f"표본 {total:,} 중 셰도우(canon) 통과 {canon_n}")
    print(f"{'패스-쌍 조합':<40}{'셰도우':>8}{'4_1층존재':>10}{'공리통과':>9}")
    for k,v in sorted(agg.items(),key=lambda kv:-kv[1][0]):
        print(f"{str(k):<42}{v[0]:>7}{v[1]:>10}{v[2]:>9}")
    print(f"\n완전 유효(4_1 + 공리) 상태: {len(hits)}개")
    json.dump(hits[:200],open("verify/n4_valid41.json","w"))
    for th,o in hits[:6]:
        print("  ell=",round(ell(th),3)," order=",o," th=",[round(v,3) for v in th])
    return hits

if __name__=="__main__":
    if "--harvest" in sys.argv:
        harvest()

# ---------------- 층 고정 조임 ---------------- — layer-fixed tightening
def tighten_fixed(th0, order, rnd, rounds=7):
    best=list(th0); ebest=ell(best); th=list(th0)
    for rd in range(rounds):
        step=0.06/(2**rd); stall=0
        for it in range(2500):
            e=ell(th)
            adv=False
            for _ in range(30):
                dv=[rnd.gauss(0,1) for _ in range(7)]
                for k in (4,5,6): dv[k]-=abs(rnd.gauss(0,1))*0.5
                nn=math.sqrt(sum(v*v for v in dv))
                cand=[th[k]+step*dv[k]/nn for k in range(7)]
                if any(cand[k]<0.05 for k in (4,5,6)): continue
                if valid_fixed(tuple(cand),order) is None: continue
                ec=ell(cand)
                if ec<e-1e-12: th=cand; adv=True; break
                if ec<e+step*0.05 and rnd.random()<0.35: th=cand; break
            if ell(th)<ebest-1e-12: ebest=ell(th); best=list(th); stall=0
            else:
                stall+=1
                if stall>300: break
    return tuple(best),ebest

if __name__=="__main__":
    if "--tighten" in sys.argv:
        rnd=random.Random(88)
        hits=json.load(open("verify/n4_valid41.json"))
        print(f"시작점 {len(hits)}개")
        res=[]
        for th0,order in hits:
            thf,ef=tighten_fixed(th0,order,rnd)
            res.append((ef,thf,order))
            print(f"  {ell(th0):.3f} -> {ef:.6f}   order={order}")
        res.sort()
        eb,tb,ob=res[0]
        print(f"\n[조임 최소] ell* = {eb:.7f}   order={ob}")
        print(f"  tau = {[round(v,6) for v in tb[:4]]}")
        print(f"  L   = {[round(v,6) for v in tb[4:]]}")
        v=valid_fixed(tb,ob)
        print(f"  교차: {[(p,round(s,4),round(u,4)) for p,s,u in v['X']]}")
        near=[t for e,t,o in res if e<eb+0.05]
        print(f"  ell*+0.05 이내: {len(near)}개")
        json.dump({"ell":eb,"theta":list(tb),"order":ob},open("verify/n4_tight41.json","w"))

# ---------------- 정확 공리 (표본 없음: affine 구간 산술) ---------------- — exact axiom (no sampling: affine interval arithmetic)
def axiom_exact_ok(c, order, tol=1e-11):
    """크리스 i의 '열린' 현이, 층에서 사이에 끼는 패스의 열린 밴드영역과
    양의 길이로 만나면 위반. 밴드=직선 슬랩×종방향 창 => f의 선형 부등식.
    It is a violation if the 'open' chord of crease i meets, in positive length, the open band
    region of a pass sandwiched between the two in the layer order. band = straight slab ×
    longitudinal window => linear inequalities in f."""
    segs=pieces(c)
    rank={p:q for q,p in enumerate(order)}
    for i in range(4):
        ra,rb=rank[i],rank[i+1]; lo,hi=min(ra,rb),max(ra,rb)
        between=[order[r] for r in range(lo+1,hi)]
        if not between: continue
        (ax,ay),(bx,by)=c['CH'][i]
        for k in between:
            P,dd,l=segs[k]
            n=(-dd[1],dd[0])
            latA=(ax-P[0])*n[0]+(ay-P[1])*n[1]; latB=(bx-P[0])*n[0]+(by-P[1])*n[1]
            lonA=(ax-P[0])*dd[0]+(ay-P[1])*dd[1]; lonB=(bx-P[0])*dd[0]+(by-P[1])*dd[1]
            # f∈[0,1]에서 |lat|<1/2 구간 — the subinterval of f∈[0,1] where |lat|<1/2
            f0,f1=0.0,1.0
            dB=latB-latA
            def clip(lo_,hi_,c0,c1,f0,f1):
                # c0 + (c1-c0) f 가 (lo_,hi_) 안에 있는 f 구간과 교차 — intersect with the f-interval on which c0 + (c1-c0) f lies inside (lo_,hi_)
                d=c1-c0
                if abs(d)<1e-15:
                    if not (lo_<c0<hi_): return 1.0,0.0
                    return f0,f1
                fa=(lo_-c0)/d; fb=(hi_-c0)/d
                if fa>fb: fa,fb=fb,fa
                return max(f0,fa),min(f1,fb)
            DEL=1e-12   # referee-H F4: 정확 flush 접촉(lat≡±1/2)의 부동소수점 knife-edge 수용 — accept the floating-point knife-edge of an exact flush contact (lat≡±1/2)
            g0,g1=clip(-0.5+DEL,0.5-DEL,latA,latB,0.0,1.0)
            if g0>=g1: continue
            g0,g1=clip(0.0+DEL,l-DEL,lonA,lonB,g0,g1)
            if g1-g0>tol: return False
    return True

def valid_fixed2(th, order):
    c=config(th)
    if c is None: return None
    X=crossings(c)
    cd=code_str(X)
    if cd is None or canon(cd)!=TARGET_CANON: return None
    if fox5_events(ev_over_of(X,order))!=25: return None
    if not axiom_exact_ok(c,order): return None
    return dict(c=c,X=X)

def tighten_fixed2(th0, order, rnd, rounds=8):
    best=list(th0); ebest=ell(best); th=list(th0)
    for rd in range(rounds):
        step=0.06/(2**rd); stall=0
        for it in range(2500):
            e=ell(th); adv=False
            for _ in range(30):
                dv=[rnd.gauss(0,1) for _ in range(7)]
                for k in (4,5,6): dv[k]-=abs(rnd.gauss(0,1))*0.5
                nn=math.sqrt(sum(v*v for v in dv))
                cand=[th[k]+step*dv[k]/nn for k in range(7)]
                if any(cand[k]<0.05 for k in (4,5,6)): continue
                if valid_fixed2(tuple(cand),order) is None: continue
                ec=ell(cand)
                if ec<e-1e-12: th=cand; adv=True; break
                if ec<e+step*0.05 and rnd.random()<0.35: th=cand; break
            if ell(th)<ebest-1e-12: ebest=ell(th); best=list(th); stall=0
            else:
                stall+=1
                if stall>300: break
    return tuple(best),ebest

if __name__=="__main__":
    if "--tighten2" in sys.argv:
        rnd=random.Random(77)
        hits=json.load(open("verify/n4_valid41.json"))
        ok=[(th,o) for th,o in hits if valid_fixed2(tuple(th),o)]
        print(f"정확 공리 재검: {len(ok)}/{len(hits)} 시작점 생존")
        res=[]
        for th0,order in ok:
            thf,ef=tighten_fixed2(th0,order,rnd)
            res.append((ef,thf,order))
            print(f"  {ell(th0):.3f} -> {ef:.6f}")
        res.sort()
        eb,tb,ob=res[0]
        print(f"\n[정확공리 조임 최소] ell* = {eb:.7f}")
        print(f"  tau = {[round(v,6) for v in tb[:4]]}")
        print(f"  L   = {[round(v,6) for v in tb[4:]]}")
        near=[t for e,t,o in res if e<eb+0.05]
        print(f"  ell*+0.05 이내: {len(near)}개")
        json.dump({"ell":eb,"theta":list(tb),"order":ob},open("verify/n4_tight41.json","w"))
