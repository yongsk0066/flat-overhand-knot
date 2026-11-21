# n4_strata2.py — n=4 층화 KKT 수확 (uniq_strata의 유사물) — n=4 stratified KKT harvest (analogue of uniq_strata)
# 인벤토리 벽 목록 위에서 |S|<=7 전 층: {g_S=0} ∪ {∇ell = Σλ∇g} 정방계 뉴턴 다중시드 — over the inventory wall list, all strata with |S|<=7: multi-seeded Newton on the square system {g_S=0} ∪ {∇ell = Σλ∇g}
import sys, math, random, json, itertools, collections
sys.path.insert(0,'verify')
from n4_tighten2 import config, pieces, crossings, code_str, canon, ev_over_of, fox5_events, valid_fixed2
from n4_walls import all_margins, cap_data, ORDER
import multiprocessing as mp

def margin_fn(label, sign_ref):
    """라벨별 매끄러운 대표 g(th). sign_ref: 부호 동결용 참조점
    Smooth representative g(th) per label. sign_ref: reference point used to freeze the sign."""
    kind=label[0]
    def lat_lon(th, i, k, e):
        c=config(th); segs=pieces(c)
        A,B=c['CH'][i]; pt=A if e=='-' else B
        P,dd,l=segs[k]
        lat=-(pt[0]-P[0])*dd[1]+(pt[1]-P[1])*dd[0]
        lon=(pt[0]-P[0])*dd[0]+(pt[1]-P[1])*dd[1]
        return lat,lon,l
    if kind=='a':
        i=int(label[3])-1; k=int(label[5]); e=label[7]
        lat0,_,_=lat_lon(sign_ref,i,k,e); s=1.0 if lat0>=0 else -1.0
        return lambda th: s*lat_lon(th,i,k,e)[0]-0.5
    if kind=='b':
        i=int(label[3])-1; k=int(label[5]); e=label[7]; hi=label.endswith('hi')
        def g(th):
            c=config(th);pt=c['CH'][i][0 if e=='-' else 1]
            form=next(form for name,form,_ in cap_data(c,k)[1] if name==('hi' if hi else 'lo'))
            return -(form[0]*pt[0]+form[1]*pt[1]+form[2])
        return g
    if kind=='d':
        i=int(label[3])-1; k=int(label[5])
        cn='hi' if 'Khi' in label else 'lo'; sn=1.0 if label.endswith('p') else -1.0
        def g(th):
            c=config(th); segs=pieces(c)
            A,B=c['CH'][i]; P,dd,l=segs[k]
            corners=next(corners for name,_,corners in cap_data(c,k)[1] if name==cn)
            cx,cy=max(corners,key=lambda pt:sn*(-(pt[0]-P[0])*dd[1]+(pt[1]-P[1])*dd[0]))
            ab=(B[0]-A[0],B[1]-A[1]); nab=math.hypot(*ab)
            return ((cx-A[0])*ab[1]-(cy-A[1])*ab[0])/nab
        s0=1.0 if g(sign_ref)>=0 else -1.0
        return lambda th: s0*g(th)
    if kind=='x':
        pr=label[3:label.index(')')]; ii,jj=map(int,pr.split(','))
        side=label[label.index(')')+1:]
        def g(th):
            c=config(th); X=crossings(c); segs=pieces(c)
            for (a,b),s,u in X:
                if (a,b)==(ii,jj):
                    _,_,li=segs[a]; _,_,lj=segs[b]
                    return {'s0':s,'u0':u,'s1':li-s,'u1':lj-u}[side]
            return 10.0
        return g
    if kind=='L':
        k=int(label[1])+3
        return lambda th: th[k]
    raise ValueError(label)

def grad(f,th,h=2e-7):
    return [ (f(tuple(th[m]+h if m==q else th[m] for m in range(7)))
             -f(tuple(th[m]-h if m==q else th[m] for m in range(7))))/(2*h) for q in range(7)]

GRAD_ELL=[0,0,0,0,1,1,1]

def solve_stratum(labels, seed_th, ref):
    fns=[margin_fn(lb,ref) for lb in labels]
    m=len(labels)
    th=list(seed_th); lam=[0.5]*m
    for it in range(60):
        try:
            G=[grad(f,th) for f in fns]
            F=[f(tuple(th)) for f in fns]+[GRAD_ELL[q]-sum(lam[j]*G[j][q] for j in range(m)) for q in range(7)]
        except Exception:
            return None
        rn=math.sqrt(sum(v*v for v in F))
        if rn<5e-8: break
        # 야코비안 (수치): 변수 = (th7, lam_m) — Jacobian (numerical): variables = (th7, lam_m)
        J=[]
        h=1e-5
        for q in range(7):
            thp=list(th); thm=list(th); thp[q]+=h; thm[q]-=h
            try:
                Gp=[grad(f,thp) for f in fns]; Gm=[grad(f,thm) for f in fns]
                Fp=[f(tuple(thp)) for f in fns]+[GRAD_ELL[r]-sum(lam[j]*Gp[j][r] for j in range(m)) for r in range(7)]
                Fm=[f(tuple(thm)) for f in fns]+[GRAD_ELL[r]-sum(lam[j]*Gm[j][r] for j in range(m)) for r in range(7)]
            except Exception:
                return None
            J.append([(Fp[r]-Fm[r])/(2*h) for r in range(m+7)])
        for j in range(m):
            col=[0.0]*(m+7)
            G=[grad(f,th) for f in fns]
            for r in range(7): col[m+r]=-G[j][r]
            J.append(col)
        # (m+7)x(m+7) 최소제곱 뉴턴 스텝 (가우스 소거) — (m+7)x(m+7) least-squares Newton step (Gaussian elimination)
        n=m+7
        A=[[J[cc][rr] for cc in range(n)] for rr in range(n)]
        b=[-F[rr] for rr in range(n)]
        # 부분 피벗 소거 — elimination with partial pivoting
        for col in range(n):
            piv=max(range(col,n),key=lambda r:abs(A[r][col]))
            if abs(A[piv][col])<1e-14: return None
            A[col],A[piv]=A[piv],A[col]; b[col],b[piv]=b[piv],b[col]
            for r in range(col+1,n):
                f2=A[r][col]/A[col][col]
                for cc in range(col,n): A[r][cc]-=f2*A[col][cc]
                b[r]-=f2*b[col]
        dx=[0.0]*n
        for r in range(n-1,-1,-1):
            ssum=b[r]-sum(A[r][cc]*dx[cc] for cc in range(r+1,n))
            dx[r]=ssum/A[r][r]
        stepn=math.sqrt(sum(v*v for v in dx))
        damp=min(1.0, 0.3/max(stepn,1e-12)) if stepn>0.3 else 1.0
        for q in range(7): th[q]+=damp*dx[q]
        for j in range(m): lam[j]+=damp*dx[m+j]
        if any(abs(v)>25 for v in th): return None
    else:
        return None
    if rn>5e-8: return None
    return tuple(th), lam, rn

def _work(args):
    labels_list, seeds, ref = args
    out=[]
    for labels in labels_list:
        for sd in seeds:
            r=solve_stratum(labels, sd, ref)
            if r is None: continue
            th,lam,rn=r
            out.append((labels, list(th), list(lam), rn))
    return out

def actual_closure_ok(th, tolerance=1e-6):
    """Numerical approximate-closure test, independent of wall inventory.

    Newton's residual tolerance is 5e-8; the geometric rejection tolerance
    is explicitly wider. This is not an exact feasibility certificate.
    """
    if any(not math.isfinite(v) for v in th) or any(v<=0 for v in th[4:]):
        return False
    return valid_fixed2(tuple(th),ORDER,distance_tolerance=tolerance,
                        parameter_tolerance=1e-8) is not None

if __name__=="__main__":
    inv=json.load(open("verify/n4_wall_inventory.json"))
    WALLS=[k for k,v in sorted(inv.items(),key=lambda kv:-kv[1]) if v>=2 and k!="(미분류)"]
    print(f"[인벤토리 벽] {len(WALLS)}종: {WALLS}")
    tstar=math.acos(-0.25)
    THSTAR=(tstar,2*math.pi-2*tstar,-(2*math.pi-2*tstar),-tstar,
            2*math.sqrt(15)/5,8*math.sqrt(15)/15,2*math.sqrt(15)/5)
    hits=[tuple(th) for th,o in json.load(open("verify/n4_valid41.json")) if th[0]>0][:6]
    seeds=[THSTAR]+hits
    print(f"[시드] {len(seeds)}개")
    ref=THSTAR
    strata=[]
    for r in range(1,8):
        strata+= list(itertools.combinations(WALLS,r))
    print(f"[층 수] {len(strata)} (|S|<=7)")
    CH=max(1,len(strata)//8)
    chunks=[strata[i:i+CH] for i in range(0,len(strata),CH)]
    with mp.get_context("fork").Pool(8) as pool:
        parts=pool.map(_work,[(ch,seeds,ref) for ch in chunks])
    sols=[s for p in parts for s in p]
    print(f"[수확] 수렴 해 {len(sols)}")
    # 필터: 유효 폐포 + λ>=0 + 코드/fox — filter: closure of the valid region + λ>=0 + code/fox
    passed=[]; nearmiss=[]
    for labels,th,lam,rn in sols:
        if not actual_closure_ok(th): continue
        if all(l>=-1e-8 for l in lam):
            passed.append((labels,th,lam))
        else:
            nearmiss.append((labels,th,lam))
    print(f"[필터] 유효+λ>=0 통과: {len(passed)},  λ<0 탈락(유효): {len(nearmiss)}")
    # 군집: θ*와 거리 — clustering: distance from θ*
    def d2(th): return sum((th[k]-THSTAR[k])**2 for k in range(7))
    clus=collections.Counter()
    others=[]
    for labels,th,lam in passed:
        if d2(th)<1e-6: clus["θ*"]+=1
        else:
            others.append((labels,[round(v,5) for v in th]))
            clus["기타"]+=1
    print(f"[군집] {dict(clus)}")
    for lb,th in others[:12]:
        print("  기타:",lb,th)
    nm_clus=collections.Counter()
    nm_others=[]
    for labels,th,lam in nearmiss:
        if d2(th)<1e-6: nm_clus["θ*(대체부분집합)"]+=1
        else:
            nm_clus["기타(안장후보)"]+=1
            nm_others.append((labels,[round(v,5) for v in th],[round(v,4) for v in lam]))
    print(f"[λ<0 근 분석] {dict(nm_clus)}")
    for lb,th,lm in nm_others[:10]:
        print("  안장후보:",lb,"\n    θ=",th,"\n    λ=",lm)
    json.dump({"walls":WALLS,"n_strata":len(strata),"n_sol":len(sols),
               "n_pass":len(passed),"others":others[:50],
               "nearmiss_clusters":{k:v for k,v in nm_clus.items()},
               "saddle_candidates":nm_others[:20]},
              open("verify/n4_strata_result.json","w"))
