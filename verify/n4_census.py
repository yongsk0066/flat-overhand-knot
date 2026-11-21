# Sample four-crossing configurations and record layer orders and fixed-order descents.
# ① 대형 블라인드 수확(세포 표 작성) ② 유효 시작점 전부에서 층고정 조임 ③ 군집 분석 vs θ* — ① large-scale blind harvest (build the cell table) ② layer-fixed tightening from every valid start point ③ cluster analysis vs θ*
import sys, math, random, json, collections
sys.path.insert(0,'verify')
from n4_tighten2 import (config, crossings, code_str, canon, orders_making_41,
                         axiom_ok, axiom_exact_ok, valid_fixed2, tighten_fixed2, ell,
                         ev_over_of, fox5_events, TARGET_CANON, general_position_ok)
import multiprocessing as mp

def _hworker(args):
    seed,n=args
    rnd=random.Random(seed)
    cells=collections.Counter()      # (pairing, tuple(order)) -> count (공리까지 통과) — passes the axioms as well
    cells_shadow=collections.Counter()
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
        if not general_position_ok(c,X): continue
        key=tuple(sorted(p for p,_,_ in X))
        cells_shadow[key]+=1
        for m,exts in orders_making_41(X):
            for o in exts:
                if axiom_exact_ok(c,o):
                    cells[(key,tuple(o))]+=1
                    hits.append((list(th),o))
                    break
    return cells_shadow,cells,hits

def _dworker(args):
    seed,th0,order=args
    rnd=random.Random(seed)
    thf,ef=tighten_fixed2(th0,order,rnd)
    # 폴리시 한 번 더 (정체 탈출) — one more polish pass (escape stagnation)
    thf2,ef2=tighten_fixed2(list(thf),order,rnd,rounds=10)
    if ef2<ef: thf,ef=thf2,ef2
    return ef,list(thf),order

if __name__=="__main__":
    TOTAL=8_000_000; P=8
    with mp.get_context("fork").Pool(P) as pool:
        parts=pool.map(_hworker,[(7000+i,TOTAL//P) for i in range(P)])
    shadow=collections.Counter(); cells=collections.Counter(); hits=[]
    for s_,c_,h_ in parts:
        shadow.update(s_); cells.update(c_); hits+=h_
    print(f"[수확] {TOTAL:,} 표본, 셰도우 통과 {sum(shadow.values())}")
    print("셰도우 세포 분포:", dict(shadow))
    print(f"\n[4_1+정확공리 세포] (쌍조합, 층순서) — 유효 상태 총 {len(hits)}개")
    for (key,o),n in sorted(cells.items(),key=lambda kv:-kv[1]):
        print(f"  {key}  층{list(o)}: {n}")
    # 하강: 최대 64 시작점 — descent: at most 64 start points
    starts=hits[:64]
    print(f"\n[하강 센서스] 시작점 {len(starts)}개 (8-way)")
    with mp.get_context("fork").Pool(P) as pool:
        res=pool.map(_dworker,[(9000+i,tuple(th),o) for i,(th,o) in enumerate(starts)])
    res.sort(key=lambda r:r[0])
    ELL_STAR=4*math.sqrt(15)/3
    clusters=collections.Counter()
    for ef,thf,o in res:
        clusters[round(ef,2)]+=1
    print("종점 ell 분포(0.01 반올림):", dict(sorted(clusters.items())))
    hit_star=sum(1 for ef,_,_ in res if abs(ef-ELL_STAR)<0.01)
    print(f"θ*(ell*={ELL_STAR:.6f}) 도달: {hit_star}/{len(res)}")
    # θ* 근접 종점의 좌표 대조 (대칭 잔차) — coordinate comparison of the endpoints near θ* (symmetry residual)
    best=res[0]
    ef,thf,o=best
    print(f"최소 종점: ell={ef:.7f}, tau={[round(v,5) for v in thf[:4]]}, L={[round(v,5) for v in thf[4:]]}")
    print(f"  대칭잔차: t4+t1={thf[3]+thf[0]:.2e}, t3+t2={thf[2]+thf[1]:.2e}, L3-L1={thf[6]-thf[4]:.2e}")
    json.dump({"clusters":{str(k):v for k,v in clusters.items()},
               "hit_star":hit_star,"n":len(res),
               "endpoints":[{"ell":ef,"theta":thf,"order":o} for ef,thf,o in res[:20]]},
              open("verify/n4_census.json","w"))
