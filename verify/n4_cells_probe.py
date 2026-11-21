# n4_cells_probe.py — 미실현 8세포의 표적 배제 프로브 [수치-강] — targeted exclusion probe for the 8 unrealized cells [numerics-strong]
# 각 세포(쌍조합): 목적함수 = 쌍조합 불일치 수 → 언덕오르기 다중시작. — per cell (pair combination): objective = number of pair-combination mismatches → multi-start hill climbing.
# 판정: (a) 셰도우 실현 가능한가 (b) 되면 4_1층+정확공리 통과 상태 존재하는가 — verdict: (a) is the shadow realizable (b) if so, does a state passing the 4_1 layers + exact axioms exist
import sys, math, random, json, collections
sys.path.insert(0,'verify')
from n4_tighten2 import config, crossings, code_str, canon, orders_making_41, axiom_exact_ok, TARGET_CANON, general_position_ok
import multiprocessing as mp

CELLS=[((0,2),(0,3),(0,4),(1,3)),((0,2),(0,3),(0,4),(1,4)),((0,2),(0,3),(1,3),(1,4)),
       ((0,2),(0,3),(1,3),(2,4)),((0,2),(1,3),(1,4),(2,4)),((0,3),(0,4),(1,4),(2,4)),
       ((0,3),(1,3),(1,4),(2,4)),((0,4),(1,3),(1,4),(2,4))]

def mismatch(th, target):
    c=config(th)
    if c is None: return 99,None
    X=crossings(c)
    pairs=sorted(p for p,_,_ in X)
    t=list(target)
    # 다중집합 대칭차 — multiset symmetric difference
    from collections import Counter
    ca,cb=Counter(pairs),Counter(t)
    d=sum((ca-cb).values())+sum((cb-ca).values())
    return d,X

def _work(args):
    cell, seed, restarts, steps = args
    rnd=random.Random(seed)
    best=99; shadow_hits=0; valid_hits=0; best_th=None
    for r in range(restarts):
        th=[rnd.choice([1,-1])*rnd.uniform(0.3,3.1) for _ in range(4)]+[rnd.uniform(0.3,4.0) for _ in range(3)]
        cur,_=mismatch(tuple(th),cell)
        stall=0; step=0.25
        for it in range(steps):
            cand=[th[k]+step*rnd.gauss(0,1) for k in range(7)]
            if any(cand[k]<0.05 for k in (4,5,6)): continue
            d,X=mismatch(tuple(cand),cell)
            if d<cur or (d==cur and rnd.random()<0.3):
                th=cand; cur=d; stall=0
            else:
                stall+=1
                if stall>60: step=max(step*0.7,0.02); stall=0
            if cur==0:
                shadow_hits+=1
                c=config(tuple(th)); X=crossings(c)
                cd=code_str(X)
                if cd and canon(cd)==TARGET_CANON and general_position_ok(c,X):
                    for m,exts in orders_making_41(X):
                        for o in exts:
                            if axiom_exact_ok(c,o):
                                valid_hits+=1; best_th=list(th)
                break
        best=min(best,cur)
    return cell,best,shadow_hits,valid_hits,best_th

if __name__=="__main__":
    jobs=[(cell, 1000+i, 60, 4000) for i,cell in enumerate(CELLS)]
    with mp.get_context("fork").Pool(6) as pool:
        res=pool.map(_work,jobs)
    print(f"{'세포':44} {'최소불일치':>6} {'셰도우':>6} {'유효':>4}")
    out={}
    for cell,best,sh,vh,bt in res:
        print(f"{str(cell):44} {best:6d} {sh:6d} {vh:4d}")
        out[str(cell)]={"best_mismatch":best,"shadow":sh,"valid":vh,"theta":bt}
    json.dump(out,open("verify/n4_cells_probe.json","w"))
    n_sh=sum(1 for c,b,s,v,t in res if s>0); n_v=sum(1 for c,b,s,v,t in res if v>0)
    print(f"\n[요약] 셰도우 실현 세포 {n_sh}/8, 유효(4_1+공리) 세포 {n_v}/8")
