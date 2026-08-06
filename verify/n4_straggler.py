# n4_straggler.py — 하강 정체점의 Farkas 배제 [referee-H F8 재현] — Farkas exclusion of descent stall points [reproducing referee-H F8]
# 각 센서스 종점: 활성 접촉(여유<1e-3) 기울기 원뿔에 ∇ell 소속 검사 — 잔차 크면 차단점 아님 — for each census endpoint: test whether ∇ell lies in the gradient cone of the active contacts (margin<1e-3); a large residual means it is not a blocked point
import sys, json, math, itertools
sys.path.insert(0,'verify')
from n4_walls import all_margins
from n4_strata2 import margin_fn, grad, GRAD_ELL
ELLSTAR=4*math.sqrt(15)/3
cen=json.load(open("verify/n4_census.json"))
print(f"{'ell-ell*':>10} {'활성':>4} {'Farkas 잔차':>12}  판정")
for e in cen["endpoints"]:
    th=tuple(e["theta"])
    m=all_margins(th)
    act=[k for k,v in m.items() if k.startswith("a:") and v<1e-3]
    fns=[margin_fn(lb,th) for lb in act]
    G=[grad(f,th) for f in fns]
    # 비음수 최소제곱: 부분집합 전수 (<=8) — nonnegative least squares: exhaustive over subsets (<=8)
    best=None
    import mpmath as mp
    for r in range(1,len(G)+1):
        for sub in itertools.combinations(range(len(G)),r):
            A=mp.matrix(7,r)
            for c,k in enumerate(sub):
                for q in range(7): A[q,c]=G[k][q]
            try: lam=mp.lu_solve(A.T*A, A.T*mp.matrix(GRAD_ELL))
            except Exception: continue
            if any(l<-1e-9 for l in lam): continue
            res=float(mp.norm(A*lam-mp.matrix(GRAD_ELL)))
            if best is None or res<best: best=res
    gap=e["ell"]-ELLSTAR
    verdict="차단(θ*)" if (best is not None and best<1e-6) else "차단 아님 — 하강 가능"
    print(f"{gap:10.5f} {len(act):4d} {best if best is not None else float('nan'):12.4g}  {verdict}")

# ---- 정체점 생성 + 판정 (짧은 예산 하강으로 정체 유도) ---- — generate stall points + adjudicate them (induce stalls via a short-budget descent)
if __name__=="__main__" and "--stalls" in sys.argv:
    import random
    from n4_tighten2 import tighten_fixed2, ell
    rnd=random.Random(4242)
    hits=json.load(open("verify/n4_valid41.json"))
    print("\n[정체 유도 하강 (rounds=4) + Farkas]")
    print(f"{'ell-ell*':>10} {'활성':>4} {'Farkas 잔차':>12}  판정")
    import mpmath as mp
    for th0,order in hits[:10]:
        thf,ef=tighten_fixed2(list(th0),order,rnd,rounds=4)
        m=all_margins(tuple(thf))
        act=[k for k,v in m.items() if k.startswith("a:") and v<1e-3]
        if not act:
            print(f"{ef-ELLSTAR:10.5f} {0:4d} {'—':>12}  활성 없음(내부) — 차단 아님")
            continue
        fns=[margin_fn(lb,tuple(thf)) for lb in act]
        G=[grad(f,tuple(thf)) for f in fns]
        best=None
        for r in range(1,len(G)+1):
            for sub in itertools.combinations(range(len(G)),r):
                A=mp.matrix(7,r)
                for c,k in enumerate(sub):
                    for q in range(7): A[q,c]=G[k][q]
                try: lam=mp.lu_solve(A.T*A, A.T*mp.matrix(GRAD_ELL))
                except Exception: continue
                if any(l<-1e-9 for l in lam): continue
                res=float(mp.norm(A*lam-mp.matrix(GRAD_ELL)))
                if best is None or res<best: best=res
        verdict="차단(θ*)" if (best is not None and best<1e-6) else "차단 아님 — 1차 하강방향 존재"
        print(f"{ef-ELLSTAR:10.5f} {len(act):4d} {best if best is not None else float('nan'):12.4g}  {verdict}")
