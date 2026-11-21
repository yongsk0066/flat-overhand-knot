# n4_kkt.py — 정확 θ*에서 KKT 차단 검사: ∇ell ∈ cone{∇g_i} (8 활성 접촉) — KKT blocked-ness check at the exact θ*: ∇ell ∈ cone{∇g_i} (8 active contacts)
import sys, math, itertools, json
sys.path.insert(0,'verify')
from n4_tighten2 import config, valid_fixed2, crossings, code_str, canon, ev_over_of, fox5_events

t=math.acos(-0.25); s=2*math.pi-2*t
L1=2*math.sqrt(15)/5; L2=8*math.sqrt(15)/15
TH=[t,s,-s,-t,L1,L2,L1]
order=[1,4,2,0,3]

def lat(pt,P,dd): return -(pt[0]-P[0])*dd[1]+(pt[1]-P[1])*dd[0]

def gvec(th):
    c=config(th)
    d=c['d']; M=c['M']; CH=c['CH']
    E1m,E1p=CH[0]; E2m,E2p=CH[1]; E3m,E3p=CH[2]; E4m,E4p=CH[3]
    P4=(M[3],d[4]); P2=(M[1],d[2]); P0=(M[0],d[0])
    return [
        -lat(E1p,*P4)-0.5,   # g1
         lat(E1m,*P2)-0.5,   # g2
         lat(E2p,*P4)-0.5,   # g3a
         lat(E2m,*P4)-0.5,   # g3b
         lat(E4m,*P0)-0.5,   # g1'
        -lat(E4p,*P2)-0.5,   # g2'
        -lat(E3p,*P0)-0.5,   # g3a'
        -lat(E3m,*P0)-0.5,   # g3b'
    ]

g0=gvec(TH)
print("활성 잔차(전부 ~0이어야):", [f"{v:.2e}" for v in g0])
# 수치 기울기 — numerical gradients
h=1e-7
G=[]
for i in range(7):
    tp=TH[:]; tm=TH[:]
    tp[i]+=h; tm[i]-=h
    gp=gvec(tp); gm=gvec(tm)
    G.append([(gp[k]-gm[k])/(2*h) for k in range(8)])
# G[i][k] = ∂g_k/∂θ_i ; 목표: ∇ell=(0,0,0,0,1,1,1) = Σ λ_k ∇g_k, λ≥0 — goal: write ∇ell as Σ λ_k ∇g_k with λ≥0
import mpmath as mp
mp.mp.dps=30
target=[0,0,0,0,1,1,1]
best=None
for r in range(4,9):
    for sub in itertools.combinations(range(8),r):
        A=mp.matrix(7,r)
        for col,k in enumerate(sub):
            for row in range(7): A[row,col]=G[row][k]
        try:
            lam=mp.lu_solve(A.T*A, A.T*mp.matrix(target))
        except Exception:
            continue
        res=A*lam-mp.matrix(target)
        rn=mp.norm(res)
        if all(l>=-mp.mpf('1e-9') for l in lam) and rn<mp.mpf('1e-6'):
            if best is None or rn<best[0]: best=(rn,sub,[float(l) for l in lam])
if best:
    rn,sub,lam=best
    names=['g1','g2','g3a','g3b',"g1'","g2'","g3a'","g3b'"]
    print("\n[KKT 차단 인증] ∇ell ∈ cone — PASS")
    print("  부분집합:", [names[k] for k in sub])
    print("  λ =", [round(l,6) for l in lam], " 잔차:", mp.nstr(rn,3))
else:
    print("\n[KKT] 원뿔 소속 실패 — 자유 하강 방향 존재 가능! Farkas 정밀검토 필요")
# 유효성 재확인 (미세 내부 섭동) — validity re-check (tiny interior perturbation)
import random
rnd=random.Random(5)
cnt=0
for _ in range(200):
    d=[rnd.gauss(0,1) for _ in range(7)]
    n=sum(x*x for x in d)**0.5
    th2=[TH[k]+1e-5*d[k]/n for k in range(7)]
    if valid_fixed2(tuple(th2),order): cnt+=1
print(f"근방 유효율(반경 1e-5, 200방향): {cnt}/200  (경계점이므로 ~절반 이하 기대)")
c=config(TH); X=crossings(c)
print("θ*에서 fox5:", fox5_events(ev_over_of(X,order)), " canon:", canon(code_str(X)))

# Test linearized accessibility with a feasible inward direction.
print("\n[접근가능성] 8접촉 원뿔의 내부 방향 탐색 + 비선형 검사")
import mpmath as mpx
G8=[[G[i][k] for i in range(7)] for k in range(8)]   # ∇g_k
best=None
rnd2=random.Random(9)
v=[rnd2.gauss(0,1) for _ in range(7)]
for it in range(4000):
    dots=[sum(G8[k][q]*v[q] for q in range(7)) for k in range(8)]
    m=min(dots); kmin=dots.index(m)
    nn=math.sqrt(sum(x*x for x in v))
    if best is None or m/nn>best: best=m/nn
    v=[v[q]+0.05*G8[kmin][q] for q in range(7)]
nn=math.sqrt(sum(x*x for x in v)); v=[x/nn for x in v]
dots=[sum(G8[k][q]*v[q] for q in range(7)) for k in range(8)]
print(f"  내부 방향 v: min ∇g·v = {min(dots):.5f} (전부 양수: {all(d>0 for d in dots)})")
from n4_walls import all_margins
for r in (1e-1,1e-2,1e-3,1e-4,1e-5,1e-6):
    th2=tuple(TH[q]+r*v[q] for q in range(7))
    ok=valid_fixed2(th2,order) is not None
    m=all_margins(th2)
    worst=min(vv for kk,vv in m.items() if kk.startswith("a:"))
    print(f"  r={r:.0e}: valid={ok}  최소 슬리버 여유={worst:.3e}  (여유/r={worst/r:.4f})")
