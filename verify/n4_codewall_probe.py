# n4_codewall_probe.py — 코드 벽(교차 소멸)·L 벽이 유효 폐포에 닿는가: 표적 최소화 — do the code walls (crossing annihilation) and the L walls touch the valid closure: targeted minimization
import sys, math, random, json
sys.path.insert(0,'verify')
from n4_tighten2 import config, pieces, crossings, valid_fixed2
ORDER=[1,4,2,0,3]
def xmargins(th):
    c=config(th)
    if c is None: return None
    segs=pieces(c); out={}
    for (i,j),s,u in crossings(c):
        _,_,li=segs[i]; _,_,lj=segs[j]
        out[f"({i},{j})s"]=s; out[f"({i},{j})u"]=u
        if li<50: out[f"({i},{j})S"]=li-s
        if lj<50: out[f"({i},{j})U"]=lj-u
    for k in (4,5,6): out[f"L{k-3}"]=th[k]
    return out
rnd=random.Random(21)
starts=[tuple(th) for th,o in json.load(open("verify/n4_valid41.json")) if th[0]>0]
keys=sorted(xmargins(starts[0]).keys())
best={k:1e9 for k in keys}
for k in keys:
    for th0 in starts[:4]:
        th=list(th0); step=0.05
        cur=xmargins(tuple(th))[k]
        for it in range(3000):
            adv=False
            for _ in range(20):
                cand=[th[q]+step*rnd.gauss(0,1) for q in range(7)]
                if not valid_fixed2(tuple(cand),ORDER): continue
                m=xmargins(tuple(cand))
                if m and m[k]<cur-1e-12:
                    th=cand; cur=m[k]; adv=True; break
            if not adv:
                step*=0.6
                if step<1e-6: break
        best[k]=min(best[k],cur)
print("[코드/L 벽 도달성 프로브] 유효 영역 내 최소값 (0 근접 = 벽이 폐포에 닿음):")
for k in keys:
    print(f"  {k:10} inf ≈ {best[k]:.4f}")
json.dump(best,open("verify/n4_codewall_inf.json","w"))
