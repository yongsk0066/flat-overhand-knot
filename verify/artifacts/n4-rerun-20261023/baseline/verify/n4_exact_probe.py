# n4_exact_probe.py — 접촉 세분 + 대칭 축소 고정밀 수렴 + PSLQ — contact refinement + symmetry-reduced high-precision convergence + PSLQ
import sys, json, math
sys.path.insert(0,'verify')
from n4_tighten2 import config, pieces, crossings, ell
import mpmath as mp
mp.mp.dps=40

st=json.load(open("verify/n4_tight41.json"))
th=st["theta"]; order=st["order"]

# ---- 1) 접촉 세분: 각 (크리스 i, 사이패스 k) 쌍의 개별 여유 + 접촉점 위치 ---- — 1) contact refinement: the individual margin of each (crease i, in-between pass k) pair + the contact point location
c=config(th); segs=pieces(c)
rank={p:q for q,p in enumerate(order)}
print("[접촉 세분]")
for i in range(4):
    ra,rb=rank[i],rank[i+1]; lo,hi=min(ra,rb),max(ra,rb)
    between=[order[r] for r in range(lo+1,hi)]
    (x1,y1),(x2,y2)=c['CH'][i]
    for k in between:
        P,dd,l=segs[k]
        best=(1e9,None)
        for t in range(0,201):
            f=t/200; pt=(x1+f*(x2-x1),y1+f*(y2-y1))
            rx,ry=pt[0]-P[0],pt[1]-P[1]
            lon=rx*dd[0]+ry*dd[1]; lat=abs(-rx*dd[1]+ry*dd[0])
            dlat=lat-0.5; dlon=max(-lon,lon-l)
            m=max(dlat,dlon) if (dlat>0 or dlon>0) else max(dlat,dlon)
            if m<best[0]: best=(m,(f,round(lon,4),round(lat,4)))
        print(f"  크리스{i+1} vs P{k}: 여유 {best[0]:+.6f}  (현비율 f={best[1][0]:.3f}, lon={best[1][1]}, lat={best[1][2]}, 밴드길이 {l:.4f})")
