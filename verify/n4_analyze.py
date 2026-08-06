# n4_analyze.py — 타이트 상태 분석: 활성 접촉, 코어 외곽 다각형, 대칭, PSLQ — tight-state analysis: active contacts, core outline polygon, symmetry, PSLQ
import sys, json, math
sys.path.insert(0,'verify')
from n4_tighten2 import config, pieces, crossings, code_str, canon, ev_over_of, fox5_events, ell, covers

def margins(th, order):
    c=config(th); X=crossings(c); segs=pieces(c)
    out=[]
    for (i,j),s,u in X:
        _,_,li=segs[i]; _,_,lj=segs[j]
        out.append((f"교차({i},{j}) 소멸여유", min(s,li-s,u,lj-u)))
    for k in (4,5,6):
        out.append((f"L{k-3} -> 0", th[k]))
    rank={p:q for q,p in enumerate(order)}
    for i in range(4):
        ra,rb=rank[i],rank[i+1]; lo,hi=min(ra,rb),max(ra,rb)
        between=[order[r] for r in range(lo+1,hi)]
        if not between: continue
        (x1,y1),(x2,y2)=c['CH'][i]
        m=1e9
        for t in range(1,161):
            f=t/161; pt=(x1+f*(x2-x1),y1+f*(y2-y1))
            for k in between:
                P,dd,l=segs[k]
                rx,ry=pt[0]-P[0],pt[1]-P[1]
                lon=rx*dd[0]+ry*dd[1]; lat=abs(-rx*dd[1]+ry*dd[0])
                # 밴드 밖이면 양수 여유: 경계까지 거리 — outside the band, positive margin: distance to the boundary
                dlat=lat-0.5; dlon=max(-lon, lon-l)
                m=min(m, max(dlat,dlon,0) if (dlat>0 or dlon>0) else 0.0)
        out.append((f"크리스{i+1} 이웃공리 (사이:{between})", m))
    return out,c,X

def outline(th):
    """코어 = 밴드 1,2,3 합집합의 외곽 (shapely)
    core = the outline of the union of bands 1,2,3 (shapely)"""
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    c=config(th)
    from n4_render import band_quads
    q=band_quads(c)
    polys=[Polygon(q[k]) for k in (1,2,3)]
    U=unary_union(polys)
    return U, c

if __name__=="__main__":
    st=json.load(open(sys.argv[1] if len(sys.argv)>1 else "verify/n4_tight41.json"))
    th,order=st["theta"],st["order"]
    print(f"ell = {ell(th):.8f}")
    print(f"tau = {[round(v,8) for v in th[:4]]}")
    print(f"L   = {[round(v,8) for v in th[4:]]}")
    ms,c,X=margins(th,order)
    print("\n[여유(margin) — 0 근처 = 활성 접촉]")
    for name,m in sorted(ms,key=lambda t:t[1]):
        flag=" <== 활성" if m<0.02 else ""
        print(f"  {m:10.5f}  {name}{flag}")
    try:
        U,c=outline(th)
        B=list(U.exterior.coords)[:-1]
        print(f"\n[코어 외곽 다각형] 꼭짓점 {len(B)}개")
        for k,(px,py) in enumerate(B):
            a=B[(k-1)%len(B)]; b=B[(k+1)%len(B)]
            v1=(a[0]-px,a[1]-py); v2=(b[0]-px,b[1]-py)
            ang=math.degrees(math.acos(max(-1,min(1,(v1[0]*v2[0]+v1[1]*v2[1])/(math.hypot(*v1)*math.hypot(*v2))))))
            e=math.hypot(b[0]-px,b[1]-py)
            print(f"  V{k}: ({px:8.4f},{py:8.4f})  내각~{ang:7.2f}°  다음변 {e:.4f}")
        print(f"  면적 = {U.area:.6f},  둘레 = {U.exterior.length:.6f}")
        holes=len(U.interiors) if U.geom_type=="Polygon" else "?"
        print(f"  구멍 수 = {holes}")
    except ImportError:
        print("\n(shapely 없음 — pip 설치 필요)")
