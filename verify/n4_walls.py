# n4_walls.py — n=4 벽 분류 (문서 06 유사물): 유효 세포 경계의 이분 표집 + 벽 동정 — n=4 wall classification (analogue of document 06): bisection sampling of the valid-cell boundary + wall identification
# 후보 벽: (a) 끝점-측면, (b) 끝점-캡, (d) 밴드코너-현 관통, (x) 교차 소멸, (L) 붕괴 — candidate walls: (a) endpoint-side, (b) endpoint-cap, (d) band-corner piercing a chord, (x) crossing vanishing, (L) collapse
import sys, math, random, json, collections
sys.path.insert(0,'verify')
from n4_tighten2 import config, pieces, crossings, valid_fixed2, ell

ORDER=[1,4,2,0,3]
PAIRS=[(0,[4,2]),(1,[4]),(2,[0]),(3,[2,0])]

def all_margins(th):
    out={}
    c=config(th)
    if c is None: return None
    segs=pieces(c)
    for i,ks in PAIRS:
        A,B=c['CH'][i]
        ab=(B[0]-A[0],B[1]-A[1]); nab2=ab[0]**2+ab[1]**2
        for k in ks:
            P,dd,l=segs[k]
            n=(-dd[1],dd[0])
            for e,pt in (('-',A),('+',B)):
                lat=-(pt[0]-P[0])*dd[1]+(pt[1]-P[1])*dd[0]
                lon=(pt[0]-P[0])*dd[0]+(pt[1]-P[1])*dd[1]
                if -0.02*l<lon<1.02*l:
                    out[f"a:c{i+1}P{k}E{e}"]=abs(lat)-0.5
                if abs(lat)<0.52:
                    if lon<0.5*l: out[f"b:c{i+1}P{k}E{e}lo"]=abs(lon)
                    elif l<50:    out[f"b:c{i+1}P{k}E{e}hi"]=abs(l-lon)
            caps=[(0.0,'lo')]+([(l,'hi')] if l<50 else [])
            for cap,cn in caps:
                for sgn,sn in ((0.5,'p'),(-0.5,'m')):
                    cx=P[0]+cap*dd[0]+sgn*n[0]; cy=P[1]+cap*dd[1]+sgn*n[1]
                    t=((cx-A[0])*ab[0]+(cy-A[1])*ab[1])/nab2
                    if 0.03<t<0.97:
                        px,py=A[0]+t*ab[0],A[1]+t*ab[1]
                        out[f"d:c{i+1}P{k}K{cn}{sn}"]=math.hypot(cx-px,cy-py)
    X=crossings(c)
    for (ii,jj),s,u in X:
        _,_,li=segs[ii]; _,_,lj=segs[jj]
        out[f"x:({ii},{jj})s0"]=s
        out[f"x:({ii},{jj})u0"]=u
        if li<50: out[f"x:({ii},{jj})s1"]=li-s
        if lj<50: out[f"x:({ii},{jj})u1"]=lj-u
    for k in (4,5,6): out[f"L{k-3}"]=th[k]
    return out

def boundary_probe(nprobe=4000, seed=11):
    rnd=random.Random(seed)
    starts=[tuple(th) for th,o in json.load(open("verify/n4_valid41.json"))]
    try:
        cen=json.load(open("verify/n4_census.json"))
        starts+= [tuple(e["theta"]) for e in cen["endpoints"]]
    except Exception: pass
    hits=collections.Counter(); co=collections.Counter()
    for p in range(nprobe):
        th0=list(rnd.choice(starts))
        # 내부 점 확보: 살짝 흔들어 유효 유지 — secure an interior point: jitter slightly while staying valid
        for _ in range(10):
            cand=[th0[k]+rnd.gauss(0,0.02) for k in range(7)]
            if valid_fixed2(tuple(cand),ORDER): th0=cand; break
        if not valid_fixed2(tuple(th0),ORDER): continue
        dv=[rnd.gauss(0,1) for _ in range(7)]
        nn=math.sqrt(sum(v*v for v in dv)); dv=[v/nn for v in dv]
        lo,hi=0.0,None
        step=0.05
        t=step
        while t<3.0:
            cand=tuple(th0[k]+t*dv[k] for k in range(7))
            if not valid_fixed2(cand,ORDER): hi=t; break
            lo=t; t+=step
        if hi is None: continue
        for _ in range(45):
            mid=(lo+hi)/2
            cand=tuple(th0[k]+mid*dv[k] for k in range(7))
            if valid_fixed2(cand,ORDER): lo=mid
            else: hi=mid
        bth=tuple(th0[k]+lo*dv[k] for k in range(7))
        m=all_margins(bth)
        if m is None: continue
        near=[(v,kk) for kk,v in m.items() if v<5e-4]
        if not near:
            # 경계인데 후보 벽 없음 → 미분류 (code/canon/fox 전환?) — on the boundary but no candidate wall → unclassified (a code/canon/fox transition?)
            hits["(미분류)"]+=1; continue
        near.sort()
        hits[near[0][1]]+=1
        if len(near)>1: co[tuple(sorted(kk for _,kk in near[:3]))]+=1
    return hits,co

if __name__=="__main__":
    hits,co=boundary_probe()
    tot=sum(hits.values())
    print(f"[경계 표집] {tot}회 경계 도달, 벽 동정:")
    for k,v in hits.most_common():
        print(f"  {k:24} {v:5}  ({100*v/tot:.1f}%)")
    print("\n[동시 근접(코너) 상위]")
    for k,v in co.most_common(10):
        print(f"  {v:4}  {k}")
    json.dump({k:v for k,v in hits.items()}, open("verify/n4_wall_inventory.json","w"))
