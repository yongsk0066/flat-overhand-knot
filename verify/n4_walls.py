# Bisect numerical figure-eight cell boundaries and classify their contacts.
# 후보 벽: (a) 끝점-측면, (b) 끝점-캡, (d) 밴드코너-현 관통, (x) 교차 소멸, (L) 붕괴 — candidate walls: (a) endpoint-side, (b) endpoint-cap, (d) band-corner piercing a chord, (x) crossing vanishing, (L) collapse
import sys, math, random, json, collections
sys.path.insert(0,'verify')
from n4_tighten2 import config, pieces, crossings, valid_fixed2, ell

ORDER=[1,4,2,0,3]
PAIRS=[(0,[4,2]),(1,[4]),(2,[0]),(3,[2,0])]

def cap_data(c, k):
    """Actual crease caps: (label, normalized inward form, two corners)."""
    from n4_geometry import halfplanes
    planes=halfplanes(c,k,sign=lambda v:int(v>0)-int(v<0))
    labels=([('lo',k-1)] if k>0 else [])+([('hi',k)] if k<4 else [])
    out=[]
    for (name,crease),form in zip(labels,planes[2:]):
        norm=math.hypot(form[0],form[1])
        out.append((name,tuple(v/norm for v in form),c['CH'][crease]))
    return planes,out

def all_margins(th):
    """Candidate contact functions, windowed by the *actual* other faces.

    These are event diagnostics, not a conjunction defining feasibility.
    In particular, a negative lateral value outside a cap is not penetration;
    callers must use valid_fixed2 for the complete layer-dependent axiom.
    """
    from n4_geometry import halfplanes
    out={}; c=config(th)
    if c is None: return None
    segs=pieces(c)
    def distance(form,pt):
        a,b,z=form
        return (a*pt[0]+b*pt[1]+z)/math.hypot(a,b)
    for i,ks in PAIRS:
        A,B=c['CH'][i]; ab=(B[0]-A[0],B[1]-A[1]); nab2=sum(v*v for v in ab)
        for k in ks:
            planes,caps=cap_data(c,k)
            P,dd,_=segs[k]
            for e,pt in (('-',A),('+',B)):
                lat=-(pt[0]-P[0])*dd[1]+(pt[1]-P[1])*dd[0]
                if all(distance(form,pt)>=-0.02 for form in planes[2:]):
                    out[f"a:c{i+1}P{k}E{e}"]=abs(lat)-0.5
                for ci,(name,form,corners) in enumerate(caps):
                    if all(distance(other,pt)>=-0.02 for j,other in enumerate(planes) if j!=ci+2):
                        out[f"b:c{i+1}P{k}E{e}{name}"]=-distance(form,pt)
            for name,form,corners in caps:
                # The real cap corners are crease endpoints, not perpendicular
                # offsets of a centreline vertex.
                for pt in corners:
                    lat=-(pt[0]-P[0])*dd[1]+(pt[1]-P[1])*dd[0]
                    side='p' if lat>=0 else 'm'
                    t=((pt[0]-A[0])*ab[0]+(pt[1]-A[1])*ab[1])/nab2
                    if 0.03<t<0.97:
                        out[f"d:c{i+1}P{k}K{name}{side}"]=abs((pt[0]-A[0])*ab[1]-(pt[1]-A[1])*ab[0])/math.sqrt(nab2)
    for (ii,jj),s,u in crossings(c):
        _,_,li=segs[ii]; _,_,lj=segs[jj]
        out[f"x:({ii},{jj})s0"]=s
        out[f"x:({ii},{jj})u0"]=u
        if ii not in (0,4): out[f"x:({ii},{jj})s1"]=li-s
        if jj not in (0,4): out[f"x:({ii},{jj})u1"]=lj-u
    for k in (4,5,6): out[f"L{k-3}"]=th[k]
    return out

def boundary_probe(nprobe=4000, seed=11, starts=None, return_records=False):
    rnd=random.Random(seed)
    if starts is None:
        starts=[tuple(th) for th,o in json.load(open("verify/n4_valid41.json"))]
        try:
            cen=json.load(open("verify/n4_census.json"))
            starts+= [tuple(e["theta"]) for e in cen["endpoints"]]
        except Exception: pass
    hits=collections.Counter(); co=collections.Counter(); records=[]; skipped=collections.Counter()
    for p in range(nprobe):
        th0=list(rnd.choice(starts))
        # 내부 점 확보: 살짝 흔들어 유효 유지 — secure an interior point: jitter slightly while staying valid
        for _ in range(10):
            cand=[th0[k]+rnd.gauss(0,0.02) for k in range(7)]
            if valid_fixed2(tuple(cand),ORDER): th0=cand; break
        if not valid_fixed2(tuple(th0),ORDER):
            skipped["invalid_start"]+=1; continue
        dv=[rnd.gauss(0,1) for _ in range(7)]
        nn=math.sqrt(sum(v*v for v in dv)); dv=[v/nn for v in dv]
        lo,hi=0.0,None
        step=0.05
        t=step
        while t<3.0:
            cand=tuple(th0[k]+t*dv[k] for k in range(7))
            if not valid_fixed2(cand,ORDER): hi=t; break
            lo=t; t+=step
        if hi is None:
            skipped["no_boundary_before_3"]+=1; continue
        for _ in range(45):
            mid=(lo+hi)/2
            cand=tuple(th0[k]+mid*dv[k] for k in range(7))
            if valid_fixed2(cand,ORDER): lo=mid
            else: hi=mid
        bth=tuple(th0[k]+lo*dv[k] for k in range(7))
        m=all_margins(bth)
        if m is None: continue
        near=[(abs(v),kk) for kk,v in m.items() if abs(v)<5e-4]
        records.append(dict(attempt=p,theta=list(bth),start=th0,direction=dv,bracket=[lo,hi],near=sorted(near)))
        if not near:
            # 경계인데 후보 벽 없음 → 미분류 (code/canon/fox 전환?) — on the boundary but no candidate wall → unclassified (a code/canon/fox transition?)
            hits["(미분류)"]+=1; continue
        near.sort()
        hits[near[0][1]]+=1
        if len(near)>1: co[tuple(sorted(kk for _,kk in near[:3]))]+=1
    if return_records:
        return hits,co,records,dict(skipped)
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
