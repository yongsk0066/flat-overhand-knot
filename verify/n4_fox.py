# n4_fox.py — 실현 배치의 매듭형 오라클: Fox 5-채색 수 (25 <=> 4_1; det: unknot 1, 3_1 3, 4_1 5) — knot-type oracle for a realized configuration: the Fox 5-coloring count (25 <=> 4_1; det: unknot 1, 3_1 3, 4_1 5)
import math

def fox5_count(X, RANK):
    """X: crossings [((i,j),s,u)] material params (piece0 param already negated upstream? here raw).
    재료 순서 이벤트: (piece, param_along_material, crossing_id, is_over)
    events in material order: (piece, param_along_material, crossing_id, is_over)"""
    ev=[]
    for cid,((i,j),s,u) in enumerate(X):
        over_i = RANK[i]<RANK[j]
        ev.append((i, -s if i==0 else s, cid, over_i))
        ev.append((j, -u if j==0 else u, cid, not over_i))
    ev.sort()
    # 호(arc): under 방문 뒤에서 절단; 폐포(장매듭 닫기)로 첫/끝 호 동일시 — arcs: cut just after each under-visit; the closure (closing the long knot) identifies the first and last arcs
    narc=sum(1 for e in ev if not e[3])
    if narc==0: return 5  # 모두 over -> unknot — all over -> unknot
    arc=0; arc_of_visit=[]
    for (p,t,cid,ov) in ev:
        arc_of_visit.append(arc)
        if not ov: arc=(arc+1)
    # 마지막 호 == 호 0 (닫기) — the last arc == arc 0 (closing)
    A=narc
    def a(k): return k % A
    # 관계: 각 교차 c: 2*over_arc - under_in - under_out = 0 (mod 5) — relation: at each crossing c: 2*over_arc - under_in - under_out = 0 (mod 5)
    rows=[]
    under_in={}; under_out={}; over_arc={}
    for idx,(p,t,cid,ov) in enumerate(ev):
        if ov: over_arc[cid]=a(arc_of_visit[idx])
        else:
            under_in[cid]=a(arc_of_visit[idx]); under_out[cid]=a(arc_of_visit[idx]+1)
    for cid in over_arc:
        r=[0]*A
        r[over_arc[cid]]+=2; r[under_in[cid]]-=1; r[under_out[cid]]-=1
        rows.append([x%5 for x in r])
    # F5 널리티 — nullity over F5
    m=[r[:] for r in rows]; rank=0; col=0
    R=len(m)
    while rank<R and col<A:
        piv=None
        for rr in range(rank,R):
            if m[rr][col]%5: piv=rr; break
        if piv is None: col+=1; continue
        m[rank],m[piv]=m[piv],m[rank]
        inv=pow(m[rank][col],3,5)  # a^-1 = a^3 mod 5
        m[rank]=[(v*inv)%5 for v in m[rank]]
        for rr in range(R):
            if rr!=rank and m[rr][col]%5:
                f=m[rr][col]
                m[rr]=[(m[rr][k]-f*m[rank][k])%5 for k in range(A)]
        rank+=1; col+=1
    nullity=A-rank
    return 5**nullity
