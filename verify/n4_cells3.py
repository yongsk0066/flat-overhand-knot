# Enumerate cyclic cuts and layer masks by integer knot determinant. Each pair of straight passes crosses at most once.
import itertools
from fractions import Fraction
from collections import defaultdict
W="ABCADCBD"

def det_knot(ev_over):
    """폐포 채색 행렬의 3x3 소행렬식 |det| (unknot 1, 3_1 3, 4_1 5)
    |det| of the 3x3 minor of the closure coloring matrix (unknot 1, 3_1 3, 4_1 5)"""
    A=4
    arc=0; visits=[]
    for cid,ov in ev_over:
        visits.append((cid,ov,arc%A))
        if not ov: arc+=1
    over_arc={}; und=defaultdict(list)
    for k,(cid,ov,a0) in enumerate(visits):
        if ov: over_arc[cid]=a0
        else: und[cid].append((a0,(a0+1)%A))
    rows=[]
    for cid in range(4):
        r=[0]*A
        ain,aout=und[cid][0]
        r[over_arc[cid]]+=2; r[ain]-=1; r[aout]-=1
        rows.append(r)
    # 3x3 소행렬 (마지막 행/열 제거) 정수 행렬식 — 3x3 minor (drop the last row/column), integer determinant
    m=[[Fraction(rows[i][j]) for j in range(3)] for i in range(3)]
    det=Fraction(1)
    for c in range(3):
        p=None
        for r_ in range(c,3):
            if m[r_][c]!=0: p=r_; break
        if p is None: return 0
        if p!=c: m[c],m[p]=m[p],m[c]; det=-det
        det*=m[c][c]
        for r_ in range(c+1,3):
            f=m[r_][c]/m[c][c]
            for cc in range(c,3): m[r_][cc]-=f*m[c][cc]
    return abs(int(det))

def pass_of(pos,cuts): return sum(1 for c in cuts if c<=pos)

cells_all=defaultdict(set); cells_mult=defaultdict(set)
count_41_masks=None
for r in range(8):
    word=W[r:]+W[:r]
    lab={}; seq=[]
    for ch in word:
        if ch not in lab: lab[ch]=len(lab)
        seq.append(lab[ch])
    first_seen={}
    for i,cid in enumerate(seq):
        if cid not in first_seen: first_seen[cid]=i
    for mask in range(16):
        ev=[]
        seen=set()
        for cid in seq:
            firstv = cid not in seen
            seen.add(cid)
            ov = bool((mask>>cid)&1) if firstv else not bool((mask>>cid)&1)
            ev.append((cid,ov))
        d=det_knot(ev)
        if r==0:
            if count_41_masks is None: count_41_masks={1:0,3:0,5:0}
            count_41_masks[d]=count_41_masks.get(d,0)+1
        if d!=5: continue
        # 컷 삽입 — insert the cuts
        for cuts in itertools.combinations_with_replacement(range(9),4):
            pos={cid:[i for i,c in enumerate(seq) if c==cid] for cid in range(4)}
            pa={cid:(pass_of(pos[cid][0],cuts),pass_of(pos[cid][1],cuts)) for cid in range(4)}
            if any(abs(p1-p2)<2 for p1,p2 in pa.values()): continue
            pairs=[tuple(sorted(v)) for v in pa.values()]
            # 층 비순환 (over 패스가 위) — layer acyclicity (the over pass goes on top)
            edges=set()
            for cid in range(4):
                p1,p2=pa[cid]
                # 첫 방문 over면 첫 방문 패스가 위 — if the first visit is over, the first-visit pass is on top
                fv_pass = pass_of(pos[cid][0],cuts)
                sv_pass = pass_of(pos[cid][1],cuts)
                ov_first = bool((mask>>cid)&1)
                top,bot=(fv_pass,sv_pass) if ov_first else (sv_pass,fv_pass)
                edges.add((top,bot))
            indeg={n:0 for n in range(5)}; g={n:set() for n in range(5)}
            for a,b in edges:
                if b not in g[a]: g[a].add(b); indeg[b]+=1
            q=[n for n in range(5) if indeg[n]==0]; seen2=0
            while q:
                n=q.pop(); seen2+=1
                for mm in g[n]:
                    indeg[mm]-=1
                    if indeg[mm]==0: q.append(mm)
            if seen2!=5: continue
            key=tuple(sorted(pairs))
            cells_all[key].add((r,cuts,mask))
            if len(set(pairs))==4:      # 쌍 중복 금지 (직선성 정리) — no duplicate pairs (straightness lemma)
                cells_mult[key].add((r,cuts,mask))
print("[마스크 판정 (회전 0)] det 분포:",count_41_masks," (기대 {1:12, 3:2, 5:2})")
print(f"[세포 전수 v3] 4_1 층 허용 쌍조합: 중복 허용 {len(cells_all)}종 / 쌍중복 금지 {len([k for k in cells_mult if cells_mult[k]])}종")
print("\n[쌍중복 금지 세포 목록]")
harvest={((0,2),(0,3),(1,4),(2,4)):"실현(수확)",((0,3),(0,4),(1,4),(2,4)):"",((0,2),(0,3),(0,4),(1,4)):""}
for k in sorted(cells_mult):
    if not cells_mult[k]: continue
    tag=harvest.get(k,"")
    print(f"  {k}  변형 {len(cells_mult[k])}개  {tag}")
print("\n[대조] 수확에서 '4_1층 존재=0'였던 쌍조합의 포함 여부:")
for k in [((0,2),(0,3),(0,4),(2,4)),((0,2),(0,4),(1,4),(2,4))]:
    print(f"  {k}: {'포함(문제!)' if k in cells_mult and cells_mult[k] else '배제 ✓'}")
print("\n[증인 쌍조합] ((0,2),(0,3),(0,4),(1,3)):",
      "포함(기하 불가 주장 철회 필요)" if ((0,2),(0,3),(0,4),(1,3)) in cells_mult and cells_mult[((0,2),(0,3),(0,4),(1,3))] else "배제")
