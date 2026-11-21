# n4_cells.py — 추상 전수: ABCADCBD에 크리스 4개 삽입(패스 분해) 전부에 대해 — abstract enumeration: over every insertion of 4 creases into ABCADCBD (pass decomposition),
# 4_1 배정 (+,-,-,+)/(-,+,+,-)가 비순환 층 그래프(=평면 접힘 층 전순서 존재)를 주는가 — does the 4_1 assignment (+,-,-,+)/(-,+,+,-) yield an acyclic layer graph (= a flat-folding layer total order exists)
import itertools, sys
WORD="ABCADCBD"
POS={ch:[i for i,c in enumerate(WORD) if c==ch] for ch in "ABCD"}
ASSIGNS=[(1,-1,-1,1),(-1,1,1,-1)]   # 문서17: V(4_1) 배정 (첫 방문 over = +1), 거울쌍 — doc 17: V(4_1) assignment (first visit over = +1), mirror pair

def pass_of(pos, cuts):
    # cuts: 오름차순 4개 (0..8), 갭 g = 위치 g 앞. pass = 잘린 구간 인덱스 — cuts: 4 values in ascending order (0..8), gap g = just before position g. pass = index of the cut interval
    return sum(1 for c in cuts if c<=pos)

ok_cells=[]
tot=0; geom=0
for cuts in itertools.combinations_with_replacement(range(9),4):
    tot+=1
    pa={ch:(pass_of(POS[ch][0],cuts), pass_of(POS[ch][1],cuts)) for ch in "ABCD"}
    # 기하 실현 필요조건: 각 교차의 두 방문은 서로 다른, 비인접 패스 — necessary condition for geometric realization: the two visits of each crossing lie on distinct, non-adjacent passes
    if any(abs(p1-p2)<2 for p1,p2 in pa.values()): continue
    geom+=1
    for asg in ASSIGNS:
        edges=[]
        for k,ch in enumerate("ABCD"):
            p1,p2=pa[ch]
            top,bot=(p1,p2) if asg[k]==1 else (p2,p1)
            edges.append((top,bot))
        # DAG 검사 (5 노드) — DAG check (5 nodes)
        from collections import defaultdict
        g=defaultdict(set)
        for a,b in edges: g[a].add(b)
        # 위상정렬 — topological sort
        import copy
        indeg={n:0 for n in range(5)}
        for a,b in edges: indeg[b]+=1
        q=[n for n in range(5) if indeg[n]==0]; seen=0
        gg={n:set() for n in range(5)}
        for a,b in edges: gg[a].add(b)
        while q:
            n=q.pop(); seen+=1
            for m2 in gg[n]:
                indeg[m2]-=1
                if indeg[m2]==0: q.append(m2)
        if seen==5:
            ok_cells.append((cuts, pa, asg))

print(f"분해 총 {tot}, 기하 필요조건(비인접·이패스) 통과 {geom}, 4_1 비순환 층 존재: {len(ok_cells)}")
seenpair={}
for cuts,pa,asg in ok_cells:
    pairs=tuple(sorted(tuple(sorted(v)) for v in pa.values()))
    seenpair.setdefault(pairs,[]).append((cuts,asg))
print(f"서로 다른 패스-쌍 조합: {len(seenpair)}")
for pairs,lst in sorted(seenpair.items()):
    cuts0,asg0=lst[0]
    print(f"  쌍 {pairs}  (예: cuts={cuts0}, asg={asg0}, 변형 {len(lst)}개)")
