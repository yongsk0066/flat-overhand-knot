# n4_cells2.py — 세포 전수 v2: 폐포점(회전) 8종 전부. 장매듭 절단 위치가 어디든, — cell enumeration v2: all 8 closure points (rotations). Wherever the long-knot cut is placed,
# 열린 재료어 = ABCADCBD의 회전. 회전 시 letter의 첫/둘째 방문이 바뀌면 배정 부호 뒤집힘. — the open material word = a rotation of ABCADCBD. Under rotation, if a letter's first/second visit swaps, the assignment sign flips.
import itertools
from collections import defaultdict
W="ABCADCBD"
ASSIGNS=[(1,-1,-1,1),(-1,1,1,-1)]
BASEPOS={ch:[i for i,c in enumerate(W) if c==ch] for ch in "ABCD"}

cells=defaultdict(list)
for r in range(8):
    word=W[r:]+W[:r]
    pos={ch:[i for i,c in enumerate(word) if c==ch] for ch in "ABCD"}
    swap={ch: ((BASEPOS[ch][0]-r)%8 > (BASEPOS[ch][1]-r)%8) for ch in "ABCD"}
    for asg0 in ASSIGNS:
        asg={ch:(a if not swap[ch] else -a) for ch,a in zip("ABCD",asg0)}
        for cuts in itertools.combinations_with_replacement(range(9),4):
            pa={ch:(sum(1 for c in cuts if c<=pos[ch][0]), sum(1 for c in cuts if c<=pos[ch][1])) for ch in "ABCD"}
            if any(abs(p1-p2)<2 for p1,p2 in pa.values()): continue
            edges=[]
            for ch in "ABCD":
                p1,p2=pa[ch]
                first,second=(p1,p2) if pos[ch][0]<pos[ch][1] else (p2,p1)
                top,bot=(first,second) if asg[ch]==1 else (second,first)
                edges.append((top,bot))
            indeg={n:0 for n in range(5)}; g={n:set() for n in range(5)}
            for a,b in edges:
                if b not in g[a]: g[a].add(b); indeg[b]+=1
            q=[n for n in range(5) if indeg[n]==0]; seen=0
            while q:
                n=q.pop(); seen+=1
                for m in g[n]:
                    indeg[m]-=1
                    if indeg[m]==0: q.append(m)
            if seen==5:
                pairs=tuple(sorted(tuple(sorted(v)) for v in pa.values()))
                cells[pairs].append((r,cuts,tuple(asg[ch] for ch in "ABCD")))

print(f"4_1 층 전순서를 허용하는 패스-쌍 조합: {len(cells)}")
for pairs,lst in sorted(cells.items()):
    revpairs=tuple(sorted(tuple(sorted((4-a,4-b))) for a,b in pairs))
    tag="자기대칭" if revpairs==pairs else f"거울={revpairs}"
    print(f"  {pairs}  변형 {len(lst)}개  [{tag}]")
