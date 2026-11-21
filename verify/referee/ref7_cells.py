"""Enumerate abstract cells for each cut of the cyclic crossing word."""
from probe_files import VERIFY

import itertools, collections, sys
from fractions import Fraction
sys.path.insert(0, str(VERIFY))

def det_events(ev):
    A=sum(1 for _,ov in ev if not ov)
    if A==0: return 1
    arc=0; vis=[]
    for cid,ov in ev:
        vis.append((cid,ov,arc))
        if not ov: arc+=1
    oa={}; und={}
    for cid,ov,a0 in vis:
        if ov: oa[cid]=a0%A
        else: und.setdefault(cid,[]).append((a0%A,(a0+1)%A))
    rows=[]
    for cid,o in oa.items():
        i_,o_=und[cid][0]; r=[0]*A; r[o]+=2; r[i_]-=1; r[o_]-=1; rows.append(r)
    m=[[Fraction(x) for x in r[1:]] for r in rows[1:]]; n=len(m)
    if n==0: return 1
    det=Fraction(1)
    for c in range(n):
        p=next((r for r in range(c,n) if m[r][c]!=0),None)
        if p is None: return 0
        if p!=c: m[c],m[p]=m[p],m[c]; det=-det
        det*=m[c][c]; inv=m[c][c]
        for r in range(c+1,n):
            f=m[r][c]/inv
            if f: m[r]=[m[r][k]-f*m[c][k] for k in range(n)]
    return abs(int(det))

W="ABCADCBD"
print("A. which of the 16 over/under masks of the CLOSED code ABCADCBD give 4_1?")
res=collections.Counter()
good=[]
for mask in range(16):
    ev=[]
    seen=set()
    for ch in W:
        k="ABCD".index(ch)
        first = ch not in seen; seen.add(ch)
        over = (((mask>>k)&1)==1) == first     # bit set <=> FIRST visit is the over-visit
        ev.append((k,over))
    D=det_events(ev); res[D]+=1
    if D==5: good.append(tuple(1 if (mask>>k)&1 else -1 for k in range(4)))
print("   determinant histogram over the 16 masks:", dict(res),
      "  (paper: 2 give 4_1, 2 give 3_1, 12 unknot)")
print("   the two 4_1 masks (as +-1 per letter, +1 = first visit on top):", good)
print("   legacy hardcoded ASSIGNS = [(1,-1,-1,1), (-1,1,1,-1)]  -> MATCH:",
      sorted(good)==sorted([(1,-1,-1,1),(-1,1,1,-1)]))

print("\nB. independent cell enumeration (rotation + cut placement + acyclicity)")
cells=collections.defaultdict(list); nvar=0
for r in range(8):
    word=W[r:]+W[:r]
    for mask in range(16):
        ev=[]; seen=set(); firstflag=[]
        for ch in word:
            k="ABCD".index(ch); first = ch not in seen; seen.add(ch)
            ev.append((k, (((mask>>k)&1)==1)==first)); firstflag.append(first)
        if det_events(ev)!=5: continue          # <- 4_1 decided by the determinant, per rotation
        for cuts in itertools.combinations_with_replacement(range(9),4):
            passof=[sum(1 for c in cuts if c<=p) for p in range(8)]
            byletter=collections.defaultdict(list)
            for p,ch in enumerate(word): byletter[ch].append(p)
            ok=True; edges=[]; pairs=[]
            for ch in "ABCD":
                p1,p2=byletter[ch]
                a,b=passof[p1],passof[p2]
                if abs(a-b)<2: ok=False; break   # straight passes cross only NON-ADJACENT passes
                over_is_first = ev[p1][1]        # is the first visit (position p1) the over one?
                edges.append((a,b) if over_is_first else (b,a))
                pairs.append(tuple(sorted((a,b))))
            if not ok: continue
            indeg={n:0 for n in range(5)}; g={n:set() for n in range(5)}
            for a,b in edges:
                if b not in g[a]: g[a].add(b); indeg[b]+=1
            q=[n for n in range(5) if indeg[n]==0]; s=0
            while q:
                n=q.pop(); s+=1
                for m2 in g[n]:
                    indeg[m2]-=1
                    if indeg[m2]==0: q.append(m2)
            if s==5:
                cells[tuple(sorted(pairs))].append((r,mask,cuts)); nvar+=1
print(f"   distinct pass-pair combinations admitting a 4_1 total layer order: {len(cells)}"
      f"   [variants {nvar}]")
sel=tuple(sorted([(0,2),(0,3),(1,4),(2,4)]))
print(f"   the realised cell {sel} present: {sel in cells}  (variants {len(cells.get(sel,[]))})")
wit=tuple(sorted([(0,2),(0,3),(0,4),(1,3)]))
print(f"   abstract witness cell {wit} present: {wit in cells}"
      f" (non-adjacency test)")
# is that witness excluded only by non-adjacency, or absent for another reason?
print("     (1,3) is non-adjacent, so the exclusion cannot be the non-adjacency rule alone;")
print("     check whether it survives WITHOUT the acyclicity test:")
cnt=0
for r in range(8):
    word=W[r:]+W[:r]
    for mask in range(16):
        ev=[]; seen=set()
        for ch in word:
            k="ABCD".index(ch); first=ch not in seen; seen.add(ch)
            ev.append((k,(((mask>>k)&1)==1)==first))
        if det_events(ev)!=5: continue
        for cuts in itertools.combinations_with_replacement(range(9),4):
            passof=[sum(1 for c in cuts if c<=p) for p in range(8)]
            byl=collections.defaultdict(list)
            for p,ch in enumerate(word): byl[ch].append(p)
            prs=[]; ok=True
            for ch in "ABCD":
                a,b=passof[byl[ch][0]],passof[byl[ch][1]]
                if abs(a-b)<2: ok=False; break
                prs.append(tuple(sorted((a,b))))
            if ok and tuple(sorted(prs))==wit: cnt+=1
print(f"     witness cell reachable ignoring acyclicity: {cnt} combinatorial variants")

print("\nC. proof that a straight pass can cross only NON-ADJACENT passes")
print("   P_i and P_{i+1} are straight and share exactly the fold point M_i.  Two distinct")
print("   straight sets sharing a point meet only there unless collinear; collinear needs")
print("   d_i = +-d_{i+1}, i.e. tau_{i+1} in {0, pi}, and config() returns None on")
print("   |d_i+d_{i+1}|<1e-9 (tau=pi) while tau=0 is not a crease.  The crossing test uses")
print("   the strict interior (1e-9 < s < l-1e-9), so the shared endpoint is never counted.")
print("   Likewise P_0,P_1 both emanate from M_0 (directions -d_0, d_1) and P_3,P_4 from M_3.")
print("   => the restriction  for i in range(5): for j in range(i+2,5)  is EXACT, not heuristic.")
