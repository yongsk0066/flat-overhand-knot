# ref9_confirm.py — REFEREE: settle the ASSIGNS question on REAL planar diagrams,
# and cross-validate the abstract cell list against the geometric harvest.
import sys, random, collections, itertools
import sympy as sp
sys.path.insert(0,'/Users/yongseok/dev/pentagon/verify')
from n4_tighten2 import (config, crossings, code_str, canon, events, ev_over_of,
                         fox5_events, orders_making_41, linear_extensions, axiom_exact_ok, TARGET_CANON)

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
    return int(abs(sp.Matrix([r[1:] for r in rows[1:]]).det()))

print("="*76)
print("A. On REAL planar diagrams: are the det=5 (=4_1) layerings exactly the ALTERNATING ones?")
print("="*76)
rnd=random.Random(2718281828)
pat=collections.Counter(); alt_stats=collections.Counter(); n=0
for _ in range(300000):
    if rnd.random()<0.5:
        rho=rnd.choice([1,-1]); sg=[rho]*4
    else: sg=[rnd.choice([1,-1]) for _ in range(4)]
    th=tuple([sg[k]*rnd.uniform(0.3,3.1) for k in range(4)]+[rnd.uniform(0.3,4.0) for _ in range(3)])
    c=config(th)
    if c is None: continue
    X=crossings(c)
    if len(X)!=4: continue
    cd=code_str(X)
    if canon(cd)!=TARGET_CANON: continue
    n+=1
    for mask in range(16):
        pairs=[]
        for cid,((i,j),_,_) in enumerate(X):
            pairs.append((i,j) if (mask>>cid)&1 else (j,i))
        if not linear_extensions(pairs): continue
        evo=[(cid,(side=='i')==((mask>>cid)&1==1)) for _,_,cid,side in events(X)]
        D=det_events(evo)
        s="".join('O' if ov else 'U' for _,ov in evo)
        alt=all(s[k]!=s[(k+1)%8] for k in range(8))
        alt_stats[(D,alt)]+=1
        if D==5:
            # express the assignment in the letters of the canonical relabelling of cd
            lab={}; word=""
            for ch in cd:
                if ch not in lab: lab[ch]=chr(ord('A')+len(lab))
                word+=lab[ch]
            seen=set(); asg={}
            for k,ch in enumerate(word):
                if ch not in seen:
                    seen.add(ch); asg[ch]=1 if evo[k][1] else -1
            pat[tuple(asg[ch] for ch in "ABCD")]+=1
print(f"  canon-passing 4-crossing diagrams scanned: {n}")
print(f"  (determinant, alternating?) histogram over acyclic layerings: {dict(alt_stats)}")
print(f"  -> every det=5 layering is alternating, and no other det is: "
      f"{all((k[1] if k[0]==5 else not k[1]) for k in alt_stats)}")
print(f"  assignment patterns (in the relabelled word, +1 = first visit is OVER) of the 4_1 layerings:")
for k,v in pat.most_common(): print(f"     {k}: {v}")
print(f"  n4_cells2.py hardcodes ASSIGNS = [(1,-1,-1,1), (-1,1,1,-1)]  -> observed: "
      f"{[k for k in pat if k in [(1,-1,-1,1),(-1,1,1,-1)]] or 'NEVER'}")
print("  NOTE: the relabelled word of a realised diagram is a ROTATION of ABCADCBD, so the")
print("  pattern is stated per-rotation; the decisive invariant is ALTERNATION, which is")
print("  rotation-independent, and the authors' pair is non-alternating in every rotation.")
for r in range(8):
    w="ABCADCBD"[r:]+"ABCADCBD"[:r]
    for asg0 in [(1,-1,-1,1),(-1,1,1,-1)]:
        seen=set(); s=""
        for ch in w:
            first=ch not in seen; seen.add(ch)
            k="ABCD".index(ch)
            s+="O" if ((asg0[k]==1)==first) else "U"
        assert not all(s[i]!=s[(i+1)%8] for i in range(8))
print("  verified: for all 8 rotations, both authors' ASSIGNS give NON-alternating sequences.")

print()
print("="*76)
print("B. abstract cell list vs the geometric harvest  (consistency test)")
print("="*76)
def enum_cells(pair):
    W="ABCADCBD"; BASE={ch:[i for i,c in enumerate(W) if c==ch] for ch in "ABCD"}
    out=collections.defaultdict(int)
    for r in range(8):
        word=W[r:]+W[:r]
        pos={ch:[i for i,c in enumerate(word) if c==ch] for ch in "ABCD"}
        swap={ch:((BASE[ch][0]-r)%8>(BASE[ch][1]-r)%8) for ch in "ABCD"}
        for a0 in pair:
            asg={ch:(a if not swap[ch] else -a) for ch,a in zip("ABCD",a0)}
            for cuts in itertools.combinations_with_replacement(range(9),4):
                pa={ch:(sum(1 for c in cuts if c<=pos[ch][0]),sum(1 for c in cuts if c<=pos[ch][1])) for ch in "ABCD"}
                if any(abs(x-y)<2 for x,y in pa.values()): continue
                edges=[((p1,p2) if asg[ch]==1 else (p2,p1)) for ch in "ABCD" for p1,p2 in [pa[ch]]]
                indeg={k:0 for k in range(5)}; g={k:set() for k in range(5)}
                for a,b in edges:
                    if b not in g[a]: g[a].add(b); indeg[b]+=1
                q=[k for k in range(5) if indeg[k]==0]; s=0
                while q:
                    k=q.pop(); s+=1
                    for m in g[k]:
                        indeg[m]-=1
                        if indeg[m]==0: q.append(m)
                if s==5: out[tuple(sorted(tuple(sorted(v)) for v in pa.values()))]+=1
    return out
cA=enum_cells([(1,-1,-1,1),(-1,1,1,-1)]); cR=enum_cells([(-1,1,-1,-1),(1,-1,1,1)])
distinctA=[k for k in cA if len(set(k))==4]; distinctR=[k for k in cR if len(set(k))==4]
print(f"  authors' (wrong) pair: {len(cA)} cells, of which {len(distinctA)} have 4 DISTINCT pass pairs")
print(f"  true 4_1 pair        : {len(cR)} cells, of which {len(distinctR)} have 4 DISTINCT pass pairs")
print("  (two straight passes meet at most once, so only 4-distinct-pair cells are realisable;")
print("   neither script applies that second geometric necessary condition)")
print(f"  the realisable-shape list under the TRUE pair: {sorted(distinctR)}")
HARVEST={((0,2),(0,3),(1,4),(2,4)):2426, ((0,3),(0,4),(1,4),(2,4)):475, ((0,2),(0,3),(0,4),(1,4)):444,
         ((0,2),(0,3),(0,4),(2,4)):0, ((0,2),(0,4),(1,4),(2,4)):0}
print("\n  doc-18 harvest table ('4_1 layer exists' column) vs the two abstract lists:")
print(f"  {'pairing':38}{'harvest':>9}{'authors 72':>12}{'referee 17':>12}")
for k,v in HARVEST.items():
    kk=tuple(sorted(k))
    print(f"  {str(k):38}{v:>9}{str(kk in cA):>12}{str(kk in cR):>12}")
ok=all((v>0)==(tuple(sorted(k)) in cR) for k,v in HARVEST.items())
okA=all((v>0)==(tuple(sorted(k)) in cA) for k,v in HARVEST.items())
print(f"\n  harvest agrees with the TRUE abstract list : {ok}")
print(f"  harvest agrees with the AUTHORS' 72 list   : {okA}")
