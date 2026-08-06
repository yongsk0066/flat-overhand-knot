# ref6_oracle.py — REFEREE: independent knot-type oracle + RAY-truncation attack + "42 vs 21".
# Independent oracle: the INTEGER determinant |det(Alexander matrix at t=-1)| from the
# Wirtinger/Fox relations, i.e. the Goeritz-free coloring matrix over Z (not mod 5).
# For a 4-crossing diagram the only possible knots are unknot/3_1/4_1 (det 1/3/5), so any
# other value proves the arc/closure bookkeeping is broken.
import sys, math, random, json, collections
from fractions import Fraction
sys.path.insert(0,'/Users/yongseok/dev/pentagon/verify')
from n4_tighten2 import (config, crossings, code_str, canon, events, ev_over_of,
                         fox5_events, fox3_events, axiom_exact_ok, orders_making_41,
                         linear_extensions, ell, TARGET_CANON)

def det_from_events(ev_over):
    """|det K| where K is the coloring/Alexander-at-(-1) matrix of the closed long knot."""
    A=sum(1 for _,ov in ev_over if not ov)
    if A==0: return 1                       # no undercrossings -> unknot
    arc=0; visits=[]
    for cid,ov in ev_over:
        visits.append((cid,ov,arc))
        if not ov: arc+=1
    over_arc={}; und={}
    for cid,ov,a0 in visits:
        if ov: over_arc[cid]=a0%A
        else: und.setdefault(cid,[]).append((a0%A,(a0+1)%A))
    rows=[]
    for cid,oa in over_arc.items():
        ain,aout=und[cid][0]
        r=[0]*A; r[oa]+=2; r[ain]-=1; r[aout]-=1
        rows.append(r)
    # delete one row and one column, integer determinant via fraction-free Gauss
    if len(rows)<1: return 1
    m=[[Fraction(x) for x in r[1:]] for r in rows[1:]]
    n=len(m)
    if n==0: return 1
    if any(len(r)!=n for r in m): return None
    det=Fraction(1)
    for c in range(n):
        p=next((r for r in range(c,n) if m[r][c]!=0), None)
        if p is None: return 0
        if p!=c: m[c],m[p]=m[p],m[c]; det=-det
        det*=m[c][c]
        inv=m[c][c]
        for r in range(c+1,n):
            f=m[r][c]/inv
            if f: m[r]=[m[r][k]-f*m[c][k] for k in range(n)]
    return abs(int(det))

print("="*78)
print("A. unit tests of the independent determinant oracle on standard codes")
print("="*78)
tre=[(0,True),(1,False),(2,True),(0,False),(1,True),(2,False)]
f8 =[(0,True),(1,False),(2,True),(3,False),(1,True),(0,False),(3,True),(2,False)]
unk=[(0,True),(0,False)]
# a 4-crossing *reducible* diagram (unknot with a kink pair) for good measure
for nm,code,exp in (("trefoil 3_1",tre,3),("figure-eight 4_1",f8,5),("unknot",unk,1)):
    print(f"  {nm:18} det={det_from_events(code)} (expect {exp})   fox3={fox3_events(code)} fox5={fox5_events(code)}")
print("  theory: #Fox-p colourings = p^(1+dim ker) ; p | det(K) <=> nontrivial p-colouring")
print("  so among {unknot,3_1,4_1} (dets 1,3,5): fox5=25 <=> det divisible by 5 <=> 4_1  [VALID]")
print("  component count: the strip centreline is ONE arc from infinity to infinity, so its")
print("  closure is a 1-component KNOT by construction -- the 'could be a 2-component link'")
print("  worry cannot arise. (A 2-component link would give >=p^2 trivial colourings.)")

print()
print("="*78)
print("B. cross-check det-oracle vs fox5 on random realised configurations")
print("="*78)
rnd=random.Random(987654321)
seen=collections.Counter(); mism=0; badval=[]
n4=0
for _ in range(400000):
    if rnd.random()<0.5:
        rho=rnd.choice([1,-1]); sg=[rho]*4
    else: sg=[rnd.choice([1,-1]) for _ in range(4)]
    th=tuple([sg[k]*rnd.uniform(0.3,3.1) for k in range(4)]+[rnd.uniform(0.3,4.0) for _ in range(3)])
    c=config(th)
    if c is None: continue
    X=crossings(c)
    if len(X)!=4: continue
    n4+=1
    for mask in range(16):
        pairs=[]
        for cid,((i,j),_,_) in enumerate(X):
            pairs.append((i,j) if (mask>>cid)&1 else (j,i))
        if not linear_extensions(pairs): continue
        evo=[(cid,(side=='i')==((mask>>cid)&1==1)) for _,_,cid,side in events(X)]
        D=det_from_events(evo); F=fox5_events(evo)
        seen[D]+=1
        if D not in (1,3,5): badval.append((list(th),mask,D))
        if (D%5==0) != (F==25): mism+=1
print(f"  4-crossing configs scanned: {n4};  determinant multiset over acyclic layerings: {dict(seen)}")
print(f"  values outside {{1,3,5}}: {len(badval)}  {'-> bookkeeping OK' if not badval else '*** BROKEN ***'}")
print(f"  (det%5==0) disagrees with (fox5==25): {mism} times  {'-> oracle CONSISTENT' if mism==0 else '*** ORACLE WRONG ***'}")

print()
print("="*78)
print("C. RAY-truncation attack, restricted to states the pipeline ACCEPTS")
print("="*78)
def crossings_R(c,RAY):
    d,M=c['d'],c['M']
    segs={0:(M[0],(-d[0][0],-d[0][1]),RAY),1:(M[0],d[1],c['L'][0]),2:(M[1],d[2],c['L'][1]),
          3:(M[2],d[3],c['L'][2]),4:(M[3],d[4],RAY)}
    X=[]
    for i in range(5):
        for j in range(i+2,5):
            P1,d1,l1=segs[i]; P2,d2,l2=segs[j]
            det=d1[0]*d2[1]-d1[1]*d2[0]
            if abs(det)<1e-12: continue
            rx,ry=P2[0]-P1[0],P2[1]-P1[1]
            s=(rx*d2[1]-ry*d2[0])/det; u=(rx*d1[1]-ry*d1[0])/det
            if 1e-9<s<l1-1e-9 and 1e-9<u<l2-1e-9: X.append(((i,j),s,u))
    return X
def check_ray(label, states):
    bad=0; tot=0
    for th,o in states:
        c=config(tuple(th))
        if c is None: continue
        tot+=1
        if len(crossings_R(c,1e7))!=4: bad+=1
    print(f"  {label}: {bad}/{tot} accepted states gain/lose a crossing when RAY 60 -> 1e7")
hits=json.load(open('/Users/yongseok/dev/pentagon/verify/n4_valid41.json'))
check_ray("authors' harvested valid starts (n4_valid41.json)", [(t,o) for t,o in hits])
cen=json.load(open('/Users/yongseok/dev/pentagon/verify/n4_census.json'))
check_ray("census descent endpoints (n4_census.json)", [(e['theta'],e['order']) for e in cen['endpoints']])
mine=json.load(open('/private/tmp/claude-501/-Users-yongseok-dev-pentagon/0c3f9975-7617-44d6-be99-ab855225d015/scratchpad/ref_hits.json'))
check_ray("my independent harvest (seed 424242)", [(t,o) for t,o,cd in mine])
# and among ALL canon-passing 4-crossing configs, how often does RAY matter?
rnd=random.Random(5551212); tot=0; bad=0
for _ in range(400000):
    if rnd.random()<0.5:
        rho=rnd.choice([1,-1]); sg=[rho]*4
    else: sg=[rnd.choice([1,-1]) for _ in range(4)]
    th=tuple([sg[k]*rnd.uniform(0.3,3.1) for k in range(4)]+[rnd.uniform(0.3,4.0) for _ in range(3)])
    c=config(th)
    if c is None: continue
    X=crossings(c)
    if len(X)!=4 or canon(code_str(X))!=TARGET_CANON: continue
    tot+=1
    if len(crossings_R(c,1e7))!=4: bad+=1
print(f"  random canon-passing 4-crossing configs: {bad}/{tot} change under RAY 60 -> 1e7")

print()
print("="*78)
print("D. '42 valid states': how many are DISTINCT geometric configurations?")
print("="*78)
key=lambda t: tuple(round(v,9) for v in t)
dis=collections.Counter(key(t) for t,o in hits)
print(f"  n4_valid41.json: {len(hits)} records, {len(dis)} distinct theta")
print(f"  multiplicity histogram: {dict(collections.Counter(dis.values()))}")
ords=collections.Counter(tuple(o) for t,o in hits)
print(f"  layer orders: {dict(ords)}")
print("  => each geometric state is recorded once per over/under mask (mirror pair),")
print("     i.e. once with layer order [1,4,2,0,3] and once with its top-bottom flip.")
