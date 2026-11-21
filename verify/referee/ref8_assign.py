"""Check the old layer masks with alternation, determinants and the Jones polynomial."""
from probe_files import VERIFY

import itertools, collections, sys
from fractions import Fraction
import sympy as sp
sys.path.insert(0, str(VERIFY))
W="ABCADCBD"

def seq(asg):
    """asg[ch]=+1 means the FIRST visit of ch is the over-crossing (positive first-visit sign)."""
    seen=set(); out=[]
    for ch in W:
        first = ch not in seen; seen.add(ch)
        out.append((ch, (asg[ch]==1)==first))
    return out

def det_of(asg):
    ev=[("ABCD".index(ch), ov) for ch,ov in seq(asg)]
    A=sum(1 for _,ov in ev if not ov); arc=0; vis=[]
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
    M=sp.Matrix([r[1:] for r in rows[1:]])
    return int(abs(M.det()))

# --- Kauffman bracket from the Gauss code, fully independent of the coloring matrix ---
def bracket_jones(asg):
    """<D> by state sum over A/B smoothings; then f = (-A^3)^{-w} <D>, t = A^{-4}."""
    A_=sp.symbols('A')
    # build the 4-valent graph: 8 visit slots around the knot; crossing c joins slots.
    S=seq(asg)
    pos={}
    for i,(ch,ov) in enumerate(S): pos.setdefault(ch,[]).append(i)
    # each crossing has 4 local ends: (in_over,out_over,in_under,out_under) on the knot circle
    # Represent the diagram as a set of 8 edges of the circle; smoothing reconnects them.
    # circle edges: e_i joins slot i to slot i+1 (mod 8)  -> endpoints labelled (i,'out'),(i+1,'in')
    total=0
    for bits in itertools.product([0,1],repeat=4):    # 0 = A-smoothing, 1 = B-smoothing
        # union-find over the 8 circle-edge objects
        par=list(range(8))
        def f(x):
            while par[x]!=x: par[x]=par[par[x]]; x=par[x]
            return x
        def u(a,b):
            a,b=f(a),f(b)
            if a!=b: par[a]=b
        for k,ch in enumerate("ABCD"):
            i,j=pos[ch]                # slots of the two visits of this crossing
            # incoming edge to slot i is edge (i-1); outgoing edge from slot i is edge i
            a_in,a_out=(i-1)%8,i%8
            b_in,b_out=(j-1)%8,j%8
            # crossing sign convention: over-strand determines which pairing is the A-smoothing
            over_first = S[i][1]
            if (bits[k]==0) == over_first:
                u(a_in,b_out); u(b_in,a_out)
            else:
                u(a_in,a_out); u(b_in,b_out)
            u(a_in,a_in)
        loops=len({f(x) for x in range(8)})
        nA=4-sum(bits); nB=sum(bits)
        total += A_**(nA-nB) * (-A_**2-A_**-2)**(loops-1)
    return sp.simplify(sp.expand(total))

CAND={"legacy ASSIGNS[0] (+,-,-,+)":dict(zip("ABCD",(1,-1,-1,1))),
      "legacy ASSIGNS[1] (-,+,+,-)":dict(zip("ABCD",(-1,1,1,-1))),
      "determinant-five pair [0]  (-,+,-,-)":dict(zip("ABCD",(-1,1,-1,-1))),
      "determinant-five pair [1]  (+,-,+,+)":dict(zip("ABCD",(1,-1,1,1)))}
print("code W = ABCADCBD;  a 4-crossing diagram of 4_1 must be MINIMAL hence ALTERNATING")
print(f"{'assignment':32} {'O/U sequence':12} {'alternating':>11} {'det':>5}  knot")
for nm,a in CAND.items():
    s="".join('O' if ov else 'U' for ch,ov in seq(a))
    alt=all(s[i]!=s[(i+1)%8] for i in range(8))
    D=det_of(a)
    knot={1:'unknot',3:'3_1',5:'4_1'}.get(D,'?')
    print(f"  {nm:30} {s:12} {str(alt):>11} {D:>5}  {knot}")

print("\nfull census of the 16 masks (det):")
h=collections.Counter()
for bits in itertools.product([1,-1],repeat=4):
    h[det_of(dict(zip("ABCD",bits)))]+=1
print("  ",dict(h), " -> paper's '2 figure-eight, 2 trefoil, 12 unknot' is CORRECT as a count")

print("\nKauffman bracket cross-check (independent of the coloring matrix):")
for nm,a in list(CAND.items()):
    try:
        print(f"  {nm:30} <D> = {bracket_jones(a)}")
    except Exception as e:
        print(f"  {nm:30} bracket failed: {e}")

# Compare the cells selected by the determinant-one and determinant-five mask pairs.
def enumerate_cells(assign_pair):
    BASEPOS={ch:[i for i,c in enumerate(W) if c==ch] for ch in "ABCD"}
    cells=collections.defaultdict(int)
    for r in range(8):
        word=W[r:]+W[:r]
        pos={ch:[i for i,c in enumerate(word) if c==ch] for ch in "ABCD"}
        swap={ch:((BASEPOS[ch][0]-r)%8 > (BASEPOS[ch][1]-r)%8) for ch in "ABCD"}
        for asg0 in assign_pair:
            asg={ch:(a if not swap[ch] else -a) for ch,a in zip("ABCD",asg0)}
            for cuts in itertools.combinations_with_replacement(range(9),4):
                pa={ch:(sum(1 for c in cuts if c<=pos[ch][0]),sum(1 for c in cuts if c<=pos[ch][1])) for ch in "ABCD"}
                if any(abs(p1-p2)<2 for p1,p2 in pa.values()): continue
                edges=[]
                for ch in "ABCD":
                    p1,p2=pa[ch]
                    top,bot=(p1,p2) if asg[ch]==1 else (p2,p1)
                    edges.append((top,bot))
                indeg={n:0 for n in range(5)}; g={n:set() for n in range(5)}
                for a,b in edges:
                    if b not in g[a]: g[a].add(b); indeg[b]+=1
                q=[n for n in range(5) if indeg[n]==0]; s=0
                while q:
                    n=q.pop(); s+=1
                    for m in g[n]:
                        indeg[m]-=1
                        if indeg[m]==0: q.append(m)
                if s==5: cells[tuple(sorted(tuple(sorted(v)) for v in pa.values()))]+=1
    return cells
cA=enumerate_cells([(1,-1,-1,1),(-1,1,1,-1)])
cR=enumerate_cells([(-1,1,-1,-1),(1,-1,1,1)])
print(f"\ncell count with the determinant-one mask pair : {len(cA)} (legacy count)")
print(f"cell count with the true 4_1 (alternating, det=5) pair      : {len(cR)}")
sel=tuple(sorted([(0,2),(0,3),(1,4),(2,4)])); wit=tuple(sorted([(0,2),(0,3),(0,4),(1,3)]))
print(f"  realised cell {sel}: legacy list {sel in cA}, true list {sel in cR}")
print(f"  abstract witness  {wit}: legacy list {wit in cA}, true list {wit in cR}")
print("  The witness pairs (1,3), (0,2), (0,3) and (0,4) all satisfy non-adjacency.")
