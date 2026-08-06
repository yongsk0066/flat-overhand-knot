# n4_witness.py — B_{4_1} 1단계: ABCADCBD의 증인(조각쌍+통과순서) 추출, 층 순서 무모순 판정 — B_{4_1} stage 1: extract witnesses for ABCADCBD (piece pairs + pass order) and test the layer order for consistency
import itertools
PAIRS=[(0,2),(0,3),(0,4),(1,3),(1,4),(2,4)]
def parity_ok(code):
    for ch in set(code):
        i,j=code.index(ch),code.rindex(ch)
        if (j-i-1)%2==1: return False
    return True
def canon(code):
    n=len(code); best=None
    for r in range(n):
        rot=code[r:]+code[:r]; lab={}; out=[]
        for ch in rot:
            if ch not in lab: lab[ch]=chr(ord('A')+len(lab))
            out.append(lab[ch])
        cand="".join(out)
        best=cand if best is None or cand<best else best
    return best
wit=[]
for sub in itertools.combinations(range(6),4):
    pairs=[PAIRS[i] for i in sub]
    marks={p:[] for p in range(5)}
    for cid,(p,q) in enumerate(pairs):
        marks[p].append(cid); marks[q].append(cid)
    for combo in itertools.product(*[list(itertools.permutations(marks[p])) for p in range(5)]):
        seq=[]; 
        for p in range(5): seq.extend(combo[p])
        code="".join(chr(ord('A')+c) for c in seq)
        if parity_ok(code) and canon(code)=="ABCADCBD":
            wit.append((pairs,combo,code))
print(f"증인 배치 수: {len(wit)}")
pairs,combo,code=wit[0]
print("대표 증인: 조각쌍 =", pairs, " 코드 =", code)
print("  조각별 통과열:", {p:[chr(ord('A')+c) for c in combo[p]] for p in range(5) if combo[p]})
# 층 순서: V(4_1) 배정 (+,-,-,+) — n4_fig8 규약: sign=+1 ⟺ '첫 통과 조각이 over' — layer order: V(4_1) assignment (+,-,-,+); n4_fig8 convention: sign=+1 ⟺ 'the first-passing piece is the over one'
for signs in [(+1,-1,-1,+1),(-1,+1,+1,-1)]:
    rel=[]
    for cid,(p,q) in enumerate(pairs):
        # 첫 통과가 속한 조각 결정 — determine which piece the first pass belongs to
        ch=chr(ord('A')+cid); i1=code.index(ch)
        cum=0; first_piece=None
        for pc in range(5):
            if i1 < cum+len(combo[pc]): first_piece=pc; break
            cum+=len(combo[pc])
        other=q if first_piece==p else p
        over,under=(first_piece,other) if signs[cid]==+1 else (other,first_piece)
        rel.append((over,under))
    # 반순서 무모순(위상 정렬) — consistency of the partial order (topological sort)
    import collections
    g=collections.defaultdict(set)
    for o,u in rel: g[o].add(u)
    order=[]; nodes=set(range(5)); ok=True
    while nodes:
        free=[v for v in nodes if not any(v in g[w] for w in nodes if w!=v)]
        if not free: ok=False; break
        v=free[0]; order.append(v); nodes.discard(v)
    print(f"배정 {signs}: 교차별 (over,under) = {rel} → 전역 층 순서 {'존재: '+str(order) if ok else '모순!'}")

# ---- 2단계: 중심선 실현 탐색 (tau1..4, L1..3) — dim 7 = 2n-1 ---- — stage 2: search for a centerline realization (tau1..4, L1..3)
import math, random
def diagram(tau,L,ray=60.0):
    d=[ (math.cos(sum(tau[:j])), math.sin(sum(tau[:j]))) for j in range(5) ]
    M=[(0.0,0.0)]
    for i in range(3): M.append((M[-1][0]+L[i]*d[i+1][0], M[-1][1]+L[i]*d[i+1][1]))
    # 조각: P0 레이(기점 M[0], 방향 -d0), P1..P3 선분, P4 레이(M[3], d4) — pieces: P0 is a ray (base M[0], direction -d0), P1..P3 are segments, P4 is a ray (M[3], d4)
    segs={0:(M[0],(-d[0][0],-d[0][1]),ray),1:(M[0],d[1],L[0]),2:(M[1],d[2],L[1]),3:(M[2],d[3],L[2]),4:(M[3],d[4],ray)}
    X=[]
    for i in range(5):
        for j in range(i+2,5):
            if (i,j)==(0,1) or (i,j)==(3,4): continue
            P1,dd1,l1=segs[i]; P2,dd2,l2=segs[j]
            det=dd1[0]*dd2[1]-dd1[1]*dd2[0]
            if abs(det)<1e-12: continue
            rx,ry=P2[0]-P1[0],P2[1]-P1[1]
            s=(rx*dd2[1]-ry*dd2[0])/det; u=(rx*dd1[1]-ry*dd1[0])/det
            if 1e-9<s<l1-1e-9 and 1e-9<u<l2-1e-9:
                X.append(((i,j),s,u))
    return X
TARGET={(0,2),(0,3),(0,4),(1,3)}
rnd=random.Random(2026)
hits=[]
for _ in range(400000):
    rho=rnd.choice([1,-1])
    tau=[rho*rnd.uniform(1.7,3.05) for _ in range(3)]+[rnd.choice([1,-1])*rnd.uniform(0.4,3.05)]
    L=[rnd.uniform(0.3,3.0) for _ in range(3)]
    X=diagram(tau,L)
    if {p for p,_,_ in X}==TARGET and len(X)==4:
        # P0 통과 순서: s 클수록 재료상 이른 것 → [B,A,C] = [(0,3),(0,2),(0,4)] — pass order along P0: larger s means earlier in the material ⟹ [B,A,C] = [(0,3),(0,2),(0,4)]
        p0=[(s,p) for p,s,u in X if p[0]==0]; p0.sort(reverse=True)
        ordP0=[p for _,p in p0]
        p3=[(u if p==(0,3) else s, p) for p,s,u in X if 3 in p and p!=(0,3) or p==(0,3)]
        if ordP0==[(0,3),(0,2),(0,4)]:
            hits.append((tau,L,X)); 
            if len(hits)>=3: break
print(f"\n[실현 탐색] 표적 교차쌍+P0 순서 일치 배치: {len(hits)}개 발견 (40만 표본)")
for tau,L,X in hits[:2]:
    print("  τ=",[round(t,3) for t in tau]," L=",[round(l,3) for l in L])
    print("   교차:",[(p,round(s,2),round(u,2)) for p,s,u in X])
