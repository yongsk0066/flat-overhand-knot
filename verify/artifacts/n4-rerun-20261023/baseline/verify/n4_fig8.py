DBG=[]
# n4_fig8.py — ABCADCBD: 그림자 강성(회전계 전수) + 회전계-인식 Kauffman 괄호 → 매듭 동정 — shadow rigidity (enumeration of rotation systems) + rotation-system-aware Kauffman bracket → knot identification
# 스무딩 = 회전계의 인접 끝 짝짓기(2택); A/B 배정·교차부호 규약은 세잎 진리표로 캘리브레이션. — smoothing = pairing adjacent ends in the rotation system (2 choices); the A/B assignment and the crossing-sign convention are calibrated against the trefoil truth table.
import itertools, sympy as sp
A_, t = sp.symbols('A t')
TRE = {sp.expand(-t**-4+t**-3+t**-1), sp.expand(-t**4+t**3+t)}
FIG8 = sp.expand(t**-2-t**-1+1-t+t**2)

def build(code):
    occ={}
    for i,ch in enumerate(code): occ.setdefault(ch,[]).append(i)
    VIS={ch:tuple(occ[ch]) for ch in occ}; n=len(code)
    def rotations(choice):
        rot={}
        for k,ch in enumerate(sorted(VIS)):
            i1,i2=VIS[ch]
            o1,b1=(i1,+1),((i1-1)%n,-1)
            o2,b2=(i2,+1),((i2-1)%n,-1)
            X,Y=(o2,b2) if choice[k]==0 else (b2,o2)
            rot[ch]=[o1,X,b1,Y]
        return rot
    def faces(rot):
        nxt={}
        for cyc in rot.values():
            for k,d in enumerate(cyc): nxt[d]=cyc[(k+1)%4]
        seen=set(); F=[]
        for d0 in [(i,s) for i in range(n) for s in (+1,-1)]:
            if d0 in seen: continue
            L=0; d=d0
            while d not in seen:
                seen.add(d); L+=1; d=nxt[(d[0],-d[1])]
            F.append(L)
        return sorted(F)
    planar=[b for b in itertools.product((0,1),repeat=len(VIS)) if len(faces(rotations(b)))==2+len(VIS)]
    return VIS,n,rotations,planar

def vertex_of(dart,VIS,n):
    i,s=dart; pos=i if s==+1 else (i+1)%n
    return next(ch for ch,(a,b) in VIS.items() if pos in (a,b))

def jones(code,rot,VIS,n,signs,CAL_A,CAL_W,bits_rot=None,dbg=False):
    chs=sorted(VIS); delta=sp.expand(-A_**2-A_**-2)
    # 끝(다트) → 정점, 짝짓기 준비 — end (dart) → vertex, prepare the pairing
    tot=sp.Integer(0)
    for bits in itertools.product((0,1),repeat=len(chs)):
        sm=dict(zip(chs,bits))
        # partner: 각 정점 cyc에서 type0=(0,1)(2,3), type1=(1,2)(3,0) — at each vertex cycle, type0=(0,1)(2,3), type1=(1,2)(3,0)
        partner={}
        for ch in chs:
            cyc=rot[ch]
            if sm[ch]==0: pairs=[(0,1),(2,3)]
            else: pairs=[(1,2),(3,0)]
            for a,b in pairs:
                partner[cyc[a]]=cyc[b]; partner[cyc[b]]=cyc[a]
        # 루프: 상태 = 다트(변의 한쪽 끝, 진행방향); 변 통과 후 도착 끝의 짝으로 전환 — loops: a state = a dart (one end of an edge, with a direction); after traversing the edge, switch to the partner of the arriving end
        seen=set(); loops=0
        for d0 in [(i,+1) for i in range(n)]:
            if d0 in seen: continue
            loops+=1; d=d0
            while d not in seen:
                seen.add(d); seen.add((d[0],-d[1]))
                arrive=(d[0],-d[1])          # 변 d 를 지나 반대쪽 끝에 도착 — traverse edge d and arrive at the opposite end
                d=partner[arrive]
        # 가중: 정점별 A/B — over 가닥 = 첫 통과(s=+1) 또는 둘째(s=-1); cyc에서 over쌍이 — weight: A/B per vertex — the over strand = the first traversal (s=+1) or the second (s=-1); if the over pair in cyc is
        # (cyc0,cyc2)이면 type0이 A (CAL_A=+1 규약) — over가 (cyc1,cyc3)이면 반대. — (cyc0,cyc2) then type0 is A (the CAL_A=+1 convention) — if the over pair is (cyc1,cyc3) it is the opposite.
        if dbg: DBG.append(loops)
        w=sp.Integer(1)
        for k,ch in enumerate(chs):
            over_first = (signs[ch]==+1)
            rbit = 0 if bits_rot is None else bits_rot[k]
            typeA = 0 if ((over_first == (rbit==0)) == (CAL_A==+1)) else 1
            w *= A_ if sm[ch]==typeA else A_**-1
        tot += w*delta**(loops-1)
    wr = CAL_W*sum(signs.values())
    return sp.expand(sp.powsimp(sp.expand((-A_**3)**(-wr)*tot).subs(A_,t**sp.Rational(-1,4))))

# ---- 캘리브레이션: 세잎 진리표 ---- — calibration: the trefoil truth table
VIS3,n3,rots3,pl3=build("ABCABC")
print(f"[검증] ABCABC 평면 회전계 {len(pl3)}개")
DBG=[]
_=jones("ABCABC",rots3(pl3[0]),VIS3,n3,dict(zip("ABC",(1,1,1))),+1,+1,pl3[0],dbg=True)
from collections import Counter as _C
print("[디버그] ABCABC 상태 루프 멀티셋:", dict(_C(DBG)), " (기대 {1:3, 2:4, 3:1})")
cal=None
for CA in (+1,-1):
    for CW in (+1,-1):
        ok=True
        for s3 in itertools.product((+1,-1),repeat=3):
            V=jones("ABCABC",rots3(pl3[0]),VIS3,n3,dict(zip("ABC",s3)),CA,CW,pl3[0])
            mixed=len(set(s3))>1
            if mixed!=(V==1) or (not mixed and V not in TRE): ok=False; break
        if ok: cal=(CA,CW); break
    if cal: break
print(f"[캘리브레이션] (A규약, writhe규약) = {cal}"); assert cal
CA,CW=cal

# ---- ABCADCBD ----
VIS4,n4,rots4,pl4=build("ABCADCBD")
print(f"\n[ABCADCBD] 평면 회전계 {len(pl4)}개; 면벡터 예: ", end="")
_,_,r,_=build("ABCADCBD")
print("(오일러 V-E+F=2 통과)")
from collections import Counter
cnt=Counter(); hits=[]
rot=rots4(pl4[0])
for s4 in itertools.product((+1,-1),repeat=4):
    V=jones("ABCADCBD",rot,VIS4,n4,dict(zip("ABCD",s4)),CA,CW,pl4[0])
    if V==1: cnt['unknot(V=1)']+=1
    elif V==FIG8: cnt['figure8']+=1; hits.append(s4)
    elif V in TRE: cnt['trefoil']+=1
    else: cnt['other']+=1
print(f"[ABCADCBD] 16 상하 배정: {dict(cnt)}")
print(f"[ABCADCBD] V(4_1) 배정: {hits}")
if len(pl4)>=2:
    V2=jones("ABCADCBD",rots4(pl4[1]),VIS4,n4,dict(zip("ABCD",hits[0] if hits else (1,1,1,1))),CA,CW,pl4[1])
    print(f"[거울 회전계 교차확인] 같은 배정의 V = {V2}")
