# n4_exact_valid.py — θ*_{4_1}의 '전' 후보 접촉 재고를 Q(√15)에서 정확 판정 [리뷰 05 #5] — exactly decide the 'entire' candidate contact inventory of θ*_{4_1} in Q(√15) [review 05 #5]
# n4_walls.all_margins의 재고(끝점-측면 a, 끝점-캡 b, 밴드코너-현 d, 교차창 x, 길이 L)를 — the inventory of n4_walls.all_margins (endpoint-side a, endpoint-cap b, band-corner-chord d, crossing window x, length L)
# 심볼릭으로 재구성: 게이팅(어떤 여유가 존재하는가)만 수치로 정하고, 각 여유의 부호는 — is rebuilt symbolically: only the gating (which margins exist) is decided numerically, and the sign of each margin
# 전부 Q(√15) 정확 원소로 판정한다. 기대: 활성 8개 = 0, 나머지 전부 > 0. — is decided entirely as an exact element of Q(√15). Expected: the 8 active ones = 0, all the rest > 0.
import sys, math, json
import sympy as sp
from sympy import Rational as R
r15=sp.sqrt(15)
RAY=None

# θ* 정확 데이터 — exact data for θ*
c1,s1=-R(1,4),r15/4
c2,s2=-R(7,8),r15/8
c3,s3=-R(7,8),-r15/8
c4,s4=-R(1,4),-r15/4
def rot(c,s,v): return sp.Matrix([c*v[0]-s*v[1], s*v[0]+c*v[1]])
d=[sp.Matrix([1,0])]
for c_,s_ in ((c1,s1),(c2,s2),(c3,s3),(c4,s4)):
    d.append(rot(c_,s_,d[-1]))
L=[2*r15/5, 8*r15/15, 2*r15/5]
M=[sp.Matrix([0,0])]
for i in range(3): M.append(sp.expand(M[-1]+L[i]*d[i+1]))
SIN=[s1,s2,s3,s4]
CH=[]
for i in range(4):
    CV=(d[i]+d[i+1])/(2*SIN[i])
    A=sp.expand(M[i]-CV); B=sp.expand(M[i]+CV)
    CH.append((A,B))
PIECES={0:(M[0],-d[0],RAY),1:(M[0],d[1],L[0]),2:(M[1],d[2],L[1]),3:(M[2],d[3],L[2]),4:(M[3],d[4],RAY)}
PAIRS=[(0,[4,2]),(1,[4]),(2,[0]),(3,[2,0])]

def simp(e):
    return sp.nsimplify(sp.simplify(sp.expand(e)),[r15])
def sgn_exact(e):
    """Q(√15) 원소의 정확 부호.
    exact sign of an element of Q(√15)."""
    e=sp.simplify(sp.expand(e))
    if e==0: return 0
    a=sp.Rational(sp.expand(e).coeff(r15,0)); b=sp.Rational(sp.expand(e).coeff(r15,1))
    # e = a + b√15: 부호 = sign 비교 a² vs 15b² — the sign is determined by comparing a² vs 15b²
    if b==0: return 1 if a>0 else -1
    if a==0: return 1 if b>0 else -1
    if a>0 and b>0: return 1
    if a<0 and b<0: return -1
    big = a*a - 15*b*b
    if a>0: return 1 if big>0 else -1
    return -1 if big>0 else 1
NUM=lambda e: float(sp.N(e,30))

zero_cnt=0; pos_cnt=0; bad=[]
def record(name, val, kind='lin'):
    global zero_cnt,pos_cnt
    s=sgn_exact(val)
    if s==0: zero_cnt+=1; tag='=0 (활성)'
    elif s>0: pos_cnt+=1; tag='>0'
    else: bad.append(name); tag='<0 !!'
    print(f"  {name:16s} {tag:10s} 값={sp.nsimplify(sp.radsimp(sp.simplify(val)),[r15])}")

print("=== (a) 끝점-측면 여유: ±lat−1/2 ===")
for i,ks in PAIRS:
    A,B=CH[i]
    for k in ks:
        P,dd,l=PIECES[k]
        for e,pt in (('-',A),('+',B)):
            lat=simp(-(pt[0]-P[0])*dd[1]+(pt[1]-P[1])*dd[0])
            lon=simp((pt[0]-P[0])*dd[0]+(pt[1]-P[1])*dd[1])
            lo_ok = True if l is None else (NUM(lon) < 1.02*NUM(l))
            if -0.02*(NUM(l) if l is not None else 1e9) < NUM(lon) and lo_ok:
                slat=sgn_exact(lat)
                record(f"a:c{i+1}P{k}E{e}", slat*lat - R(1,2))
print("=== (b) 끝점-캡 종방향 여유 ===")
for i,ks in PAIRS:
    A,B=CH[i]
    for k in ks:
        P,dd,l=PIECES[k]
        for e,pt in (('-',A),('+',B)):
            lat=simp(-(pt[0]-P[0])*dd[1]+(pt[1]-P[1])*dd[0])
            lon=simp((pt[0]-P[0])*dd[0]+(pt[1]-P[1])*dd[1])
            if abs(NUM(lat))<0.52:
                if NUM(lon) < 0.5*(NUM(l) if l is not None else 1e9):
                    record(f"b:c{i+1}P{k}E{e}lo", sgn_exact(lon)*lon)
                elif l is not None:
                    record(f"b:c{i+1}P{k}E{e}hi", sgn_exact(l-lon)*(l-lon))
print("=== (d) 밴드코너-현 거리² ===")
for i,ks in PAIRS:
    A,B=CH[i]
    ab=sp.expand(B-A); nab2=simp(ab.dot(ab))
    for k in ks:
        P,dd,l=PIECES[k]
        n_=sp.Matrix([-dd[1],dd[0]])
        caps=[(sp.Integer(0),'lo')]+([(l,'hi')] if l is not None else [])
        for cap,cn in caps:
            for sgn_,sn in ((R(1,2),'p'),(-R(1,2),'m')):
                Cpt=sp.expand(P+cap*dd+sgn_*n_)
                t=simp((Cpt-A).dot(ab))/nab2
                if 0.03<NUM(t)<0.97:
                    proj=sp.expand(A+t*ab)
                    record(f"d:c{i+1}P{k}K{cn}{sn}", simp((Cpt-proj).dot(Cpt-proj)))
print("=== (x) 교차창 여유 ===")
for i in range(5):
    for j in range(i+2,5):
        P1,d1,l1=PIECES[i]; P2,d2,l2=PIECES[j]
        det=simp(d1[0]*d2[1]-d1[1]*d2[0])
        if det==0: continue
        r_=sp.expand(P2-P1)
        s=simp((r_[0]*d2[1]-r_[1]*d2[0])/det); u=simp((r_[0]*d1[1]-r_[1]*d1[0])/det)
        inside = (NUM(s) > 1e-9 and (l1 is None or NUM(s)<NUM(l1)-1e-9) and
                  NUM(u) > 1e-9 and (l2 is None or NUM(u)<NUM(l2)-1e-9))
        # 각주: i=0 조각은 후방 레이 — s>0 쪽만 재료 — footnote: piece i=0 is the backward ray — only the s>0 side is material
        if inside:
            record(f"x:({i},{j})s0", s)
            record(f"x:({i},{j})u0", u)
            if l1 is not None: record(f"x:({i},{j})s1", l1-s)
            if l2 is not None: record(f"x:({i},{j})u1", l2-u)
print("=== (L) 패널 길이 ===")
for k,Lv in enumerate(L): record(f"L{k+1}", Lv)

print(f"\n[요약] 활성(=0) {zero_cnt}개 (기대 8), 양수 {pos_cnt}개, 음수 {len(bad)}개 {bad if bad else ''}")
assert zero_cnt==8 and not bad
print("결론: θ*_{4_1}의 전 후보 접촉 재고가 Q(√15)에서 정확 판정 — 활성 8개 = 0, 나머지 전부 강양수")
