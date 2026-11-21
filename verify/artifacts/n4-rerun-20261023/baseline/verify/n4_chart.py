# n4_chart.py — 유리 차트 u_i=tan(tau_i/2) 에서 8개 벽의 다항화 + θ* 정확 소멸 검증 — polynomialization of the 8 walls in the rational chart u_i=tan(tau_i/2) + verification of exact vanishing at θ*
# cos t=(1-u^2)/(1+u^2), sin t=2u/(1+u^2); 현 벡터 C_i=(d_i+d_{i+1})(1+u_{i+1}^2)/(4 u_{i+1}) — chord vector C_i
import sympy as sp, itertools, random, math

u1,u2,u3,u4,l1,l2,l3=sp.symbols('u1 u2 u3 u4 l1 l2 l3')
U=[u1,u2,u3,u4]; L=[l1,l2,l3]

def rot(u):  # 회전행렬 (cos,sin) 유리형 — rotation matrix, rational form in (cos,sin)
    c=(1-u**2)/(1+u**2); s=2*u/(1+u**2)
    return sp.Matrix([[c,-s],[s,c]])

d=[sp.Matrix([1,0])]
for k in range(4): d.append(sp.simplify(rot(U[k])*d[k]))
M=[sp.Matrix([0,0])]
for i in range(3): M.append(M[-1]+L[i]*d[i+1])
# 크리스 i (piece i - piece i+1): 현 벡터 CV_i = (d_i+d_{i+1})/(2 sin τ_{i+1}) = (d_i+d_{i+1})(1+u^2)/(4u) — crease i; chord vector CV_i
CV=[sp.simplify((d[i]+d[i+1])*(1+U[i]**2)/(4*U[i])) for i in range(4)]
E={}
for i in range(4):
    E[(i,'-')]=M[i]-CV[i]; E[(i,'+')]=M[i]+CV[i]
def lat(pt,base,dd): return -(pt[0]-base[0])*dd[1]+(pt[1]-base[1])*dd[0]

# 벽 8종 (인벤토리): a:c{i}P{k}E{e} — 기준(밴드) = P4:(M3,d4), P2:(M1,d2), P0:(M0,d0) — 8 wall kinds (inventory); reference (band)
BASE={4:(M[3],d[4]), 2:(M[1],d[2]), 0:(M[0],d[0])}
WALL_DEFS={
 "a:c1P4E+": (0,'+',4), "a:c1P2E-": (0,'-',2),
 "a:c2P4E-": (1,'-',4), "a:c2P4E+": (1,'+',4),
 "a:c3P0E-": (2,'-',0), "a:c3P0E+": (2,'+',0),
 "a:c4P0E-": (3,'-',0), "a:c4P2E+": (3,'+',2),
}
# θ* 정확값 (Q(√15)) — exact value of θ* (Q(√15))
r15=sp.sqrt(15)
SUB={u1:r15/3, u2:r15, u3:-r15, u4:-r15/3, l1:2*r15/5, l2:8*r15/15, l3:2*r15/5}

# 수치 대조용 — for numerical cross-checking
import sys
sys.path.insert(0,'verify')
from n4_walls import all_margins

def to_theta(uv,lv):
    taus=[2*math.atan(x) for x in uv]
    return tuple(taus)+tuple(lv)

print("[다항화] 벽별: 분모 소거 → Z[u,L] 다항 / θ* 소멸 / 표본 부호·배율")
rnd=random.Random(3)
import json
samples=[]
hits=json.load(open("verify/n4_valid41.json"))
for th,o in hits:
    if th[0]<0: continue           # 거울 배향 제외 (u 부호 고정: +,+,-,-) — exclude the mirror orientation (u signs fixed: +,+,-,-)
    uv=[math.tan(t/2) for t in th[:4]]; lv=list(th[4:])
    samples.append((uv,lv,th))
samples=samples[:12]
print(f"[표본] 정배향 {len(samples)}개")
results={}
def swap_e(i,e):
    # 셀에서 u3,u4<0 ⟹ 크리스 3,4(i=2,3)의 CV 방향 반전: config의 E± 와 라벨 스왑 — in this cell u3,u4<0 ⟹ the CV direction of creases 3,4 (i=2,3) is reversed: labels swap relative to config's E±
    if i in (2,3): return '-' if e=='+' else '+'
    return e
for name,(i,e0,k) in WALL_DEFS.items():
    e=swap_e(i,e0)
    base,dd=BASE[k]
    g=sp.together(lat(E[(i,e)],base,dd))
    num,den=sp.fraction(g)
    num=sp.expand(num); den=sp.factor(den)
    # 벽 다항 = num - (±1/2) den  형태가 아니라: 벽함수 s*lat-1/2 ⟹ 다항 대표 = s*num - den/2 → ×2: 2 s num - den — the wall polynomial is not of the form num - (±1/2) den: the wall function is s*lat-1/2 ⟹ polynomial representative = s*num - den/2 → ×2: 2 s num - den
    # 부호 s와 den 부호는 표본으로 캘리브레이션 — the sign s and the sign of den are calibrated from the samples
    latf=sp.lambdify((u1,u2,u3,u4,l1,l2,l3), g, 'math')
    denf=sp.lambdify((u1,u2,u3,u4,l1,l2,l3), den, 'math')
    svals=[]; dvals=[]
    for uv,lv,th in samples:
        try:
            svals.append(latf(*uv,*lv)); dvals.append(denf(*uv,*lv))
        except Exception: pass
    s=1 if sum(1 for v in svals if v>0)>=len(svals)/2 else -1
    dsign="+" if all(v>0 for v in dvals) else ("-" if all(v<0 for v in dvals) else "혼합!")
    W=sp.expand(2*s*num - den)   # 항등식은 den 부호 무관; 유효측은 sign(W)=sign(den) — the identity does not depend on den's sign; the valid side is sign(W)=sign(den)
    # θ* 소멸 검증 (다항 대표) — verify vanishing at θ* (polynomial representative)
    z=sp.simplify(W.subs(SUB))
    deg=sp.total_degree(W)
    nterm=len(W.as_ordered_terms())
    # 수치 벽 대조: s*lat-1/2 vs all_margins의 |lat|-1/2 — 표본 3개 — numerical wall cross-check: s*lat-1/2 vs all_margins' |lat|-1/2 — 3 samples
    ok=[]
    for uv,lv,th in samples[:3]:
        m=all_margins(th)
        if m and name in m:
            v1=s*latf(*uv,*lv)-0.5
            ok.append(abs(v1-m[name])<1e-9)
    results[name]=(z,deg,nterm,dsign,all(ok) if ok else None)
    print(f"  {name}: θ*소멸={'0 ✓' if z==0 else z} 차수={deg} 항={nterm} 분모부호={dsign} 수치대조={results[name][4]}")
print("\n[u* 정확값]  u=(√15/3, √15, -√15, -√15/3), L=(2√15/5, 8√15/15, 2√15/5)")
