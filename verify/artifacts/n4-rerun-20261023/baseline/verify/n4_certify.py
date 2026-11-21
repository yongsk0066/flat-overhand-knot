# n4_certify.py — B_{4_1} 타이트 상태의 정확 인증 (Q(√15) 정확대수 + KKT 수치) — exact certification of the B_{4_1} tight state (Q(√15) exact algebra + KKT numerics)
# θ* : cos t = -1/4, s = 2π-2t, L = (6,8,6)/√15;  ell* = 20/√15 = 4√15/3
import sympy as sp

r15=sp.sqrt(15)
ct=-sp.Rational(1,4); st_=r15/4              # cos t, sin t (t∈(π/2,π))
c2t=2*ct**2-1; s2t=2*st_*ct                  # cos 2t, sin 2t
# t+s = 2π-t ⟹ cos=ct, sin=-st ; t+s/2 = π ⟹ e_s=(-1,0)
L1=2*r15/5; L2=8*r15/15
d0=sp.Matrix([1,0]); d1=sp.Matrix([ct,st_]); d2=sp.Matrix([ct,-st_])  # d2=e(t+s)
M0=sp.Matrix([0,0]); M1=L1*d1; M2=M1+L2*d2; M3=M2+L1*d1
# 반각: cos(t/2)=√(3/8)=√6/4, sin(t/2)=√(10)/4 — half-angle
ch2=sp.sqrt(6)/4; sh2=sp.sqrt(10)/4
rt=1/(2*sh2); et=sp.Matrix([ch2,sh2])
# s/2 = π-t ⟹ sin(s/2)=sin t, e(t+s/2)=(-1,0)
rs=1/(2*st_); es=sp.Matrix([-1,0])
E1p=M0+rt*et; E1m=M0-rt*et
E2p=M1+rs*es; E2m=M1-rs*es
E3p=M2+rs*es; E3m=M2-rs*es
E4p=M3+rt*et; E4m=M3-rt*et

def ok(name, expr):
    residual = sp.simplify(expr)
    assert residual == 0, (name, residual)
    print(f"  {name:34} PASS")
print("[접촉 등식 — Q(√15,√10,√6) 정확]")
ok("M3.y - 1", M3[1]-1)
ok("g1: (M3.y-E1p.y) - 1/2", (M3[1]-E1p[1])-sp.Rational(1,2))
lat2=lambda P: -(P[0]-M1[0])*d2[1]+(P[1]-M1[1])*d2[0]
ok("g2: lat_P2(E1m) - 1/2", lat2(E1m)-sp.Rational(1,2))
ok("g3a: (E2p.y-M3.y) - 1/2", (E2p[1]-M3[1])-sp.Rational(1,2))
ok("g3b: (E2m.y-M3.y) - 1/2", (E2m[1]-M3[1])-sp.Rational(1,2))
print("[대칭 이미지]")
ok("g1': (E4m.y) - 1/2", E4m[1]-sp.Rational(1,2))            # P0 상모서리 y=1/2 — upper edge
lat2b=lambda P: -(P[0]-M2[0])*d2[1]+(P[1]-M2[1])*d2[0]        # P2 기준 (M2 쪽) — referenced to P2 (M2 side)
ok("g2': -lat_P2(E4p) - 1/2", -lat2b(E4p)-sp.Rational(1,2))
ok("g3a': -(E3p.y) - 1/2", -E3p[1]-sp.Rational(1,2))          # P0 하모서리 y=-1/2 — lower edge
ok("g3b': -(E3m.y) - 1/2", -E3m[1]-sp.Rational(1,2))
print("[코어 육각형 꼭짓점 (정확)]")
V=[E4m if False else None]
hexv={"V0=E4-": (E4m, sp.Matrix([-13*r15/30, sp.Rational(1,2)])),
      "W1=E4+": (E4p, sp.Matrix([-7*r15/30, sp.Rational(3,2)])),
      "V4=E2-": (E2m, sp.Matrix([ r15/30, sp.Rational(3,2)])),
      "V5=E1+": (E1p, sp.Matrix([ 3*r15/30, sp.Rational(1,2)])),
      "W2=E1-": (E1m, sp.Matrix([-3*r15/30,-sp.Rational(1,2)])),
      "V9=E3+": (E3p, sp.Matrix([-11*r15/30,-sp.Rational(1,2)]))}
for name,(E,T) in hexv.items():
    ok(name, sp.simplify((E-T).norm()))
print("[크리스 끝점 일치 항등식 (신규 구조)]")
ok("E2+ = E4+ (이중점 상)", sp.simplify((E2p-E4p).norm()))
ok("E3- = E1- (이중점 하)", sp.simplify((E3m-E1m).norm()))
print("[변 길이]")
seq=["V0=E4-","W1=E4+","V4=E2-","V5=E1+","W2=E1-","V9=E3+"]
tgt=[2*sp.sqrt(10)/5, 4/r15, 4/r15, 2*sp.sqrt(10)/5, 4/r15, 4/r15]
for k in range(6):
    a=hexv[seq[k]][1]; b=hexv[seq[(k+1)%6]][1]
    ok(f"{seq[k]}->{seq[(k+1)%6]} = {tgt[k]}", sp.simplify((b-a).norm()-tgt[k]))
print("[내각]")
import sympy as sp2
angs=[]
for k in range(6):
    p=hexv[seq[k]][1]; a=hexv[seq[(k-1)%6]][1]; b=hexv[seq[(k+1)%6]][1]
    v1=a-p; v2=b-p
    cosang=sp.simplify((v1.dot(v2))/(v1.norm()*v2.norm()))
    angs.append(cosang)
    print(f"  {seq[k]}: cos(내각) = {cosang}  (= {sp.deg(sp.acos(cosang)).evalf(8)}°)")
print("[총합 검산]")
ok("ell* - 4√15/3", 2*L1+L2-4*r15/3)
ok("둘레 - (4√10/5+16/√15)", sum(((hexv[seq[(k+1)%6]][1]-hexv[seq[k]][1]).norm() for k in range(6)), sp.S(0))-(4*sp.sqrt(10)/5+16/r15))
A=sp.S(0)  # 신발끈 — shoelace
for k in range(6):
    a=hexv[seq[k]][1]; b=hexv[seq[(k+1)%6]][1]
    A+=a[0]*b[1]-b[0]*a[1]
print("  면적 =", sp.simplify(abs(A)/2), "=", sp.simplify(abs(A)/2).evalf(12))

# ---- [referee-H F7] 타일링의 정확 증명: 스케일 x_s=2√15·x 에서 전 좌표 유리 ---- — exact proof of the tiling: at the scale x_s=2√15·x all coordinates are rational
print("[타일링 정확 증명 (스케일 좌표, 유리 포함-배제)]")
S=2*r15
def sc(P): return (sp.simplify(sp.expand(S*P[0])), sp.simplify(sp.expand(P[1])))   # x축만 스케일 (y는 이미 유리) — scale the x-axis only (y is already rational)
# 코어 패널: P1=quad(E1-,E1+,E2+,E2-)? 패널 k = 크리스 k, k+1 현 사이 밴드 사각형 — core panels: P1=quad(E1-,E1+,E2+,E2-)? panel k = the band quadrilateral between the chords of creases k and k+1
Pan1=[sc(E1m),sc(E1p),sc(E2p),sc(E2m)]
Pan2=[sc(E2m),sc(E2p),sc(E3p),sc(E3m)]
Pan3=[sc(E3m),sc(E3p),sc(E4p),sc(E4m)]
HEXs=[sc(v) for v in [sp.Matrix([-13*r15/30,sp.Rational(1,2)]),sp.Matrix([-7*r15/30,sp.Rational(3,2)]),
      sp.Matrix([r15/30,sp.Rational(3,2)]),sp.Matrix([3*r15/30,sp.Rational(1,2)]),
      sp.Matrix([-3*r15/30,-sp.Rational(1,2)]),sp.Matrix([-11*r15/30,-sp.Rational(1,2)])]]
def rat(pts):
    out=[]
    for x_,y_ in pts:
        assert sp.simplify(x_).is_rational and sp.simplify(y_).is_rational, (x_,y_)
        out.append((sp.Rational(x_),sp.Rational(y_)))
    return out
Pan1,Pan2,Pan3,HEXs=rat(Pan1),rat(Pan2),rat(Pan3),rat(HEXs)
def convex_order(pts):
    cx=sum(p[0] for p in pts)/len(pts); cy=sum(p[1] for p in pts)/len(pts)
    from functools import cmp_to_key
    def compare(p, q):
        px, py = p[0]-cx, p[1]-cy
        qx, qy = q[0]-cx, q[1]-cy
        hp = 0 if py > 0 or (py == 0 and px >= 0) else 1
        hq = 0 if qy > 0 or (qy == 0 and qx >= 0) else 1
        if hp != hq:
            return -1 if hp < hq else 1
        turn = px*qy-py*qx
        return -1 if turn > 0 else (1 if turn < 0 else 0)
    result = sorted(pts, key=cmp_to_key(compare))
    for i in range(len(result)):
        a, b, c = result[i], result[(i+1)%len(result)], result[(i+2)%len(result)]
        assert (b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0]) > 0
    return result
Pan1,Pan2,Pan3,HEXs=convex_order(Pan1),convex_order(Pan2),convex_order(Pan3),convex_order(HEXs)
def area(poly):
    A=sp.S(0)
    for i in range(len(poly)):
        x1,y1=poly[i]; x2,y2=poly[(i+1)%len(poly)]
        A+=x1*y2-x2*y1
    return abs(A)/2
def clip(subj, cl):
    """Sutherland-Hodgman: subj ∩ (볼록 cl) — 전부 유리
    Sutherland-Hodgman: subj ∩ (convex cl) — everything rational"""
    def side(p,a,b): return (b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])
    # cl 방향 통일 (반시계) — normalize the orientation of cl (counter-clockwise)
    if area(cl)>0 and sum((cl[(i+1)%len(cl)][0]-cl[i][0])*(cl[(i+1)%len(cl)][1]+cl[i][1]) for i in range(len(cl)))>0:
        cl=cl[::-1]
    out=subj[:]
    for i in range(len(cl)):
        a,b=cl[i],cl[(i+1)%len(cl)]
        inp=out; out=[]
        if not inp: break
        for j in range(len(inp)):
            P=inp[j]; Q=inp[(j+1)%len(inp)]
            sp_=side(P,a,b); sq=side(Q,a,b)
            if sp_>=0: out.append(P)
            if (sp_>0 and sq<0) or (sp_<0 and sq>0):
                t=sp.Rational(sp_, sp_-sq)
                out.append((P[0]+t*(Q[0]-P[0]), P[1]+t*(Q[1]-P[1])))
    return out
a1,a2,a3=area(Pan1),area(Pan2),area(Pan3)
a12=area(clip(Pan1,Pan2)) if clip(Pan1,Pan2) else 0
a23=area(clip(Pan2,Pan3)) if clip(Pan2,Pan3) else 0
a13=area(clip(Pan1,Pan3)) if clip(Pan1,Pan3) else 0
tri=clip(clip(Pan1,Pan2),Pan3)
a123=area(tri) if tri else 0
ahex=area(HEXs)
union=a1+a2+a3-a12-a23-a13+a123
print(f"  |P1|={a1} |P2|={a2} |P3|={a3}  ∩12={a12} ∩23={a23} ∩13={a13} ∩123={a123}")
print(f"  합집합 = {union},  육각형 = {ahex},  일치: {sp.simplify(union-ahex)==0}")
inside=all(v in HEXs for v in Pan1+Pan2+Pan3)
assert inside
verts_on_hex=sum(1 for v in Pan1+Pan2+Pan3 if v in HEXs)
print(f"  패널 꼭짓점 12개 중 육각형 꼭짓점과 일치: {verts_on_hex}/12 (볼록성 ⟹ 합집합 ⊆ 육각형)")
print(f"  평균 겹침 = {(a1+a2+a3)}/{ahex} = {sp.Rational(a1+a2+a3,ahex)}")
assert sp.simplify(union-ahex)==0
print("  ⟹ 노치 없음 [정확대수] — shapely 1e-15 검증을 대체")
