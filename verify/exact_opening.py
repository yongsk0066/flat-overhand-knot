"""개방 보조정리 9.8의 정확대수 승격 (완결성 정리의 마지막 조건부 제거).
Promotion of the opening lemma 9.8 to exact algebra (removing the last proviso of the completeness theorem).

A. theta†: 활성 벽 {w1,w4,w5,w6} 기울기 4x5 행렬의 랭크 = 4 를
A. theta†: the rank of the 4x5 gradient matrix of the active walls {w1,w4,w5,w6} equals 4,
   몫체 Q[x]/(p_x)에서 4x4 소행렬식 != 0 으로 정확 판정.
   decided exactly in the quotient field Q[x]/(p_x) via a 4x4 minor != 0.
   => 전 행 랭크 => 임의 목표(전부 1)의 해 delta 존재 => 1계 동시 개방.
   => full row rank => a solution delta exists for any target (here all ones) => simultaneous first-order opening.
B. 오각형 P*: 명시적 유리 방향 delta*에 대해 8벽 방향도함수 전부 > 0 을
B. pentagon P*: for an explicit rational direction delta*, all 8 wall directional derivatives are > 0,
   Q(sqrt5, x0) (4차 체)에서 구간산술로 엄밀 판정.
   decided rigorously by interval arithmetic in Q(sqrt5, x0) (a degree-4 field).
C. 연속성 논증(사람 증명, 문서에 수록): 1계 동시 개방 + 유리식 C^1
C. continuity argument (human proof, recorded in the document): simultaneous first-order opening + C^1 rational functions
   => 작은 t>0에서 전 벽 양수·재료 정합 유지 => 열린 유효 상태의 극한.
   => for small t>0 all walls stay positive and the material stays consistent => a limit of open valid states.
"""
from sympy import symbols, Rational, S, cancel, diff, Poly, sqrt
import mpmath
from mpmath import mp, iv

x, g, h, a, b = symbols('x g h a b', real=True)
R = Rational
mp.dps = 120

def refl_lin(d):
    dx, dy = d; n = dx**2 + dy**2
    return lambda v: (cancel(((dx**2-dy**2)*v[0] + 2*dx*dy*v[1])/n),
                      cancel((2*dx*dy*v[0] + (dy**2-dx**2)*v[1])/n))
def cross(u,v): return u[0]*v[1]-u[1]*v[0]
def sub2(u,v): return (u[0]-v[0], u[1]-v[1])
def add2(u,v): return (u[0]+v[0], u[1]+v[1])

E=(S(0),S(0)); F=(-x,S(1)); G=(g,S(1)); H=(h,S(0))
d2 = refl_lin(sub2(F,E))((S(1),S(0)))
w0 = refl_lin(sub2(G,H))((S(1),S(0)))
Q1 = (g + a*w0[0], 1 + a*w0[1]); Q0 = (h + b*w0[0], b*w0[1])
rl3 = refl_lin(sub2(Q0,Q1))
rho3 = lambda P: add2(Q1, rl3(sub2(P,Q1)))
dD = rl3(w0); e1p = rho3(G)
M2s = ((g+h)/2, R(1,2))
s23 = cancel(cross(sub2(Q1,M2s), sub2(Q0,Q1)) / cross(w0, sub2(Q0,Q1)))
x0n = 0.3221274973112272
gn=(1-x0n**2)/(2*x0n); hn=1.0363863004307971; bn=1.7132441870641958
XN = {x:x0n, g:gn, h:hn, a:hn, b:bn}
sA = 1 if float(cross(sub2(F,E), d2).subs(XN)) > 0 else -1
pD = (M2s[0] + s23*w0[0] + dD[0], M2s[1] + s23*w0[1] + dD[1])
sD = 1 if float(cross(sub2(pD,e1p), dD).subs(XN)) > 0 else -1
w1 = -sA*cross(sub2(G,E), d2)
w7 =  sA*cross(sub2(Q0,E), d2) - 1
w4 = -sD*cross(sub2(H,e1p), dD)
w8 =  sD*cross(sub2(F,e1p), dD) - 1
sI = cancel(cross(sub2(F,G), d2)/cross(w0, d2))
w5 = a - sI
Xl = {x:0.25, g:(1-0.0625)/0.5, h:(1-3*0.0625)/0.5, a:(1-3*0.0625)/0.5, b:1.60}
w6raw = cross(sub2(E,rho3(H)), dD)
s6 = 1 if float(w6raw.subs(Xl)) > 0 else -1
w6 = s6*w6raw
den9 = cross(sub2(F,E), sub2(Q0,Q1))
w9  = cancel(cross(sub2(Q1,E), sub2(Q0,Q1))/den9) - 1
w10 = cancel(cross(sub2(Q1,E), sub2(F,E))/den9) - 1
VS = (x, g, h, a, b)

# ---------- A. theta†: 몫체 랭크 판정 ---------- — A. theta†: rank decision in the quotient field ----------
import sys
print("=== A. theta†에서 활성 벽 기울기의 랭크 [몫체 정확] ===")
p_x = x**12 + 4*x**10 - 7*x**4 + 20*x**2 - 2
PX = Poly(p_x, x, domain='QQ')
h_expr = (2*x**11 + 8*x**9 - x**7 - 4*x**5 - 13*x**3 + 40*x)/12
b_expr = (x**11 + 4*x**9 - 7*x**3 + 22*x)/4
g_expr = (1 - x**2)/(2*x)

class Fe:
    __slots__=('p',)
    def __init__(s, val): s.p = val if isinstance(val, Poly) else Poly(val, x, domain='QQ').rem(PX)
    def add(s,o): return Fe((s.p+o.p).rem(PX))
    def sub(s,o): return Fe((s.p-o.p).rem(PX))
    def mul(s,o): return Fe((s.p*o.p).rem(PX))
    def inv(s): return Fe(s.p.invert(PX))
    def is_zero(s): return s.p.is_zero

def fe_eval(e, env):
    if e.is_Symbol: return env[e]
    if e.is_Integer or e.is_Rational: return Fe(Poly(R(e), x, domain='QQ'))
    if e.is_Add:
        out = fe_eval(e.args[0], env)
        for t in e.args[1:]: out = out.add(fe_eval(t, env))
        return out
    if e.is_Mul:
        out = fe_eval(e.args[0], env)
        for t in e.args[1:]: out = out.mul(fe_eval(t, env))
        return out
    if e.is_Pow:
        base = fe_eval(e.base, env); n = int(e.exp)
        if n < 0: base = base.inv(); n = -n
        out = Fe(Poly(1, x, domain='QQ'))
        for _ in range(n): out = out.mul(base)
        return out
    raise TypeError(str(e.func))
def to_fe(expr):
    n_, d_ = expr.as_numer_denom()
    return Fe(Poly(n_, x, domain='QQ')).mul(Fe(Poly(d_, x, domain='QQ')).inv())
ENV = {x: Fe(Poly(x, x, domain='QQ')), g: to_fe(g_expr), h: to_fe(h_expr), a: to_fe(h_expr), b: to_fe(b_expr)}

ACT = [('w1',w1),('w4',w4),('w5',w5),('w6',w6)]
Mfe = [[fe_eval(diff(wf, v), ENV) for v in VS] for _, wf in ACT]
# 4x4 소행렬 (열 1..4 = g,h,a,b) 행렬식 — 여인수 전개 — determinant of the 4x4 minor (columns 1..4 = g,h,a,b) — cofactor expansion
def det_fe(M):
    n = len(M)
    if n == 1: return M[0][0]
    tot = Fe(Poly(0, x, domain='QQ'))
    for j in range(n):
        minor = [row[:j]+row[j+1:] for row in M[1:]]
        term = M[0][j].mul(det_fe(minor))
        tot = tot.add(term) if j % 2 == 0 else tot.sub(term)
    return tot
sub4 = [[Mfe[i][j] for j in (1,2,3,4)] for i in range(4)]
D = det_fe(sub4)
print(f"  det(열 g,h,a,b) = 0 인가? {D.is_zero()}")
assert not D.is_zero()
# 수치 참고값 — numerical reference value
xd = mpmath.findroot(lambda t: t**12+4*t**10-7*t**4+20*t**2-2, mpmath.mpf('0.3221274973112272'))
dv = sum(mpmath.mpf(str(R(c)))*xd**int(m[0]) for m, c in zip(Poly(D.p, x).monoms(), Poly(D.p, x).coeffs()))
print(f"  det 수치값 = {mpmath.nstr(dv, 10)}  => 랭크 4 [정확]")
print("  => 전 행 랭크: 모든 목표벡터(예: 전부 1)에 해 delta 존재 — 1계 동시 개방 [정확]")

# ---------- B. 오각형: 유리 방향의 8벽 동시 개방 ---------- — B. pentagon: simultaneous opening of the 8 walls along a rational direction ----------
print("\n=== B. P*에서 유리 방향 delta*의 8벽 방향도함수 > 0 [구간산술 엄밀] ===")
iv.dps = 60
s5iv = iv.sqrt(iv.mpf(5))
x0iv = iv.sqrt(1 - 2/s5iv)
g0iv = (1 - x0iv**2)/(2*x0iv); h0iv = (1 - 3*x0iv**2)/(2*x0iv)
# b0 = dot(F-H, w0) — 구간으로 직접 평가 — evaluated directly in interval arithmetic
def iv_eval(e, env):
    if e.is_Symbol: return env[e]
    if e.is_Integer or e.is_Rational:
        q = R(e); return iv.mpf(int(q.p))/iv.mpf(int(q.q))
    if e.is_Add:
        out = iv.mpf(0)
        for t in e.args: out = out + iv_eval(t, env)
        return out
    if e.is_Mul:
        out = iv.mpf(1)
        for t in e.args: out = out * iv_eval(t, env)
        return out
    if e.is_Pow:
        base = iv_eval(e.base, env); n = int(e.exp)
        if n < 0: return iv.mpf(1)/(base**(-n))
        return base**n
    raise TypeError(str(e.func))
from sympy import Symbol
envP = {x: x0iv, g: g0iv, h: h0iv, a: h0iv}
b0iv = iv_eval((w0[0]*(F[0]-H[0]) + w0[1]*(F[1]-H[1])), {x:x0iv, g:g0iv, h:h0iv})
envP[b] = b0iv
# 방향 delta*: R2/R3의 개방 방향을 단순 유리수로 (수치로 먼저 근사 해 확인) — direction delta*: the opening direction of R2/R3 expressed in simple rationals (an approximate solution is checked numerically first)
# 8벽 전부의 도함수가 ~1이 되는 최소자승 해를 배정밀로 구해 유리 근사 — compute in double precision the least-squares solution making all 8 wall derivatives ~1, then approximate it rationally
import math
XNP = {x: float(mpmath.mpf(x0iv.mid)), g: float(mpmath.mpf(g0iv.mid)), h: float(mpmath.mpf(h0iv.mid)),
       a: float(mpmath.mpf(h0iv.mid)), b: float(mpmath.mpf(b0iv.mid))}
W8 = [('w1',w1),('w4',w4),('w5',w5),('w6',w6),('w7',w7),('w8',w8),('w9',w9),('w10',w10)]
Gn = [[float(diff(wf, v).subs(XNP)) for v in VS] for _, wf in W8]
AtA = [[sum(Gn[k][i]*Gn[k][j] for k in range(8)) for j in range(5)] for i in range(5)]
Atb = [sum(Gn[k][i] for k in range(8)) for i in range(5)]
Msys = [AtA[i][:]+[Atb[i]] for i in range(5)]
for c in range(5):
    p_ = max(range(c,5), key=lambda r0: abs(Msys[r0][c])); Msys[c],Msys[p_]=Msys[p_],Msys[c]
    for r0 in range(5):
        if r0!=c:
            f0=Msys[r0][c]/Msys[c][c]; Msys[r0]=[Msys[r0][k]-f0*Msys[c][k] for k in range(6)]
dsol = [Msys[i][5]/Msys[i][i] for i in range(5)]
dstar = [R(round(v*1000), 1000) for v in dsol]
print(f"  delta* (유리 근사) = {dstar}")
allpos = True
for nm, wf in W8:
    dd = sum(dstar[i]*diff(wf, VS[i]) for i in range(5))
    val = iv_eval(dd, envP)
    pos = val.a > 0
    allpos &= bool(pos)
    print(f"  D_delta {nm} ⊂ [{mpmath.nstr(mpmath.mpf(val.a.a), 8)}, {mpmath.nstr(mpmath.mpf(val.b.a), 8)}]  양수={pos}")
print("\n====", "보조정리 9.8 승격: 두 타이트 상태의 1계 동시 개방 [정확대수/구간산술]" if (not D.is_zero()) and allpos else "실패", "====")
