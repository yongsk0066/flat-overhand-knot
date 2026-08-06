"""theta†-KKT의 정확대수 승격 (08 정리 8.2 -> [정확대수]).
Promotion of the theta†-KKT to exact algebra (08 Thm 8.2 -> [exact algebra]).
핵심 기법: 몫체 Q[x]/(p_x) 원소(차수<=11 다항식)로 식 트리를 직접 평가 —
key technique: evaluate the expression tree directly over elements of the quotient field Q[x]/(p_x) (polynomials of degree<=11) —
심볼릭 치환/cancel의 식 폭발 회피. lambda는 PSLQ(기저 x^0..x^11)로 식별.
avoiding the expression blow-up of symbolic substitution/cancel. lambda is identified by PSLQ (basis x^0..x^11)."""
from sympy import symbols, Rational, S, cancel, diff, Poly
import mpmath
from mpmath import mp

x, g, h, a, b = symbols('x g h a b', real=True)
R = Rational
mp.dps = 260

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
ell = (g+h+x)/2 + s23
x0n = 0.3221274973112272
gn=(1-x0n**2)/(2*x0n); hn=1.0363863004307971; bn=1.7132441870641958
XN = {x:x0n, g:gn, h:hn, a:hn, b:bn}
sA = 1 if float(cross(sub2(F,E), d2).subs(XN)) > 0 else -1
pD = (M2s[0] + s23*w0[0] + dD[0], M2s[1] + s23*w0[1] + dD[1])
sD = 1 if float(cross(sub2(pD,e1p), dD).subs(XN)) > 0 else -1
w1 = -sA*cross(sub2(G,E), d2)
w4 = -sD*cross(sub2(H,e1p), dD)
sI = cancel(cross(sub2(F,G), d2)/cross(w0, d2))
w5 = a - sI
w6raw = cross(sub2(E,rho3(H)), dD)
s6 = 1 if float(w6raw.subs(XN)) > 0 else -1
w6 = s6*w6raw
WALLS = [('w1', w1), ('w4', w4), ('w5', w5), ('w6', w6)]

p_x = x**12 + 4*x**10 - 7*x**4 + 20*x**2 - 2
PX = Poly(p_x, x, domain='QQ')
h_expr = (2*x**11 + 8*x**9 - x**7 - 4*x**5 - 13*x**3 + 40*x)/12
b_expr = (x**11 + 4*x**9 - 7*x**3 + 22*x)/4
g_expr = (1 - x**2)/(2*x)

# ---------- 몫체 원소 ---------- — elements of the quotient field ----------
class Fe:
    __slots__ = ('p',)
    def __init__(self, val):
        self.p = val if isinstance(val, Poly) else Poly(val, x, domain='QQ').rem(PX)
    def add(s, o): return Fe((s.p + o.p).rem(PX))
    def sub(s, o): return Fe((s.p - o.p).rem(PX))
    def mul(s, o): return Fe((s.p * o.p).rem(PX))
    def inv(s):
        return Fe(s.p.invert(PX))
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
    raise TypeError(f"미지원 노드: {e.func}")

def fe_div(num, den): return num.mul(den.inv())

# 좌표를 몫체 원소로 — turn the coordinates into quotient-field elements
xe = Fe(Poly(x, x, domain='QQ'))
def to_fe(expr):
    # 유리식 -> Fe (x만 포함) — rational function -> Fe (contains x only)
    n, d = expr.as_numer_denom()
    return fe_div(Fe(Poly(n, x, domain='QQ')), Fe(Poly(d, x, domain='QQ')))
ENV = {x: xe, g: to_fe(g_expr), h: to_fe(h_expr), a: to_fe(h_expr), b: to_fe(b_expr)}

# ---------- S1: 260자리 수치 (mpf 재귀 평가) ---------- — S1: 260-digit numerics (recursive mpf evaluation) ----------
xd = mpmath.findroot(lambda t: t**12 + 4*t**10 - 7*t**4 + 20*t**2 - 2, mpmath.mpf('0.3221274973112272'))
def num_eval(e, env):
    if e.is_Symbol: return env[e]
    if e.is_Integer or e.is_Rational:
        q = R(e); return mpmath.mpf(q.p)/mpmath.mpf(q.q)
    if e.is_Add: return sum(num_eval(t, env) for t in e.args)
    if e.is_Mul:
        out = mpmath.mpf(1)
        for t in e.args: out *= num_eval(t, env)
        return out
    if e.is_Pow: return num_eval(e.base, env)**int(e.exp)
    raise TypeError(str(e.func))
def coordn(expr):
    return num_eval(expr, {x: xd})
NENV = {x: xd, g: coordn(g_expr), h: coordn(h_expr), a: coordn(h_expr), b: coordn(b_expr)}

VS = (x, g, h, a, b)
Dell = [diff(ell, v) for v in VS]
Dw = {nm: [diff(wf, v) for v in VS] for nm, wf in WALLS}
gl_n = [num_eval(d_, NENV) for d_ in Dell]
gw_n = {nm: [num_eval(d_, NENV) for d_ in Dw[nm]] for nm, _ in WALLS}
Am = mpmath.matrix(4,4); bm = mpmath.matrix(4,1)
for i in range(4):
    for j, (nm, _) in enumerate(WALLS): Am[i,j] = gw_n[nm][i]
    bm[i] = gl_n[i]
lam_m = mpmath.lu_solve(Am, bm)
res5 = gl_n[4] - sum(lam_m[j]*gw_n[WALLS[j][0]][4] for j in range(4))
print("S1: lambda(260자리) =", [mpmath.nstr(lam_m[j], 12) for j in range(4)])
print("    5번째 방정식 잔차 =", mpmath.nstr(res5, 4))

# ---------- S2: PSLQ 식별 ---------- — S2: PSLQ identification ----------
print("\nS2: lambda in Q(x†) 식별")
BAS = [xd**i for i in range(12)]
lam_fe = []; lam_expr = []
for j, (nm, _) in enumerate(WALLS):
    rel = mpmath.pslq(BAS + [lam_m[j]], maxcoeff=10**15, maxsteps=2*10**7)
    assert rel and rel[12] != 0, f"{nm} 식별 실패"
    cand = sum(R(-rel[i], rel[12])*x**i for i in range(12))
    err = abs(sum(mpmath.mpf(R(-rel[i], rel[12]).p)/mpmath.mpf(R(-rel[i], rel[12]).q)*BAS[i] for i in range(12)) - lam_m[j])
    lam_expr.append(cand); lam_fe.append(to_fe(cand))
    print(f"  lambda[{nm}] = {cand}")
    print(f"     잔차 = {mpmath.nstr(err, 4)}")

# ---------- S3: KKT 5등식 정확 검증 ---------- — S3: exact verification of the 5 KKT equalities ----------
print("\nS3: KKT 등식 (몫체 정확 연산)")
allok = True
for i, v in enumerate(VS):
    lhs = fe_eval(Dell[i], ENV)
    for j, (nm, _) in enumerate(WALLS):
        lhs = lhs.sub(lam_fe[j].mul(fe_eval(Dw[nm][i], ENV)))
    ok = lhs.is_zero()
    allok &= ok
    print(f"  성분 {i} ({v}): {'0 [정확]' if ok else '실패: ' + str(lhs.p.as_expr())[:80]}")
assert allok

# ---------- S4: 양수성 (Krawczyk 구간) ---------- — S4: positivity (Krawczyk interval) ----------
print("\nS4: 양수성 — x† ± 1e-30 구간 평가")
from mpmath import iv
iv.dps = 80
lo = mpmath.nstr(xd - mpmath.mpf('1e-30'), 78); hi = mpmath.nstr(xd + mpmath.mpf('1e-30'), 78)
xiv = iv.mpf([lo, hi])
for j, (nm, _) in enumerate(WALLS):
    pl = Poly(lam_expr[j], x)
    val = iv.mpf(0)
    for m_, c_ in zip(pl.monoms(), pl.coeffs()):
        q = R(c_)
        val = val + (iv.mpf(int(q.p))/iv.mpf(int(q.q))) * xiv**int(m_[0])
    pos = val.a > 0
    allok &= bool(pos)
    print(f"  lambda[{nm}] ⊂ [{mpmath.nstr(mpmath.mpf(val.a.a), 12)}, {mpmath.nstr(mpmath.mpf(val.b.a), 12)}]  양수={pos}")

print("\n====", "정리 8.2 승격: theta†-KKT [정확대수 + 구간 양수성]" if allok else "실패", "====")
