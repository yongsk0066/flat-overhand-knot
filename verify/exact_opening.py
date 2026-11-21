"""Exact opening certificates for theta-dagger and the regular pentagon P*.

A. At theta-dagger, the 4x4 minor of the active wall gradients {w1, w4, w5, w6}
   in columns (g, h, a, b) is certified nonzero in Q[x]/(p_x), so a simultaneous
   first-order opening direction delta exists with grad w_j . delta = 1.
B. At P*, the explicit rational direction
   delta_op = (-93/1000, 2197/1000, 1021/500, 509/250, 419/200)
   is certified by outward-rounded interval arithmetic on an interval-Newton
   enclosure of x0 to have all eight wall directional derivatives in [0.984, 1.018],
   and exact core directional derivative 8277/2000.

All band orientations are fixed by exact interior probes, never by floating-point signs.
"""
import sympy as sp
from sympy import symbols, Rational as R, S, cancel, diff, Poly
from mpmath import mp, iv, mpf, findroot, nstr

def _require(cond, msg="verification failed"):
    if not cond:
        raise RuntimeError(msg)


def certify():
    x, g, h, a, b = symbols("x g h a b", real=True)

    def refl_lin(d):
        dx, dy = d
        n = dx ** 2 + dy ** 2
        return lambda v: (
            cancel(((dx ** 2 - dy ** 2) * v[0] + 2 * dx * dy * v[1]) / n),
            cancel((2 * dx * dy * v[0] + (dy ** 2 - dx ** 2) * v[1]) / n),
        )

    def cross(u, v):
        return u[0] * v[1] - u[1] * v[0]

    def sub2(u, v):
        return (u[0] - v[0], u[1] - v[1])

    def add2(u, v):
        return (u[0] + v[0], u[1] + v[1])

    E = (S(0), S(0))
    F = (-x, S(1))
    G = (g, S(1))
    H = (h, S(0))
    d2 = refl_lin(sub2(F, E))((S(1), S(0)))
    w0 = refl_lin(sub2(G, H))((S(1), S(0)))
    Q1 = (g + a * w0[0], 1 + a * w0[1])
    Q0 = (h + b * w0[0], b * w0[1])
    rl3 = refl_lin(sub2(Q0, Q1))
    rho3 = lambda P: add2(Q1, rl3(sub2(P, Q1)))
    dD = rl3(w0)
    e1p = rho3(G)

    w1 = -cross(sub2(G, E), d2)
    w4 = -cross(sub2(H, e1p), dD)
    sI = cancel(cross(sub2(F, G), d2) / cross(w0, d2))
    w5 = a - sI
    w6 = cross(sub2(E, rho3(H)), dD)
    w7 = cross(sub2(Q0, E), d2) - 1
    w8 = cross(sub2(F, e1p), dD) - 1
    den9 = cross(sub2(F, E), sub2(Q0, Q1))
    w9 = cancel(cross(sub2(Q1, E), sub2(Q0, Q1)) / den9) - 1
    w10 = cancel(cross(sub2(Q1, E), sub2(F, E)) / den9) - 1
    VS = (x, g, h, a, b)

    # Exact w6 orientation check at rational slack x = 1/4.
    xr = R(1, 4)
    slack_env = {
        x: xr,
        g: (1 - xr ** 2) / (2 * xr),
        h: (1 - 3 * xr ** 2) / (2 * xr),
        a: (1 - 3 * xr ** 2) / (2 * xr),
        b: (1 + xr ** 2) / (2 * xr),
    }
    _require(sp.simplify(w6.subs(slack_env)) > 0)

    # --- Part A: rank-4 opening minor at theta-dagger in Q[x]/(p_x) ---
    p_x = x ** 12 + 4 * x ** 10 - 7 * x ** 4 + 20 * x ** 2 - 2
    PX = Poly(p_x, x, domain="QQ")
    h_expr = (2 * x ** 11 + 8 * x ** 9 - x ** 7 - 4 * x ** 5 - 13 * x ** 3 + 40 * x) / 12
    b_expr = (x ** 11 + 4 * x ** 9 - 7 * x ** 3 + 22 * x) / 4
    g_expr = (1 - x ** 2) / (2 * x)

    class Fe:
        __slots__ = ("p",)

        def __init__(self, val):
            self.p = val if isinstance(val, Poly) else Poly(val, x, domain="QQ").rem(PX)

        def add(self, o):
            return Fe((self.p + o.p).rem(PX))

        def sub(self, o):
            return Fe((self.p - o.p).rem(PX))

        def mul(self, o):
            return Fe((self.p * o.p).rem(PX))

        def inv(self):
            return Fe(self.p.invert(PX))

        def is_zero(self):
            return self.p.is_zero

    def fe_eval(expr, env):
        if expr.is_Symbol:
            return env[expr]
        if expr.is_Integer or expr.is_Rational:
            return Fe(Poly(R(expr), x, domain="QQ"))
        if expr.is_Add:
            out = fe_eval(expr.args[0], env)
            for t in expr.args[1:]:
                out = out.add(fe_eval(t, env))
            return out
        if expr.is_Mul:
            out = fe_eval(expr.args[0], env)
            for t in expr.args[1:]:
                out = out.mul(fe_eval(t, env))
            return out
        if expr.is_Pow:
            base = fe_eval(expr.base, env)
            n = int(expr.exp)
            if n < 0:
                base = base.inv()
                n = -n
            out = Fe(Poly(1, x, domain="QQ"))
            for _ in range(n):
                out = out.mul(base)
            return out
        raise TypeError(str(expr.func))

    def to_fe(expr):
        n_, d_ = expr.as_numer_denom()
        return Fe(Poly(n_, x, domain="QQ")).mul(Fe(Poly(d_, x, domain="QQ")).inv())

    ENV = {
        x: Fe(Poly(x, x, domain="QQ")),
        g: to_fe(g_expr),
        h: to_fe(h_expr),
        a: to_fe(h_expr),
        b: to_fe(b_expr),
    }

    # Exact A- and D-band interior probes at theta-dagger.
    A_in = ((E[0] + F[0]) / 2 + d2[0], (E[1] + F[1]) / 2 + d2[1])
    D_in = ((Q1[0] + Q0[0]) / 2 + dD[0], (Q1[1] + Q0[1]) / 2 + dD[1])
    _require(fe_eval(cross(sub2(A_in, E), d2) - R(1, 2), ENV).is_zero())
    _require(fe_eval(cross(sub2(D_in, e1p), dD) - R(1, 2), ENV).is_zero())

    ACT = [("w1", w1), ("w4", w4), ("w5", w5), ("w6", w6)]
    Mfe = [[fe_eval(diff(wf, v), ENV) for v in VS] for _, wf in ACT]

    def det_fe(M):
        n = len(M)
        if n == 1:
            return M[0][0]
        tot = Fe(Poly(0, x, domain="QQ"))
        for j in range(n):
            minor = [row[:j] + row[j + 1:] for row in M[1:]]
            term = M[0][j].mul(det_fe(minor))
            tot = tot.add(term) if j % 2 == 0 else tot.sub(term)
        return tot

    D = det_fe([[Mfe[i][j] for j in (1, 2, 3, 4)] for i in range(4)])
    _require(not D.is_zero(), "Active-wall 4x4 minor at theta-dagger vanished")
    print("[exact] theta-dagger: 4x4 minor in (g,h,a,b) is nonzero in Q[x]/(p_x) => rank 4 opening")

    # --- Part B: simultaneous 8-wall opening at P* along delta_op ---
    dstar = (R(-93, 1000), R(2197, 1000), R(1021, 500), R(509, 250), R(419, 200))
    _require(sum(dstar) / 2 == R(8277, 2000), "Pentagon opening direction has the wrong core slope")

    mp.dps = 70
    iv.dps = 70
    x0_approx = findroot(lambda tt: 5 * tt ** 4 - 10 * tt ** 2 + 1, mpf("0.3249196962"))
    rad = mpf("1e-50")
    x0iv = iv.mpf([x0_approx - rad, x0_approx + rad])
    m_iv = iv.mpf(x0_approx)
    N_iv = m_iv - (5 * m_iv ** 4 - 10 * m_iv ** 2 + 1) / (20 * x0iv ** 3 - 20 * x0iv)
    _require(x0iv.a < N_iv.a and N_iv.b < x0iv.b, "Interval Newton contraction failed for x0")

    g0iv = (1 - x0iv ** 2) / (2 * x0iv)
    h0iv = (1 - 3 * x0iv ** 2) / (2 * x0iv)
    b0iv = (1 + x0iv ** 2) / (2 * x0iv)
    envP = {x: x0iv, g: g0iv, h: h0iv, a: h0iv, b: b0iv}

    def iv_eval(expr, env):
        if expr.is_Symbol:
            return env[expr]
        if expr.is_Integer or expr.is_Rational:
            q = R(expr)
            return iv.mpf(int(q.p)) / iv.mpf(int(q.q))
        if expr.is_Add:
            out = iv.mpf(0)
            for t in expr.args:
                out = out + iv_eval(t, env)
            return out
        if expr.is_Mul:
            out = iv.mpf(1)
            for t in expr.args:
                out = out * iv_eval(t, env)
            return out
        if expr.is_Pow:
            base = iv_eval(expr.base, env)
            n = int(expr.exp)
            if n < 0:
                return iv.mpf(1) / (base ** (-n))
            return base ** n
        raise TypeError(str(expr.func))

    W8 = [
        ("w1", w1),
        ("w4", w4),
        ("w5", w5),
        ("w6", w6),
        ("w7", w7),
        ("w8", w8),
        ("w9", w9),
        ("w10", w10),
    ]
    for nm, wf in W8:
        dd = sum(dstar[i] * diff(wf, VS[i]) for i in range(5))
        val = iv_eval(dd, envP)
        _require(val.a > 0, f"Directional derivative of {nm} not strictly positive: {val}")
        _require(val.a > mpf("0.984") and val.b < mpf("1.018"), (nm, val))
        print(f"  D_delta {nm} in [{nstr(mpf(val.a.a), 10)}, {nstr(mpf(val.b.a), 10)}] > 0")

    print("PASS: first-order simultaneous opening certified at theta-dagger and P*")


if __name__ == "__main__":
    certify()
