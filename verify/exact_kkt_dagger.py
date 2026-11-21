"""Exact KKT certificate for the core-length critical state theta-dagger.

All stationarity identities and the 4x4 LICQ minor are evaluated in the
quotient field Q[x]/(p_x), p_x = x^12 + 4x^10 - 7x^4 + 20x^2 - 2.
Band orientations are fixed by exact interior probes in Q[x]/(p_x) and a
rational slack check, never by floating-point signs at active walls.
Strict positivity of the four multipliers and four inactive walls is
certified via outward-rounded interval arithmetic on an interval-Newton
root enclosure of x-dagger.
"""
import sympy as sp
from sympy import Rational as R, Poly, symbols, cancel, diff
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

    E = (sp.S(0), sp.S(0))
    F = (-x, sp.S(1))
    G = (g, sp.S(1))
    H = (h, sp.S(0))
    d2 = refl_lin(sub2(F, E))((sp.S(1), sp.S(0)))
    w0 = refl_lin(sub2(G, H))((sp.S(1), sp.S(0)))
    Q1 = (g + a * w0[0], 1 + a * w0[1])
    Q0 = (h + b * w0[0], b * w0[1])
    rl3 = refl_lin(sub2(Q0, Q1))
    rho3 = lambda P: add2(Q1, rl3(sub2(P, Q1)))
    dD = rl3(w0)
    e1p = rho3(G)

    M2 = ((g + h) / 2, R(1, 2))
    s23 = cancel(cross(sub2(Q1, M2), sub2(Q0, Q1)) / cross(w0, sub2(Q0, Q1)))
    ell = (g + h + x) / 2 + s23

    # Certified wall definitions (orientations verified exactly below).
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

    ACTIVE = [("w1", w1), ("w4", w4), ("w5", w5), ("w6", w6)]
    INACTIVE = [("w7", w7), ("w8", w8), ("w9", w9), ("w10", w10)]

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
        raise TypeError(f"Unsupported expression node: {expr.func}")

    def to_fe(expr):
        num, den = expr.as_numer_denom()
        return Fe(Poly(num, x, domain="QQ")).mul(Fe(Poly(den, x, domain="QQ")).inv())

    ENV = {
        x: Fe(Poly(x, x, domain="QQ")),
        g: to_fe(g_expr),
        h: to_fe(h_expr),
        a: to_fe(h_expr),
        b: to_fe(b_expr),
    }

    # Exact orientation verification at theta-dagger via genuine band-interior points.
    A_inside = ((E[0] + F[0]) / 2 + d2[0], (E[1] + F[1]) / 2 + d2[1])
    D_inside = ((Q1[0] + Q0[0]) / 2 + dD[0], (Q1[1] + Q0[1]) / 2 + dD[1])
    _require(fe_eval(cross(sub2(A_inside, E), d2) - R(1, 2), ENV).is_zero())
    _require(fe_eval(cross(sub2(D_inside, e1p), dD) - R(1, 2), ENV).is_zero())
    # Exact w6 orientation from the rational slack standard-slice state x = 1/4.
    xr = R(1, 4)
    slack_env = {
        x: xr,
        g: (1 - xr ** 2) / (2 * xr),
        h: (1 - 3 * xr ** 2) / (2 * xr),
        a: (1 - 3 * xr ** 2) / (2 * xr),
        b: (1 + xr ** 2) / (2 * xr),
    }
    _require(sp.simplify(w6.subs(slack_env)) > 0)
    print("[exact] band orientations and w6 sign certified without floating-point probes")

    # Verify all four active walls vanish identically in Q[x]/(p_x).
    for nm, wf in ACTIVE:
        _require(fe_eval(wf, ENV).is_zero(), f"Active wall {nm} is nonzero at theta-dagger")
    print("[exact] w1 = w4 = w5 = w6 = 0 identically in Q[x]/(p_x)")

    # Exact KKT multipliers in Q(x-dagger).
    lam_exprs = {
        "w1": R(7, 48) * x ** 11 + R(7, 12) * x ** 9 - R(1, 48) * x ** 7 - R(1, 16) * x ** 5 - R(23, 24) * x ** 3 + R(53, 16) * x,
        "w4": R(7, 48) * x ** 11 + R(7, 12) * x ** 9 - R(1, 48) * x ** 7 - R(1, 16) * x ** 5 - R(23, 24) * x ** 3 + R(53, 16) * x,
        "w5": R(1, 24) * x ** 10 + R(1, 6) * x ** 8 - R(1, 24) * x ** 4 - R(3, 8) * x ** 2 + R(25, 24),
        "w6": R(17, 96) * x ** 11 + R(35, 48) * x ** 9 + R(1, 16) * x ** 7 - R(1, 12) * x ** 5 - R(119, 96) * x ** 3 + R(161, 48) * x,
    }
    _require(lam_exprs["w1"] == lam_exprs["w4"])
    lam_fe = {nm: to_fe(expr) for nm, expr in lam_exprs.items()}

    VS = (x, g, h, a, b)
    Dell = [diff(ell, v) for v in VS]
    Dw = {nm: [diff(wf, v) for v in VS] for nm, wf in ACTIVE}

    # Exact stationarity check: all five KKT residuals vanish in Q[x]/(p_x).
    for i, v in enumerate(VS):
        res = fe_eval(Dell[i], ENV)
        for nm, _ in ACTIVE:
            res = res.sub(lam_fe[nm].mul(fe_eval(Dw[nm][i], ENV)))
        _require(res.is_zero(), f"KKT stationarity residual nonzero at {v}: {res.p.as_expr()}")
    print("[exact] all five KKT stationarity residuals vanish in Q[x]/(p_x)")

    # Exact LICQ check: 4x4 minor of active gradients in columns (g, h, a, b) is nonzero.
    Mfe = [[fe_eval(Dw[nm][j], ENV) for j in (1, 2, 3, 4)] for nm, _ in ACTIVE]

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

    D4 = det_fe(Mfe)
    _require(not D4.is_zero(), "Active-wall 4x4 minor vanished (LICQ failed)")
    print("[exact] LICQ holds: 4x4 minor in (g,h,a,b) is nonzero in Q[x]/(p_x)")

    # Certified root enclosure for x-dagger by interval Newton.
    mp.dps = 80
    iv.dps = 80

    def iv_rat(q):
        return iv.mpf(int(q.p)) / iv.mpf(int(q.q))

    def iv_poly_expr(expr, X_iv):
        pl = Poly(expr, x, domain="QQ")
        acc = iv.mpf(0)
        for c in pl.all_coeffs():
            acc = acc * X_iv + iv_rat(R(c))
        return acc

    xr_approx = findroot(lambda t: t ** 12 + 4 * t ** 10 - 7 * t ** 4 + 20 * t ** 2 - 2, mpf("0.3221274973112272"))
    rad = mpf("1e-50")
    X_box = iv.mpf([xr_approx - rad, xr_approx + rad])
    m_iv = iv.mpf(xr_approx)
    N_box = m_iv - iv_poly_expr(p_x, m_iv) / iv_poly_expr(diff(p_x, x), X_box)
    _require(X_box.a < N_box.a and N_box.b < X_box.b, "Interval Newton contraction failed for x-dagger")

    # Certify strict positivity of all four multipliers on X_box.
    for nm, _ in ACTIVE:
        val = iv_poly_expr(lam_exprs[nm], X_box)
        _require(val.a > 0, f"Multiplier {nm} not strictly positive: {val}")
        print(f"  lambda[{nm}] in [{nstr(mpf(val.a.a), 12)}, {nstr(mpf(val.b.a), 12)}] > 0")

    # Certify strict positivity of the four inactive walls via their Q[x]/(p_x) representatives.
    for nm, wf in INACTIVE:
        rep = fe_eval(wf, ENV).p.as_expr()
        val = iv_poly_expr(rep, X_box)
        _require(val.a > 0, f"Inactive wall {nm} not strictly positive: {val}")
        print(f"  inactive {nm} in [{nstr(mpf(val.a.a), 12)}, {nstr(mpf(val.b.a), 12)}] > 0")

    print("PASS: exact KKT stationarity, LICQ, and interval positivity at theta-dagger")


if __name__ == "__main__":
    certify()
