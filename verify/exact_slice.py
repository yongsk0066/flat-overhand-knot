"""Exact algebraic derivation of the standard symmetric slice and regular pentagon.

Derives G, H, k, and t from the four slice contact equations, proves the exact
sum-of-squares monotonicity of Delta(x), verifies all five equal sides and 108-deg
interior angles of the regular pentagon at x0 = sqrt(1 - 2/sqrt(5)), and checks
the polynomial linear_combination certificates used in formal/Formal/Pentagon.lean.
"""
from sympy import (
    symbols,
    Rational,
    sqrt,
    simplify,
    cancel,
    expand,
    together,
    solve,
    diff,
    S,
)

def _require(cond, msg="verification failed"):
    if not cond:
        raise RuntimeError(msg)

x = symbols("x", positive=True)
R = Rational

ok_all = True


def check(name, expr):
    global ok_all
    z = simplify(expr)
    good = z == 0
    print(f"[{'PASS' if good else 'FAIL'}] {name}" + ("" if good else f"   residual = {z}"))
    ok_all = ok_all and good


def report(name, val):
    print(f"       {name} = {val}")


def refl_lin(d):
    dx, dy = d
    n = dx ** 2 + dy ** 2
    return lambda v: (
        cancel(((dx ** 2 - dy ** 2) * v[0] + 2 * dx * dy * v[1]) / n),
        cancel((2 * dx * dy * v[0] + (dy ** 2 - dx ** 2) * v[1]) / n),
    )


def refl_pt(P, d):
    L = refl_lin(d)

    def rho(Q):
        w = L((Q[0] - P[0], Q[1] - P[1]))
        return (cancel(P[0] + w[0]), cancel(P[1] + w[1]))

    return rho


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def addm(a, d, s):
    return (a[0] + s * d[0], a[1] + s * d[1])


def line_meet_param(P1, d1, P2, d2):
    return cancel(cross(sub(P2, P1), d2) / cross(d1, d2))


def certify():
    global ok_all
    ok_all = True

    print("======== A. Derivation of the standard slice (C1-C4) ========")
    E = (S(0), S(0))
    F = (-x, S(1))
    L1 = refl_lin(sub(F, E))
    d2 = L1((S(1), S(0)))
    check("A1a: d2_x = (x^2-1)/(1+x^2)", d2[0] - (x ** 2 - 1) / (1 + x ** 2))
    check("A1b: d2_y = -2x/(1+x^2)", d2[1] + 2 * x / (1 + x ** 2))

    gs = symbols("gs")
    g_roots = solve(cross((gs, S(1)), d2), gs)
    report("C1 root g", g_roots)
    g = [r0 for r0 in g_roots if simplify(r0 - (1 - x ** 2) / (2 * x)) == 0]
    _require(g and len(g_roots) == 1, f"Unexpected C1 roots: {g_roots}")
    g = g[0]
    check("A2: C1 => g = (1-x^2)/(2x) uniquely", S(0))

    hs_ = symbols("hs_")
    w0_h = refl_lin((g - hs_, S(1)))((S(1), S(0)))
    h_roots = solve(together(cross(sub(F, (hs_, S(0))), w0_h)), hs_)
    report("C2 roots h", h_roots)
    h_expect = (1 - 3 * x ** 2) / (2 * x)
    h_good = [r0 for r0 in h_roots if simplify(r0 - h_expect) == 0]
    _require(h_good, f"Expected h root not found: {h_roots}")
    h = h_good[0]

    H = (h, S(0))
    G = (g, S(1))
    w0 = refl_lin(sub(G, H))((S(1), S(0)))
    check("A3a: r_x = (x^2-1)/(1+x^2)", w0[0] - (x ** 2 - 1) / (1 + x ** 2))
    check("A3b: r_y = 2x/(1+x^2)", w0[1] - 2 * x / (1 + x ** 2))
    check("A4: |r|^2 = 1", dot(w0, w0) - 1)

    lam = line_meet_param(G, w0, F, d2)
    check("A5: mu = (1+x^2)^2/(4x(1-x^2))", lam - (1 + x ** 2) ** 2 / (4 * x * (1 - x ** 2)))

    kap = symbols("kap")
    Q1 = addm(G, w0, kap)
    rho3 = refl_pt(F, sub(Q1, F))
    kap_roots = solve(together(cross(sub(rho3(H), G), w0)), kap)
    report("C4 roots k", [simplify(r0) for r0 in kap_roots])
    k_expect = (1 - 3 * x ** 2) / (2 * x)
    k_good = [r0 for r0 in kap_roots if simplify(r0 - k_expect) == 0]
    _require(k_good, f"Expected k root not found: {kap_roots}")
    k = k_good[0]
    K = addm(G, w0, k)
    rho3 = refl_pt(F, sub(K, F))
    t_val = cancel(dot(sub(rho3(H), G), w0))
    check("A6: t = (1-x^2)/x on the forward slice branch", t_val - (1 - x ** 2) / x)

    print("\n======== B. Gap function and SOS monotonicity ========")
    Delta = cancel(k - lam)
    check("B1: Delta = (5x^4-10x^2+1)/(4x(1-x^2))", Delta - (5 * x ** 4 - 10 * x ** 2 + 1) / (4 * x * (1 - x ** 2)))
    dDelta = cancel(diff(Delta, x))
    check(
        "B2: Delta' = -(5x^6-5x^4+7x^2+1)/(4x^2(1-x^2)^2)",
        dDelta + (5 * x ** 6 - 5 * x ** 4 + 7 * x ** 2 + 1) / (4 * x ** 2 * (1 - x ** 2) ** 2),
    )
    check(
        "B3: SOS 5x^6-5x^4+7x^2+1 = 5x^2(x^2-1/2)^2 + (23/4)x^2 + 1",
        expand(5 * x ** 6 - 5 * x ** 4 + 7 * x ** 2 + 1 - (5 * x ** 2 * (x ** 2 - R(1, 2)) ** 2 + R(23, 4) * x ** 2 + 1)),
    )

    print("\n======== C. Core length ========")
    M1 = (-x / 2, R(1, 2))
    M2 = ((h + g) / 2, R(1, 2))
    s12 = cancel(M2[0] - M1[0])
    check("C1: |M1M2| = (1-x^2)/(2x)", s12 - (1 - x ** 2) / (2 * x))
    s23 = line_meet_param(M2, w0, F, sub(K, F))
    check("C2: |M2M3| = (1-x^2)/(2x)", cancel(s23) - (1 - x ** 2) / (2 * x))
    check("C3: core = 1/x - x", cancel(s12 + s23) - (1 / x - x))

    print("\n======== D. Regular pentagon at x0 = sqrt(1 - 2/sqrt(5)) ========")
    x0 = sqrt(1 - 2 / sqrt(5))
    at0 = lambda e: e.subs(x, x0)
    P5 = [E, (at0(h), S(0)), (at0(g), S(1)), (at0(K[0]), at0(K[1])), (-x0, S(1))]
    names = ["E", "H", "G", "K", "F"]
    L2ref = simplify((P5[1][0] - P5[0][0]) ** 2 + (P5[1][1] - P5[0][1]) ** 2)
    check("D0: |EH|^2 = 2 - 2*sqrt(5)/5", L2ref - (2 - 2 * sqrt(5) / 5))
    for i in range(5):
        a_, b_ = P5[i], P5[(i + 1) % 5]
        check(f"D1: |{names[i]}{names[(i + 1) % 5]}|^2 = |EH|^2", (b_[0] - a_[0]) ** 2 + (b_[1] - a_[1]) ** 2 - L2ref)
    cos108 = -(sqrt(5) - 1) / 4
    for i in range(5):
        a_, v_, c_ = P5[(i - 1) % 5], P5[i], P5[(i + 1) % 5]
        dp = (a_[0] - v_[0]) * (c_[0] - v_[0]) + (a_[1] - v_[1]) * (c_[1] - v_[1])
        check(f"D2: interior angle at {names[i]} = 108 deg", dp - L2ref * cos108)
    exitd = refl_lin(sub((at0(K[0]), at0(K[1])), (-x0, S(1))))((at0(w0[0]), at0(w0[1])))
    d2_0 = (at0(d2[0]), at0(d2[1]))
    check("D3: tail angle = 72 deg", dot(exitd, d2_0) - (sqrt(5) - 1) / 4)

    print("\n======== E. Lean linear_combination certificates ========")
    s_, X_ = symbols("s_ X_")
    Gp = 5 * X_ ** 4 - 10 * X_ ** 2 + 1
    Ap = s_ * X_ ** 2 - (s_ - 2)
    Bp = s_ ** 2 - 5
    check("E1: x0_root certificate", expand(Gp - ((s_ * X_ ** 2 - s_ - 2) * Ap - (X_ ** 2 - 1) ** 2 * Bp)))
    check("E2: side_eq_of_contact certificate", expand((1 - 3 * X_ ** 2) ** 2 - (1 + X_ ** 2) * (2 * X_) ** 2 - Gp))
    z = symbols("z")
    check("E3: closure_root_of_unity certificate", expand((z - 1) * (1 + z + z ** 2 + z ** 3 + z ** 4) - (z ** 5 - 1)))
    check("E4: p_pos SOS certificate", expand(5 * X_ ** 3 - 5 * X_ ** 2 + 7 * X_ + 1 - (5 * X_ * (X_ - R(1, 2)) ** 2 + R(23, 4) * X_ + 1)))

    _require(ok_all, "One or more standard-slice exact identities failed")
    print("\nPASS: all standard-slice and regular-pentagon identities verified")


if __name__ == "__main__":
    certify()
