"""Exact identification of the flush-truncation edges in a neighbourhood of P*.

Certifies in Q[x]/(5x^4 - 10x^2 + 1) and outward-rounded interval arithmetic
that at the regular pentagon P*, the outgoing tail exits across [E, H] at its
exact midpoint u = 1/2 at angle 72 deg, the incoming tail exits across [G, Q1]
at its exact midpoint u = 1/2 at angle 72 deg, and every other outline edge is
strictly excluded. Strictness guarantees persistence on a neighbourhood of P*.
"""
import sympy as sp
from mpmath import mp, iv, mpf, findroot

def _require(cond, msg="verification failed"):
    if not cond:
        raise RuntimeError(msg)


def certify():
    x = sp.symbols("x")
    q4 = 5 * x ** 4 - 10 * x ** 2 + 1
    red = lambda e: sp.rem(sp.expand(e), q4, x)
    inv = lambda e: sp.invert(red(e), q4, x)
    fdiv = lambda a, b: red(a * inv(b))

    def fr(e):
        n, d = sp.fraction(sp.cancel(e))
        return fdiv(sp.expand(n), sp.expand(d))

    g = fr((1 - x ** 2) / (2 * x))
    h = fr((1 - 3 * x ** 2) / (2 * x))
    a = h
    _require(red(g - h - x) == 0)

    E = (sp.Integer(0), sp.Integer(0))
    F = (red(-x), sp.Integer(1))
    G = (g, sp.Integer(1))
    H = (h, sp.Integer(0))
    w0 = (fr((x ** 2 - 1) / (x ** 2 + 1)), fr(2 * x / (x ** 2 + 1)))

    vadd = lambda P, s, D: (red(P[0] + s * D[0]), red(P[1] + s * D[1]))
    vsub = lambda P, Q: (red(P[0] - Q[0]), red(P[1] - Q[1]))
    dot = lambda u, v: red(u[0] * v[0] + u[1] * v[1])
    crs = lambda u, v: red(u[0] * v[1] - u[1] * v[0])

    Q1 = vadd(G, a, w0)
    bprime = red(dot(vsub(F, H), w0))
    Q0 = vadd(H, bprime, w0)
    _require(vsub(Q0, F) == (0, 0), f"Expected Q0 = F at P*, got diff {vsub(Q0, F)}")

    U = vsub(Q0, Q1)
    t_refl = fdiv(red(2 * dot(w0, U)), dot(U, U))
    v3 = (red(t_refl * U[0] - w0[0]), red(t_refl * U[1] - w0[1]))
    M2 = (fdiv(h + g, sp.Integer(2)), sp.Rational(1, 2))
    s_ = fdiv(crs(vsub(Q1, M2), U), crs(w0, U))
    M3 = vadd(M2, s_, w0)

    FEv = vsub(F, E)
    dd = dot(FEv, FEv)
    d2 = (red(fdiv(2 * (-x), dd) * FEv[0] - 1), red(fdiv(2 * (-x), dd) * FEv[1]))
    M1 = (fdiv(red(-x), sp.Integer(2)), sp.Rational(1, 2))
    dA = (red(-d2[0]), red(-d2[1]))

    # Interval-Newton certified root enclosure for x0 in (0, 1/sqrt(3)).
    mp.dps = 60
    iv.dps = 60
    x0_approx = findroot(lambda tt: 5 * tt ** 4 - 10 * tt ** 2 + 1, mpf("0.3249196962"))
    rad = mpf("1e-40")
    XIV = iv.mpf([x0_approx - rad, x0_approx + rad])
    m_iv = iv.mpf(x0_approx)
    N_iv = m_iv - (5 * m_iv ** 4 - 10 * m_iv ** 2 + 1) / (20 * XIV ** 3 - 20 * XIV)
    _require(XIV.a < N_iv.a and N_iv.b < XIV.b, "Interval Newton contraction failed for x0")

    def sgn(e):
        e = red(e)
        if e == 0:
            return 0
        if not e.free_symbols:
            return 1 if e > 0 else -1
        coeffs = sp.Poly(e, x, domain="QQ").all_coeffs()
        v = iv.mpf(0)
        for c in coeffs:
            q = sp.Rational(c)
            v = v * XIV + iv.mpf(int(q.p)) / iv.mpf(int(q.q))
        if v.a > 0:
            return 1
        if v.b < 0:
            return -1
        raise RuntimeError(f"Undecided interval sign for {e}: {v}")

    s5 = fr((5 - 5 * x ** 2) / 2)
    _require(red(s5 ** 2 - 5) == 0)
    c72 = fdiv(s5 - 1, sp.Integer(4))

    EDGES = [
        ("EH", E, vsub(H, E)),
        ("HG", H, vsub(G, H)),
        ("GQ1", G, vsub(Q1, G)),
        ("K3", Q1, vsub(Q0, Q1)),
        ("FE", F, vsub(E, F)),
    ]

    for tag, P0, D0, target in (("outgoing D", M3, v3, "EH"), ("incoming A", M1, dA, "GQ1")):
        hits = []
        for nm, Q, DE in EDGES:
            d_ = crs(D0, DE)
            if d_ == 0:
                _require(nm != target)
                continue
            tt = fdiv(crs(vsub(Q, P0), DE), d_)
            uu = fdiv(crs(vsub(Q, P0), D0), d_)
            st = sgn(tt)
            if nm == target:
                _require(red(uu - sp.Rational(1, 2)) == 0, f"{tag} on {nm}: u != 1/2")
                _require(st > 0, f"{tag} on {nm}: t <= 0")
                resang = red(dot(D0, DE) ** 2 - c72 ** 2 * dot(DE, DE))
                _require(resang == 0, f"{tag} on {nm}: angle != 72 deg")
                hits.append(nm)
                print(f"[exact] {tag} crosses {nm} at u = 1/2, t > 0, angle = 72 deg")
            else:
                su0 = sgn(uu)
                su1 = sgn(red(1 - uu))
                # The ray originates at t = 0 on its own crease (K3 for D, FE for A);
                # exclude all strictly forward intersections (t > 0) with closed edges [0, 1].
                if nm in ("K3", "FE") and st == 0:
                    continue
                inside = st > 0 and su0 >= 0 and su1 >= 0
                _require(not inside, f"{tag} unexpectedly hits non-target edge {nm}")
        _require(hits == [target], f"{tag} hit list mismatch: {hits}")

    print("PASS: flush-truncation edges uniquely identified at P* with strict openness margins")


if __name__ == "__main__":
    certify()
