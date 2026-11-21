"""Exact certificates for local vs. global layer order on the three-crease branch.

1. Cyclic local-order counterexample (exact rational arithmetic):
   theta = (1/3, 301/3, 100, 100, 302/3) realizes the alternating trefoil
   crossing data XYZXYZ with all eight walls strictly positive, no triple
   overlap of the infinite width-one bands (minimum Euclidean separation 851/40),
   and actual panel overlap intervals separated by at least 21275/936 widths.
   Assigning A > B near crease 1 yields a locally consistent cyclic order
   A > B > D > A that admits an explicit strong epsilon-realization for all
   0 < epsilon <= 1/100 via a smooth height ramp on the free incoming tail A.

2. Common four-pass overlap lemma (exact algebra + Sturm + outward-rounded intervals):
   - On the entire standard slice 0 < x < 1/sqrt(3), the point (h(x)/2, 3/4)
     lies strictly in the interior of all four actual panels A, B, C, D.
   - At P*, at theta-dagger, throughout a 10^-3 coordinate box around each,
     and along the entire certified second-family arc [x-dagger, x0] (plus a
     10^-4 coordinate neighbourhood), the fixed point (1/2, 3/4) lies strictly
     in the interior of all four actual panels A, B, C, D.
   Consequently, local layer order together with connectedness of convex
   pairwise overlaps forces the global total order B > D > A > C throughout
   these loci.
"""
from fractions import Fraction as F
import sympy as sp
from sympy import Rational as R
from mpmath import mp, iv, mpf, findroot


def _require(cond, msg="verification failed"):
    if not cond:
        raise RuntimeError(msg)


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _mul(a, s):
    return (a[0] * s, a[1] * s)


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def _cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def _refl(axis, v):
    s = 2 * _dot(axis, v) / _dot(axis, axis)
    return _sub(_mul(axis, s), v)


def _panel_halfplanes(theta):
    """Return the 14 inward linear half-plane functionals (origin, dir) -> value > 0
    for the four actual obliquely crease-bounded passes A, B, C, D.

    Each entry is (panel_name, label, value_at_P(P), scale_norm_sq) where value > 0
    means strictly on the pass-interior side of that boundary line.
    """
    x, g, h, a, b = theta
    E = (0, 0)
    Fp = (-x, 1)
    G = (g, 1)
    H = (h, 0)
    dA = _refl(_sub(Fp, E), (1, 0))
    dB = (1, 0)
    dC = _refl(_sub(G, H), (1, 0))
    Q1 = _add(G, _mul(dC, a))
    Q0 = _add(H, _mul(dC, b))
    U = _sub(Q0, Q1)
    dD = _refl(U, dC)
    r3G = _add(Q1, _refl(U, _sub(G, Q1)))
    # Inward signed distances (up to positive chord-length scaling on crease lines):
    # A: band 0 < cross(P-E, dA) < 1, behind crease EF (cross(F-E, P-E) > 0 since M1-dA gives > 0).
    # B: 0 < y < 1, between EF (cross(F-E, P-E) < 0) and HG (cross(G-H, P-H) > 0).
    # C: band 0 < cross(dC, P-H) < 1? Check: G-H=(q,1), cross(dC, G-H) = ((q^2-1) - 2q^2)/(q^2+1) = -1,
    #    so band of C is 0 < cross(P-H, dC) < 1, between HG (cross(G-H, P-H) < 0) and Q0Q1 (cross(Q1-Q0, P-Q0) < 0).
    # D: band 0 < cross(P-r3G, dD) < 1, ahead of Q0Q1 (cross(Q1-Q0, P-Q0) > 0).
    FE = _sub(Fp, E)
    HG = _sub(G, H)
    Q0Q1 = _sub(Q1, Q0)

    def eval_all(P):
        cA = _cross(_sub(P, E), dA)
        cC = _cross(_sub(P, H), dC)
        cD = _cross(_sub(P, r3G), dD)
        s1 = _cross(FE, _sub(P, E))
        s2 = _cross(HG, _sub(P, H))
        s3 = _cross(Q0Q1, _sub(P, Q0))
        return [
            ("A", "e_A", cA),
            ("A", "e_B", 1 - cA),
            ("A", "k1", -s1),
            ("B", "y=0", P[1]),
            ("B", "y=1", 1 - P[1]),
            ("B", "k1", -s1),
            ("B", "k2", s2),
            ("C", "eC0", cC),
            ("C", "l2", 1 - cC),
            ("C", "k2", s2),
            ("C", "k3", -s3),
            ("D", "eD1", cD),
            ("D", "eD0", 1 - cD),
            ("D", "k3", -s3),
        ]

    return eval_all, dict(E=E, F=Fp, G=G, H=H, Q1=Q1, Q0=Q0, dA=dA, dB=dB, dC=dC, dD=dD, r3G=r3G)


def certify_counterexample():
    theta = (F(1, 3), F(301, 3), F(100), F(100), F(302, 3))
    x, g, h, a, b = theta
    q = g - h
    p = q + a - b
    _require(x == F(1, 3) and q == F(1, 3) and p == F(-1, 3))

    eval_all, geom = _panel_halfplanes(theta)
    E, Fp, G, H = geom["E"], geom["F"], geom["G"], geom["H"]
    Q1, Q0 = geom["Q1"], geom["Q0"]
    dA, dB, dC, dD = geom["dA"], geom["dB"], geom["dC"], geom["dD"]
    r3G = geom["r3G"]

    _require(dA == (F(-4, 5), F(-3, 5)))
    _require(dB == (F(1), F(0)))
    _require(dC == (F(-4, 5), F(3, 5)))
    _require(dD == (F(7, 25), F(-24, 25)))
    _require(_cross(dA, dB) == F(3, 5))
    _require(_cross(dB, dC) == F(3, 5))
    _require(_cross(dC, dD) == F(3, 5))

    # Verify the three non-adjacent centreline crossings and passage order XYZXYZ.
    M1 = (-x / 2, F(1, 2))
    M2 = ((g + h) / 2, F(1, 2))
    M3 = ((Q0[0] + Q1[0]) / 2, (Q0[1] + Q1[1]) / 2)
    L1 = _dot(_sub(M2, M1), dB)
    L2 = _dot(_sub(M3, M2), dC)
    _require(L1 == F(301, 3) and L2 == F(301, 3))

    def intersect_lines(P1, v1, P2, v2):
        det = _cross(v1, v2)
        dp = _sub(P2, P1)
        return _cross(dp, v2) / det, _cross(dp, v1) / det

    tX_A, tX_C = intersect_lines(M1, dA, M2, dC)
    tY_A, tY_D = intersect_lines(M1, dA, M3, dD)
    tZ_B, tZ_D = intersect_lines(M1, dB, M3, dD)
    _require((tX_A, tX_C) == (F(-1505, 24), F(1505, 24)))
    _require((tY_A, tY_D) == (F(-1505, 39), F(1505, 39)))
    _require((tZ_B, tZ_D) == (F(301, 8), F(1505, 24)))
    # Passage order along A: X before Y (tX_A < tY_A < 0); along D: Y before Z (0 < tY_D < tZ_D).
    _require(tX_A < tY_A < 0)
    _require(0 < tZ_B < L1)
    _require(0 < tX_C < L2)
    _require(0 < tY_D < tZ_D)

    # Verify all eight chart walls are strictly positive.
    U = _sub(Q0, Q1)
    r3H = _add(Q1, _refl(U, _sub(H, Q1)))
    den9 = _cross(_sub(Fp, E), U)
    walls = [
        -_cross(G, dA),
        -_cross(_sub(H, r3G), dD),
        a - _cross(_sub(Fp, G), dA) / _cross(dC, dA),
        _cross(_sub(E, r3H), dD),
        _cross(Q0, dA) - 1,
        _cross(_sub(Fp, r3G), dD) - 1,
        _cross(Q1, U) / den9 - 1,
        _cross(Q1, Fp) / den9 - 1,
    ]
    _require(walls == [F(297, 5), F(297, 5), F(445, 12), F(178, 5), F(891, 25), F(891, 25), F(297, 8), F(297, 8)])
    _require(all(w > 0 for w in walls))

    # Check empty triple intersections of the infinite width-one support bands.
    bands = {"A": (M1, dA), "B": (M1, dB), "C": (M2, dC), "D": (M3, dD)}

    def band_corners(k1, k2):
        O1, v1 = bands[k1]
        O2, v2 = bands[k2]
        det = _cross(v1, v2)
        pts = []
        for s1 in (F(-1, 2), F(1, 2)):
            for s2 in (F(-1, 2), F(1, 2)):
                r1 = s1 + _cross(v1, O1)
                r2 = s2 + _cross(v2, O2)
                pts.append(((r1 * v2[0] - v1[0] * r2) / det, (-v1[1] * r2 + r1 * v2[1]) / det))
        return pts

    min_triple_gap = None
    for k1, k2, k3 in [("A", "B", "C"), ("A", "B", "D"), ("A", "C", "D"), ("B", "C", "D")]:
        O3, v3 = bands[k3]
        vals = [_cross(v3, _sub(pt, O3)) for pt in band_corners(k1, k2)]
        lo, hi = min(vals), max(vals)
        gap = (lo - F(1, 2)) if lo > F(1, 2) else (F(-1, 2) - hi)
        _require(gap > 0, (k1, k2, k3, lo, hi))
        min_triple_gap = gap if min_triple_gap is None else min(min_triple_gap, gap)
    _require(min_triple_gap == F(851, 40))

    # Verify longitudinal separation of pairwise overlaps along each support band
    # (which lower-bounds the separation of the clipped actual-panel overlaps).
    min_long_gap = None
    for k in ("A", "B", "C", "D"):
        Ok, vk = bands[k]
        intervals = []
        for other in ("A", "B", "C", "D"):
            if other == k:
                continue
            ts = [_dot(_sub(pt, Ok), vk) for pt in band_corners(k, other)]
            intervals.append((min(ts), max(ts), other))
        intervals.sort()
        for (a0, b0, _), (a1, b1, _) in zip(intervals, intervals[1:]):
            gap = a1 - b0
            _require(gap > 22)
            min_long_gap = gap if min_long_gap is None else min(min_long_gap, gap)
    _require(min_long_gap == F(21275, 936))

    # Verify free-tail ramp bounds (A transitions on s in [-25, -15] between eps/4 and eps):
    # A's overlap with D ends at -8861/234 < -37.86 and with B starts at -3/2 = -1.5,
    # leaving > 12 strip widths on both sides of [-25, -15].
    O_A, v_A = bands["A"]
    ad_max = max(_dot(_sub(pt, O_A), v_A) for pt in band_corners("A", "D"))
    ab_min = min(_dot(_sub(pt, O_A), v_A) for pt in band_corners("A", "B"))
    _require(ad_max < F(-37) and ab_min == F(-3, 2), f"Free-tail clearance failed: {ad_max}, {ab_min}")

    # Verify the smooth C^2 isometric height ramp p(t) = 10t^3 - 15t^4 + 6t^5 on [-25, -15]:
    # p(0)=0, p(1)=1, p'(0)=p'(1)=p''(0)=p''(1)=0, max_{[0,1]} p'(t) = p'(1/2) = 15/8,
    # so |z'(s)| = (3*eps/40)*p'(t) <= (9/64)*eps < 1 for 0 < eps <= 1/100, and the total
    # horizontal foreshortening int_{-25}^{-15} (1 - sqrt(1 - z'^2)) ds <= int z'^2 ds <= (405/2048)*eps^2.
    t_sym = sp.symbols("t", real=True)
    p_ramp = 10 * t_sym ** 3 - 15 * t_sym ** 4 + 6 * t_sym ** 5
    dp = sp.diff(p_ramp, t_sym)
    d2p = sp.diff(dp, t_sym)
    _require(p_ramp.subs(t_sym, 0) == 0 and p_ramp.subs(t_sym, 1) == 1, "Ramp endpoint values wrong")
    _require(dp.subs(t_sym, 0) == 0 and dp.subs(t_sym, 1) == 0, "Ramp C^1 endpoint derivatives wrong")
    _require(d2p.subs(t_sym, 0) == 0 and d2p.subs(t_sym, 1) == 0, "Ramp C^2 endpoint derivatives wrong")
    _require(sp.simplify(dp - 30 * t_sym ** 2 * (1 - t_sym) ** 2) == 0, "Ramp derivative factorization wrong")
    _require(dp.subs(t_sym, R(1, 2)) == R(15, 8), "Ramp peak derivative wrong")
    # 10 * (3/40)^2 * (15/8)^2 = 405/2048; exact int_0^1 10*(3/40*p')^2 dt = 9/112 < 405/2048.
    exact_int = sp.integrate(10 * (R(3, 40) * dp) ** 2, (t_sym, 0, 1))
    _require(exact_int == R(9, 112) and exact_int < R(405, 2048), f"Ramp foreshortening bound wrong: {exact_int}")
    print(f"[exact] counterexample (1/3,301/3,100,100,302/3): triple gap = {min_triple_gap}, longitudinal gap = {min_long_gap}, ramp int(z'^2)/eps^2 = {exact_int} < 405/2048")


def certify_standard_slice_overlap():
    x = sp.symbols("x", positive=True)
    theta = (x, (1 - x ** 2) / (2 * x), (1 - 3 * x ** 2) / (2 * x), (1 - 3 * x ** 2) / (2 * x), (1 + x ** 2) / (2 * x))
    eval_all, _ = _panel_halfplanes(theta)
    P = (theta[2] / 2, R(3, 4))
    # Every factorized functional must be strictly positive on 0 < x^2 < 1/3.
    expected = {
        ("A", "e_A"): (3 * x ** 2 + 1) / (4 * (x ** 2 + 1)),
        ("A", "e_B"): (x ** 2 + 3) / (4 * (x ** 2 + 1)),
        ("A", "k1"): 1 / (4 * x),
        ("B", "y=0"): R(3, 4),
        ("B", "y=1"): R(1, 4),
        ("B", "k1"): 1 / (4 * x),
        ("B", "k2"): 1 / (4 * x),
        ("C", "eC0"): (3 * x ** 2 + 1) / (4 * (x ** 2 + 1)),
        ("C", "l2"): (x ** 2 + 3) / (4 * (x ** 2 + 1)),
        ("C", "k2"): 1 / (4 * x),
        ("C", "k3"): (1 + x ** 2 - 4 * x ** 4) / (4 * x * (x ** 2 + 1)),
        ("D", "eD1"): (3 * x ** 2 - 2 * x + 1) * (3 * x ** 2 + 2 * x + 1) / (4 * (x ** 2 + 1) ** 2),
        ("D", "eD0"): (3 + 6 * x ** 2 - 5 * x ** 4) / (4 * (x ** 2 + 1) ** 2),
        ("D", "k3"): (1 + x ** 2 - 4 * x ** 4) / (4 * x * (x ** 2 + 1)),
    }
    for panel, lbl, val in eval_all(P):
        _require(sp.simplify(val - expected[(panel, lbl)]) == 0, f"Slice functional mismatch {(panel, lbl, val)}")
    # Exact closed-interval certification on y = x^2 in [0, 1/3] with zero endpoint gap:
    # 1 - 3y vanishes only at y = 1/3 and is strictly positive on [0, 1/3);
    # 3x^2 +- 2x + 1 = 3(x +- 1/3)^2 + 2/3 > 0 identically;
    # 1 + y - 4y^2 and 3 + 6y - 5y^2 are strictly positive at 0 and 1/3 with zero Sturm roots on [0, 1/3].
    _require(sp.expand(3 * x ** 2 - 2 * x + 1 - (3 * (x - R(1, 3)) ** 2 + R(2, 3))) == 0, "SOS identity 1 failed")
    _require(sp.expand(3 * x ** 2 + 2 * x + 1 - (3 * (x + R(1, 3)) ** 2 + R(2, 3))) == 0, "SOS identity 2 failed")
    y = sp.symbols("y")
    for poly_y in (1 + y - 4 * y ** 2, 3 + 6 * y - 5 * y ** 2):
        _require(poly_y.subs(y, R(0)) > 0 and poly_y.subs(y, R(1, 3)) > 0, f"Endpoint sign failed for {poly_y}")
        p_ = sp.Poly(poly_y, y, domain="QQ")
        chain = sp.sturm(p_)

        def sc(t):
            vals = [c.eval(t) for c in chain if c.eval(t) != 0]
            return sum(1 for a_, b_ in zip(vals, vals[1:]) if (a_ > 0) != (b_ > 0))

        _require(sc(R(0)) - sc(R(1, 3)) == 0, f"Sturm root found on [0, 1/3] for {poly_y}")
    print("[exact] standard slice: (h(x)/2, 3/4) lies strictly inside all four panels for all 0 < x < 1/sqrt(3)")


def certify_second_family_and_nbhds():
    # 1. Exact Sturm certificate that the branch h(x) of N5(x,h)=0 stays in (1.035, 1.053)
    #    over [x-dagger, x0] subset [0.322, 0.325].
    x, h = sp.symbols("x h")
    N5 = (
        8 * h ** 3 * x ** 3
        + 4 * h ** 2 * x ** 4
        + 4 * h ** 2 * x ** 2
        + 2 * h * x ** 5
        - 4 * h * x ** 3
        - 6 * h * x
        + x ** 6
        + 3 * x ** 4
        + 3 * x ** 2
        + 1
    )

    def sturm_count(poly_x, lo, hi):
        chain = sp.sturm(sp.Poly(poly_x, x, domain="QQ"))
        def sgnchg(t):
            vals = [c.eval(t) for c in chain]
            vals = [v for v in vals if v != 0]
            return sum(1 for a_, b_ in zip(vals, vals[1:]) if (a_ > 0) != (b_ > 0))
        return sgnchg(lo) - sgnchg(hi)

    lo_x, hi_x = R(322, 1000), R(325, 1000)
    for h_bound in (R(1035, 1000), R(1053, 1000)):
        p_b = sp.expand(N5.subs(h, h_bound))
        _require(p_b.subs(x, lo_x) != 0 and p_b.subs(x, hi_x) != 0)
        _require(sturm_count(p_b, lo_x, hi_x) == 0, f"N5 crossed h={h_bound} on [0.322, 0.325]")

    # Confirm endpoints lie strictly inside (1.035, 1.053) and certify root enclosures.
    mp.dps = 60
    iv.dps = 60
    x0_approx = findroot(lambda t: 5 * t ** 4 - 10 * t ** 2 + 1, mpf("0.3249196962"))
    xd_approx = findroot(lambda t: t ** 12 + 4 * t ** 10 - 7 * t ** 4 + 20 * t ** 2 - 2, mpf("0.3221274973"))
    rad = mpf("1e-40")
    X0_iv = iv.mpf([x0_approx - rad, x0_approx + rad])
    XD_iv = iv.mpf([xd_approx - rad, xd_approx + rad])
    # Interval Newton verification for both roots.
    m0 = iv.mpf(x0_approx)
    N0 = m0 - (5 * m0 ** 4 - 10 * m0 ** 2 + 1) / (20 * X0_iv ** 3 - 20 * X0_iv)
    _require(X0_iv.a < N0.a and N0.b < X0_iv.b)
    md = iv.mpf(xd_approx)
    Nd = md - (md ** 12 + 4 * md ** 10 - 7 * md ** 4 + 20 * md ** 2 - 2) / (
        12 * XD_iv ** 11 + 40 * XD_iv ** 9 - 28 * XD_iv ** 3 + 40 * XD_iv
    )
    _require(XD_iv.a < Nd.a and Nd.b < XD_iv.b)
    _require(float(XD_iv.a) > 0.322 and float(X0_iv.b) < 0.325)

    h0_iv = (1 - 3 * X0_iv ** 2) / (2 * X0_iv)
    hd_iv = (2 * XD_iv ** 11 + 8 * XD_iv ** 9 - XD_iv ** 7 - 4 * XD_iv ** 5 - 13 * XD_iv ** 3 + 40 * XD_iv) / 12
    _require(float(hd_iv.a) > 1.035 and float(h0_iv.b) < 1.053)

    P_iv = (iv.mpf(1) / 2, iv.mpf(3) / 4)

    def check_box(theta_iv, min_required=0):
        eval_all, _ = _panel_halfplanes(theta_iv)
        worst = None
        for panel, lbl, val in eval_all(P_iv):
            _require(val.a > min_required, f"Box failed {panel}/{lbl}: {val}")
            worst = float(val.a) if worst is None else min(worst, float(val.a))
        return worst

    # 2. Certified 10^-3 neighbourhoods of P* and theta-dagger in all five coordinates.
    eps_pt = iv.mpf(["-0.001", "0.001"])
    P_star_iv = (
        X0_iv + eps_pt,
        (1 - X0_iv ** 2) / (2 * X0_iv) + eps_pt,
        h0_iv + eps_pt,
        h0_iv + eps_pt,
        (1 + X0_iv ** 2) / (2 * X0_iv) + eps_pt,
    )
    D_star_iv = (
        XD_iv + eps_pt,
        (1 - XD_iv ** 2) / (2 * XD_iv) + eps_pt,
        hd_iv + eps_pt,
        hd_iv + eps_pt,
        (1 + XD_iv ** 2) / (2 * XD_iv) + eps_pt,
    )
    m_p = check_box(P_star_iv, min_required=0.05)
    m_d = check_box(D_star_iv, min_required=0.05)
    print(f"[interval] (1/2, 3/4) strictly inside all four panels on [P* +- 10^-3] (min functional {m_p:.4f}) and [theta-dagger +- 10^-3] (min {m_d:.4f})")

    # 3. Entire second-family rectangle [0.322, 0.325] x [1.035, 1.053] plus 10^-4 in all 5 coords.
    eps_arc = iv.mpf(["-0.0001", "0.0001"])
    worst_arc = None
    nx, nh = 16, 12
    for ix in range(nx):
        xa = iv.mpf(322) / 1000 + iv.mpf(3 * ix) / (1000 * nx)
        xb = iv.mpf(322) / 1000 + iv.mpf(3 * (ix + 1)) / (1000 * nx)
        X_sub = iv.mpf([xa.a, xb.b])
        for ih in range(nh):
            ha = iv.mpf(1035) / 1000 + iv.mpf(18 * ih) / (1000 * nh)
            hb = iv.mpf(1035) / 1000 + iv.mpf(18 * (ih + 1)) / (1000 * nh)
            H_sub = iv.mpf([ha.a, hb.b])
            th_sub = (
                X_sub + eps_arc,
                (1 - X_sub ** 2) / (2 * X_sub) + eps_arc,
                H_sub + eps_arc,
                H_sub + eps_arc,
                (1 + X_sub ** 2) / (2 * X_sub) + eps_arc,
            )
            w = check_box(th_sub, min_required=0.05)
            worst_arc = w if worst_arc is None else min(worst_arc, w)
    print(f"[interval] (1/2, 3/4) strictly inside all four panels on entire second family + 10^-4 tube ({nx*nh} boxes, min functional {worst_arc:.4f} > 0.05)")
    print("PASS: cyclic counterexample and common four-pass overlap lemma certified")


def certify():
    certify_counterexample()
    certify_standard_slice_overlap()
    certify_second_family_and_nbhds()


if __name__ == "__main__":
    certify()
