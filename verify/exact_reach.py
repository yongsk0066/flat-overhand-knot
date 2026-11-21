"""Certify monotone descent on the second family (Theorem C).

The two wall numerators vanish at the pentagon. The core derivative is
Aell/(2*x**2*N5_h), with Aell and N5_h both negative in the interior.
"""

import sympy as sp
from mpmath import iv

from ribbon_geometry import family_polynomial

R = sp.Rational
x, h = sp.symbols("x h")
X_INTERVAL = (R(322, 1000), R(325, 1000))
HEIGHT_INTERVAL = (R(1035, 1000), R(1053, 1000))
SAMPLE = R(1293, 4000)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def polynomial(expr):
    return sp.Poly(expr, x, domain="QQ")


def root_count(expr, bounds):
    poly = polynomial(expr)
    require(not poly.is_zero and all(poly.eval(t) != 0 for t in bounds),
            "Sturm interval has a zero polynomial or a root at its boundary")
    return poly.count_roots(*bounds)


def enclosure(value):
    value = R(value)
    return iv.mpf(int(value.p)) / int(value.q)


def horner(poly, argument):
    value = iv.mpf(0)
    for coefficient in poly.all_coeffs():
        value = value * argument + enclosure(coefficient)
    return value


def isolate_root(poly, bounds):
    roots = poly.intervals(eps=R(1, 10**35), inf=bounds[0], sup=bounds[1])
    require(len(roots) == 1 and roots[0][1] == 1, "expected one simple root in the stated interval")
    return roots[0][0]


def check_endpoint(name, modulus, x_box, height, vanishing, nonzero):
    require(root_count(modulus, X_INTERVAL) == root_count(modulus, x_box) == 1,
            f"{name}: endpoint root isolation failed")
    X = iv.mpf([enclosure(x_box[0]).a, enclosure(x_box[1]).b])
    numerator, denominator = sp.fraction(sp.cancel(height))
    D = horner(polynomial(denominator), X)
    require(0 not in D, f"{name}: height denominator contains zero")
    H = horner(polynomial(numerator), X) / D
    require(H.a > enclosure(HEIGHT_INTERVAL[0]) and H.b < enclosure(HEIGHT_INTERVAL[1]),
            f"{name}: wrong height branch")
    for label, expr in {**vanishing, **nonzero}.items():
        num, den = sp.fraction(sp.cancel(expr.subs(h, height)))
        require(sp.gcd(den, modulus, x) == 1, f"{name}: singular {label}")
        zero = sp.rem(num, modulus, x) == 0
        require(zero == (label in vanishing), f"{name}: endpoint identity for {label}")
    print(f"[endpoint] {name}: isolated in {x_box}, exact endpoint identities verified")


def certify():
    iv.dps = 60
    g = (1 - x**2) / (2*x)
    b = (1 + x**2) / (2*x)
    q = g - h
    N5 = family_polynomial(x, h)
    N5h = sp.diff(N5, h)
    Aell = sp.expand(-2*x**2*sp.diff(N5, x) + (x**2 - 1)*N5h)

    W7 = 2*b*(q + x)*(1 - q*x)/(q**2 + 1) - 2*h*x - (x**2 + 1)
    P7n, P7d = sp.fraction(sp.cancel(W7))
    require(sp.cancel(P7d - 4*x**2*(1 + q**2)) == 0, "w7 denominator identity failed")
    # The reflected unit vector is rational. With F=(-x,1),
    # t9 - 1 = det(Q1-F,Q0-Q1) / det(F,Q0-Q1).
    direction = sp.Matrix([(q**2 - 1)/(1 + q**2), 2*q/(1 + q**2)])
    F = sp.Matrix([-x, 1])
    Q1 = sp.Matrix([g, 1]) + h*direction
    Q0 = sp.Matrix([h, 0]) + b*direction
    U = Q0 - Q1
    W9 = (Q1 - F).row_join(U).det() / F.row_join(U).det()
    P9n, P9d = sp.fraction(sp.cancel(W9))

    dagger = x**12 + 4*x**10 - 7*x**4 + 20*x**2 - 2
    pentagon = 5*x**4 - 10*x**2 + 1
    hd = (2*x**11 + 8*x**9 - x**7 - 4*x**5 - 13*x**3 + 40*x)/12
    hp = (1 - 3*x**2)/(2*x)
    check_endpoint("dagger", dagger, (R("0.32212749"), R("0.32212750")), hd,
                   {"N5": N5, "Aell": Aell}, {"N5_h": N5h, "w7": P7n, "w9": P9n, "w9 denominator": P9d})
    check_endpoint("pentagon", pentagon, (R("0.32491969"), R("0.32491970")), hp,
                   {"N5": N5, "w7": P7n, "w9": P9n}, {"Aell": Aell, "N5_h": N5h, "w9 denominator": P9d})

    for bound in HEIGHT_INTERVAL:
        require(root_count(N5.subs(h, bound), X_INTERVAL) == 0, "height-strip boundary crossed")
    signs = {"P7n": (P7n, 1), "P9n": (P9n, 1), "P9d": (P9d, 1),
             "Aell": (Aell, -1), "N5_h": (N5h, -1)}
    endpoint_factors = {"P7n": pentagon, "P9n": pentagon, "Aell": dagger}
    for name, (expr, _) in signs.items():
        resultant = polynomial(sp.resultant(N5, expr, h))
        if name in endpoint_factors:
            resultant, rem = resultant.div(polynomial(endpoint_factors[name]))
            require(rem.is_zero, f"{name}: expected endpoint factor missing")
        require(root_count(resultant, X_INTERVAL) == 0, f"{name}: unexpected resultant root")
        print(f"[Sturm] {name}: zero residual roots on [0.322, 0.325]")

    # Constant root count in the strip follows from its root-free boundary,
    # nonzero leading coefficient and the nonvanishing discriminant above.
    sample_poly = sp.Poly(N5.subs(x, SAMPLE), h, domain="QQ")
    require(sample_poly.count_roots(*HEIGHT_INTERVAL) == 1, "sample is not a unique strip root")
    low, high = isolate_root(sample_poly, HEIGHT_INTERVAL)
    center = (low + high)/2
    radius = R(1, 10**20)
    H = iv.mpf([enclosure(center - radius).a, enclosure(center + radius).b])
    mid = enclosure(center)
    derivative = horner(sample_poly.diff(), H)
    require(0 not in derivative, "sample derivative contains zero")
    image = mid - horner(sample_poly, mid)/derivative
    require(H.a < image.a and image.b < H.b, "sample interval Newton inclusion failed")

    for name, (expr, sign) in signs.items():
        value = horner(sp.Poly(expr.subs(x, SAMPLE), h, domain="QQ"), H)
        require((sign*value).a > 0, f"{name}: wrong interior sign")
        print(f"[interval] {name}: interior sign {sign:+d}")
    print("PASS: the two wall pairs are positive away from the pentagon and core length strictly decreases from the pentagon to dagger")


if __name__ == "__main__":
    certify()
