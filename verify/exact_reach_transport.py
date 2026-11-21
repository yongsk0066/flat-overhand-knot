"""Finite-step algebra and exact real inequalities for transport to the trefoil minimum."""

import sympy as sp
import z3


def identity(label, actual, expected):
    if sp.cancel(actual - expected) != 0:
        raise RuntimeError(label)


def finite_steps():
    q, a, b, m, k, x, z, s = sp.symbols("q a b m k x z s")
    n = 1 + q * q
    C = n / q

    def T(y):
        return (q + y) * (1 - q * y)

    def A(y):
        return n * y / T(y)

    def L(y):
        return y - 2 * q - 1 / y

    def F(y):
        return n * (1 / (y - q) + y / (q * y + 1))

    def R(y):
        return 2 * n * y / (1 + q * y)

    def offset(y):
        return y + 1 / y - C

    def walls(m, k, x, z):
        return k - 2 * q - z - A(x) * (m + x), m - 2 * q - x - A(z) * (k + z)

    identity("angular secant", (A(b) - A(a)) / (b - a), q * n * (1 + a * b) / (T(a) * T(b)))
    identity("angular secant majorant", n * (1 + a * b) - (q + a) * (q + b), (1 - q * a) * (1 - q * b))
    identity("corner angle square", n * (1 + a * a) - 2 * T(a), (q * a + q + a - 1) ** 2)
    identity("reciprocal floor", offset(a), (a - q) * (a - 1 / q) / a)
    identity("coordinatewise angle floor", A(1 / a) * (a + 1 / a), F(a))
    identity("mixed face", C + A(1 / m) * offset(m), R(m))
    identity("mixed switch", R(1 / q), C)
    identity("positive mixed secant", (R(b) - R(a)) / (b - a), 2 * n / ((1 + q * a) * (1 + q * b)))
    identity("mixed lower bound", q * R(a) - n, n * (q * a - 1) / (1 + q * a))
    identity("reciprocal length secant", (L(b) - L(a)) / (b - a), 1 + 1 / (a * b))
    identity(
        "opposite length secant",
        (F(b) - F(a)) / (b - a),
        n * (1 / ((1 + q * a) * (1 + q * b)) - 1 / ((a - q) * (b - q))),
    )
    mixed_left, mixed_right = walls(m, k, 1 / m, R(m) - k)
    identity("mixed left wall", mixed_left, 2 * k - 2 * q - R(m) - F(m))
    identity("mixed right wall", mixed_right, L(m) - A(R(m) - k) * R(m))
    corner_left, corner_right = walls(C - x, C - z, x, z)
    identity("corner left wall", corner_left, 2 * (C / 2 - q - z - (C / 2) * A(x)))
    identity("corner right wall", corner_right, 2 * (C / 2 - q - x - (C / 2) * A(z)))
    identity("symmetric corner product", A(1 / q - s), C / (s * (1 + q * q / (1 - s * q))))
    identity(
        "angle at the symmetric switch", (1 / q - s) - 1 / (s + q), (1 - q * (s + q)) * s / (q * (s + q))
    )
    length_poly = (
        -(q**5)
        - 4 * s * q**4
        - (4 * s * s + 2) * q**3
        - 6 * s * q * q
        + (s**4 - 4 * s * s - 1) * q
        - 2 * s
    )
    d = s + q
    identity("symmetric feasible interval", (L(d) - F(d)) * d * (q * d + 1) * (d - q), length_poly)
    identity(
        "symmetric chord concavity",
        -length_poly.diff(q, 2),
        20 * q**3 + 48 * s * q * q + (24 * s * s + 12) * q + 12 * s,
    )
    print("PASS: transport faces, finite differences, switch values and symmetric length envelope")


def universal_signs():
    count = 0

    def prove(label, assumptions, conclusion):
        nonlocal count
        solver = z3.SolverFor("QF_NRA")
        solver.set(timeout=60000)
        solver.add(*assumptions, z3.Not(conclusion))
        result = solver.check()
        if result != z3.unsat:
            detail = solver.model() if result == z3.sat else solver.reason_unknown()
            raise RuntimeError(f"{label}: {result}, {detail}")
        count += 1

    q, a, b, x, z, m, k, s = z3.Reals("q a b x z m k s")
    u, v, ax, az, px, pz, t = z3.Reals("u v ax az px pz t")
    h, g, rate, lower, upper, gap, step, shift = z3.Reals("h g rate lower upper gap step shift")
    first, second, new_boundary, new_ratio, old_ratio = z3.Reals(
        "first second new_boundary new_ratio old_ratio"
    )
    n = 1 + q * q

    def T(y):
        return (q + y) * (1 - q * y)

    def domain(y):
        return [q > 0, y - q > 0, (y - q) ** 2 > n]

    prove(
        "length domain",
        [q > 0, x > 0, m > 2 * q + x, m * x >= 1],
        z3.And(m > 1, m - q > 0, (m - q) ** 2 > n),
    )
    prove(
        "opposite contraction",
        [
            q > 0,
            x > 0,
            z > 0,
            h > 0,
            g > 0,
            ax > 0,
            az > 0,
            g >= ax * (h + q + x),
            h >= az * (g + q + z),
        ],
        ax * az < 1,
    )
    prove(
        "positive offsets",
        [ax > 0, az > 0, ax * az < 1, v >= ax * u, u >= az * v],
        z3.And(u >= 0, v >= 0),
    )
    prove("strict offsets", [ax > 0, az > 0, ax * az < 1, v > ax * u, u > az * v], z3.And(u > 0, v > 0))
    prove(
        "pole-free floor angle",
        domain(m) + [m * x == 1, ax * T(x) == n * x],
        z3.And(q * x < 1, 0 < ax, ax < 1),
    )
    prove(
        "radial reduction",
        [u >= 0, v >= 0, t >= 0, t <= 1, ax >= px, az >= pz, v >= ax * u, u >= az * v],
        z3.And(t * v >= px * t * u, t * u >= pz * t * v),
    )
    prove(
        "length walls under angle reduction",
        [
            q > 0,
            m > 0,
            k > 0,
            x >= a,
            a > 0,
            z >= b,
            b > 0,
            ax >= px,
            px > 0,
            az >= pz,
            pz > 0,
            k - 2 * q - z >= ax * (m + x),
            m - 2 * q - x >= az * (k + z),
        ],
        z3.And(k - 2 * q - b >= px * (m + a), m - 2 * q - a >= pz * (k + b)),
    )
    prove("ordered floor offsets", [m >= k, k > 1, u * m == m * m + 1, v * k == k * k + 1], u >= v)
    prove(
        "floor and corner thresholds",
        [q > 0, m > q, m > 0, u * q * m == (m - q) * (q * m - 1)],
        z3.And((u >= 0) == (q * m >= 1), (u <= 0) == (m * (n - q * m) >= q)),
    )
    prove("right-first floor", [lower >= v, v >= 0, u >= lower, az > 0, az < 1], u >= az * v)
    prove(
        "left-first floor",
        [ax > 0, u >= 0, v >= ax * u, 0 < az, az * ax < 1],
        z3.And(ax * u >= 0, az * (ax * u) <= u),
    )
    prove("mixed reverse offset", [u >= 0, ax > 0, az > 0, az * ax < 1, v == ax * u], u >= az * v)
    # Here lower=(q+a)(q+b), upper=n(1+ab), t=(1-qa)(1-qb).
    # Their gap is the positive product t by the checked polynomial identity.
    prove(
        "lower angular secant",
        [q > 0, lower > 0, 0 < t, t < 1, upper == lower + t, rate * lower * t == q * upper],
        rate > q,
    )
    prove(
        "corner secant inside unit interval",
        [a > 0, b > 0, a < 1, b < 1],
        (1 + a * a) * (1 + b * b) < 2 * (1 + a * b),
    )
    # Normalising each denominator by n(1+y^2)/2 gives factors in (0,1].
    # The checked square and the preceding inequality bound their product.
    prove(
        "corner scaled secant",
        [
            0 < first,
            first <= 1,
            0 < second,
            second <= 1,
            lower > 0,
            upper > lower,
            rate * lower * first * second == upper,
        ],
        rate > 1,
    )
    prove(
        "one-sided floor secants",
        domain(a)
        + domain(b)
        + [a <= b, first * (1 + q * a) * (1 + q * b) == n, second * (a - q) * (b - q) == n],
        z3.And(first - second > -1, first - second < 1),
    )
    prove(
        "floor contraction and inherited opposite wall",
        [step > 0, lower > 1, upper > -1, upper < 1],
        z3.And(step * (lower + upper) > 0, step * (lower - upper) > 0),
    )
    prove(
        "mixed bounds for the shorter new length",
        [
            q > 0,
            0 < a,
            a < b,
            q * a >= 1,
            (1 + q * a) * new_boundary == 2 * n * a,
            lower * a * b == a * b + 1,
        ],
        z3.And(q * new_boundary >= n, lower <= n),
    )
    prove(
        "two finite mixed improvements",
        [
            q > 0,
            step > 0,
            shift > 0,
            upper > -1,
            rate > q,
            new_boundary > 0,
            q * new_boundary >= lower,
            old_ratio > 0,
            new_ratio > 0,
            old_ratio - new_ratio == rate * (step + shift),
        ],
        z3.And(
            2 * step + shift + upper * step > 0,
            -lower * step + new_boundary * (old_ratio - new_ratio) + old_ratio * shift > 0,
        ),
    )
    prove(
        "mixed stops before a second switch",
        [b > a, a > 1, shift > 0, gap * a * b == (b - a) * (a * b - 1)],
        shift + gap > 0,
    )
    prove("corner angle upper bound", [q > 0, a > 0, q * a < 1, ax * T(a) == n * a, ax < 1], a < 1)
    prove("corner equalisation", [rate > 1, step > 0], step - rate * step < 0)
    prove("symmetric corner domain", [s > 2, q > 0, q * (s + q) <= 1], z3.And(q < 1, s * q < 1))
    prove(
        "symmetric corner factor monotonicity",
        [s > 0, 0 < a, a < b, b < 1, s * b < 1],
        z3.And(b * (1 + a * a) > a * (1 + b * b), b * b * (1 - s * a) > a * a * (1 - s * b)),
    )
    prove(
        "symmetric floor is a single interval in q",
        [s > 0, 0 < a, a <= q, q <= b, a * (s + a) >= 1, s * s > 1 + b * b],
        z3.And(q * (s + q) >= 1, s * s > 1 + q * q),
    )
    length = 8 * s**12 - 88 * s**10 + 267 * s**8 - 456 * s**6 + 88 * s**4 - 80 * s**2 - 16
    low, high = z3.RealVal("2.7496"), z3.RealVal("2.7497")
    qlo, qhi = z3.RealVal("0.35"), z3.RealVal("0.36")
    derivative = -5 * q**4 - 16 * s * q**3 - (12 * s * s + 6) * q * q - 12 * s * q + s**4 - 4 * s * s - 1
    prove("independent length isolation", [s > 0, length == 0], z3.And(s > low, s < high))
    prove(
        "independent middle-tilt bracket",
        [low <= s, s <= high],
        z3.And(z3.substitute(derivative, (q, qlo)) > 0, z3.substitute(derivative, (q, qhi)) < 0),
    )
    prove("terminal region", [s >= low, q >= qlo, q <= qhi], z3.And(q * (s + q) > 1, s * s > 1 + q * q))
    prove(
        "strict crossing orientation",
        [q > 0, x > 0, z > 0, q * x < 1, q * z < 1, T(x) * T(z) > n * n * x * z],
        z3.And(q + x + z - q * x * z > 0, 1 - q * x - q * z - x * z > 0),
    )
    prove(
        "signs at a possible clipped far corner",
        [q > 0, g > 0, z > 0, 2 * q * (g + q + z) < n],
        q * q + 2 * q * z - 1 < 0,
    )
    prove(
        "a clipped negative far wall contradicts the opposite material inequality",
        [
            v < 0,
            -1 < pz,
            pz < 0,
            v + pz * u >= 0,
            az > 0,
            u >= az * v,
            lower + v > 0,
            h > 0,
            h < lower + u,
            h >= az * (lower + v),
        ],
        z3.BoolVal(False),
    )
    prove(
        "opening gives strictly positive side margins",
        [
            q > 0,
            x > 0,
            z > 0,
            h > 0,
            g > 0,
            ax > 0,
            az > 0,
            g >= ax * (h + q + x),
            h >= az * (g + q + z),
        ],
        z3.And(g - ax * h > 0, h - az * g > 0),
    )
    prove("closed crease gains a positive vertical gap", [v >= 0, t > 0, g > 0], v + 2 * t * g > 0)
    print(f"PASS: {count} exact real counterexample queries are unsatisfiable")


if __name__ == "__main__":
    finite_steps()
    universal_signs()
