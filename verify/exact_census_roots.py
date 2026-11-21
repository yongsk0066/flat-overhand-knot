"""Arb ball-arithmetic certification of all 72 isolated census roots and 32 inconclusive centres.

For every one of the 72 regular census candidates in verify/uniq_certified.json:
  1. Re-isolates the stratum KKT root in 100-digit Arb ball arithmetic via a
     Banach/Krawczyk contraction box of half-width 2^-66 (plus a wider enclosing
     contraction box containing the saved decimal seed).
  2. Certifies that the recorded disposition (W_k < 0 or lambda(W_8) < 0) holds
     on the isolated root box.
  3. Evaluates the true geometric sliver walls directly on the isolated root box
     and certifies that all 72 regular roots violate a geometric sliver wall with
     upper bound strictly below -0.315 (including row 96, whose archived verdict
     uses lambda(W_8) < 0 and whose true wall W_7 is approximately -2).

For the 32 Krawczyk-inconclusive saved centres, verifies the exact 4 / 6 / 18 / 4
classification at tolerance 10^-6, and certifies that all 104 saved decimal
centres lie in radius-10^-12 boxes violating a true geometric sliver wall below
-0.315. All checks raise RuntimeError unconditionally (safe under python -O).
"""
from collections import Counter
from decimal import Decimal
import itertools
import json
from pathlib import Path
import re
from flint import arb, arb_mat, ctx

_HERE = Path(__file__).resolve().parent
VERIFY_DIR = _HERE if (_HERE / "uniq_certified.json").exists() else Path.cwd() / "verify"
N = 5
UPPER = list(itertools.combinations_with_replacement(range(N), 2))
INDEX = {(i, j): k for k, (i, j) in enumerate(UPPER)}


def _require(cond, msg="verification failed"):
    if not cond:
        raise RuntimeError(msg)


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _scale(a, s):
    return (a[0] * s, a[1] * s)


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def _cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


class Jet:
    """Scalar, five first partials, and upper triangle of the Hessian."""

    def __init__(self, value, gradient=None, hessian=None):
        self.value = value
        self.grad = [0] * N if gradient is None else gradient
        self.hess = [0] * len(UPPER) if hessian is None else hessian

    def __add__(self, o):
        if not isinstance(o, Jet):
            return Jet(self.value + o, self.grad, self.hess)
        return Jet(
            self.value + o.value,
            [a + b for a, b in zip(self.grad, o.grad)],
            [a + b for a, b in zip(self.hess, o.hess)],
        )

    __radd__ = __add__

    def __neg__(self):
        return self * -1

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return -self + o

    def __mul__(self, o):
        if not isinstance(o, Jet):
            return Jet(self.value * o, [a * o for a in self.grad], [a * o for a in self.hess])
        g = [a * o.value + self.value * b for a, b in zip(self.grad, o.grad)]
        h = [
            self.hess[k] * o.value
            + self.value * o.hess[k]
            + self.grad[i] * o.grad[j]
            + self.grad[j] * o.grad[i]
            for k, (i, j) in enumerate(UPPER)
        ]
        return Jet(self.value * o.value, g, h)

    __rmul__ = __mul__


def _poly_walls(theta):
    x, g, h, a, b = [Jet(val, [int(i == j) for j in range(N)]) for i, val in enumerate(theta)]
    q = g - h
    n = q * q + 1
    W = (q * q - 1, 2 * q)
    F = (-x, 1)
    G = (g, 1)
    H = (h, 0)
    U = _add(_scale(_sub(H, G), n), _scale(W, b - a))
    normU = _dot(U, U)
    R = _sub(_scale(U, 2 * _dot(W, U)), _scale(W, normU))
    Q = _add(_scale(G, n), _scale(W, a))
    V4 = _add(_sub(_scale(_sub(H, G), n * normU), _scale(W, a * normU)), _scale(R, a))
    V8 = _add(_scale(_sub(_scale(_sub(F, G), n), _scale(W, a)), normU), _scale(R, a))
    B = _sub(U, _scale(W, b))
    V6 = _sub(_scale(Q, -normU), _sub(_scale(U, 2 * _dot(B, U)), _scale(B, normU)))
    return {
        "W1": 2 * g * x + x * x - 1,
        "W4": -_cross(V4, R),
        "W5": a * (q + x) * (1 - q * x) - x * (g + x) * n,
        "W6": _cross(V6, R),
        "W7": 2 * b * (q + x) * (1 - q * x) - (2 * h * x + x * x + 1) * n,
        "W8": _cross(V8, R) - n * n * normU * normU,
        "W9": _cross(Q, U) - n * _cross(F, U),
        "W10": _cross(Q, F) - _cross(F, U),
    }


def _kkt_system(S, z):
    n = len(z)
    curves = _poly_walls(z[:N])
    walls = [curves[name] for name in S]
    lam = z[N:]
    half = arb(1) / 2
    F = [w.value for w in walls]
    F += [half - sum(l * w.grad[i] for w, l in zip(walls, lam)) for i in range(N)]
    J = [w.grad + [0] * len(S) for w in walls]
    for i in range(N):
        deriv = [
            -sum(l * w.hess[INDEX[min(i, j), max(i, j)]] for w, l in zip(walls, lam))
            for j in range(N)
        ]
        J.append(deriv + [-w.grad[i] for w in walls])
    _require(len(F) == len(J) == n, f"KKT system dimension mismatch for {S}")
    return arb_mat([[v] for v in F]), arb_mat(J)


def _true_geometric_slivers(theta):
    """Evaluate true unit-width geometric sliver walls directly in Arb ball arithmetic."""
    x, g, h, a, b = theta
    q = g - h
    p = q + a - b
    nq = 1 + q * q
    np_ = 1 + p * p
    nx = 1 + x * x
    c5 = (q + x) * (1 - q * x)
    u = ((q * q - 1) / nq, 2 * q / nq)
    d = (((q * q - 1) * (p * p - 1) + 4 * p * q) / (nq * np_), 2 * (p - q) * (1 + p * q) / (nq * np_))
    K = (g + a * u[0], 1 + a * u[1])
    X = (h + b * u[0], b * u[1])
    walls = {
        "W1": (2 * g * x + x * x - 1) / nx,
        "W4": -_cross(_sub((h, arb(0)), K), d),
        "W6": -_cross(X, d),
        "W7": (2 * b * c5 - (2 * h * x + x * x + 1) * nq) / (nq * nx),
        "W8": _cross(_sub((-x, arb(1)), K), d) - 1,
    }
    if c5 > 0 or c5 < 0:
        walls["W5"] = a - x * (g + x) * nq / c5
    return walls


def _isolate_root(S, seed, exponent=66):
    n = len(seed)
    initial = [arb(str(v)).mid() for v in seed]
    z = list(initial)
    for _ in range(5):
        F, J = _kkt_system(S, z)
        dz = J.solve(F)
        z = [(v - dz[i, 0]).mid() for i, v in enumerate(z)]
    F, J = _kkt_system(S, z)
    Y = J.inv().mid()
    I = arb_mat([[int(i == j) for j in range(n)] for i in range(n)])
    A = Y * F

    def certify_radius(rad):
        X = [v + arb(0, rad) for v in z]
        _, JX = _kkt_system(S, X)
        E = I - Y * JX
        rows = [sum(abs(E[i, j]).abs_upper() for j in range(n)) for i in range(n)]
        bounds = [rows[i] + abs(A[i, 0]).abs_upper() / rad for i in range(n)]
        _require(all(row < 1 for row in bounds), f"Krawczyk contraction failed for {S}")
        return X

    small = arb(1) / (2 ** exponent)
    X = certify_radius(small)
    delta = [abs(z[i] - initial[i]) for i in range(n)]
    large = 2 * max([small] + [v.abs_upper() for v in delta])
    wide = certify_radius(large)
    _require(all(wide[i].contains(initial[i]) for i in range(n)), f"Wide box missed seed for {S}")
    return X


def certify():
    ctx.dps = 100
    records = json.loads((VERIFY_DIR / "uniq_certified.json").read_text())
    census = json.loads((VERIFY_DIR / "uniq_census.json").read_text())["census"]
    candidates = [row for row in census if row["kind"] == "other"]
    def key(row):
        return tuple(row["S"]), tuple(row["theta"])
    lookup = {key(row): row for row in candidates}
    record_keys = {key(row) for row in records}
    _require(len(records) == len(record_keys) == len(candidates) == len(lookup) == 104
             and record_keys == lookup.keys(),
             "Census rows must match 104 distinct raw candidate keys exactly")

    threshold = arb("-0.315")
    reg_count = 0
    reg_wall_counts = Counter()
    worst_reg_upper = None
    inc_buckets = Counter()
    saved_centres_ok = 0

    for idx, rec in enumerate(records):
        S = rec["S"]
        theta_f = rec["theta"]
        # Check saved approximate centre at radius 10^-12 for all 104 rows.
        saved_box = [arb(str(v), "1e-12") for v in theta_f]
        saved_walls = _true_geometric_slivers(saved_box)
        _require(
            any(val < threshold for val in saved_walls.values()),
            f"Saved centre row {idx} ({S}) failed geometric sliver < -0.315",
        )
        saved_centres_ok += 1

        if rec["verdict"] == "KRAWCZYK-FAIL":
            xd, gd, hd, ad, bd = [Decimal(repr(v)) for v in theta_f]
            qd = gd - hd
            _require(qd + xd != 0 and Decimal(1) - qd * xd != 0, f"Exact zero denominator at row {idx}")
            tol = Decimal("1e-6")
            if xd < 0:
                inc_buckets["x<0"] += 1
            elif abs(qd + xd) < tol:
                inc_buckets["|q+x|<1e-6"] += 1
            elif abs(Decimal(1) - qd * xd) < tol:
                inc_buckets["|1-qx|<1e-6"] += 1
            else:
                _require(hd < 0 and ad < 0, f"Unexpected inconclusive centre at row {idx}: {theta_f}")
                inc_buckets["other(h,a<0)"] += 1
            continue

        ref = lookup[(tuple(S), tuple(theta_f))]
        X = _isolate_root(S, ref["theta"] + ref["lam"])
        geom_walls = _true_geometric_slivers(X[:N])

        # Verify recorded disposition on the isolated root box.
        m_wall = re.match(r"(W\d+)<0", rec["verdict"])
        if m_wall:
            wname = m_wall.group(1)
            _require(wname in geom_walls and geom_walls[wname] < 0, f"Row {idx} verdict {wname}<0 failed on isolated root")
        else:
            m_lam = re.match(r"λ\((W\d+)\)<0", rec["verdict"])
            _require(m_lam is not None and m_lam.group(1) == "W8", f"Unexpected verdict on row {idx}: {rec['verdict']}")
            _require(X[N + S.index("W8")] < 0, f"Row {idx} lambda(W8)<0 failed on isolated root")

        # Certify that EVERY regular root violates a true geometric sliver wall < -0.315.
        neg_walls = [(nm, val) for nm, val in geom_walls.items() if val < threshold]
        _require(bool(neg_walls), f"Row {idx} ({S}) has no geometric sliver < -0.315 on isolated root box")
        best_nm, best_val = min(neg_walls, key=lambda kv: float(kv[1].mid()))
        upper = float(best_val.mid() + best_val.rad())
        worst_reg_upper = upper if worst_reg_upper is None else max(worst_reg_upper, upper)
        reg_wall_counts[best_nm] += 1
        reg_count += 1

    _require(reg_count == 72, f"Expected 72 certified regular roots, got {reg_count}")
    _require(saved_centres_ok == 104, f"Expected 104 saved centres checked, got {saved_centres_ok}")
    expected_inc = {"x<0": 4, "|q+x|<1e-6": 6, "|1-qx|<1e-6": 18, "other(h,a<0)": 4}
    _require(dict(inc_buckets) == expected_inc, f"Inconclusive bucket mismatch: {dict(inc_buckets)}")

    print(f"[Arb] 72/72 regular census roots freshly isolated (half-width 2^-66) and excluded by true geometric slivers: {dict(reg_wall_counts)}")
    print(f"[Arb] weakest isolated-root geometric-sliver upper bound across all 72 roots: {worst_reg_upper:.10f} < -0.315")
    print(f"[exact] 32 Krawczyk-inconclusive saved centres classified: {dict(inc_buckets)}")
    print(f"[Arb] 104/104 saved decimal centres violate a true geometric sliver (< -0.315) on radius-10^-12 boxes")
    print("PASS: 72 isolated roots excluded, 32 inconclusive centres classified (root existence remains open)")


if __name__ == "__main__":
    certify()
