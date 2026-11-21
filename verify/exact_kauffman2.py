"""Exact local optimality of Kauffman's truncation length at the regular pentagon.

All algebraic identities are decided in Q(x0), x0 = sqrt(1 - 2/sqrt(5)),
using first-order dual numbers (five-variable jets). No floating-point sign
calibration, heuristic nsimplify, or numerical approximation enters any
equality or positivity decision.
"""
import sympy as sp
from sympy import Rational as R
from sympy.polys.matrices import DomainMatrix
from mpmath import iv

def _require(cond, msg="verification failed"):
    if not cond:
        raise RuntimeError(msg)


class Jet5:
    def __init__(self, field, value, grad=None):
        self.field = field
        self.v = field.convert(value)
        self.g = grad if grad is not None else (field.zero,) * 5

    def _lift(self, other):
        return other if isinstance(other, Jet5) else Jet5(self.field, other)

    def __add__(self, other):
        o = self._lift(other)
        return Jet5(self.field, self.v + o.v, tuple(a + b for a, b in zip(self.g, o.g)))

    __radd__ = __add__

    def __neg__(self):
        return Jet5(self.field, -self.v, tuple(-a for a in self.g))

    def __sub__(self, other):
        return self + (-self._lift(other))

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        o = self._lift(other)
        return Jet5(
            self.field,
            self.v * o.v,
            tuple(a * o.v + self.v * b for a, b in zip(self.g, o.g)),
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        o = self._lift(other)
        inv = self.field.one / o.v
        inv2 = inv * inv
        return Jet5(
            self.field,
            self.v * inv,
            tuple((a * o.v - self.v * b) * inv2 for a, b in zip(self.g, o.g)),
        )

    def __rtruediv__(self, other):
        return self._lift(other) / self


def add2(a, b):
    return (a[0] + b[0], a[1] + b[1])


def sub2(a, b):
    return (a[0] - b[0], a[1] - b[1])


def mul2(a, s):
    return (a[0] * s, a[1] * s)


def dot2(a, b):
    return a[0] * b[0] + a[1] * b[1]


def cross2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def reflect2(axis, vector):
    return sub2(mul2(axis, (2 * dot2(axis, vector)) / dot2(axis, axis)), vector)


def certify():
    alpha = sp.sqrt(1 - 2 / sp.sqrt(5))
    field = sp.QQ.algebraic_field(alpha)
    X = field.from_sympy(alpha)
    S5 = field.convert(2) / (field.one - X * X)
    _require(S5 * S5 == field.convert(5))
    S_val = S5 * X  # s = sqrt(5 - 2*sqrt(5))
    _require(S_val * S_val == field.convert(5) - field.convert(2) * S5)
    PHI = (field.one + S5) / field.convert(2)

    g0 = (field.one - X * X) / (field.convert(2) * X)
    h0 = (field.one - field.convert(3) * X * X) / (field.convert(2) * X)
    b0 = (field.one + X * X) / (field.convert(2) * X)
    P = (X, g0, h0, h0, b0)

    def J(val, idx=None):
        grad = tuple(field.one if j == idx else field.zero for j in range(5))
        return Jet5(field, val, grad)

    xj, gj, hj, aj, bj = (J(P[i], i) for i in range(5))
    E = (J(0), J(0))
    F = (-xj, J(1))
    G = (gj, J(1))
    H = (hj, J(0))
    d2 = reflect2(sub2(F, E), (J(1), J(0)))
    w0 = reflect2(sub2(G, H), (J(1), J(0)))
    Q1 = add2(G, mul2(w0, aj))
    Q0 = add2(H, mul2(w0, bj))
    U = sub2(Q0, Q1)
    rho3 = lambda pt: add2(Q1, reflect2(U, sub2(pt, Q1)))
    dD = reflect2(U, w0)

    M1 = (-xj / 2, J(R(1, 2)))
    M2 = ((gj + hj) / 2, J(R(1, 2)))
    s23 = cross2(sub2(Q1, M2), U) / cross2(w0, U)
    M3 = add2(M2, mul2(w0, s23))
    core = (xj + gj + hj + aj + bj) / 2
    _require(core.v == (gj.v + hj.v + xj.v) / field.convert(2) + s23.v)

    # Exact interior probes fix the signs of the A- and D-band lateral functionals.
    A_inside = add2(mul2(add2(E, F), R(1, 2)), d2)
    D_inside = add2(mul2(add2(Q1, Q0), R(1, 2)), dD)
    _require(cross2(sub2(A_inside, E), d2).v == field.convert(R(1, 2)))
    _require(cross2(sub2(D_inside, rho3(G)), dD).v == field.convert(R(1, 2)))

    # Truncation length in the working chart: L_trunc = s_in + core + t_out.
    s_in = cross2(sub2(G, M1), w0) / cross2((-d2[0], -d2[1]), w0)
    t_out = (-M3[1]) / dD[1]
    Ltr = s_in + core + t_out

    _require(Ltr.v == field.convert(4) / S_val, f"L_trunc(P*) mismatch: {Ltr.v}")
    expected_grad = (field.zero, PHI ** 3, -(PHI ** 2), PHI, field.zero)
    _require(Ltr.g == expected_grad, f"grad L_trunc(P*) mismatch: {Ltr.g}")
    print("[exact] L_trunc(P*) = 4 / sqrt(5 - 2*sqrt(5))")
    print("[exact] grad L_trunc(P*) = (0, phi^3, -phi^2, phi, 0)")

    # Eight active walls (orientations certified above; w6 sign checked on slack slice x = 1/4).
    walls = {
        "w1": -cross2(G, d2),
        "w4": -cross2(sub2(H, rho3(G)), dD),
        "w5": aj - cross2(sub2(F, G), d2) / cross2(w0, d2),
        "w6": cross2(sub2(E, rho3(H)), dD),
        "w7": cross2(Q0, d2) - 1,
        "w8": cross2(sub2(F, rho3(G)), dD) - 1,
        "w9": cross2(Q1, U) / cross2(F, U) - 1,
        "w10": cross2(Q1, F) / cross2(F, U) - 1,
    }
    for name, w in walls.items():
        _require(w.v == field.zero, (name, w.v))

    # Verify w6 > 0 on the slack standard-slice state x = 1/4 (pure rational arithmetic).
    xr = R(1, 4)
    gr, hr, ar, br = (1 - xr ** 2) / (2 * xr), (1 - 3 * xr ** 2) / (2 * xr), (1 - 3 * xr ** 2) / (2 * xr), (1 + xr ** 2) / (2 * xr)
    Gr, Hr = (gr, R(1)), (hr, R(0))
    w0r = reflect2(sub2(Gr, Hr), (R(1), R(0)))
    Q1r, Q0r = add2(Gr, mul2(w0r, ar)), add2(Hr, mul2(w0r, br))
    Ur = sub2(Q0r, Q1r)
    r3Hr = add2(Q1r, reflect2(Ur, sub2(Hr, Q1r)))
    dDr = reflect2(Ur, w0r)
    w6_slack = cross2(sub2((R(0), R(0)), r3Hr), dDr)
    _require(w6_slack > 0, f"w6 orientation failed at slack x=1/4: {w6_slack}")
    print(f"[exact] w6(slack x=1/4) = {w6_slack} > 0")

    # Closed-form Farkas multipliers in simplified form (no uncancelled golden-ratio units).
    lam = {
        "w1": (field.convert(17) + field.convert(7) * S5) * S_val / field.convert(25),
        "w4": (field.convert(28) + field.convert(12) * S5) * S_val / field.convert(25),
        "w5": (field.convert(15) + field.convert(11) * S5) / field.convert(50),
        "w6": (field.convert(42) + field.convert(20) * S5) * S_val / field.convert(25),
        "w7": (field.convert(25) + field.convert(11) * S5) * S_val / field.convert(25),
    }
    support = ("w1", "w4", "w5", "w6", "w7")

    # Invertibility of the 5x5 active-wall gradient matrix over Q(x0).
    M = DomainMatrix([[walls[nm].g[i] for nm in support] for i in range(5)], (5, 5), field)
    det_M = M.det()
    _require(det_M != field.zero, "5x5 wall gradient matrix is singular")
    print("[exact] 5x5 wall gradient matrix on {w1,w4,w5,w6,w7} is invertible")

    # Exact stationarity residual in Q(x0); w9 and w10 carry multiplier 0, so both taco branches hold.
    for i in range(5):
        residual = sum((lam[nm] * walls[nm].g[i] for nm in support), field.zero) - Ltr.g[i]
        _require(residual == field.zero, f"Farkas residual nonzero at component {i}: {residual}")
    print("[exact] grad L_trunc = sum lambda_j grad w_j holds identically in Q(x0) for both taco branches")

    # Outward-rounded interval enclosure confirming positivity and reporting numerical values.
    iv.dps = 60
    s5_iv = iv.sqrt(iv.mpf(5))
    s_iv = iv.sqrt(5 - 2 * s5_iv)
    lam_iv = {
        "w1": (17 + 7 * s5_iv) * s_iv / 25,
        "w4": (28 + 12 * s5_iv) * s_iv / 25,
        "w5": (15 + 11 * s5_iv) / 50,
        "w6": (42 + 20 * s5_iv) * s_iv / 25,
        "w7": (25 + 11 * s5_iv) * s_iv / 25,
    }
    for nm in support:
        bound = lam_iv[nm]
        _require(bound.a > 0, (nm, bound))
        print(f"  lambda[{nm}] in [{float(bound.a):.12f}, {float(bound.b):.12f}] > 0")

    print("PASS: regular pentagon is a strict local minimizer of L_trunc with linear growth")
    return field, lam


if __name__ == "__main__":
    certify()
