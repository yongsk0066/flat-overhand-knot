"""Recover the pulling chart and crossing order from crease midpoints and reflections."""

import sympy as sp

from ribbon_geometry import DaggerField


def det(left, right):
    return left[0] * right[1] - left[1] * right[0]


def reflect(axis, direction):
    return (2 * axis.dot(direction) * axis / axis.dot(axis) - direction).applyfunc(sp.factor)


def equal(label, actual, expected):
    if sp.cancel(actual - expected) != 0:
        raise RuntimeError(label)


def certify():
    x, q, z, m, k = sp.symbols("x q z m k")
    n = 1 + q * q
    height_scale = n / q
    lx, lz = (m - q) / 2, (k - q) / 2
    h, a = (m - 2 * q - x) / 2, (k - 2 * q - z) / 2
    u, v = m + x - height_scale, k + z - height_scale

    def T(y):
        return (q + y) * (1 - q * y)

    def angle(y):
        return n * y / T(y)

    horizontal = sp.Matrix([1, 0])
    fold1 = sp.Matrix([-x, 1])
    fold2 = sp.Matrix([q, 1])
    fold3 = reflect(fold2, sp.Matrix([-z, 1]))  # Q0 to Q1
    incoming = reflect(fold1, horizontal)
    middle = reflect(fold2, horizontal)
    outgoing = reflect(fold3, middle)

    mid1 = -lx * horizontal
    mid2 = sp.zeros(2, 1)
    mid3 = lz * middle
    E, F = mid1 - fold1 / 2, mid1 + fold1 / 2
    H, G = mid2 - fold2 / 2, mid2 + fold2 / 2
    Q0, Q1 = mid3 - fold3 / 2, mid3 + fold3 / 2
    side_a = lambda point: det(point - E, incoming)
    side_d = lambda point: det(point - Q1, outgoing)

    left_scale = T(x) / (n * (1 + x * x))
    right_scale = T(z) / (n * (1 + z * z))
    equal("incoming side normalization", side_a(F), 1)
    equal("outgoing side normalization", side_d(Q0), 1)
    equal("middle crease outside A", -side_a(G), (m * x - 1) / (1 + x * x))
    equal("middle crease outside D", -side_d(H), (k * z - 1) / (1 + z * z))
    equal("strict other endpoint H", side_a(H), -2 * h * x / (1 + x * x))
    equal("strict other endpoint G", side_d(G), -2 * a * z / (1 + z * z))
    equal("far crease Q1", side_a(Q1) - 1, left_scale * (k - 2 * q - z - angle(x) * (m + x)))
    equal("far crease Q0", side_a(Q0) - 1, left_scale * (v - angle(x) * u))
    equal("first crease E", side_d(E) - 1, right_scale * (m - 2 * q - x - angle(z) * (k + z)))
    equal("first crease F", side_d(F) - 1, right_scale * (u - angle(z) * v))
    equal("corner vertical displacement", Q0[1] - F[1], v / height_scale)
    equal("strict opposite end height", Q1[1] - F[1], 2 * q * a / n)
    equal("corner horizontal displacement", Q0[0] - F[0], (u + (q * q - 1) * v / n) / 2)

    major = q + x + z - q * x * z
    minor = 1 - q * x - q * z - x * z
    jac = det(incoming, outgoing)
    equal("angle product gap", T(x) * T(z) - x * z * n * n, q * major * minor)
    equal("nonparallel distant creases", n * det(fold1, -fold3), major + q * minor)
    equal("AD orientation", jac, 2 * major * minor / (n * (1 + x * x) * (1 + z * z)))
    carrier = det(fold1, -fold3)
    tail_cap = lambda point: det(point - E, fold1)
    equal("incoming cap increases towards Q1", tail_cap(Q1) - tail_cap(Q0), carrier)
    taco_t = det(Q1 - E, -fold3) - carrier
    taco_u = det(Q1 - E, fold1) - carrier
    equal("outgoing carrier end", 2 * taco_t, v + (q * q + 2 * q * z - 1) * u / n)
    equal("incoming carrier end", 2 * taco_u, u + (q * q + 2 * q * x - 1) * v / n)
    equal("cap at Q0", tail_cap(Q0), taco_u)
    equal("AC orientation", det(middle, incoming), 2 * left_scale)
    equal("BD orientation", -outgoing[1], 2 * right_scale)
    equal("last fold turn", det(middle, outgoing), 2 * z / (1 + z * z))
    equal("first fold turn", det(incoming, horizontal), 2 * x / (1 + x * x))
    equal("middle fold turn", det(horizontal, middle), 2 * q / n)

    # X=AC, Y=AD, Z=BD. On A these parameters point back from M1.
    c_x = det(mid1 - mid2, incoming) / det(middle, incoming)
    back_b_z = det(mid3 - mid2, outgoing) / det(-horizontal, outgoing)
    a_x = det(mid2 - mid1, middle) / det(-incoming, middle)
    a_y = det(mid3 - mid1, outgoing) / det(-incoming, outgoing)
    d_y = det(mid3 - mid1, incoming) / jac
    d_z = det(mid2 - mid3, horizontal) / det(outgoing, horizontal)
    equal("AC on C", c_x, angle(x) * lx)
    equal("BD on B", back_b_z, angle(z) * lz)
    equal("AD on A", a_y, (lx - back_b_z) * (-outgoing[1]) / jac)
    equal("AD on D", d_y, (lz - c_x) * det(middle, incoming) / jac)
    equal("material order on A", a_x - a_y, (lz - c_x) * det(middle, outgoing) / jac)
    equal("material order on D", d_z - d_y, (lx - back_b_z) * det(incoming, horizontal) / jac)
    print("PASS: six contact functions, crease heights and nonparallel carriers from midpoints")
    print("PASS: interior crossing parameters and both strict passage-order gaps for XYZXYZ")

    time = sp.symbols("time")
    dilation = {m: m + 2 * time * h, k: k + 2 * time * a}
    for margin, slope in [
        (m * x - 1, 2 * x * h),
        (k * z - 1, 2 * z * a),
        (k - 2 * q - z - angle(x) * (m + x), 2 * (a - angle(x) * h)),
        (m - 2 * q - x - angle(z) * (k + z), 2 * (h - angle(z) * a)),
        (v - angle(x) * u, 2 * (a - angle(x) * h)),
        (u - angle(z) * v, 2 * (h - angle(z) * a)),
        (Q0[1] - F[1], 2 * a / height_scale),
    ]:
        equal(
            "finite opening of a clearance",
            margin.subs(dilation, simultaneous=True),
            margin + time * slope,
        )
    print("PASS: fixed-angle dilation opens all six side margins and the distant crease gap")

    field = DaggerField()
    endpoint = field.state()
    x0, g0, h0, a0, b0 = endpoint["theta"]
    q0 = field.red(g0 - h0)
    s0 = endpoint["lengths"]["core"]
    d0 = field.red(s0 + q0)
    floor_angle = field.div(1, d0)
    floor_height = field.red((d0 - 2 * q0 - floor_angle) / 2)
    reconstructed = (
        floor_angle,
        floor_height + q0,
        floor_height,
        floor_height,
        floor_height + q0 + floor_angle,
    )
    for actual, expected in zip(reconstructed, (x0, g0, h0, a0, b0), strict=True):
        if field.red(actual - expected) != 0:
            raise RuntimeError("symmetric path endpoint")
    for margin in (q0 * d0 - 1, s0 * s0 - 1 - q0 * q0, s0 - 2, h0, x0, q0):
        if field.sign(margin) != 1:
            raise RuntimeError("endpoint chart margin")
    print("PASS: all five terminal coordinates and strict chart margins in the dagger field")


if __name__ == "__main__":
    certify()
