"""Independent rational panel-clipping and Gauss-word checks along representative paths."""

from fractions import Fraction as Q
from collections import Counter

one = Q(1)
grid = list(range(65))
counts = Counter()


def difference(a, b):
    return (a[0] - b[0], a[1] - b[1])


def plus(a, b):
    return (a[0] + b[0], a[1] + b[1])


def times(t, v):
    return (t * v[0], t * v[1])


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def reflect(axis, v):
    return difference(times(2 * dot(axis, v) / dot(axis, axis), axis), v)


def mid(a, b):
    return times(one / 2, plus(a, b))


def meet(p, v, q, w):
    determinant = det(v, w)
    if not determinant:
        raise ValueError("parallel carriers")
    return det(difference(q, p), w) / determinant, det(difference(q, p), v) / determinant


def band(origin, other, direction, interior):
    """Three affine half-planes obtained from side and cap carriers with an interior probe."""
    planes = []
    for pt, d in [(origin, direction), (other, direction), (origin, difference(other, origin))]:
        sense = det(d, difference(interior, pt))
        if not sense:
            raise ValueError("interior probe on a cap or a side")
        sign = 1 if sense > 0 else -1
        planes.append(lambda v, p=pt, d=d, sign=sign: sign * det(d, difference(v, p)))
    return planes


def penetrates(chord, halfplanes):
    low, high = Q(0), one
    for f in halfplanes:
        start, end = map(f, chord)
        slope = end - start
        if not slope:
            if start <= 0:
                return False
        elif slope > 0:
            low = max(low, -start / slope)
        else:
            high = min(high, -start / slope)
    return low < high


def checked(state, label, strict=False):
    q, m, k, x, z = state
    h = (m - 2 * q - x) / 2
    a = (k - 2 * q - z) / 2
    if min(q, x, z, h, a) <= 0:
        raise ValueError((label, "negative material", state))
    E, F, H, G = (Q(0), Q(0)), (-x, one), (h, Q(0)), (h + q, one)
    e1 = (one, Q(0))
    incoming = reflect(difference(F, E), e1)
    middle = reflect(difference(G, H), e1)
    Q1 = plus(G, times(a, middle))
    Q0 = plus(H, times(a + q + z, middle))
    outgoing = reflect(difference(Q0, Q1), middle)
    M1, M2, M3 = mid(E, F), mid(H, G), mid(Q1, Q0)
    L1, L2 = dot(difference(M2, M1), e1), dot(difference(M3, M2), middle)
    if L1 + L2 != (m + k) / 2 - q:
        raise ValueError((label, "wrong core"))
    A = band(E, F, incoming, difference(M1, incoming))
    D = band(Q0, Q1, outgoing, plus(M3, outgoing))
    creases = [(E, F), (H, G), (Q1, Q0)]
    for chord, plane in [(creases[0], D), (creases[1], A), (creases[1], D), (creases[2], A)]:
        if penetrates(chord, plane):
            raise ValueError((label, "penetrating crease", state))
    if strict:
        for pt, panel in [(G, A), (H, D), (Q0, A), (Q1, A), (E, D), (F, D)]:
            if min(f(pt) for f in panel[:2]) >= 0:
                raise ValueError((label, "non-slack side wall"))
    t, u = meet(E, difference(F, E), Q1, difference(Q0, Q1))
    if (0 < t < 1 and 0 < u < 1) or (strict and min(t, u) <= 1):
        raise ValueError((label, "crease intersection", t, u))
    data = [
        (M1, incoming, None, Q(0)),
        (M1, e1, Q(0), L1),
        (M2, middle, Q(0), L2),
        (M3, outgoing, Q(0), None),
    ]
    passages = [[] for _ in data]
    points = []
    for i, j in [(0, 2), (0, 3), (1, 3)]:
        p, ip, il, ih = data[i]
        b, jp, jl, jh = data[j]
        a1, b1 = meet(p, ip, b, jp)
        for value, lo, hi in [(a1, il, ih), (b1, jl, jh)]:
            if (lo is not None and value <= lo) or (hi is not None and value >= hi):
                raise ValueError((label, "missing crossing", i, j, state))
        pt = plus(p, times(a1, ip))
        for left, right in creases:
            if (
                det(difference(right, left), difference(pt, left)) == 0
                and dot(difference(pt, left), difference(pt, right)) <= 0
            ):
                raise ValueError((label, "crossing on a crease"))
        points.append(pt)
        passages[i].append((a1, (i, j)))
        passages[j].append((b1, (i, j)))
    raw = [pair for part in passages for time, pair in sorted(part)]
    names = {v: chr(ord("X") + i) for i, v in enumerate(dict.fromkeys(raw))}
    if len(set(points)) != 3 or "".join(names[i] for i in raw) != "XYZXYZ":
        raise ValueError((label, "wrong passage word", state))
    counts["panel checks"] += 1


def least(q, m, k):
    C = q + 1 / q
    if m <= 1 / q:
        return C - m, C - k, "corner"
    R = 2 * (1 + q * q) * m / (q * m + 1)
    if R - k > 1 / k:
        return 1 / m, R - k, "mixed"
    return 1 / m, 1 / k, "both floors"


def interp(a, b, t):
    return a + (b - a) * t


def check_example(name, q, m, k, x, z, reversed_material=False):
    def original(q, m, k, x, z):
        return (q, k, m, z, x) if reversed_material else (q, m, k, x, z)

    checked(original(q, m, k, x, z), name, strict=True)
    C = q + 1 / q
    U, V = m + x - C, k + z - C
    ratio = max(Q(0), (m + 1 / m - C) / U, (k + 1 / k - C) / V)
    x1, z1 = C - m + ratio * U, C - k + ratio * V
    if x1 == 1 / m and z1 == 1 / k:
        expected = "corner and both floors" if ratio == 0 else "both floors"
    elif ratio == 0:
        expected = "corner"
    else:
        expected = "left floor" if x1 == 1 / m else "right floor"
    counts["ray hits " + expected] += 1
    xl, zl, _ = least(q, m, k)
    for xa, za, xb, zb in [(x, z, x1, z1), (x1, z1, xl, zl)]:
        for j in grid:
            t = Q(j, 64)
            checked(original(q, m, k, interp(xa, xb, t), interp(za, zb, t)), name)
    d = (m + k) / 2
    mids = {interp(m, d, Q(j, 64)) for j in grid}
    if d <= 1 / q <= m:
        mids.add(1 / q)
    for mt in mids:
        kt = 2 * d - mt
        yt, zt, regime = least(q, mt, kt)
        counts[regime] += 1
        checked(original(q, mt, kt, yt, zt), name)
    s = d - q
    q_proxy = Q(354730, 1000000)
    s_proxy = Q(2749631, 1000000)
    for j in grid:
        qt = interp(q, q_proxy, Q(j, 64))
        dt = s + qt
        y = max(1 / dt, 1 / qt - s)
        checked(original(qt, dt, dt, y, y), name)
    for j in grid:
        st = interp(s, s_proxy, Q(j, 64))
        dt = st + q_proxy
        checked(original(q_proxy, dt, dt, 1 / dt, 1 / dt), name)
    print(
        "PASS",
        name,
        "reversed" if reversed_material else "forward",
        "initial ray hits",
        expected,
        "initial balanced regime",
        least(q, m, k)[2],
    )


def main():
    counts.clear()
    cases = [
        ("corner to symmetric corner", Q(1, 5), Q(99, 20), Q(97, 20), Q(27, 100), Q(37, 100)),
        ("mixed to both floors", Q(1, 5), Q(21, 4), Q(5), Q(21, 100), Q(43, 100)),
        ("mixed to corner", Q(1, 5), Q(126, 25), Q(97, 20), Q(23, 100), Q(41, 100)),
        ("both floors", Q(1, 2), Q(6), Q(5), Q(19, 100), Q(21, 100)),
        ("large middle tilt", Q(2), Q(20), Q(18), Q(6, 100), Q(61, 1000)),
        ("two floor contacts together", Q(7, 20), Q(7, 2), Q(7, 2), Q(3, 10), Q(3, 10)),
        ("corner meeting both floors", Q(1, 5), Q(5), Q(5), Q(21, 100), Q(21, 100)),
    ]
    for case in cases:
        for rev in (False, True):
            check_example(*case, reversed_material=rev)
    for key in (
        "corner",
        "mixed",
        "both floors",
        "ray hits left floor",
        "ray hits right floor",
        "ray hits corner",
        "ray hits both floors",
        "ray hits corner and both floors",
    ):
        if not counts[key]:
            raise ValueError(f"missed {key}")
    for key in sorted(counts):
        print(key, counts[key])
    print("PASS: independent exact rational clipping, crease exclusion, all three crossings and XYZXYZ")
    print("The rational endpoint is an interior proxy. The algebraic endpoint is checked separately.")


if __name__ == "__main__":
    main()
