"""Actual crease-bounded pass geometry, shared by probes and exact certificates.

Each pass is the intersection of its width-one band and the appropriate
crease half-planes. In particular, a cap is a crease, not the perpendicular
line through the corresponding centreline vertex. No ray truncation is used.
The same affine clipping code accepts floats and exact ordered field elements.
"""
import math


def cross(u, v):
    return u[0] * v[1] - u[1] * v[0]


def sub(u, v):
    return (u[0] - v[0], u[1] - v[1])


def halfplanes(c, k, *, sign, half=0.5):
    """Return (a,b,c) for the inequalities a*x+b*y+c >= 0 of pass k."""
    directions, vertices, chords = c['d'], c['M'], c['CH']
    origin = vertices[0] if k == 0 else vertices[k - 1]
    direction = directions[k]
    nx, ny = -direction[1], direction[0]
    offset = nx * origin[0] + ny * origin[1]
    result = [(nx, ny, half - offset), (-nx, -ny, half + offset)]
    # Inside lies forward of the entrance crease and backward of the exit.
    caps = ([(k - 1, direction)] if k > 0 else [])
    if k < 4:
        caps.append((k, (-direction[0], -direction[1])))
    for i, inward in caps:
        edge = sub(chords[i][1], chords[i][0])
        orientation = sign(cross(edge, inward))
        if orientation == 0:
            raise ValueError('degenerate crease: no inward half-plane')
        a, b = -orientation * edge[1], orientation * edge[0]
        result.append((a, b, -(a * vertices[i][0] + b * vertices[i][1])))
    return result


def open_chord_interval(chord, planes, *, sign, zero=0, one=1,
                        distance_tolerance=0.0):
    """Parameter interval where the open chord is in every open half-plane.

    Return None when empty, otherwise (lo,hi). With tolerance zero and an
    exact sign function every decision is exact; numeric probes may shrink
    each half-plane by a stated geometric tolerance.
    """
    start, end = chord
    lo, hi = zero, one
    for a, b, offset in planes:
        threshold = (distance_tolerance * math.hypot(float(a), float(b))
                     if distance_tolerance else zero)
        constant = a * start[0] + b * start[1] + offset - threshold
        slope = a * (end[0] - start[0]) + b * (end[1] - start[1])
        direction = sign(slope)
        if direction == 0:
            if sign(constant) <= 0:
                return None
            continue
        root = -constant / slope
        if direction > 0:
            if sign(root - lo) > 0:
                lo = root
        elif sign(root - hi) < 0:
            hi = root
        if sign(hi - lo) <= 0:
            return None
    return lo, hi


def forbidden_pairs(order):
    rank = {p: i for i, p in enumerate(order)}
    for i in range(4):
        lo, hi = sorted((rank[i], rank[i + 1]))
        for k in order[lo + 1:hi]:
            yield i, k


def numeric_axiom_ok(c, order, *, distance_tolerance=1e-12,
                     parameter_tolerance=1e-11):
    """Floating-point probe, not an exact certificate."""
    sign = lambda v: int(v > 0) - int(v < 0)
    try:
        passes = [halfplanes(c, k, sign=sign) for k in range(5)]
        # Complete adjacent creases must not cross inside the material band.
        for i, chord in enumerate(c['CH']):
            for k in (i, i + 1):
                for x, y in chord:
                    if any(a * x + b * y + z < -distance_tolerance * math.hypot(a, b)
                           for a, b, z in passes[k]):
                        return False
        for i, k in forbidden_pairs(order):
            interval = open_chord_interval(
                c['CH'][i], passes[k], sign=sign,
                distance_tolerance=distance_tolerance)
            if interval is not None and interval[1] - interval[0] > parameter_tolerance:
                return False
    except (ValueError, ZeroDivisionError):
        return False
    return True
