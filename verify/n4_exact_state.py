"""The four-crease figure-eight state in Q(sqrt(15)), with signed crease ends."""

import sympy as sp


def tight_state():
    root = sp.sqrt(15)
    turns = [
        (-sp.Rational(1, 4), root / 4),
        (-sp.Rational(7, 8), root / 8),
        (-sp.Rational(7, 8), -root / 8),
        (-sp.Rational(1, 4), -root / 4),
    ]
    lengths = [2 * root / 5, 8 * root / 15, 2 * root / 5]
    directions = [sp.Matrix([1, 0])]
    for cosine, sine in turns:
        x, y = directions[-1]
        directions.append(sp.Matrix([cosine * x - sine * y, sine * x + cosine * y]))

    midpoints = [sp.Matrix([0, 0])]
    for length, direction in zip(lengths, directions[1:4], strict=True):
        midpoints.append(sp.expand(midpoints[-1] + length * direction))

    chords = []
    for i, (_, sine) in enumerate(turns):
        half_chord = (directions[i] + directions[i + 1]) / (2 * sine)
        start = sp.expand(midpoints[i] - half_chord)
        end = sp.expand(midpoints[i] + half_chord)
        chords.append((start, end))
    return dict(d=directions, M=midpoints, CH=chords, L=lengths)
