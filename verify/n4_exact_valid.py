"""Exact figure-eight validity against every actual crease-bounded pass.

All half-plane and crossing decisions are made in Q(sqrt(15)). There is no
floating-point gating of the candidate inventory and no perpendicular cap.
The eight lateral contacts are checked separately from full nonpenetration.
"""
import sympy as sp
from sympy import Rational as R
from n4_exact_state import tight_state
r15=sp.sqrt(15)
RAY=None

config_exact = tight_state()
d, M, CH, L = (config_exact[k] for k in ("d", "M", "CH", "L"))
PIECES={0:(M[0],-d[0],RAY),1:(M[0],d[1],L[0]),2:(M[1],d[2],L[1]),3:(M[2],d[3],L[2]),4:(M[3],d[4],RAY)}

def simp(e):
    return sp.simplify(sp.expand(e))
def sgn_exact(e):
    """Q(√15) 원소의 정확 부호.
    exact sign of an element of Q(√15)."""
    e=sp.simplify(sp.expand(e))
    if e==0: return 0
    a=sp.Rational(sp.expand(e).coeff(r15,0)); b=sp.Rational(sp.expand(e).coeff(r15,1))
    assert sp.expand(e-a-b*r15) == 0, e
    # e = a + b√15: 부호 = sign 비교 a² vs 15b² — the sign is determined by comparing a² vs 15b²
    if b==0: return 1 if a>0 else -1
    if a==0: return 1 if b>0 else -1
    if a>0 and b>0: return 1
    if a<0 and b<0: return -1
    big = a*a - 15*b*b
    if a>0: return 1 if big>0 else -1
    return -1 if big>0 else 1

from n4_geometry import halfplanes, open_chord_interval, forbidden_pairs

ORDER = [1, 4, 2, 0, 3]
planes = [halfplanes(config_exact, k, sign=sgn_exact, half=R(1, 2))
          for k in range(5)]

def evaluate(form, point):
    a, b, c = form
    return simp(a * point[0] + b * point[1] + c)

# Every complete crease belongs to both adjoining closed passes. In
# particular, the two caps of a core panel do not cross inside its band.
for i, chord in enumerate(CH):
    for k in (i, i + 1):
        for endpoint in chord:
            assert all(sgn_exact(evaluate(f, endpoint)) >= 0 for f in planes[k])
for k in range(5):
    if k == 0:
        interior = M[0] - d[0]
    elif k == 4:
        interior = M[3] + d[4]
    else:
        interior = (M[k - 1] + M[k]) / 2
    assert all(sgn_exact(evaluate(f, interior)) > 0 for f in planes[k])
print('[exact panels] all complete creases and five nonempty pass interiors verified')

pairs = list(forbidden_pairs(ORDER))
assert len(pairs) == 6
for i, k in pairs:
    overlap = open_chord_interval(CH[i], planes[k], sign=sgn_exact,
                                  zero=R(0), one=R(1))
    assert overlap is None, (i + 1, k, overlap)
    print(f'[exact nonpenetration] crease {i + 1} / pass {k}: empty open intersection')

# Check all non-adjacent centreline pairs, including the parallel pairs.
actual_crossings = []
for i in range(5):
    for j in range(i + 2, 5):
        P1, v1, l1 = PIECES[i]
        P2, v2, l2 = PIECES[j]
        determinant = simp(v1[0] * v2[1] - v1[1] * v2[0])
        delta = P2 - P1
        if sgn_exact(determinant) == 0:
            assert sgn_exact(simp(delta[0] * v1[1] - delta[1] * v1[0])) != 0
            continue
        t = simp((delta[0] * v2[1] - delta[1] * v2[0]) / determinant)
        u = simp((delta[0] * v1[1] - delta[1] * v1[0]) / determinant)
        inside = (sgn_exact(t) > 0 and sgn_exact(u) > 0
                  and (l1 is None or sgn_exact(l1 - t) > 0)
                  and (l2 is None or sgn_exact(l2 - u) > 0))
        if inside:
            actual_crossings.append((i, j))
            print(f'[exact crossing] {(i, j)}: parameters {t}, {u}')
assert actual_crossings == [(0, 2), (0, 3), (1, 4), (2, 4)]

# CV-oriented endpoint labels; negative turns reverse the config labels.
contacts = [(0, 1, 4), (0, 0, 2), (1, 0, 4), (1, 1, 4),
            (2, 0, 0), (2, 1, 0), (3, 1, 0), (3, 0, 2)]
for i, e, k in contacts:
    point = CH[i][e]
    lateral = simp((point - (M[0] if k == 0 else M[k - 1])).dot(
        sp.Matrix([-d[k][1], d[k][0]])))
    assert simp(lateral ** 2 - R(1, 4)) == 0
print('[exact active contacts] all eight lateral equalities verified')
print('PASS: full actual-panel nonpenetration and the four-crossing cell, exactly in Q(sqrt(15))')
