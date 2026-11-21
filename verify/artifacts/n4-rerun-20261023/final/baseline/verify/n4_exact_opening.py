"""Exact actual-panel opening path, with a certified whole interval.

Angles stay fixed and L(s)=L(0)+s*(8sqrt(15)/15,4sqrt(15)/5,
8sqrt(15)/15). For 0<s<=1/1000 every forbidden chord lies strictly
outside the corresponding infinite band. The finite/half-infinite actual
passes are subsets of those bands, so their caps need no approximation.
Panel validity and crossing windows are checked on the entire interval by
exact affine endpoint signs.
"""
import contextlib
import io
import sympy as sp
from sympy import Rational as R
with contextlib.redirect_stdout(io.StringIO()):
    import n4_exact_valid as base
from n4_geometry import halfplanes, forbidden_pairs

s = sp.Symbol('s', real=True)
S = R(1, 1000)
r15 = sp.sqrt(15)
d = base.d
L = [old + s*step for old, step in zip(base.L, [8*r15/15, 4*r15/5, 8*r15/15])]
M = [sp.Matrix([0, 0])]
for i in range(3):
    M.append(sp.expand(M[-1] + L[i]*d[i+1]))
CH = [tuple(sp.expand(M[i] + endpoint - base.M[i]) for endpoint in base.CH[i])
      for i in range(4)]
config = dict(d=d, M=M, CH=CH, L=L)
planes = [halfplanes(config, k, sign=base.sgn_exact, half=R(1, 2)) for k in range(5)]


def affine_positive(expr, *, allow_zero_at_start=False, allow_identically_zero=False):
    expr = sp.expand(expr)
    assert sp.Poly(expr, s).degree() <= 1 or expr == 0, expr
    lo, hi = sp.simplify(expr.subs(s, 0)), sp.simplify(expr.subs(s, S))
    if allow_identically_zero and expr == 0:
        return
    assert base.sgn_exact(lo) >= (0 if allow_zero_at_start else 1), (expr, lo)
    assert base.sgn_exact(hi) > 0, (expr, hi)


def at(form, point):
    a, b, c = form
    return sp.expand(a*point[0]+b*point[1]+c)


# Adjacent complete creases remain within their actual panels.
for i, chord in enumerate(CH):
    for k in (i, i+1):
        for endpoint in chord:
            for form in planes[k]:
                affine_positive(at(form, endpoint), allow_zero_at_start=True,
                                allow_identically_zero=True)
for i, k in forbidden_pairs(base.ORDER):
    origin = M[0] if k == 0 else M[k-1]
    lateral = [sp.expand((p-origin).dot(sp.Matrix([-d[k][1], d[k][0]]))) for p in CH[i]]
    side = base.sgn_exact(sp.simplify(sum(lateral).subs(s, 0)))
    assert side != 0
    for e in lateral:
        affine_positive(side*e-R(1, 2), allow_zero_at_start=True)
    print(f'[exact opening] crease {i+1} / pass {k}: both endpoints outside the whole band for 0<s<=1/1000')

pieces = {0:(M[0],-d[0],None), 1:(M[0],d[1],L[0]),
          2:(M[1],d[2],L[1]), 3:(M[2],d[3],L[2]), 4:(M[3],d[4],None)}
found = []
events = {i: [] for i in range(5)}
for i in range(5):
    for j in range(i+2,5):
        p, v, li = pieces[i]; q, w, lj = pieces[j]
        det = sp.simplify(v[0]*w[1]-v[1]*w[0])
        delta = q-p
        if det == 0:
            offset = sp.expand(delta[0]*v[1]-delta[1]*v[0])
            sign = base.sgn_exact(offset.subs(s, 0))
            affine_positive(sign*offset)
            continue
        t = sp.expand((delta[0]*w[1]-delta[1]*w[0])/det)
        u = sp.expand((delta[0]*v[1]-delta[1]*v[0])/det)
        affine_positive(t); affine_positive(u)
        if li is not None: affine_positive(li-t)
        if lj is not None: affine_positive(lj-u)
        found.append((i,j))
        events[i].append(t)
        events[j].append(u)
        crossing_point = p+t*v
        for ci, (a, b) in enumerate(CH):
            edge = b-a
            gap = sp.expand(edge[0]*(crossing_point[1]-a[1])
                            -edge[1]*(crossing_point[0]-a[0]))
            side = base.sgn_exact(gap.subs(s, 0))
            assert side != 0, ('crossing on crease carrier', i, j, ci)
            affine_positive(side*gap)
assert found == [(0,2),(0,3),(1,4),(2,4)]
for parameters in events.values():
    for i, t in enumerate(parameters):
        for u in parameters[i+1:]:
            gap = sp.expand(t-u)
            side = base.sgn_exact(gap.subs(s, 0))
            assert side != 0
            affine_positive(side*gap)
assert sp.simplify(sp.diff(sum(L),s)-28*r15/15) == 0
print('[exact opening] all panel and crossing-window inequalities hold throughout [0,1/1000]')
print('[exact opening] core=4sqrt(15)/3+(28sqrt(15)/15)s; reversal is a strictly decreasing pulling path')
print('PASS: reachable from actual open states; no numerical gating or finite sampling')
