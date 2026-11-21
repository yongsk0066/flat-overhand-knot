"""Exact field checks, independent of numerical fitting.

Use the finite-field coefficient API to avoid deprecated comparisons of
modular integers in SymPy 1.14's high-level factor-list sorting.
"""
import os
os.environ["SYMPY_GROUND_TYPES"] = "python"
import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_factor

x, h, y = sp.symbols("x h y")
p = x**12 + 4*x**10 - 7*x**4 + 20*x**2 - 2
q = y**6 + 4*y**5 - 7*y**2 + 20*y - 2
n5 = (8*h**3*x**3 + 4*h**2*x**4 + 4*h**2*x**2 + 2*h*x**5
      - 4*h*x**3 - 6*h*x + x**6 + 3*x**4 + 3*x**2 + 1)
k = 2*x**2*sp.diff(n5, x) + (1-x**2)*sp.diff(n5, h)
assert sp.Poly(p, x).is_irreducible
assert sp.Poly(q, y).is_irreducible
assert sp.Poly(p, x, extension=sp.sqrt(5)).is_irreducible
assert sp.Poly(p, x).count_roots(-sp.oo, sp.oo) == 2
assert sp.Poly(n5.subs(x, sp.Rational(1, 5)), h).count_roots(-sp.oo, sp.oo) == 3
assert sp.expand(sp.resultant(n5, k, h) + 2**18*x**12*(x*x+1)**2*p) == 0
expected = {11: [6], 23: [1, 5], 37: [1, 1, 1, 3]}
disc = sp.discriminant(q, y)
assert disc == 2**19 * 3**6 * 2531
for prime, degrees in expected.items():
    factors = gf_factor([1, 4, 0, 0, -7, 20, -2], prime, ZZ)[1]
    actual = sorted(len(f) - 1 for f, multiplicity in factors
                    for _ in range(multiplicity))
    assert disc % prime != 0 and actual == degrees
    print(f"mod {prime}: factor degrees {actual}")
print("PASS: irreducibility over Q and Q(sqrt(5)); two real roots; resultant.")
print("PASS: three real ordinates at x=1/5; real-branch selection is necessary.")
print(f"PASS: discriminant = {disc}; cycle types give S6 by Jordan's theorem.")
