"""Exact first-order sufficiency on the stated smooth contact chamber.

The chamber is {g_i >= 0}, with the eight functions in n4_exact_kkt.py.
Seven independent active gradients with strictly positive multipliers imply
that the objective grows strictly along every nonzero vector of its tangent
cone, hence strict local minimality in that chamber. This certificate does
not assume that the chamber exhausts the actual geometric cell. In the
CV endpoint convention, g2, g3a, g3ap and g2p contact a foreign crease cap
as well as a lateral edge. Their necessity needs a corner analysis; it
does not follow from the usual interior-edge nonpenetration argument.
"""
import contextlib
import io
import runpy
import sympy as sp
from pathlib import Path

with contextlib.redirect_stdout(io.StringIO()):
    data = runpy.run_path(str(Path(__file__).with_name('n4_exact_kkt.py')))
J, A, lam = data['J'], data['A'], data['lam']
assert all(data['X'](g) == 0 for _, g in data['G'])
assert sp.simplify(A.det() + sp.sqrt(15)/2) == 0
assert all(value > 0 for value in lam)
assert A*lam == data['grad_ell']
assert J*data['v_int'] == sp.ones(8,1)
assert sp.simplify(data['Lt']-32*sp.sqrt(15)/15) == 0
print('[exact KKT] eight primal equalities; seven positive multipliers; residual exactly zero')
print('[exact rank] determinant = -sqrt(15)/2')
print('[exact opening] J v = (1,1,1,1,1,1,1,1)')
print('PASS: strict local minimum in the explicitly stated eight-contact chamber')
print('Scope: full-cell local optimality additionally requires necessity of these contact inequalities')
