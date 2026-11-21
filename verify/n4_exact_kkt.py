"""Exact KKT multipliers, opening ray, and truncation length for the figure-eight state.

All contact equalities, Jacobian entries, KKT multipliers, opening rates, and
truncation stubs are verified in Q(sqrt(15)) without heuristic nsimplify or
floating-point calibration.
"""
import sympy as sp
from sympy import Rational as R

def _require(cond, msg="verification failed"):
    if not cond:
        raise RuntimeError(msg)

t1, t2, t3, t4, l1, l2, l3 = sp.symbols("t1 t2 t3 t4 l1 l2 l3", real=True)
TH = [t1, t2, t3, t4, l1, l2, l3]
phi = [0, t1, t1 + t2, t1 + t2 + t3, t1 + t2 + t3 + t4]
d = [sp.Matrix([sp.cos(p), sp.sin(p)]) for p in phi]
M = [sp.Matrix([0, 0])]
for i in range(3):
    M.append(M[-1] + TH[4 + i] * d[i + 1])

# Signed chord half-vectors V_i^+- = M_i +- (d_{i-1} + d_i) / (2 sin tau_i), matching Section 11.
CV = [(d[i] + d[i + 1]) / (2 * sp.sin(TH[i])) for i in range(4)]
Vm = [M[i] - CV[i] for i in range(4)]
Vp = [M[i] + CV[i] for i in range(4)]

lat = lambda pt, P, dd: -(pt[0] - P[0]) * dd[1] + (pt[1] - P[1]) * dd[0]

r15 = sp.sqrt(15)
SUB = {
    sp.cos(t1): -R(1, 4),
    sp.sin(t1): r15 / 4,
    sp.cos(t2): -R(7, 8),
    sp.sin(t2): r15 / 8,
    sp.cos(t3): -R(7, 8),
    sp.sin(t3): -r15 / 8,
    sp.cos(t4): -R(1, 4),
    sp.sin(t4): -r15 / 4,
    l1: 2 * r15 / 5,
    l2: 8 * r15 / 15,
    l3: 2 * r15 / 5,
}


def X(e):
    return sp.simplify(sp.expand_trig(e).subs(SUB))


# Eight smooth contact chamber functions defining C in Section 11.
G = [
    ("g1", -lat(Vp[0], M[3], d[4]) - R(1, 2)),
    ("g2", lat(Vm[0], M[1], d[2]) - R(1, 2)),
    ("g3a", lat(Vp[1], M[3], d[4]) - R(1, 2)),
    ("g3b", lat(Vm[1], M[3], d[4]) - R(1, 2)),
    ("g1p", lat(Vp[3], M[0], d[0]) - R(1, 2)),
    ("g2p", -lat(Vm[3], M[1], d[2]) - R(1, 2)),
    ("g3ap", -lat(Vp[2], M[0], d[0]) - R(1, 2)),
    ("g3bp", -lat(Vm[2], M[0], d[0]) - R(1, 2)),
]

for nm, g_expr in G:
    val = X(g_expr)
    _require(val == 0, f"Active contact {nm} nonzero: {val}")

names = [nm for nm, _ in G]
J = sp.Matrix([[X(sp.diff(g_expr, v)) for v in TH] for _, g_expr in G])
grad_ell = sp.Matrix([0, 0, 0, 0, 1, 1, 1])

# Invertible 7-generator subsystem (excluding g2).
sel = [0, 2, 3, 4, 5, 6, 7]
A = J[sel, :].T
_require(sp.simplify(A.det() + r15 / 2) == 0, f"Unexpected det(A): {A.det()}")

lam = sp.simplify(A.solve(grad_ell))
expected_lam = [sp.Rational(k, 45) * r15 for k in (21, 3, 21, 13, 4, 17, 5)]
for nm, lv, exp_lv in zip([names[i] for i in sel], lam, expected_lam):
    _require(sp.simplify(lv - exp_lv) == 0, f"Multiplier mismatch for {nm}: {lv} != {exp_lv}")
    _require(lv > 0, f"Multiplier not positive for {nm}: {lv}")

res = sp.simplify(A * lam - grad_ell)
_require(all(entry == 0 for entry in res), f"KKT residual nonzero: {list(res)}")

# Exact truncation length via the flush-outline exit rule on the certified hexagon.
V = [
    sp.Matrix([-13 * r15 / 30, R(1, 2)]),
    sp.Matrix([-7 * r15 / 30, R(3, 2)]),
    sp.Matrix([r15 / 30, R(3, 2)]),
    sp.Matrix([3 * r15 / 30, R(1, 2)]),
    sp.Matrix([-3 * r15 / 30, -R(1, 2)]),
    sp.Matrix([-11 * r15 / 30, -R(1, 2)]),
]


def ray_exit(P0, D0):
    ts = []
    for k in range(6):
        A_, B_ = V[k], V[(k + 1) % 6]
        e = B_ - A_
        den = sp.simplify(D0[0] * (-e[1]) + D0[1] * e[0])
        if den == 0:
            continue
        t = sp.simplify((e[0] * (P0[1] - A_[1]) - e[1] * (P0[0] - A_[0])) / (D0[0] * e[1] - D0[1] * e[0]))
        u_pt = P0 + t * D0
        uu = sp.simplify(e.dot(u_pt - A_) / e.dot(e))
        if t.is_positive and bool(uu >= 0) and bool(uu <= 1):
            ts.append(sp.radsimp(t))
    _require(ts, "Ray failed to exit hexagon")
    return sp.Max(*ts)


s_in = ray_exit(sp.Matrix([0, 0]), sp.Matrix([-1, 0]))
s_out = ray_exit(sp.Matrix([-r15 / 3, 1]), sp.Matrix([1, 0]))
_require(sp.simplify(s_in - 2 * r15 / 5) == 0, f"s_in mismatch: {s_in}")
_require(sp.simplify(s_out - 2 * r15 / 5) == 0, f"s_out mismatch: {s_out}")

Lc = 4 * r15 / 3
Lt = sp.simplify(s_in + Lc + s_out)
_require(sp.simplify(Lt - 32 * r15 / 15) == 0, f"L_trunc mismatch: {Lt}")
_require(sp.simplify(Lt / Lc - R(8, 5)) == 0, f"Ratio mismatch: {Lt / Lc}")

# Exact simultaneous opening direction: J * v_int = (1,...,1).
A7 = J[sel, :]
v_int = sp.simplify(A7.solve(sp.ones(7, 1)))
expected_v = sp.Matrix([0, 0, 0, 0, 8 * r15 / 15, 4 * r15 / 5, 8 * r15 / 15])
_require(all(sp.simplify(a - b) == 0 for a, b in zip(v_int, expected_v)), f"v_int mismatch: {v_int}")
_require(sp.simplify(grad_ell.dot(v_int) - 28 * r15 / 15) == 0)
_require(all(sp.simplify(entry - 1) == 0 for entry in J * v_int), f"J*v_int mismatch: {J * v_int}")

if __name__ == "__main__":
    print("[exact] eight primal contacts vanish in Q(sqrt(15))")
    print("[exact] det(A) = -sqrt(15)/2; seven multipliers = (21,3,21,13,4,17,5)*sqrt(15)/45 > 0")
    print("[exact] stubs s_in = s_out = 2*sqrt(15)/5; core = 4*sqrt(15)/3; L_trunc = 32*sqrt(15)/15 (ratio 8/5)")
    print("[exact] opening direction v = (0,0,0,0, 8,12,8)*sqrt(15)/15 satisfies J*v = 1 and grad(core).v = 28*sqrt(15)/15")
    print("PASS: exact KKT, opening ray, and truncation length for the figure-eight state")
