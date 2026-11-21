"""Endpoint and cap data for local equivalence in the actual figure-eight cell."""
import sympy as sp

from n4_exact_state import tight_state
from n4_geometry import forbidden_pairs, halfplanes, open_chord_interval


def equal(label, actual, expected):
    if sp.simplify(actual - expected) != 0:
        raise RuntimeError(label)


def sign(value):
    value = sp.simplify(value)
    if value.is_zero:
        return 0
    if value.is_positive:
        return 1
    if value.is_negative:
        return -1
    raise ArithmeticError(f"undecided exact sign for {value}")


def certify():
    state = tight_state()
    chords = state["CH"]
    planes = [halfplanes(state,j,sign=sign,half=sp.Rational(1,2)) for j in range(5)]
    pairs = sorted((i+1,j) for i,j in forbidden_pairs([1,4,2,0,3]))
    if pairs != [(1,2),(1,4),(2,4),(3,0),(4,0),(4,2)]:
        raise RuntimeError("wrong six forbidden pairs")
    for i,j in pairs:
        if open_chord_interval(chords[i-1],planes[j],sign=sign,zero=sp.S.Zero,one=sp.S.One) is not None:
            raise RuntimeError(f"penetration at crease {i}, pass {j}")

    endpoint_signs = [
        (1,4,(-1,1,1),(0,1,1)), (4,0,(1,-1,1),(1,0,1)),
        (2,4,(1,0,1),(1,0,0)), (3,0,(0,1,1),(0,1,0)),
        (1,2,(1,0,1,0),(1,-1,1,1)), (4,2,(0,1,0,1),(-1,1,1,1)),
    ]
    for crease, panel, minus, plus in endpoint_signs:
        for point, expected in zip(chords[crease-1], (minus,plus), strict=True):
            actual = tuple(sign(a*point[0]+b*point[1]+c) for a,b,c in planes[panel])
            if actual != expected:
                raise RuntimeError(f"wrong cap or side sign at crease {crease}, pass {panel}")

    root = sp.sqrt(15)
    dx, dy = sp.symbols("dx dy")
    delta = sp.Matrix([dx,dy])
    for label, k, u, tail, tail_side, middle_side, middle_cap in [
        ("upper",2,4,4,1,0,2), ("lower",3,1,0,0,1,3),
    ]:
        T,K = chords[k-1]
        U,V = chords[u-1]
        if (K-U).applyfunc(sp.simplify) != sp.zeros(2,1):
            raise RuntimeError(f"{label} endpoints do not coincide")
        u_normal = -sp.Matrix(planes[tail][tail_side][:2])
        v_normal = sp.Matrix(planes[2][middle_side][:2])
        C = sp.Matrix(planes[tail][-1][:2])
        C = C / abs(C[0])
        B = sp.Matrix(planes[2][middle_cap][:2])
        B = B / abs(B[1])
        u_value, v_value = u_normal.dot(delta), v_normal.dot(delta)
        equal(f"{label} tail cap", C.dot(delta), (4*v_value-2*u_value)/root)
        equal(f"{label} pass 2 cap", -B.dot(delta), u_value)
        equal(f"{label} tail side slope", u_normal.dot(T-K), 0)
        equal(f"{label} tail cap slope", C.dot(T-K), 4*root/15)
        equal(f"{label} pass 2 side slope", -v_normal.dot(V-U), sp.Rational(1,2))
        equal(f"{label} pass 2 cap slope", B.dot(V-U), 1)
        if sign(sp.Matrix.hstack(u_normal, v_normal).det()) == 0:
            raise RuntimeError(f"{label} lateral coordinates are dependent")
        print(f"PASS: {label} coupled corner, cap normals and all four directional slopes")
    print("PASS: six actual panel pairs and lateral/cap signs for full-cell local minimality")


if __name__ == "__main__":
    certify()
