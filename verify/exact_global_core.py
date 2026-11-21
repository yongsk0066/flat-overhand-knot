"""Four necessary trefoil walls, the sharp length bound and its unique equality case."""
import sympy as sp

from ribbon_geometry import DaggerField, geometry


def equal(label, actual, expected):
    if sp.cancel(actual - expected) != 0:
        raise RuntimeError(label)


def det(u, v):
    return u[0] * v[1] - u[1] * v[0]


def certify():
    x, q, z, h, a, m, k, s = sp.symbols("x q z h a m k s")
    n = 1 + q**2

    def T(y):
        return (q + y) * (1 - q*y)

    def A(y):
        return n*y / T(y)

    def side(y):
        return (q**2 + 2*q*y - 1) / n

    chart = geometry((x, h+q, h, a, a+q+z), sp.factor)
    vertices = {key: sp.Matrix(pt) for key, pt in chart["vertices"].items()}
    E, F, G, H, Q0, Q1 = (vertices[key] for key in ("E", "F", "G", "H", "Q0", "Q1"))
    dA, dC, dD = (sp.Matrix(chart[key]) for key in ("d2", "w0", "dD"))
    e0, e1 = (sp.Matrix(chart[key]) for key in ("r3H", "r3G"))
    U = Q0 - Q1
    D = det(F-E, U)
    w5 = a - det(F-G, dA)/det(dC, dA)
    left = a - A(x)*(h+q+x)
    right = h - A(z)*(a+q+z)
    w1, w4 = -det(G-E, dA), -det(H-e1, dD)
    w6, w7, w8 = det(E-e0, dD), det(Q0-E, dA)-1, det(F-e1, dD)-1
    equal("incoming side", w1, (2*x*(h+q)+x*x-1)/(1+x*x))
    equal("outgoing side", w4, (2*z*(a+q)+z*z-1)/(1+z*z))
    equal("incoming sliver", w5, left)
    equal("outgoing sliver", w6, 2*T(z)*right/(n*(1+z*z)))
    equal("centreline intersection", det(sp.Matrix([h+(q+x)/2,0]),-dA)/det(-dA,dC),
          x*n*(h+(q+x)/2)/T(x))
    equal("far side of incoming strip", det(Q1-E, dA)-1, 2*T(x)*left/(n*(1+x*x)))
    equal("positive determinant identity", T(x)*T(z)-x*z*n*n,
          q*((1-x*z)*n*D-q*(1+x*x)*(1+z*z)))
    equal("sum of wall ratios", T(x)/x+T(z)/z,
          2*(1-q*q)+q*(x+z)*(1/(x*z)-1))

    t_num, u_num = det(Q1-E,U), det(Q1-E,F-E)
    equal("t numerator", t_num, a+q+z+side(z)*h)
    equal("u numerator", u_num, h+q+x+side(x)*a)
    equal("first taco wall", t_num-D, a+side(z)*(h+q+x))
    equal("second taco wall", u_num-D, h+side(x)*(a+q+z))
    for name, wall, rate in [
        ("w1", w1, 2*x*h/(1+x*x)), ("w4", w4, 2*z*a/(1+z*z)),
        ("w5", w5, a-A(x)*h), ("w6", w6, 2*T(z)*(h-A(z)*a)/(n*(1+z*z))),
        ("w7", w7, 2*T(x)*(a-A(x)*h)/(n*(1+x*x))),
        ("w8", w8, 2*T(z)*(h-A(z)*a)/(n*(1+z*z))),
        ("t", t_num, a+side(z)*h), ("u", u_num, h+side(x)*a),
    ]:
        equal(f"length dilation of {name}", h*wall.diff(h)+a*wall.diff(a), rate)

    change = {h: (m-2*q-x)/2, a: (k-2*q-z)/2}
    equal("core", chart["lengths"]["core"].subs(change), (m+k)/2-q)
    equal("left angle bound", w1.subs(change), (m*x-1)/(1+x*x))
    equal("right angle bound", w4.subs(change), (k*z-1)/(1+z*z))
    equal("coupled length", 2*left.subs(change), k-2*q-z-A(x)*(m+x))
    equal("reversed coupled length", 2*right.subs(change), m-2*q-x-A(z)*(k+z))
    equal("angle derivative", A(x).diff(x), q*n*(1+x*x)/T(x)**2)
    Fm = n*(1/(m-q)+m/(q*m+1))
    Lm = m-2*q-1/m
    equal("angle relaxation", A(1/m)*(m+1/m), Fm)
    equal("length derivative", Fm.diff(m), n*((q*m+1)**-2-(m-q)**-2))
    equal("length divided difference", (Fm.subs(m,k)-Fm)/(k-m),
          n*(1/((q*k+1)*(q*m+1))-1/((k-q)*(m-q))))

    P = -q**5-4*s*q**4-(4*s**2+2)*q**3-6*s*q**2+(s**4-4*s**2-1)*q-2*s
    equal("reduced length polynomial", ((Lm-Fm)*m*(q*m+1)*(m-q)).subs(m,s+q), P)
    equal("strict concavity", -P.diff(q,2), 20*q**3+48*s*q**2+(24*s**2+12)*q+12*s)
    length_poly = 8*s**12-88*s**10+267*s**8-456*s**6+88*s**4-80*s**2-16
    equal("length discriminant", sp.discriminant(P,q), -16*s**4*length_poly)
    sturm = sp.Poly(length_poly, s)
    if sturm.count_roots(0, sp.oo) != 1 or sturm.count_roots(2, 3) != 1:
        raise RuntimeError("wrong positive length root count")
    if any(c >= 0 for c in sp.Poly(P.subs(s,1), q).all_coeffs()):
        raise RuntimeError("wrong sign anchor for the discriminant proof")

    field = DaggerField()
    if field.poly.count_roots(field.lo, field.hi) != 1:
        raise RuntimeError("dagger abscissa was not isolated")
    state = field.state()
    x0, g0, h0, a0, b0 = state["theta"]
    q0, s0 = field.red(g0-h0), state["lengths"]["core"]
    if any(field.sign(v) <= 0 for v in (q0,h0,s0,1-q0*x0)):
        raise RuntimeError("dagger chart signs")
    for label, value in [
        ("corner slope", b0-a0-q0-x0), ("symmetric lengths", a0-h0),
        ("angle reconstruction", (s0+q0)*x0-1),
        ("length reconstruction", 2*h0+2*q0+x0-(s0+q0)),
        ("length root", length_poly.subs(s,s0)),
        ("attainment", P.subs({s:s0,q:q0})),
        ("supporting tangent", P.diff(q).subs({s:s0,q:q0})),
    ]:
        if field.red(value) != 0:
            raise RuntimeError(label)
    print("PASS: four geometric walls, positive carrier orientation and length dilation")
    print("PASS: global core bound, unique equality case and independent discriminant/Sturm check")


if __name__ == "__main__":
    certify()
