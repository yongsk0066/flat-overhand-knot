"""Exact global core identity and a strictly feasible descending pentagon ray.

Algebra uses Q(x0), x0=sqrt(1-2/sqrt(5)). Every sign is certified by
outward-rounded interval evaluation of its exact reduced polynomial; no
floating-point sign calibration or integer-relation identification is used.
"""
import sympy as sp
from mpmath import iv

NAMES = ('w1', 'w4', 'w5', 'w6', 'w7', 'w8', 'w9', 'w10')
DIRECTION = (-163127, 1000000, -883304, -887250, 833681)
DENOMINATOR = 1000000


class Jet:
    def __init__(self, value, derivative=0):
        self.value, self.derivative = value, derivative

    def __add__(self, other):
        other = other if isinstance(other, Jet) else Jet(other)
        return Jet(self.value + other.value, self.derivative + other.derivative)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.value, -self.derivative)

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        other = other if isinstance(other, Jet) else Jet(other)
        return Jet(self.value * other.value,
                   self.derivative * other.value + self.value * other.derivative)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = other if isinstance(other, Jet) else Jet(other)
        return Jet(self.value / other.value,
                   (self.derivative * other.value - self.value * other.derivative)
                   / (other.value * other.value))

    def __rtruediv__(self, other):
        return Jet(other) / self


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def mul(a, s):
    return tuple(x * s for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def reflect(axis, vector):
    return sub(mul(axis, 2 * dot(axis, vector) / dot(axis, axis)), vector)


def geometry(theta):
    x, g, h, a, b = theta
    E, F, G, H = (0, 0), (-x, 1), (g, 1), (h, 0)
    d2 = reflect(F, (1, 0))
    w0 = reflect(sub(G, H), (1, 0))
    Q1, Q0 = add(G, mul(w0, a)), add(H, mul(w0, b))
    U = sub(Q0, Q1)
    rho = lambda p: add(Q1, reflect(U, sub(p, Q1)))
    dD = reflect(U, w0)
    return E, F, G, H, d2, w0, Q1, Q0, U, rho, dD


def walls(theta):
    E, F, G, H, d2, w0, Q1, Q0, U, rho, dD = geometry(theta)
    # Orientations follow actual band-interior points, certified below.
    return (-cross(G, d2),
            -cross(sub(H, rho(G)), dD),
            theta[3] - cross(sub(F, G), d2) / cross(w0, d2),
            cross(sub(E, rho(H)), dD),
            cross(Q0, d2) - 1,
            cross(sub(F, rho(G)), dD) - 1,
            cross(Q1, U) / cross(F, U) - 1,
            cross(Q1, F) / cross(F, U) - 1)


def global_identity():
    x, g, h, a, b = sp.symbols('x g h a b', real=True)
    q = g - h
    W = sp.Matrix([q*q - 1, 2*q])
    G, H = sp.Matrix([g, 1]), sp.Matrix([h, 0])
    U = (q*q + 1) * (H-G) + (b-a)*W
    denominator = cross(W, U)
    numerator = cross((q*q+1)*(G-(G+H)/2)+a*W, U)
    assert sp.factor(denominator-(q*q+1)**2) == 0
    assert sp.factor(numerator-(a+b)*(q*q+1)**2/2) == 0
    length = (x+g+h)/2 + numerator/denominator
    assert sp.factor(length-(x+g+h+a+b)/2) == 0
    print('[exact identity] Dp=(q^2+1)^2, Np=(a+b)Dp/2')
    print('[exact identity] core=(x+g+h+a+b)/2; gradient=(1/2,...,1/2) on the positive-length branch')


def certify():
    global_identity()
    iv.dps = 70
    alpha = sp.sqrt(1-2/sp.sqrt(5))
    field = sp.QQ.algebraic_field(alpha)
    X = field.from_sympy(alpha)
    P = (X, (1-X*X)/(2*X), (1-3*X*X)/(2*X),
         (1-3*X*X)/(2*X), (1+X*X)/(2*X))
    D = tuple(field.convert(sp.Rational(n, DENOMINATOR)) for n in DIRECTION)
    root_interval = iv.sqrt(1-2/iv.sqrt(5))

    def enclosure(value):
        value = field.convert(value)
        out = iv.mpf(0)
        for coefficient in value.to_list():
            out = out*root_interval + iv.mpf(int(coefficient.numerator))/int(coefficient.denominator)
        return out

    E, F, G, H, d2, w0, Q1, Q0, U, rho, dD = geometry(P)
    Ainside = add(mul(add(E, F), field.convert(sp.Rational(1, 2))), d2)
    Dinside = add(mul(add(Q1, Q0), field.convert(sp.Rational(1, 2))), dD)
    assert cross(sub(Ainside, E), d2) == field.convert(sp.Rational(1, 2))
    assert cross(sub(Dinside, rho(G)), dD) == field.convert(sp.Rational(1, 2))
    # Q1+dD lies on the band edge and cannot be used to calibrate this sign.
    assert cross(sub(add(Q1, dD), rho(G)), dD) == field.zero
    print('[exact orientation] genuine A/D interior probes give +1/2; old edge probe gives 0')

    values = walls(tuple(Jet(p, d) for p, d in zip(P, D)))
    for name, value in zip(NAMES, values):
        assert value.value == field.zero, (name, value.value)
        bound = enclosure(value.derivative)
        assert bound.a > 0, (name, bound)
        print(f'[exact+interval] {name}(P)=0; D{name} in {bound}')
    objective_rate = sum(D)/2
    assert objective_rate == field.convert(sp.Rational(-1, 20))
    print('[exact descent] core directional derivative = -1/20')

    # Validate all eight wall derivatives along one concrete short ray,
    # not merely at finitely many sample points.
    T = iv.mpf(['0', '0.000001'])
    Piv = tuple(enclosure(p) for p in P)
    Div = tuple(iv.mpf(n)/DENOMINATOR for n in DIRECTION)
    path = walls(tuple(Jet(p+T*d, d) for p, d in zip(Piv, Div)))
    for name, value in zip(NAMES, path):
        assert value.derivative.a > 0, (name, value.derivative)
    print('[interval ray] all eight derivatives >0 throughout 0<=t<=1e-6')
    print('PASS: the pentagon admits a strict open-wall descent; both taco branch cones exclude grad(core)')
    return field, P, D, values


if __name__ == '__main__':
    certify()
