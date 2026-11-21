"""Shared width-one trefoil geometry for exact certificates and illustrations.

The A tail is cut on ell_2 = line(G, w0), and the D tail on y=0.
The far edge of A is e_B = line(F, d2), NOT the horizontal line y=1.
`geometry` uses only rational operations, so the identical construction works
with floats, mpmath numbers, and reduction in a specified number field.
"""
from functools import lru_cache
import math

LAYER_ORDER = ["B", "D", "A", "C"]


def geometry(theta, reduce=lambda value: value, divide=lambda a, b: a / b):
    R, div = reduce, divide
    x, g, h, a, b = map(R, theta)
    zero, one = R(0), R(1)
    add = lambda u, v: (R(u[0]+v[0]), R(u[1]+v[1]))
    sub = lambda u, v: (R(u[0]-v[0]), R(u[1]-v[1]))
    mul = lambda u, t: (R(u[0]*t), R(u[1]*t))
    dot = lambda u, v: R(u[0]*v[0]+u[1]*v[1])
    cross = lambda u, v: R(u[0]*v[1]-u[1]*v[0])

    def meet(p, u, q, v):
        return add(p, mul(u, div(cross(sub(q, p), v), cross(u, v))))

    def reflect(u, v):
        return sub(mul(u, div(R(2)*dot(u, v), dot(u, u))), v)

    E, F, G, H = (zero, zero), (R(-x), one), (g, one), (h, zero)
    q = R(g-h)
    d2 = (div(R(x*x-one), R(x*x+one)), div(R(-2*x), R(x*x+one)))
    w0 = (div(R(q*q-one), R(q*q+one)), div(R(2*q), R(q*q+one)))
    dA = mul(d2, R(-1))
    K, X = add(G, mul(w0, a)), add(H, mul(w0, b))
    u = sub(X, K)
    dD = reflect(u, w0)
    r3H = add(K, reflect(u, sub(H, K)))
    r3G = add(K, reflect(u, sub(G, K)))
    Y = meet(F, d2, r3H, dD)
    M1 = (div(R(-x), R(2)), div(one, R(2)))
    M2 = (div(R(h+g), R(2)), div(one, R(2)))
    M3 = meet(M2, w0, K, u)
    Anegative, Apositive = meet(E, dA, G, w0), meet(F, dA, G, w0)
    Dnegative = meet(X, dD, E, (one, zero))
    Dpositive = meet(K, dD, E, (one, zero))
    midpoint = lambda p, q: mul(add(p, q), div(one, R(2)))
    Aend, Dend = midpoint(Anegative, Apositive), midpoint(Dnegative, Dpositive)
    core = R(M2[0]-M1[0]+dot(sub(M3, M2), w0))
    truncation = R(core+dot(sub(Aend, M1), dA)+dot(sub(Dend, M3), dD))
    return {
        "theta": (x, g, h, a, b),
        "vertices": dict(E=E, F=F, G=G, H=H, K=K, X=X, Y=Y, Q0=X, Q1=K),
        "panels": {"A": [E,F,Apositive,Anegative], "B": [E,H,G,F],
                   "C": [H,G,K,X], "D": [X,K,Dpositive,Dnegative]},
        "truncatedEnds": {"A": [Anegative,Apositive], "D": [Dnegative,Dpositive]},
        "centerline": {"core": [M1,M2,M3], "truncated": [Aend,M1,M2,M3,Dend]},
        "lengths": {"core": core, "truncation": truncation},
        "d2": d2, "dA": dA, "w0": w0, "dD": dD, "r3H": r3H, "r3G": r3G,
    }


def standard_theta(x):
    g = (1-x*x)/(2*x)
    h = (1-3*x*x)/(2*x)
    return (x,g,h,h,(1+x*x)/(2*x))


def family_polynomial(x, h):
    """N5 from exact_reach.py; choose its branch through the pentagon."""
    return (8*h**3*x**3 + 4*h**2*x**4 + 4*h**2*x**2 + 2*h*x**5
            - 4*h*x**3 - 6*h*x + x**6 + 3*x**4 + 3*x**2 + 1)


def numeric_state(theta, *, pentagon=False):
    """JSON-ready state. Extra display tails are distinct from flush cuts."""
    state = geometry(tuple(float(v) for v in theta))
    v = state["vertices"]
    add = lambda p,d,t: [p[0]+t*d[0], p[1]+t*d[1]]
    # Rendering extents, not truncation lengths. Retain the historical view.
    Aend = [add(v["E"],state["dA"],1.9), add(v["F"],state["dA"],1.9)]
    Dend = [add(v["X"],state["dD"],1.5), add(v["K"],state["dD"],1.5)]
    state["displayEnds"] = {"A": Aend, "D": Dend}
    state["displayPanels"] = {"A": [v["E"],v["F"],Aend[1],Aend[0]],
                              "B": state["panels"]["B"], "C": state["panels"]["C"],
                              "D": [v["X"],v["K"],Dend[1],Dend[0]]}
    mid = lambda p,q: [(p[0]+q[0])/2,(p[1]+q[1])/2]
    state["centerline"]["display"] = [mid(*Aend),*state["centerline"]["core"],mid(*Dend)]
    state["outline"] = ["E","H","G","K","F"] if pentagon else ["E","H","G","K","X","Y","F"]
    state["notch"] = None if pentagon else notch_metrics(v["X"],v["Y"],v["F"])
    for internal in ("d2","dA","w0","dD","r3H","r3G"):
        del state[internal]
    return state


def notch_metrics(X, Y, F):
    xf = (X[0]-F[0],X[1]-F[1])
    mouth2 = xf[0]**2+xf[1]**2
    t = ((Y[0]-F[0])*xf[0]+(Y[1]-F[1])*xf[1])/mouth2
    foot = (F[0]+t*xf[0],F[1]+t*xf[1])
    xy, yf = math.dist(X,Y), math.dist(Y,F)
    cosine = ((X[0]-Y[0])*(F[0]-Y[0])+(X[1]-Y[1])*(F[1]-Y[1]))/(xy*yf)
    return dict(mouth=math.sqrt(mouth2),legXY=xy,legYF=yf,depth=math.dist(Y,foot),
                reflexDegrees=360-math.degrees(math.acos(max(-1,min(1,cosine)))),foot=foot)


class DaggerField:
    """Q[x]/p with exact rational root isolation and rational interval signs.

    No floating-point root estimate is used to certify signs. Sturm root count
    identifies one real root in the stated rational interval, and exact
    bisection refines that interval. Horner interval arithmetic is rational.
    """
    def __init__(self):
        import sympy as sp
        self.sp = sp
        self.x = sp.Symbol("x")
        x = self.x
        self.p = x**12+4*x**10-7*x**4+20*x**2-2
        self.poly = sp.Poly(self.p,x)
        self.lo = sp.Rational(32212749731122717,10**17)
        self.hi = sp.Rational(32212749731122719,10**17)
        assert self.poly.count_roots(self.lo,self.hi) == 1
        assert self.poly.eval(self.lo)*self.poly.eval(self.hi) < 0
        while self.hi-self.lo > sp.Rational(1,10**65):
            self.refine()

    def refine(self):
        mid = (self.lo+self.hi)/2
        if self.poly.eval(self.lo)*self.poly.eval(mid) < 0:
            self.hi = mid
        else:
            self.lo = mid

    @lru_cache(maxsize=None)
    def red(self, expr):
        return self.sp.rem(self.sp.expand(expr),self.p,self.x)

    def div(self,a,b):
        return self.red(self.red(a)*self.sp.invert(self.red(b),self.p,self.x))

    @lru_cache(maxsize=None)
    def sign(self, expr):
        expr = self.red(expr)
        if expr == 0:
            return 0
        coeffs = self.sp.Poly(expr,self.x).all_coeffs()
        for _ in range(256):
            low = high = self.sp.Rational(0)
            for c in coeffs:
                products = [low*self.lo,low*self.hi,high*self.lo,high*self.hi]
                low,high = min(products)+c,max(products)+c
            if low > 0: return 1
            if high < 0: return -1
            self.refine()
        raise ArithmeticError("Could not separate nonzero field element from zero")

    def number(self, expr, digits=50):
        return self.red(expr).subs(self.x,(self.lo+self.hi)/2).evalf(digits)

    def state(self):
        x,R = self.x,self.red
        h = R((2*x**11+8*x**9-x**7-4*x**5-13*x**3+40*x)/12)
        return geometry((x,self.div(1-x*x,2*x),h,h,self.div(1+x*x,2*x)),R,self.div)
