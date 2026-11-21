"""Exact notch measurements from the same geometry used by the figures."""
from ribbon_geometry import DaggerField


def certify():
    field=DaggerField()
    x,R,div=field.x,field.red,field.div
    V=field.state()["vertices"]
    X,Y,F=V["X"],V["Y"],V["F"]
    sub=lambda a,b:(R(a[0]-b[0]),R(a[1]-b[1]))
    dot=lambda a,b:R(a[0]*b[0]+a[1]*b[1])
    cross=lambda a,b:R(a[0]*b[1]-a[1]*b[0])
    XF,XY,FY=sub(X,F),sub(X,Y),sub(F,Y)
    mouth2=dot(XF,XF)
    leg2=dot(XY,XY)
    assert R(leg2-dot(FY,FY)) == 0, "The notch must be isosceles"
    depth2=div(R(cross(XF,sub(Y,F))**2),mouth2)
    cosine=div(dot(XY,FY),leg2)
    assert R(depth2-leg2+mouth2/4) == 0, "Isosceles depth identity"
    Q=field.sp.Rational
    expected={
        "mouth_squared":-x**10/6-Q(11,12)*x**8-Q(5,4)*x**6-Q(7,12)*x**4+Q(7,4)*x**2-Q(1,6),
        "leg_squared":-x**10/48-x**8/6-Q(5,24)*x**6-x**4/3+x**2/16,
        "depth_squared":x**10/48+x**8/16+Q(5,48)*x**6-Q(3,16)*x**4-Q(3,8)*x**2+Q(1,24),
        "cosine_reflex":-x**10/24-x**8/8+x**6/8-Q(11,24)*x**4-Q(11,12)*x**2-Q(1,4),
    }
    actual=dict(mouth_squared=mouth2,leg_squared=leg2,depth_squared=depth2,cosine_reflex=cosine)
    for name,value in actual.items():
        assert R(value-expected[name]) == 0, name
        assert all(degree[0]%2 == 0 for degree,_ in field.sp.Poly(value,x).terms()), name
        print(f"[exact] {name} = {value}")
    assert all(field.sign(v)>0 for v in [mouth2,leg2,depth2])
    assert field.sign(cosine+1)>0 and field.sign(1-cosine)>0
    import math
    print("mouth =",math.sqrt(float(field.number(mouth2))))
    print("leg =",math.sqrt(float(field.number(leg2))))
    print("depth =",math.sqrt(float(field.number(depth2))))
    print("reflex angle (degrees) =",360-math.degrees(math.acos(float(field.number(cosine)))))
    print("PASS: exact formulas, equal legs, depth identity, positivity and sextic subfield")
    return True


if __name__ == "__main__":
    certify()
