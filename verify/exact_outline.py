"""Exact full-silhouette certificate for the correctly flush-cut dagger state.

All signs use a Sturm-isolated root and rational interval arithmetic. No
midpoint sampling or floating geometry is used. Convex clipping and exact
area identities prove full polygon unions, and boundary cancellation proves
the stated simple heptagon. A linear separating half-plane for each infinite
pass excludes EVERY point of the open notch triangle, including the tails.
"""
from ribbon_geometry import DaggerField


def certify():
    field = DaggerField()
    R, div, sign = field.red, field.div, field.sign
    state = field.state()
    V = state["vertices"]
    add = lambda p,q: (R(p[0]+q[0]),R(p[1]+q[1]))
    sub = lambda p,q: (R(p[0]-q[0]),R(p[1]-q[1]))
    mul = lambda p,t: (R(p[0]*t),R(p[1]*t))
    cross = lambda p,q: R(p[0]*q[1]-p[1]*q[0])
    dot = lambda p,q: R(p[0]*q[0]+p[1]*q[1])
    orient = lambda a,b,p: cross(sub(b,a),sub(p,a))
    equal = lambda p,q: R(p[0]-q[0]) == 0 and R(p[1]-q[1]) == 0

    def clean(poly):
        result=[]
        for p in poly:
            if not result or not equal(result[-1],p): result.append(p)
        if len(result)>1 and equal(result[0],result[-1]): result.pop()
        while len(result)>=3:
            for i in range(len(result)):
                if orient(result[i-1],result[i],result[(i+1)%len(result)]) == 0:
                    result.pop(i)
                    break
            else: break
        return result

    def area2(poly):
        if len(poly)<3: return R(0)
        return R(sum(cross(a,b) for a,b in zip(poly,poly[1:]+poly[:1])))

    def convex(poly):
        poly=clean(poly)
        if sign(area2(poly))<0: poly.reverse()
        assert len(poly)>=3 and all(sign(orient(poly[i-1],poly[i],poly[(i+1)%len(poly)]))>0
                                    for i in range(len(poly))), "Nonconvex or degenerate polygon"
        assert all(all(sign(orient(a,b,p))>=0 for p in poly)
                   for a,b in zip(poly,poly[1:]+poly[:1])), "Polygon is not its convex hull"
        return poly

    def contains(poly,p):
        return all(sign(orient(a,b,p))>=0 for a,b in zip(poly,poly[1:]+poly[:1]))

    def clip(subject,clipping):
        """Sutherland-Hodgman clipping, entirely in the exact number field."""
        result=subject[:]
        for a,b in zip(clipping,clipping[1:]+clipping[:1]):
            if not result: break
            output=[]
            previous=result[-1]
            fp=orient(a,b,previous)
            for current in result:
                fc=orient(a,b,current)
                if (sign(fp)>=0)!=(sign(fc)>=0):
                    t=div(fp,R(fp-fc))
                    output.append(add(previous,mul(sub(current,previous),t)))
                if sign(fc)>=0: output.append(current)
                previous,fp=current,fc
            result=clean(output)
        return clean(result)

    def same_polygon(p,q):
        p,q=clean(p),clean(q)
        return len(p)==len(q) and any(all(equal(p[i],q[(i+j)%len(q)]) for i in range(len(p)))
                                     for j in range(len(q)))

    def union_is(a,b,target):
        # Convexity makes vertex containment sufficient for whole-polygon
        # containment. Exact area equality then excludes holes: finite unions
        # of nondegenerate closed polygons are closures of their interiors.
        assert all(contains(target,p) for p in a+b), "Union escapes proposed polygon"
        assert R(area2(a)+area2(b)-area2(clip(a,b))-area2(target)) == 0, "Union leaves a hole"

    # The computed flush intersections must agree with these exact vertices.
    expected={"A":["E","F","K","G"],"B":["E","H","G","F"],
              "C":["H","G","K","X"],"D":["X","K","H","E"]}
    panels={}
    for name,names in expected.items():
        assert all(equal(actual,V[label]) for actual,label in zip(state["panels"][name],names)), name
        panels[name]=convex(state["panels"][name])
    print("[exact] Four flush-cut panels: A=EFKG, B=EHGF, C=HGKX, D=XKHE")

    P1=convex([V[n] for n in ["E","H","G","K","F"]])
    P2=convex([V[n] for n in ["E","H","G","K","X"]])
    I=convex([V[n] for n in ["E","H","G","K","Y"]])
    T1=convex([V[n] for n in ["E","Y","F"]])
    T2=convex([V[n] for n in ["K","X","Y"]])
    union_is(panels["A"],panels["B"],P1)
    union_is(panels["C"],panels["D"],P2)
    assert same_polygon(clip(P1,P2),I), "Wrong intersection of convex pentagons"
    union_is(I,T1,P1)
    union_is(I,T2,P2)
    assert all(area2(clip(a,b)) == 0 for a,b in [(I,T1),(I,T2),(T1,T2)])
    print("[exact] A union B = EHGKF; C union D = EHGKX; their intersection = EHGKY")

    # Cancel every internal edge of I + EYF + KXY. This is a full planar
    # subdivision, rather than a collection of probes on proposed edges.
    boundary=set()
    for names in [["E","H","G","K","Y"],["E","Y","F"],["K","X","Y"]]:
        for a,b in zip(names,names[1:]+names[:1]):
            if (b,a) in boundary: boundary.remove((b,a))
            else: boundary.add((a,b))
    names=["E","H","G","K","X","Y","F"]
    assert boundary == set(zip(names,names[1:]+names[:1]))
    outline=[V[n] for n in names]
    assert all(not equal(a,b) for i,a in enumerate(outline) for b in outline[i+1:])

    def on_segment(a,b,p):
        return orient(a,b,p) == 0 and sign(dot(sub(p,a),sub(p,b)))<=0

    def segments_meet(a,b,c,d):
        signs=[sign(orient(a,b,c)),sign(orient(a,b,d)),sign(orient(c,d,a)),sign(orient(c,d,b))]
        return ((signs[0]*signs[1]<0 and signs[2]*signs[3]<0)
                or on_segment(a,b,c) or on_segment(a,b,d) or on_segment(c,d,a) or on_segment(c,d,b))

    for i in range(7):
        for j in range(i+1,7):
            if j==i+1 or (i==0 and j==6): continue
            assert not segments_meet(outline[i],outline[(i+1)%7],outline[j],outline[(j+1)%7])
    turns={n:sign(orient(outline[i-1],outline[i],outline[(i+1)%7])) for i,n in enumerate(names)}
    assert turns == {n:(-1 if n=="Y" else 1) for n in names}
    assert R(area2(outline)-area2(I)-area2(T1)-area2(T2)) == 0
    print("[exact] Complete simple boundary E-H-G-K-X-Y-F; sole reflex vertex Y")

    # For each untruncated pass, omit only its artificial flush-cut edge.
    # These directed edges have their pass on the left. Linearity then
    # separates the entire open triangle, not just its barycentre.
    untruncated={
        "A":[("E","G"),("K","F"),("F","E")],
        "B":[("E","H"),("H","G"),("G","F"),("F","E")],
        "C":[("H","G"),("G","K"),("K","X"),("X","H")],
        "D":[("H","K"),("K","X"),("X","E")],
    }
    triangle=[V[n] for n in ["X","Y","F"]]
    assert sign(area2(triangle)) != 0
    for name,edges in untruncated.items():
        assert all(all(sign(orient(V[a],V[b],p))>=0 for p in panels[name]) for a,b in edges)
        certificates=[]
        for a,b in edges:
            values=[sign(orient(V[a],V[b],p)) for p in triangle]
            if max(values)<=0 and min(values)<0: certificates.append((a+b,values))
        assert certificates, f"No complete open-notch separation from pass {name}"
        print(f"[exact] Open triangle XYF disjoint from untruncated {name}: {certificates[0]}")
    print("PASS: full four-panel silhouette and whole open notch, with exact rational signs")
    return True


if __name__ == "__main__":
    certify()
