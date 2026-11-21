"""Independent rational-coordinate check of the two cap contacts and the cell KKT system.

Plane coordinates are (physical x/sqrt(15), physical y). A unit direction is
(dx, sqrt(15)*dy), and the length coordinates are physical lengths/sqrt(15).
"""
import sympy as sp

R = sp.Rational
STAR = (R(1,5), R(1,15), -R(1,15), -R(1,5), R(2,5), R(8,15), R(2,5))


def equal(label, actual, expected):
    if sp.cancel(actual - expected) != 0:
        raise RuntimeError(label)


def shape(coords):
    """Build rotations and signed chords from cot(turn/2)/sqrt(15)."""
    dirs = [sp.Matrix([1,0])]
    mid = [sp.zeros(2,1)]
    chords = []
    for i, cot in enumerate(coords[:4]):
        dx, dy = dirs[-1]
        den = 15*cot*cot+1
        cosine, sine = (15*cot*cot-1)/den, 2*cot/den
        half = sp.Matrix([cot*dx-dy, 15*cot*dy+dx])/2
        chords.append((mid[i]-half, mid[i]+half))
        rotated = sp.Matrix([cosine*dx-15*sine*dy, sine*dx+cosine*dy]).applyfunc(sp.cancel)
        dirs.append(rotated)
        if i<3:
            mid.append(mid[-1]+coords[4+i]*sp.Matrix([rotated[0],15*rotated[1]]))
    normals = [sp.Matrix([-15*dy,dx]) for dx,dy in dirs]

    def lateral(pass_id, point):
        return normals[pass_id].dot(point-mid[max(0,pass_id-1)])

    half = R(1,2)
    g = [-lateral(4,chords[0][1])-half, lateral(2,chords[0][0])-half,
         lateral(4,chords[1][1])-half, lateral(4,chords[1][0])-half,
         lateral(0,chords[3][1])-half, -lateral(2,chords[3][0])-half,
         -lateral(0,chords[2][1])-half, -lateral(0,chords[2][0])-half]
    return chords, normals, g


def certify():
    chords, normals, contacts = shape(STAR)
    if any(contacts):
        raise RuntimeError("nonzero lateral contact at the hexagon")
    for label, k, u, tail, orientation in [("upper",2,4,4,1),("lower",3,1,0,-1)]:
        T,K = chords[k-1]
        U,V = chords[u-1]
        if K != U:
            raise RuntimeError(f"{label} double point")
        # Tail interior is forward (+X) for panel 4 and backward (-X) for panel 0.
        tail_edge = V-U
        C = orientation*sp.Matrix([1, -tail_edge[0]/tail_edge[1]])
        # The middle panel is below the upper horizontal cap and above the lower one.
        cap_edge = T-K
        B = orientation*sp.Matrix([cap_edge[1]/cap_edge[0], -1])
        un, vn = orientation*normals[tail], orientation*normals[2]
        alpha, beta = sp.Matrix.hstack(vn,un).inv()*C
        equal(label+" tail cap v coefficient", alpha, R(4,15))
        equal(label+" tail cap u coefficient", beta, -R(2,15))
        if -B != un:
            raise RuntimeError(label+" pass 2 cap coordinates")
        equal(label+" lateral slope along KT", un.dot(T-K), 0)
        equal(label+" cap slope along KT", C.dot(T-K), R(4,15))
        equal(label+" lateral slope along UV", -vn.dot(V-U), R(1,2))
        equal(label+" cap slope along UV", B.dot(V-U), 1)

    t = sp.symbols("t")
    jacobian = sp.zeros(8,7)
    for column in range(7):
        coordinates = list(STAR)
        coordinates[column] += t
        for row, g in enumerate(shape(coordinates)[2]):
            jacobian[row,column] = sp.cancel(g/t).subs(t,0)
    selected = jacobian[[0,2,3,4,5,6,7],:]
    equal("rational gradient determinant", selected.det(), -R(284765625,2048))
    gradient = sp.Matrix([0,0,0,0,1,1,1])
    multipliers = selected.T.inv()*gradient
    if list(multipliers) != [R(i,45) for i in (21,3,21,13,4,17,5)]:
        raise RuntimeError("wrong rational-coordinate multipliers")
    print("PASS: both cap contacts, seven independent gradients and positive multipliers over Q")


if __name__ == "__main__":
    certify()
