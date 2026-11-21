"""Export shared, width-normalized figure geometry and computed path samples.

Run from any directory. Use --check to detect a stale committed artifact.
The paths are numerical samples of the branch certified in exact_reach.py;
they are not independent continuum certificates. No curve-control points or
rounded screen coordinates are used as geometry.
"""
import argparse
import json
from pathlib import Path

from mpmath import mp
from ribbon_geometry import DaggerField, LAYER_ORDER, family_polynomial, geometry, numeric_state, standard_theta

OUTPUT = Path(__file__).resolve().parent / "artifacts" / "figure_geometry.json"


def check_sample(theta):
    s=geometry(theta)
    v=s["vertices"]
    sub=lambda p,q:(p[0]-q[0],p[1]-q[1])
    cross=lambda p,q:p[0]*q[1]-p[1]*q[0]
    dot=lambda p,q:p[0]*q[0]+p[1]*q[1]
    d2,dD,w0=s["d2"],s["dD"],s["w0"]
    lamI=cross(sub(v["F"],v["G"]),d2)/cross(w0,d2)
    U=sub(v["X"],v["K"])
    den=cross(sub(v["F"],v["E"]),U)
    active=[-cross(v["G"],d2),-cross(sub(v["H"],s["r3G"]),dD),
            theta[3]-lamI,cross(sub(v["E"],s["r3H"]),dD)]
    remaining=[cross(v["X"],d2)-1,cross(sub(v["F"],s["r3G"]),dD)-1,
               cross(v["K"],U)/den-1,cross(v["K"],sub(v["F"],v["E"]))/den-1]
    assert max(abs(z) for z in active)<mp.mpf("1e-45")
    assert min(remaining)>-mp.mpf("1e-45")
    M1,M2,M3=s["centerline"]["core"]
    assert M2[0]>M1[0] and dot(sub(M3,M2),w0)>0
    assert abs(s["lengths"]["core"]-(theta[2]+theta[4]))<mp.mpf("1e-45")
    return max(abs(z) for z in active)


def build_data():
    with mp.workdps(65):
        field=DaggerField()
        exact=field.state()
        dagger_theta=[mp.mpf(str(field.number(t,65))) for t in exact["theta"]]
        x0=mp.sqrt(1-2/mp.sqrt(5))
        pentagon_theta=standard_theta(x0)
        pentagon=numeric_state(pentagon_theta,pentagon=True)
        dagger=numeric_state(dagger_theta)
        assert abs(pentagon["lengths"]["truncation"]-4/(5-2*5**.5)**.5)<1e-13
        assert abs(dagger["notch"]["legXY"]-dagger["notch"]["legYF"])<1e-13

        standard=[]
        for i in range(301):
            p=mp.mpf(i)/300
            x=mp.mpf("0.25")+(x0-mp.mpf("0.25"))*p
            standard.append(dict(progress=float(p),x=float(x),ell=float(1/x-x)))

        second=[]
        h=pentagon_theta[2]
        previous=mp.inf
        residual=mp.mpf(0)
        for i in range(401):
            p=mp.mpf(i)/400
            x=x0+(dagger_theta[0]-x0)*p
            h=mp.findroot(lambda hh:family_polynomial(x,hh),h,tol=mp.mpf("1e-60"))
            assert abs(family_polynomial(x,h))<mp.mpf("1e-55")
            theta=(x,(1-x*x)/(2*x),h,h,(1+x*x)/(2*x))
            residual=max(residual,check_sample(theta))
            ell=h+theta[4]
            assert ell<previous
            previous=ell
            second.append(dict(progress=float(p),x=float(x),h=float(h),ell=float(ell),
                               theta=list(map(float,theta))))
        assert abs(h-dagger_theta[2])<mp.mpf("1e-55")
        assert abs(second[0]["ell"]-pentagon["lengths"]["core"])<1e-13
        assert abs(second[-1]["ell"]-dagger["lengths"]["core"])<1e-13
        return {
            "schemaVersion":1,"width":1,"coordinateConvention":"Cartesian x,y; y increases upward",
            "layerOrder":LAYER_ORDER,
            "source":{"generator":"verify/export_figure_geometry.py", "geometry":"verify/ribbon_geometry.py",
                      "exactOutline":"verify/exact_outline.py", "exactNotch":"verify/exact_notch.py",
                      "pathCertificate":"verify/exact_reach.py",
                      "pathStatus":"Numerical samples of N5=0 on the certified branch; linear interpolation is illustrative",
                      "panelOrder":"A=[E,F,Apositive,Anegative]; B=[E,H,G,F]; C=[H,G,Q1,Q0]; D=[Q0,Q1,Dpositive,Dnegative]",
                      "endOrder":"negative material edge, positive material edge; E/H/Q0 are negative",
                      "cuts":"A is cut at ell_2=line(G,w0); D at y=0. displayPanels extend beyond these cuts."},
            "states":{"pentagon":pentagon,"dagger":dagger},
            "paths":{"standard":standard,"second":second},
            "validation":{"secondFamilySampleCount":len(second),"maxActiveWallResidual":float(residual)},
        }


def serialized_data():
    return json.dumps(build_data(),ensure_ascii=False,indent=2,allow_nan=False)+"\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="compare with an existing export")
    parser.add_argument("--output", type=Path, default=OUTPUT, help="geometry JSON destination")
    args = parser.parse_args()
    content = serialized_data()
    if args.check:
        if args.output.read_text(encoding="utf-8") != content:
            raise SystemExit(f"Stale geometry export: {args.output}")
        print("PASS: geometry export matches current source")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(content, encoding="utf-8")
        print("Wrote", args.output)


if __name__=="__main__":
    main()
