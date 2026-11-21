"""Draw the loose and tight figure-eight states in top and hidden-line views."""
import argparse
import json
import math
from pathlib import Path
from svg_labels import math_text, save_svg
from n4_exact_state import tight_state
from n4_tighten2 import config
from n4_render import band_quads
from shapely.geometry import Polygon, LineString
from shapely.ops import unary_union
HATCH_ANG={0:0,1:45,2:90,3:135,4:20}
def hatch_segs(T, quad, ang, gap):
    P=Polygon([T(p) for p in quad])
    if not P.is_valid: P=P.buffer(0)
    minx,miny,maxx,maxy=P.bounds
    cx,cy=(minx+maxx)/2,(miny+maxy)/2
    L=math.hypot(maxx-minx,maxy-miny)
    dx,dy=math.cos(math.radians(ang)),math.sin(math.radians(ang))
    nx,ny=-dy,dx
    out=[]
    k0=int(-L/gap)-1; k1=int(L/gap)+1
    for k in range(k0,k1+1):
        ox,oy=cx+k*gap*nx, cy+k*gap*ny
        seg=LineString([(ox-L*dx,oy-L*dy),(ox+L*dx,oy+L*dy)]).intersection(P)
        gs=getattr(seg,'geoms',[seg]) if not seg.is_empty else []
        for g in gs:
            if g.geom_type=='LineString' and g.length>1e-6:
                (a,b)=(g.coords[0],g.coords[-1]); out.append((a[0],a[1],b[0],b[1]))
    return out

cols={0:"#8ecae6",1:"#219ebc",2:"#ffb703",3:"#fb8500",4:"#c1121f"}

def render_pair(th, order, out, hexoverlay=None, doubles=None, tails=False, labels_b=None, check=False, pdf=False):
    c = th if isinstance(th, dict) else config(tuple(th))
    q = band_quads(c)
    PW,PH=560,520; GAP=36; W=PW*2+GAP+30; H=PH+56
    pts=[p for k in q for p in q[k]]
    if hexoverlay: pts=pts+hexoverlay
    xs=[p[0] for p in pts]; ys=[p[1] for p in pts]
    x0,x1,y0,y1=min(xs),max(xs),min(ys),max(ys)
    pad=0.35; sc=min((PW-26)/(x1-x0+2*pad),(PH-40)/(y1-y0+2*pad))
    def mkT(ox): return lambda p:(ox+13+(p[0]-x0+pad)*sc, 18+PH-20-(p[1]-y0+pad)*sc)
    def ps(T,poly): return " ".join(f"{T(p)[0]:.1f},{T(p)[1]:.1f}" for p in poly)
    S=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">','<rect width="100%" height="100%" fill="white"/>',
       '<defs><marker id="ar" markerWidth="9" markerHeight="7" refX="8" refY="3.5" orient="auto"><polygon points="0 0, 9 3.5, 0 7" fill="#155"/></marker></defs>']
    b2t=list(reversed(order))
    def chords(T):
        for a,b in c['CH']:
            A,B=T(a),T(b)
            S.append(f'<line x1="{A[0]:.1f}" y1="{A[1]:.1f}" x2="{B[0]:.1f}" y2="{B[1]:.1f}" stroke="#111" stroke-width="2.4"/>')
    def mdots(T):
        for i,M in enumerate(c['M']):
            A=T(M); S.append(f'<circle cx="{A[0]:.1f}" cy="{A[1]:.1f}" r="4" fill="#d00"/>')
            if hexoverlay:
                dx,dy,anchor=[(13,7,'start'),(0,-15,'middle'),(0,29,'middle'),(-14,4,'end')][i]
                S.append(math_text(A[0]+dx,A[1]+dy,f'M_{i+1}',24,anchor,halo=2.5))
            else:
                S.append(f'<text x="{A[0]+7:.0f}" y="{A[1]-7:.0f}" font-family="Georgia" font-style="italic" font-size="24">M{i+1}</text>')
    def hexdraw(T, with_labels):
        if not hexoverlay: return
        S.append(f'<polygon points="{ps(T,hexoverlay)}" fill="none" stroke="#0a0a0a" stroke-width="3.4"/>')
        for p in hexoverlay:
            dbl = doubles and any(abs(p[0]-d[0])<1e-9 and abs(p[1]-d[1])<1e-9 for d in doubles)
            A=T(p)
            S.append(f'<circle cx="{A[0]:.1f}" cy="{A[1]:.1f}" r="{6.5 if dbl else 5}" fill="{"#d00000" if dbl else "#111"}" stroke="white" stroke-width="1.5"/>')
        if with_labels and labels_b:
            for lbl,p,dx,dy,anchor in labels_b:
                A=T(p)
                S.append(math_text(A[0]+dx,A[1]+dy,lbl,27,anchor,halo=3))
    def tarrows(T):
        if not tails: return
        A=T((x0-0.1,0.0)); B=T((x0+0.7,0.0))
        S.append(f'<line x1="{B[0]:.0f}" y1="{B[1]:.0f}" x2="{A[0]+24:.0f}" y2="{A[1]:.0f}" stroke="#155" stroke-width="2" marker-end="url(#ar)"/>')
        A2=T((x1+0.1,1.0)); B2=T((x1-0.7,1.0))
        S.append(f'<line x1="{B2[0]:.0f}" y1="{B2[1]:.0f}" x2="{A2[0]-24:.0f}" y2="{A2[1]:.0f}" stroke="#155" stroke-width="2" marker-end="url(#ar)"/>')
    # (a) opaque
    Ta=mkT(10)
    for k in b2t:
        S.append(f'<polygon points="{ps(Ta,q[k])}" fill="none" stroke="white" stroke-width="9"/>')
        S.append(f'<polygon points="{ps(Ta,q[k])}" fill="{cols[k]}" fill-opacity="0.95" stroke="{cols[k]}" stroke-width="2.2"/>')
    chords(Ta); hexdraw(Ta, False); tarrows(Ta); mdots(Ta)
    S.append(f'<text x="{10+PW/2}" y="{H-10}" font-family="Georgia" font-size="18" text-anchor="middle">(a)</text>')
    # (b) hidden-line
    Tb=mkT(10+PW+GAP)
    for k in b2t:
        S.append(f'<polygon points="{ps(Tb,q[k])}" fill="{cols[k]}" fill-opacity="0.22" stroke="none"/>')
        for hx1,hy1,hx2,hy2 in hatch_segs(Tb,q[k],HATCH_ANG[k],9.0):
            S.append(f'<line x1="{hx1:.1f}" y1="{hy1:.1f}" x2="{hx2:.1f}" y2="{hy2:.1f}" stroke="#333" stroke-width="0.8" opacity="0.32"/>')
    for depth,k in enumerate(order):
        above=[Polygon(q[j]) for j in order[:depth]]
        U=unary_union(above) if above else None
        quad=q[k]
        for i in range(4):
            e=LineString([quad[i],quad[(i+1)%4]])
            vis=e.difference(U) if U else e
            hid=e.intersection(U) if U else None
            def emit(geom,style):
                if geom is None or geom.is_empty: return
                for g in getattr(geom, "geoms", (geom,)):
                    if g.geom_type != "LineString" or g.length < 1e-9: continue
                    d=" ".join(f"{Tb(p)[0]:.1f},{Tb(p)[1]:.1f}" for p in g.coords)
                    S.append(f'<polyline points="{d}" fill="none" {style}/>')
            emit(vis, f'stroke="{cols[k]}" stroke-width="2.4"')
            emit(hid, f'stroke="{cols[k]}" stroke-width="1.5" stroke-dasharray="6,4" stroke-opacity="0.85"')
    chords(Tb); hexdraw(Tb, True); tarrows(Tb); mdots(Tb)
    S.append(f'<text x="{10+PW+GAP+PW/2}" y="{H-10}" font-family="Georgia" font-size="18" text-anchor="middle">(b)</text>')
    S.append('</svg>')
    save_svg('\n'.join(S),out,check=check,pdf=pdf)


def tight_drawing():
    exact = tight_state()

    def point(vector):
        return tuple(map(float, vector))

    # Convert each exact endpoint once so coincident boundaries stay coincident
    # during clipping. The plotting convention reverses the two negative turns.
    chords = [(point(a), point(b)) if i < 2 else (point(b), point(a))
              for i, (a, b) in enumerate(exact["CH"])]
    (e1m, e1p), (e2m, e2p), (_, e3p), (e4m, _) = chords
    corners = [e4m, e2p, e2m, e1p, e1m, e3p]
    state = {
        "d": [point(v) for v in exact["d"]],
        "M": [point(v) for v in exact["M"]],
        "L": list(map(float, exact["L"])),
        "CH": chords,
    }
    return state, corners

if __name__=="__main__":
    parser=argparse.ArgumentParser(description='Render paired figure-eight diagrams, retaining the stored geometry.')
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true')
    group.add_argument('--pdf',action='store_true')
    parser.add_argument('--only',choices=('loose','tight','all'),default='all')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    if args.only in ('loose','all'):
        th,order=json.loads((root/'verify/n4_valid41.json').read_text())[1]
        render_pair(th, order, root/'paper/figs/fig41-loose.svg',check=args.check,pdf=args.pdf)
    if args.only == 'loose':
        raise SystemExit(0)
    state, HEX = tight_drawing()
    DBL = [HEX[1], HEX[4]]
    LB=[('E_4^{−}',HEX[0],-13,0,'end'),('E_4^{+} = E_2^{+}',HEX[1],-14,-28,'end'),
        ('E_2^{−}',HEX[2],12,-17,'start'),('E_1^{+}',HEX[3],14,19,'start'),
        ('E_1^{−} = E_3^{−}',HEX[4],8,35,'start'),('E_3^{+}',HEX[5],-13,30,'end')]
    render_pair(state, [1,4,2,0,3], root/'paper/figs/fig41-tight.svg', hexoverlay=HEX, doubles=DBL, tails=True, labels_b=LB,check=args.check,pdf=args.pdf)
