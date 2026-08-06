# tref_render.py — 세잎 본판: 은선(hidden-line) 단일 스타일 fig1(오각형)·fig2(칠각형+인셋) — trefoil master version: single hidden-line style, fig1 (pentagon) and fig2 (heptagon + inset)
import sys, math
sys.path.insert(0,'verify')
from uniq_setup import _config, _slice_theta, _add, _mul
from shapely.geometry import Polygon, LineString
from shapely.ops import unary_union
import math as _math
HATCH_ANG={'A':0,'B':45,'C':90,'D':135}
def hatch_segs(P, ang, gap=9.0):
    if not P.is_valid: P=P.buffer(0)
    minx,miny,maxx,maxy=P.bounds
    cx,cy=(minx+maxx)/2,(miny+maxy)/2
    L=_math.hypot(maxx-minx,maxy-miny)
    dx,dy=_math.cos(_math.radians(ang)),_math.sin(_math.radians(ang))
    nx,ny=-dy,dx
    out=[]
    for k in range(int(-L/gap)-1,int(L/gap)+2):
        ox,oy=cx+k*gap*nx, cy+k*gap*ny
        seg=LineString([(ox-L*dx,oy-L*dy),(ox+L*dx,oy+L*dy)]).intersection(P)
        gs=getattr(seg,'geoms',[seg]) if not seg.is_empty else []
        for g in gs:
            if g.geom_type=='LineString' and g.length>1e-6:
                a,b=g.coords[0],g.coords[-1]; out.append((a[0],a[1],b[0],b[1]))
    return out

cols={'A':"#1f77b4",'B':"#ff7f0e",'C':"#9467bd",'D':"#2ca02c"}
ORDER=['B','D','A','C']

def quads_of(th, ex=1.5):
    c=_config(th)
    E,F,G,H,Q0,Q1,v3=c['E'],c['F'],c['G'],c['H'],c['Q0'],c['Q1'],c['v3']
    dA=(-c['d2'][0],-c['d2'][1])   # 접힌 A 방향 = κ1-반사(−x̂) — folded pass-A direction = the κ1-reflection of −x̂ (crosses the body, exits through [G,Q1])
    exA=ex+0.4
    return {'A':[E,F,_add(F,_mul(dA,exA)),_add(E,_mul(dA,exA))],
            'B':[E,H,G,F],'C':[H,G,Q1,Q0],
            'D':[Q0,Q1,_add(Q1,_mul(v3,ex)),_add(Q0,_mul(v3,ex))]}, c

def render(th, out, overlay, vlabels, inset=None, W=1060, H=640):
    q,c=quads_of(th)
    pts=[p for k in q for p in q[k]]+overlay
    xs=[p[0] for p in pts]; ys=[p[1] for p in pts]
    x0,x1,y0,y1=min(xs),max(xs),min(ys),max(ys)
    pad=0.3
    DW = W-330 if inset else W-60
    sc=min(DW/(x1-x0+2*pad),(H-80)/(y1-y0+2*pad))
    def T(p): return (30+(p[0]-x0+pad)*sc, H-42-(p[1]-y0+pad)*sc)
    def ps(poly): return " ".join(f"{T(p)[0]:.1f},{T(p)[1]:.1f}" for p in poly)
    S=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">','<rect width="100%" height="100%" fill="white"/>']
    for k in reversed(ORDER):
        S.append(f'<polygon points="{ps(q[k])}" fill="{cols[k]}" fill-opacity="0.22" stroke="none"/>')
        for hx1,hy1,hx2,hy2 in hatch_segs(Polygon([T(p) for p in q[k]]),HATCH_ANG.get(k,20)):
            S.append(f'<line x1="{hx1:.1f}" y1="{hy1:.1f}" x2="{hx2:.1f}" y2="{hy2:.1f}" stroke="#333" stroke-width="0.8" opacity="0.30"/>')
    for depth,k in enumerate(ORDER):
        above=[Polygon(q[j]) for j in ORDER[:depth]]
        U=unary_union(above) if above else None
        quad=q[k]
        for i in range(4):
            e=LineString([tuple(quad[i]),tuple(quad[(i+1)%4])])
            vis=e.difference(U) if U else e
            hid=e.intersection(U) if U else None
            def emit(geom,style):
                if geom is None or geom.is_empty: return
                if geom.geom_type=="LineString": gs=[geom]
                elif hasattr(geom,"geoms"): gs=[g for g in geom.geoms if g.geom_type=="LineString"]
                else: return
                for g in gs:
                    if g.length<1e-9: continue
                    d=" ".join(f"{T(p)[0]:.1f},{T(p)[1]:.1f}" for p in g.coords)
                    S.append(f'<polyline points="{d}" fill="none" {style}/>')
            emit(vis, f'stroke="{cols[k]}" stroke-width="2.4"')
            emit(hid, f'stroke="{cols[k]}" stroke-width="1.5" stroke-dasharray="6,4" stroke-opacity="0.85"')
    for a,b in [(c['E'],c['F']),(c['H'],c['G']),(c['Q1'],c['Q0'])]:
        A_,B_=T(a),T(b)
        S.append(f'<line x1="{A_[0]:.1f}" y1="{A_[1]:.1f}" x2="{B_[0]:.1f}" y2="{B_[1]:.1f}" stroke="#c1121f" stroke-width="3.2"/>')
    S.append(f'<polygon points="{ps(overlay)}" fill="none" stroke="#0a0a0a" stroke-width="3.4"/>')
    for p in overlay:
        A_=T(p); S.append(f'<circle cx="{A_[0]:.1f}" cy="{A_[1]:.1f}" r="4.6" fill="#111" stroke="white" stroke-width="1.3"/>')
    for lbl,p,dx,dy in vlabels:
        A_=T(p); S.append(f'<text x="{A_[0]+dx:.0f}" y="{A_[1]+dy:.0f}" font-family="Georgia" font-style="italic" font-size="25">{lbl}</text>')
    E,H_,G,Q0,Q1,F_=c['E'],c['H'],c['G'],c['Q0'],c['Q1'],c['F']
    for lbl,pos,fill in (("pass A (ray-in)",(( (E[0]+F_[0])/2-c['d2'][0]*2.35),(0.5-c['d2'][1]*2.35)),cols['A']),("B",((E[0]+H_[0])/2+0.12,0.30),"#b45f06"),
                         ("C",((H_[0]+Q1[0])/2,(1+Q1[1])/2),"#6a3d9a"),
                         ("pass D (ray-out)",(Q0[0]+c['v3'][0]*1.55,Q0[1]+c['v3'][1]*1.55),"#1a7a1a")):
        A_=T(pos); S.append(f'<text x="{A_[0]:.0f}" y="{A_[1]:.0f}" font-family="Georgia" font-style="italic" font-size="26" fill="{fill}" text-anchor="middle">{lbl}</text>')
    if inset:
        Xp,Yp,Fp=inset
        cx,cy,r=W-170,150,105
        S.append(f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="#fafafa" stroke="#888" stroke-width="1.4"/>')
        gx=sum(p[0] for p in (Xp,Yp,Fp))/3; gy=sum(p[1] for p in (Xp,Yp,Fp))/3
        scale=r*1.15/max(math.hypot(p[0]-gx,p[1]-gy) for p in (Xp,Yp,Fp))
        def TI(p): return (cx+(p[0]-gx)*scale, cy-(p[1]-gy)*scale)
        tri=" ".join(f"{TI(p)[0]:.1f},{TI(p)[1]:.1f}" for p in (Xp,Yp,Fp))
        S.append(f'<polygon points="{tri}" fill="#fdd" stroke="#d62728" stroke-width="2"/>')
        for lbl,p,dx,dy in (("X",Xp,0,-12),("Y",Yp,16,5),("F",Fp,-4,22)):
            A_=TI(p)
            dot="#d62728" if lbl=="Y" else "#111"
            S.append(f'<circle cx="{A_[0]:.1f}" cy="{A_[1]:.1f}" r="5" fill="{dot}"/>')
            S.append(f'<text x="{A_[0]+dx:.0f}" y="{A_[1]+dy:.0f}" font-family="Georgia" font-style="italic" font-size="20" text-anchor="middle">{lbl}</text>')
        A_=T(Fp)
        S.append(f'<line x1="{A_[0]:.1f}" y1="{A_[1]:.1f}" x2="{cx-r*0.85:.0f}" y2="{cy+r*0.55:.0f}" stroke="#999" stroke-dasharray="3,3"/>')
        S.append(f'<circle cx="{A_[0]:.1f}" cy="{A_[1]:.1f}" r="13" fill="none" stroke="#999" stroke-width="1.2"/>')
        S.append(f'<text x="{cx}" y="{cy+r+48}" font-family="Helvetica" font-size="17" text-anchor="middle">depth 0.0293w &#183; reflex angle at Y = 249.517&#176;</text>')
    S.append('</svg>')
    open(out,"w").write("\n".join(S))
    print("saved",out)

x0v=math.sqrt((5-2*math.sqrt(5))/5)
th0=_slice_theta(x0v); c0=_config(th0)
PENT=[c0['E'],c0['H'],c0['G'],c0['Q1'],c0['F']]
VL=[("E",c0['E'],-16,20),("H",c0['H'],12,20),("G",c0['G'],14,4),("K = Q1",c0['Q1'],0,-12),("F = Q0",c0['F'],-52,-4)]
render(th0,"paper/figs/fig1-pentagon.svg",PENT,VL,W=980)

xd=0.32212749731122718
gd=(1-xd*xd)/(2*xd); hd=1.0363863004307971; bd=1.7132441870641958
thd=(xd,gd,hd,hd,bd); cd=_config(thd)
r3H=cd['rho3'](cd['H']); dD=cd['rl3'](cd['w0'])
tpar=(1.0-r3H[1])/dD[1]
Y=(r3H[0]+tpar*dD[0],1.0)
HEP=[cd['E'],cd['H'],cd['G'],cd['Q1'],cd['Q0'],Y,cd['F']]
VL2=[("E",cd['E'],-16,20),("H",cd['H'],12,20),("G",cd['G'],14,4),("K = Q1",cd['Q1'],0,-12),("X = Q0",cd['Q0'],-56,-8)]
render(thd,"paper/figs/theta-dagger-heptagon.svg",HEP,VL2,inset=(cd['Q0'],Y,cd['F']))
