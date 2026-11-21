# fig0-notation: 좌(재료 좌표: 스트립+크리스 현+패스) / 우(접힌 상태: 느슨한 슬라이스 x=0.25) — left (material coordinates: strip + crease chords + passes) / right (folded state: loose slice x=0.25)
import argparse
from pathlib import Path
from svg_labels import math_text, save_svg
from ribbon_geometry import geometry, standard_theta
# Same width-one slice coordinates, without importing the census computations.
state=geometry(standard_theta(0.25))
c=dict(state['vertices'],d2=state['d2'])
E,F,G,H,Q0,Q1=c['E'],c['F'],c['G'],c['H'],c['Q0'],c['Q1']
M1,M2,M3=state['centerline']['core']
v3=state['dD']
W,Hh=1180,470
S=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{Hh}" viewBox="0 0 {W} {Hh}">','<rect width="100%" height="100%" fill="white"/>']
# ---- 좌: 재료 좌표 ---- — left: material coordinates
x0,y0,sc=40,120,44   # 스트립 [0,10]x[0,1] — strip [0,10]x[0,1]
def TL(p): return (x0+p[0]*sc, y0+ (1-p[1])*sc*2.2)
S.append(f'<rect x="{x0}" y="{y0}" width="{10*sc}" height="{2.2*sc}" fill="#f2f6fa" stroke="#333" stroke-width="1.5"/>')
S.append(f'<line x1="{x0}" y1="{y0+1.1*sc}" x2="{x0+10*sc}" y2="{y0+1.1*sc}" stroke="#666" stroke-width="1" stroke-dasharray="6,5"/>')
cre=[(2.6,+1),(5.2,-1),(7.8,+1)]
for i,(cx,eps) in enumerate(cre):
    x1,y1=TL((cx-0.45*eps,0)); x2,y2=TL((cx+0.45*eps,1))
    S.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#c1121f" stroke-width="4"/>')
    S.append(math_text((x1+x2)/2,y0-14,f'κ_{i+1} (ε_{i+1}={"+" if eps>0 else "−"})',21,'middle'))
for lbl,px in (("A",1.3),("B",3.9),("C",6.5),("D",8.9)):
    X,Y=TL((px,0.5))
    S.append(f'<text x="{X:.0f}" y="{Y+7:.0f}" font-family="Georgia" font-style="italic" font-size="24" text-anchor="middle">{lbl}</text>')
S.append(f'<line x1="{x0-14}" y1="{y0}" x2="{x0-14}" y2="{y0+2.2*sc}" stroke="#333" stroke-width="1"/>')
S.append(f'<text x="{x0-26}" y="{y0+1.1*sc+7}" font-family="Georgia" font-size="20" text-anchor="middle">1</text>')
S.append(f'<text x="{x0+5*sc}" y="{y0+2.2*sc+34}" font-family="Georgia" font-size="20" text-anchor="middle">material strip R&#8202;w with complete creases</text>')
# 화살표 — arrow
S.append('<defs><marker id="ar" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto"><polygon points="0 0,10 4,0 8" fill="#333"/></marker></defs>')
S.append(f'<line x1="{x0+10*sc+18}" y1="{y0+1.1*sc}" x2="{x0+10*sc+64}" y2="{y0+1.1*sc}" stroke="#333" stroke-width="2.4" marker-end="url(#ar)"/>')
S.append(f'<text x="{x0+10*sc+41}" y="{y0+1.1*sc-14}" font-family="Georgia" font-size="18" text-anchor="middle">fold flat</text>')
# ---- 우: 접힌 상태 ---- — right: the folded state
pts=[E,F,G,H,Q0,Q1,M1,M2,M3]
ex=2.0
dA=(-c['d2'][0],-c['d2'][1])   # 접힌 A 방향 = κ1-반사 — folded pass-A direction (the κ1-reflection)
Atail=[(E[0]+ex*dA[0], E[1]+ex*dA[1]),(F[0]+ex*dA[0], F[1]+ex*dA[1])]
Dt=[(Q0[0]+ex*v3[0], Q0[1]+ex*v3[1]),(Q1[0]+ex*v3[0], Q1[1]+ex*v3[1])]
allp=pts+Atail+Dt
xs=[p[0] for p in allp]; ys=[p[1] for p in allp]
bx0,bx1,by0,by1=min(xs),max(xs),min(ys),max(ys)
rx0,ry0=560,34; rw,rh=596,400
sc2=min(rw/(bx1-bx0+0.6), rh/(by1-by0+0.6))
def TR(p): return (rx0+(p[0]-bx0+0.3)*sc2, ry0+rh-(p[1]-by0+0.3)*sc2)
def quad(ps,fill,op):
    s=" ".join(f"{TR(p)[0]:.1f},{TR(p)[1]:.1f}" for p in ps)
    S.append(f'<polygon points="{s}" fill="none" stroke="white" stroke-width="8"/>')
    S.append(f'<polygon points="{s}" fill="{fill}" fill-opacity="0.16" stroke="{fill}" stroke-width="2.2"/>')
# 화가 순서 = 물리 층 순서 (B>D>A>C ⟹ 아래부터 C,A,D,B) — painter's order = physical layer order (B>D>A>C ⟹ from the bottom: C,A,D,B)
quad([H,G,Q1,Q0],"#9467bd",0.42)                                     # C (최하) — C (bottom-most)
quad([Atail[0],E,F,Atail[1]],"#8ecae6",0.5)                          # A
quad([Q0,Q1,Dt[1],Dt[0]],"#2ca02c",0.42)                             # D
quad([E,H,G,F],"#ffb703",0.5)                                        # B (최상) — B (top-most)
for a,b,lbl,dx,dy,anchor in ((E,F,'κ_1',-18,-14,'end'),(H,G,'κ_2',14,-14,'start'),(Q1,Q0,'κ_3',-22,-15,'end')):
    A_,B_=TR(a),TR(b)
    S.append(f'<line x1="{A_[0]:.1f}" y1="{A_[1]:.1f}" x2="{B_[0]:.1f}" y2="{B_[1]:.1f}" stroke="#c1121f" stroke-width="4"/>')
    S.append(math_text((A_[0]+B_[0])/2+dx,(A_[1]+B_[1])/2+dy,lbl,20,anchor,halo=2.5))
for m,lbl,dx,anchor in ((M1,'M_1',-15,'end'),(M2,'M_2',10,'start'),(M3,'M_3',-20,'end')):
    A_=TR(m)
    S.append(f'<circle cx="{A_[0]:.1f}" cy="{A_[1]:.1f}" r="5" fill="#d00"/>')
    S.append(math_text(A_[0]+dx,A_[1]+18,lbl,19,anchor,halo=2.5))
for lbl,pos in (("A",((E[0]+F[0])/2+dA[0]*(ex+0.32),0.5+dA[1]*(ex+0.32))),("B",(E[0]+0.40,0.11)),("C",((H[0]+Q1[0])/2+0.2,(1+Q1[1])/2)),("D",(Dt[0][0]*0.72+Dt[1][0]*0.28,(Dt[0][1]*0.72+Dt[1][1]*0.28)))):
    A_=TR(pos)
    S.append(f'<text x="{A_[0]:.0f}" y="{A_[1]:.0f}" font-family="Georgia" font-style="italic" font-size="24" text-anchor="middle">{lbl}</text>')
S.append(f'<text x="{rx0+rw/2}" y="{ry0+rh+30}" font-family="Georgia" font-size="20" text-anchor="middle">a loose flat folded state (layer order B &gt; D &gt; A &gt; C)</text>')
S.append('</svg>')
if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Render Figure 1 with mathematical labels outside the crease edges.')
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true')
    group.add_argument('--pdf',action='store_true')
    args=parser.parse_args()
    save_svg('\n'.join(S),Path(__file__).resolve().parents[1]/'paper/figs/fig0-notation.svg',check=args.check,pdf=args.pdf)
