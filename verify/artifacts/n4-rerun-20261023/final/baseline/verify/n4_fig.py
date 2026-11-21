# n4_fig.py — 최종 그림: B_{4_1} 타이트 상태 (정확 θ*) + 육각형 강조 — final figure: the B_{4_1} tight state (exact θ*) + hexagon highlighted
import sys, math
sys.path.insert(0,'verify')
from n4_tighten2 import config, pieces
from n4_render import band_quads

t=math.acos(-0.25); s=2*math.pi-2*t
L1=2*math.sqrt(15)/5; L2=8*math.sqrt(15)/15
TH=[t,s,-s,-t,L1,L2,L1]
order=[1,4,2,0,3]
c=config(TH); q=band_quads(c)
r15=math.sqrt(15)
HEX=[(-13*r15/30,0.5),(-7*r15/30,1.5),(r15/30,1.5),(3*r15/30,0.5),(-3*r15/30,-0.5),(-11*r15/30,-0.5)]
DOUBLE=[(-7*r15/30,1.5),(-3*r15/30,-0.5)]   # E4+=E2+, E1-=E3-
W=1040; H=760; pad=0.75
xs=[p[0] for k in q for p in q[k]]; ys=[p[1] for k in q for p in q[k]]
x0,x1=min(xs),max(xs); y0,y1=min(ys),max(ys)
sc=min((W-90)/(x1-x0+2*pad),(H-110)/(y1-y0+2*pad))
def T(p): return (45+(p[0]-x0+pad)*sc, H-55-(p[1]-y0+pad)*sc)
cols={0:"#8ecae6",1:"#219ebc",2:"#ffb703",3:"#fb8500",4:"#c1121f"}
S=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
   f'<rect width="{W}" height="{H}" fill="white"/>']
for k in reversed(order):
    quad=q[k]; pts=" ".join(f"{T(p)[0]:.1f},{T(p)[1]:.1f}" for p in quad)
    S.append(f'<polygon points="{pts}" fill="none" stroke="white" stroke-width="9"/>')
    S.append(f'<polygon points="{pts}" fill="{cols[k]}" fill-opacity="0.16" stroke="{cols[k]}" stroke-width="2.4"/>')
for a,b in c['CH']:
    A,B=T(a),T(b)
    S.append(f'<line x1="{A[0]:.1f}" y1="{A[1]:.1f}" x2="{B[0]:.1f}" y2="{B[1]:.1f}" stroke="#111" stroke-width="2.4"/>')
hexpts=" ".join(f"{T(p)[0]:.1f},{T(p)[1]:.1f}" for p in HEX)
S.append(f'<polygon points="{hexpts}" fill="none" stroke="#0a0a0a" stroke-width="3.6"/>')
for p in HEX:
    A=T(p); dbl = any(abs(p[0]-dp[0])<1e-9 and abs(p[1]-dp[1])<1e-9 for dp in DOUBLE)
    S.append(f'<circle cx="{A[0]:.1f}" cy="{A[1]:.1f}" r="{7 if dbl else 5.5}" fill="{"#d00000" if dbl else "#111"}" stroke="white" stroke-width="1.6"/>')
# 꼬리 화살표 — tail arrows
A=T((x0-0.15,0.0)); B=T((x0+0.75,0.0))
S.append(f'<line x1="{B[0]:.0f}" y1="{B[1]:.0f}" x2="{A[0]+26:.0f}" y2="{A[1]:.0f}" stroke="#155" stroke-width="2" marker-end="url(#ar)"/>')
A2=T((x1+0.15,1.0)); B2=T((x1-0.75,1.0))
S.append(f'<line x1="{B2[0]:.0f}" y1="{B2[1]:.0f}" x2="{A2[0]-26:.0f}" y2="{A2[1]:.0f}" stroke="#155" stroke-width="2" marker-end="url(#ar)"/>')
S.insert(1,'<defs><marker id="ar" markerWidth="9" markerHeight="7" refX="8" refY="3.5" orient="auto"><polygon points="0 0, 9 3.5, 0 7" fill="#155"/></marker></defs>')
lg=" &gt; ".join(f"P{p}" for p in order)
import os
if not os.environ.get("PAPER"):
    S.append(f'<text x="45" y="30" font-family="Georgia" font-size="18" font-style="italic">Tight figure-eight state: hexagon, ell* = 4√15/3 = 5.16398  (layers {lg})</text>')
for lbl,p,dx,dy in [("E4-",HEX[0],-46,4),("E4+=E2+",HEX[1],-32,-12),("E2-",HEX[2],12,-8),("E1+",HEX[3],14,20),("E1-=E3-",HEX[4],10,22),("E3+",HEX[5],-44,22)]:
    A=T(p); S.append(f'<text x="{A[0]+dx:.0f}" y="{A[1]+dy:.0f}" font-family="Georgia" font-size="15" font-style="italic">{lbl}</text>')
S.append('</svg>')
import os
out = "paper/figs/fig41-tight.svg" if os.environ.get("PAPER") else "paper/figs/fig6-fig8-tight.svg"
open(out,"w").write("\n".join(S))
print("->",out)
print("saved paper/figs/fig6-fig8-tight.svg")
