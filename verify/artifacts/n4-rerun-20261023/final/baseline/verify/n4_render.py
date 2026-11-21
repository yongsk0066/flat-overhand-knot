# n4_render.py — B_{4_1} 폭-1 접힘 상태 SVG 렌더 (밴드 사각형 + 크리스 + 중심선) — SVG render of a width-1 folded state of B_{4_1} (band quads + creases + centerline)
import sys, json, math
sys.path.insert(0,'verify')
from n4_tighten2 import config, pieces, crossings, valid_fixed, ell

def band_quads(c):
    CH=c['CH']; d=c['d']; M=c['M']
    def side(pt,k):
        P,dd,_=pieces(c)[k]
        return 1 if (-(pt[0]-P[0])*dd[1]+(pt[1]-P[1])*dd[0])>=0 else -1
    quads={}
    for k in (1,2,3):
        a0,b0=CH[k-1]; a1,b1=CH[k]
        # 측 부호로 짝 맞춤 (piece k 기준) — pair up by lateral sign (relative to piece k)
        p0=(a0,b0) if side(a0,k)>0 else (b0,a0)
        p1=(a1,b1) if side(a1,k)>0 else (b1,a1)
        quads[k]=[p0[0],p1[0],p1[1],p0[1]]
    EXT=2.2
    for k,ch,sgn in ((0,CH[0],-1),(4,CH[3],1)):
        P,dd,_=pieces(c)[k]
        a,b=ch
        if k==0: ex=(c['d'][0][0]*-EXT, c['d'][0][1]*-EXT)
        else:    ex=(c['d'][4][0]*EXT,  c['d'][4][1]*EXT)
        quads[k]=[a,b,(b[0]+ex[0],b[1]+ex[1]),(a[0]+ex[0],a[1]+ex[1])]
    return quads

def render(state_file, out_svg, title=""):
    st=json.load(open(state_file))
    th=st["theta"]; order=st.get("order",[1,4,2,0,3])
    c=config(th); q=band_quads(c)
    pts=[p for quad in q.values() for p in quad]+[p for ch in c['CH'] for p in ch]
    xs=[p[0] for p in pts]; ys=[p[1] for p in pts]
    x0,x1,y0,y1=min(xs),max(xs),min(ys),max(ys)
    W=980; pad=0.55
    sc=(W-80)/max(x1-x0+2*pad, (y1-y0+2*pad))
    H=int((y1-y0+2*pad)*sc)+80
    def T(p): return (40+(p[0]-x0+pad)*sc, H-40-(p[1]-y0+pad)*sc)
    cols={0:"#8ecae6",1:"#219ebc",2:"#ffb703",3:"#fb8500",4:"#c1121f"}
    S=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
       f'<rect width="{W}" height="{H}" fill="white"/>']
    # 층 순서 아래->위 — layer order, bottom -> top
    for k in reversed(order):
        quad=q[k]; ptsS=" ".join(f"{T(p)[0]:.1f},{T(p)[1]:.1f}" for p in quad)
        S.append(f'<polygon points="{ptsS}" fill="none" stroke="white" stroke-width="9"/>')
        S.append(f'<polygon points="{ptsS}" fill="{cols[k]}" fill-opacity="0.16" stroke="{cols[k]}" stroke-width="2.4"/>')
    for i,(a,b) in enumerate(c['CH']):
        A,B=T(a),T(b)
        S.append(f'<line x1="{A[0]:.1f}" y1="{A[1]:.1f}" x2="{B[0]:.1f}" y2="{B[1]:.1f}" stroke="#111" stroke-width="3"/>')
    segs=pieces(c)
    for k in range(5):
        P,dd,l=segs[k]; l=min(l,2.2)
        A=T(P); B=T((P[0]+dd[0]*l,P[1]+dd[1]*l))
        S.append(f'<line x1="{A[0]:.1f}" y1="{A[1]:.1f}" x2="{B[0]:.1f}" y2="{B[1]:.1f}" stroke="#000" stroke-width="1" stroke-dasharray="5,4" opacity="0.6"/>')
    for i,M in enumerate(c['M']):
        A=T(M); S.append(f'<circle cx="{A[0]:.1f}" cy="{A[1]:.1f}" r="4" fill="#d00"/>')
        S.append(f'<text x="{A[0]+7:.1f}" y="{A[1]-7:.1f}" font-family="Georgia" font-size="15" font-style="italic">M{i+1}</text>')
    if title:
        lg=" &gt; ".join(f"P{p}" for p in order)
        S.append(f'<text x="40" y="28" font-family="Georgia" font-size="17">{title}  (layers {lg},  ell = {ell(th):.4f})</text>')
    S.append('</svg>')
    open(out_svg,"w").write("\n".join(S))
    print("saved", out_svg)

if __name__=="__main__":
    render(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv)>3 else "")
