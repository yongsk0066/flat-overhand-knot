"""Reproduce Figure 3's schematic, with labels clear of both family curves.

Curve control points are schematic illustration coordinates, not measurements.
Figure 4 carries the computed family graph. This generator preserves Figure 3's
curve and contour geometry and only separates its annotations from that geometry.
"""
import argparse
from pathlib import Path

from svg_labels import math_text, save_svg


OUTPUT = Path(__file__).resolve().parents[1] / 'paper/figs/fig3-mountain-pass.svg'


def svg_content():
    s = ['<svg xmlns="http://www.w3.org/2000/svg" width="880" height="440" viewBox="0 0 880 440">',
         '<title>Schematic pentagonal contact-mode junction</title>',
         '<rect width="880" height="440" fill="white"/>']
    for rx, ry, opacity in ((215,133,.12),(170,105,.18),(125,78,.26),(82,51,.36)):
        s.append(f'<ellipse cx="310" cy="230" rx="{rx}" ry="{ry}" fill="none" stroke="#d62728" stroke-width="1.4" opacity="{opacity}"/>')
    s += ['<path id="standard-slice" d="M 800,70 C 700,100 620,155 560,230" fill="none" stroke="#1f77b4" stroke-width="3.2"/>',
          '<path id="second-family" d="M 560,230 C 490,290 400,285 310,230 C 250,196 235,140 232,95" fill="none" stroke="#d62728" stroke-width="3.2"/>',
          '<defs><marker id="ar" markerWidth="10" markerHeight="8" refX="8" refY="4" orient="auto"><path d="M0,0 L10,4 L0,8 z" fill="#444"/></marker></defs>',
          '<path d="M 745,95 C 675,120 625,160 592,198" fill="none" stroke="#444" stroke-width="1.6" marker-end="url(#ar)"/>',
          '<path d="M 536,252 C 480,288 420,288 366,262" fill="none" stroke="#444" stroke-width="1.6" marker-end="url(#ar)"/>',
          '<circle cx="560" cy="230" r="6" fill="#000"/>',
          '<path d="M 568,238 L 596,273" fill="none" stroke="#666" stroke-width="1.2"/>',
          math_text(607,281,'P^{*}',22),
          '<text x="607" y="307" font-family="Helvetica,Arial,sans-serif" font-size="16">contact-mode junction (8 walls)</text>',
          '<circle cx="310" cy="230" r="6" fill="#d62728"/>',
          '<rect x="245" y="289" width="185" height="31" rx="4" fill="white" fill-opacity="0.94"/>',
          math_text(252,311,'θ^{†}',21,color='#d62728'),
          '<text x="284" y="311" font-family="Helvetica,Arial,sans-serif" font-size="16" fill="#d62728">strict minimum</text>',
          '<text x="640" y="52" font-family="Helvetica,Arial,sans-serif" font-size="15" fill="#1f77b4" text-anchor="end">standard slice</text>',
          '<text x="225" y="72" font-family="Helvetica,Arial,sans-serif" font-size="15" fill="#d62728" text-anchor="middle">second family</text>',
          '<text x="440" y="402" font-family="Helvetica,Arial,sans-serif" font-size="13" fill="#777" text-anchor="middle">Schematic contact-mode junction (not a smooth saddle); contours indicate core length on the wall stratum.</text>',
          '</svg>']
    return '\n'.join(s)+'\n'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true')
    group.add_argument('--pdf',action='store_true')
    args=parser.parse_args()
    save_svg(svg_content(),OUTPUT,check=args.check,pdf=args.pdf)


if __name__=='__main__':
    main()
