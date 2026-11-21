"""Render Figure 4 from the common, computed geometry source.

Only the certified second-family interval [x_dagger, x_pentagon] is plotted.
The red polyline samples 401 roots of N5 on this interval. The blue polyline
evaluates the exact standard-slice formula 1/x - x at those abscissas.
The right panel is a labelled y-axis zoom, not a different parametrization.
Use --check to check the SVG, or --pdf to also regenerate its PDF.
"""
import argparse
from html import escape
from pathlib import Path
import subprocess

from export_figure_geometry import build_data


OUTPUT = Path(__file__).resolve().parents[1] / "paper/figs/fig4-two-families.svg"
BLUE, RED = "#1476b8", "#c62332"


def svg_content(data):
    second = list(reversed(data["paths"]["second"]))
    standard = [(p["x"], 1 / p["x"] - p["x"]) for p in second]
    red = [(p["x"], p["ell"]) for p in second]
    xd, ld = red[0]
    xp, lp = red[-1]
    assert len(red) == 401 and xd < xp and ld < lp
    assert all(a[0] < b[0] and a[1] < b[1] for a, b in zip(red, red[1:]))
    assert abs(standard[-1][1] - lp) < 1e-13
    assert abs(ld - data["states"]["dagger"]["lengths"]["core"]) < 1e-13
    dx = xp - xd
    xmin, xmax = xd - dx * .035, xp + dx * .035
    s = ['<svg xmlns="http://www.w3.org/2000/svg" width="1160" height="560" viewBox="0 0 1160 560">',
         '<title>Standard slice and certified second family, with a core-length zoom</title>',
         '<desc>Blue is 1/x minus x. Red is computed from N5=0 on the certified branch from the pentagon to theta dagger. The right panel changes only the vertical scale.</desc>',
         '<rect width="1160" height="560" fill="white"/>',
         '<style>text{font-family:Helvetica,Arial,sans-serif;fill:#242424;font-size:20px} .tick{font-size:18px;fill:#555} .label{font-size:21px}</style>',
         '<defs><marker id="red-arrow" markerWidth="7" markerHeight="7" refX="5" refY="3.5" orient="auto"><path d="M0,0 L7,3.5 L0,7" fill="none" stroke="#c62332" stroke-width="1.3"/></marker></defs>']

    def text(x, y, value, anchor="start", color=None, cls=None, extra=""):
        color_attr = f' style="fill:{color}"' if color else ""
        class_attr = f' class="{cls}"' if cls else ""
        s.append(f'<text x="{x:.3f}" y="{y:.3f}" text-anchor="{anchor}"{color_attr}{class_attr}{extra}>{escape(value)}</text>')

    def line(a, b, color, width=1, extra=""):
        s.append(f'<line x1="{a[0]:.3f}" y1="{a[1]:.3f}" x2="{b[0]:.3f}" y2="{b[1]:.3f}" stroke="{color}" stroke-width="{width}" {extra}/>')

    panels = [(106, "Full comparison", 2.7485, 2.785, [2.75, 2.76, 2.77, 2.78]),
              (696, "Same x-range, vertical zoom", 2.7492, 2.7535, [2.75, 2.751, 2.752, 2.753])]
    for i, (left, title, ymin, ymax, yticks) in enumerate(panels):
        top, width, height = 110, 400, 340
        def point(x, y):
            return left + (x - xmin) / (xmax - xmin) * width, top + (ymax - y) / (ymax - ymin) * height
        text(left + width / 2, 84, title, anchor="middle", cls="label")
        s.append(f'<clipPath id="plot-{i}"><rect x="{left}" y="{top}" width="{width}" height="{height}"/></clipPath>')
        for y in yticks:
            sy = point(xmin, y)[1]
            line((left, sy), (left + width, sy), "#dedede")
            text(left - 10, sy + 6, f"{y:.3f}" if i else f"{y:.2f}", anchor="end", cls="tick")
        for x in (.3225, .3235, .3245):
            sx = point(x, ymin)[0]
            line((sx, top), (sx, top + height), "#eeeeee")
            text(sx, top + height + 28, f"{x:.4f}", anchor="middle", cls="tick")
        s.append(f'<rect x="{left}" y="{top}" width="{width}" height="{height}" fill="none" stroke="#777" stroke-width="1.2"/>')
        text(left + width / 2, 508, "x", anchor="middle")
        text(left - 72, top + height / 2, "core length ℓ", anchor="middle", extra=f' transform="rotate(-90 {left-72} {top+height/2})"')
        s.append(f'<g clip-path="url(#plot-{i})">')
        for name, values, color in (("standard", standard, BLUE), ("second", red, RED)):
            points = " ".join(f"{sx:.4f},{sy:.4f}" for sx, sy in (point(x, y) for x, y in values))
            s.append(f'<polyline id="{name}-{i}" points="{points}" fill="none" stroke="{color}" stroke-width="3" stroke-linejoin="round"/>')
        # Arrow points in the certified pulling direction: P* -> theta dagger.
        line(point(*red[255]), point(*red[240]), RED, 2.2, 'marker-end="url(#red-arrow)"')
        s.append('</g>')
        for x, y, color in ((xd, ld, RED), (xp, lp, "#111111")):
            sx, sy = point(x, y)
            s.append(f'<circle cx="{sx:.3f}" cy="{sy:.3f}" r="5" fill="{color}"/>')
        px, py = point(xp, lp)
        tx, ty = point(xd, ld)
        text(px - 12, py - 15 if i == 0 else py + 26, "P*", anchor="end")
        text(tx + 12, ty - 12, "θ†", color=RED)
        if i == 1:
            text(left + 12, top + 29, f"Δℓ = {lp-ld:.7f}", cls="label")
            text(left + 12, top + 56, f"({100*(lp-ld)/lp:.3f}% of ℓ(P*))", cls="tick")

    line((140, 29), (184, 29), BLUE, 3)
    text(196, 36, "standard slice: ℓ = 1/x − x", color=BLUE)
    line((634, 29), (678, 29), RED, 3)
    text(690, 36, "second family: N₅(x,h) = 0", color=RED)
    text(580, 548, "Computed branch samples; arrow shows decreasing core length from P* to θ†.", anchor="middle", cls="tick")
    s.append('</svg>')
    return "\n".join(s) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--pdf", action="store_true")
    args = parser.parse_args()
    content = svg_content(build_data())
    if args.check:
        assert OUTPUT.read_text() == content, "Figure 4 is stale: run render_families.py --pdf"
        print("PASS: Figure 4 matches the computed branch geometry")
        return
    OUTPUT.write_text(content)
    print("Wrote", OUTPUT)
    if args.pdf:
        subprocess.run(["rsvg-convert", "-f", "pdf", str(OUTPUT), "-o", str(OUTPUT.with_suffix(".pdf"))], check=True)
        print("Wrote", OUTPUT.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
