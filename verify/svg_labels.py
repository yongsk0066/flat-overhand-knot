"""Small SVG helpers for the paper's mathematical labels and artifact checks.

Scripts are native SVG tspans, not baseline digits or font-specific Unicode
substitutes. This keeps subscripts and superscripts legible in librsvg PDFs.
"""
from html import escape
from pathlib import Path
import re
import subprocess


SCRIPT = re.compile(r"([_^])(?:\{([^}]+)\}|([A-Za-z0-9+*†−-]+))")


def math_text(x, y, expression, size=25, anchor="start", color="#111", halo=0):
    """A small TeX-like subset: Q_0, E_4^{+}, P^{*}, and plain math text."""
    spans, end, previous_sub = [], 0, None
    for match in SCRIPT.finditer(expression):
        plain = expression[end:match.start()]
        if plain:
            spans.append(escape(plain))
            previous_sub = None
        kind, braced, simple = match.groups()
        script = braced if braced is not None else simple
        shift = -.28 * size if kind == "_" else .40 * size
        dx = -.52 * .68 * size * len(previous_sub) if kind == "^" and previous_sub else 0
        spans.append(f'<tspan baseline-shift="{shift:.2f}px" dx="{dx:.2f}" font-size="{.68*size:.2f}" font-style="normal">{escape(script)}</tspan>')
        previous_sub = script if kind == "_" else None
        end = match.end()
    if end < len(expression):
        spans.append(escape(expression[end:]))
    body = "".join(spans)
    attrs = f'x="{x:.2f}" y="{y:.2f}" font-family="Georgia,Times New Roman,serif" font-style="italic" font-size="{size}" text-anchor="{anchor}" xml:space="preserve"'
    under = f'<text {attrs} fill="white" stroke="white" stroke-width="{halo}" stroke-linejoin="round">{body}</text>' if halo else ""
    return under + f'<text {attrs} fill="{color}">{body}</text>'


def save_svg(content, path, check=False, pdf=False, pdf_path=None):
    path = Path(path)
    if check:
        assert not pdf, "--check and --pdf cannot be combined"
        assert path.read_text() == content, f"Stale figure: {path}"
        print("verified", path)
        return
    path.write_text(content)
    print("saved", path)
    if pdf:
        output = Path(pdf_path) if pdf_path else path.with_suffix(".pdf")
        subprocess.run(["rsvg-convert", "-f", "pdf", str(path), "-o", str(output)], check=True)
        print("saved", output)
