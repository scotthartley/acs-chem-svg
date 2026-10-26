"""Rasterizer acceptance test (design.md §4.3, §7 [VERIFY]).

Renders kerned <text> runs with resvg (resvg-py), with the font database restricted to
one exact font file, and compares the ink box of each run in the PNG with the
fontTools prediction. Passes when every edge lands within TOL_PT, for system
Helvetica Regular/Bold/Oblique and for the TeX Gyre Heros fallback. Oblique is also
checked with the full default font list, as renders use it: font-style="italic" must pick
Helvetica Oblique, not the upright face or Heros Italic. A negative
control (no fonts loaded) must render no ink, which proves that resvg is not
silently substituting a font.

Usage: uv run python scripts/check_rasterizer.py
"""

from __future__ import annotations

import io
import sys

import numpy as np
import resvg_py
from PIL import Image

from acs_chem_svg.preview import default_font_files, rasterize

from acs_chem_svg.fontmetrics import SYSTEM_HELVETICA, FontFace, load_fallback, load_system_helvetica

SIZE = 10.0  # pt
ZOOM = 20  # px per pt -> 0.05 pt per px
TOL_PT = 0.1
# CH3 exercises Heros GPOS kerning (C-H); 11 exercises Helvetica legacy kern; n for n-BuLi
RUNS = ["H", "O", "NH", "OH", "Cl", "Br", "CH3", "H2N", "(1)", "12", "11", "n", "t-Bu"]
DEFAULT_FAMILY = "Helvetica, 'TeX Gyre Heros', 'Nimbus Sans', Arial, sans-serif"  # svg_writer.FONT_FAMILY


def _svg(family: str, weight: str, runs: list[str], style: str = "normal") -> tuple[str, list[tuple[float, float]]]:
    origins = [(5.0, 15.0 + 15.0 * i) for i in range(len(runs))]
    h = 15.0 * len(runs) + 10.0
    body = "".join(
        f'<text x="{x}" y="{y}" font-family="{family}" font-weight="{weight}" font-style="{style}" '
        f'font-size="{SIZE}">{t}</text>'
        for t, (x, y) in zip(runs, origins)
    )
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="60pt" height="{h}pt" '
           f'viewBox="0 0 60 {h}">{body}</svg>')
    return svg, origins


def _render(svg: str, font_files: list[str]) -> np.ndarray:
    img = Image.open(io.BytesIO(rasterize(svg, ZOOM, font_files=font_files))).convert("RGBA")
    return np.asarray(img)[:, :, 3].astype(float) / 255.0


def _ink_box(alpha: np.ndarray, y0: float, y1: float) -> tuple[float, float, float, float] | None:
    """Edges in pt at the 50% coverage level, restricted to rows y0..y1 (pt)."""
    a = alpha[int(y0 * ZOOM):int(y1 * ZOOM)]
    cols, rows = a.sum(axis=0), a.sum(axis=1)
    if cols.max() == 0:
        return None
    xs, ys = np.nonzero(cols)[0], np.nonzero(rows)[0]
    # sub-pixel edges: partial coverage of the first/last column or row
    xmin = xs[0] + 1 - min(1.0, a[:, xs[0]].max())
    xmax = xs[-1] + min(1.0, a[:, xs[-1]].max())
    ymin = ys[0] + 1 - min(1.0, a[ys[0]].max())
    ymax = ys[-1] + min(1.0, a[ys[-1]].max())
    return xmin / ZOOM, ymin / ZOOM + y0, xmax / ZOOM, ymax / ZOOM + y0


def _predict(face: FontFace, text: str, origin: tuple[float, float]):
    xmin, ymin, xmax, ymax = face.text_ink_bounds(text, SIZE)
    ox, oy = origin
    return ox + xmin, oy - ymax, ox + xmax, oy - ymin  # y down


def check(face: FontFace, family: str, weight: str, style: str = "normal",
          font_files: list[str] | None = None) -> bool:
    svg, origins = _svg(family, weight, RUNS, style)
    alpha = _render(svg, font_files or [str(face.path)])
    ok = True
    worst = 0.0
    for text, origin in zip(RUNS, origins):
        got = _ink_box(alpha, origin[1] - 12, origin[1] + 3)
        want = _predict(face, text, origin)
        if got is None:
            print(f"  {text!r}: NO INK rendered")
            ok = False
            continue
        d = max(abs(g - w) for g, w in zip(got, want))
        worst = max(worst, d)
        flag = "ok" if d <= TOL_PT else "FAIL"
        ok &= d <= TOL_PT
        print(f"  {text!r:7} predicted {tuple(round(v, 3) for v in want)} "
              f"rendered {tuple(round(v, 3) for v in got)}  max|Δ| {d:.3f} pt {flag}")
    print(f"  => {'PASS' if ok else 'FAIL'} (worst {worst:.3f} pt, tol {TOL_PT} pt)")
    return ok


def negative_control() -> bool:
    svg, _ = _svg("Helvetica", "normal", ["H"])
    try:
        alpha = _render(svg, [])
        ok = alpha.max() == 0
    except ValueError:  # resvg-py raises on an SVG with nothing to draw
        ok = True
    print(f"negative control (no fonts loaded): {'PASS, no ink' if ok else 'FAIL, ink rendered'}")
    return ok


def main() -> int:
    print(f"resvg {resvg_py.__resvg_version__} via resvg-py {resvg_py.__version__}; "
          f"{SIZE:g} pt text at {ZOOM} px/pt")
    results = [negative_control()]
    for sub, weight, style in (("Regular", "normal", "normal"), ("Bold", "bold", "normal"),
                               ("Oblique", "normal", "italic")):
        h = load_system_helvetica(sub)
        print(f"\nHelvetica {sub}: font db = [{SYSTEM_HELVETICA}]")
        results.append(check(h, "Helvetica", weight, style))
        fb = load_fallback(sub)
        print(f"\nFallback {sub}: font db = [{fb.path.name}]")
        results.append(check(fb, "TeX Gyre Heros", weight, style))
    print("\nHelvetica Oblique via the default font list and family:")
    results.append(check(load_system_helvetica("Oblique"), DEFAULT_FAMILY, "normal", "italic", default_font_files()))
    print("\nHelvetica Bold via the default font list and family (svg_writer's default output):")
    results.append(check(load_system_helvetica("Bold"), DEFAULT_FAMILY, "bold", "normal", default_font_files()))
    # svg_writer --affinity names the face; rasterize() must map it back (not drop the text)
    print("\nHelvetica Bold as written for Affinity (font-family=\"Helvetica Bold\"):")
    results.append(check(load_system_helvetica("Bold"), "Helvetica Bold", "bold", "normal", default_font_files()))
    print("\nHelvetica Oblique as written for Affinity (font-family=\"Helvetica Oblique\"):")
    results.append(check(load_system_helvetica("Oblique"), "Helvetica Oblique", "normal", "italic", default_font_files()))
    passed = all(results)
    print(f"\nOVERALL: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
