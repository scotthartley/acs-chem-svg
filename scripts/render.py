"""Render a spec to SVG + PNG: an M2/M3 scheme (steps), or M1 species side by side.

Layout is acs_chem_svg.layout.render_spec, shared with scripts/gallery.py. Writes
<stem>.svg and <stem>.png (4 px/pt) next to the spec, or in --out, and prints the drawing width.
--crops also writes zoomed crops (design.md §6.5, 16 px/pt) of every curved arrow, every
labeled atom with a charge or dots, every resonance arrow and every bracket corner:
<stem>_crop_arrow-<step>-<index>.png, <stem>_crop_<species>-atom<id>.png,
<stem>_crop_resonance-<step>-<k>.png and <stem>_crop_bracket-<step>-<left|right>-<top|bottom>.png.
--maps also writes <stem>_maps.png: the drawing with every mapped atom's map number in small red
type beside it, for choosing arrows and orientations. The overlay is preview only: it goes into
a copy of the SVG text, never into <stem>.svg.
--affinity writes each arrowhead as a separate filled shape instead of an SVG <marker>: Affinity
drops markers on import (the shafts arrive with no heads). Markers follow edits in Inkscape.

Usage: uv run python scripts/render.py spec.yaml [--out DIR] [--px-per-pt 4] [--crops] [--maps] [--affinity]
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from acs_chem_svg.geometry import union_box
from acs_chem_svg.layout import SchemeLayout, build_scheme, render_spec
from acs_chem_svg.preview import crops, rasterize, rasterize_file
from acs_chem_svg.spec import load_spec
from acs_chem_svg.svg_writer import FONT_FAMILY, PADDING, faces_used, num, to_svg

MAP_COLOR = "#d01c5b"


def crop_boxes(lay: SchemeLayout, offset: tuple[float, float]) -> dict[str, tuple]:
    """name -> box in the SVG frame: each curved arrow, each labeled atom with a charge or dots
    (the union of its label runs, charge bars and dots), each resonance arrow, and each bracket
    corner (6 pt around it)."""
    dx, dy = offset

    def moved(b):
        return (b[0] + dx, b[1] + dy, b[2] + dx, b[3] + dy)
    out = {f"arrow-{a.step}-{a.index}": moved(a.group.bbox()) for a in lay.curved}
    for name, p in lay.species.items():
        dots = p.drawer.dot_sets()
        for i in p.drawer.labels:
            a = p.mol.atoms[i]
            if not (a.charge or dots.get(i)):
                continue
            src = f"atom:{a.id}"
            out[f"{name}-atom{a.id}"] = moved(union_box(leaf.bbox() for leaf in p.group.leaves()
                                                        if getattr(leaf, "source", "") == src))
    for res in lay.resonance:
        for k, (_, _, g) in enumerate(res.arrows, 1):
            out[f"resonance-{res.n}-{k}"] = moved(g.bbox())
        for g in res.group.children:
            pts = g.children[0].points
            for corner, (x, y) in (("top", pts[1]), ("bottom", pts[2])):
                out[f"bracket-{res.n}-{g.id.rsplit('-', 1)[-1]}-{corner}"] = moved((x - 6, y - 6, x + 6, y + 6))
    return out


def map_labels(lay: SchemeLayout, offset: tuple[float, float]) -> list[tuple[str, float, float]]:
    """(map number, x, y) in the SVG frame for every mapped atom drawn. An unlabeled carbon's
    number sits 0.45 bond lengths out, away from its bonds; a labeled atom's 0.6 out on the
    diagonal farthest from its bonds (its lone pairs sit on the sides, a charge top right)."""
    dx, dy = offset
    bl = lay.preset.lengths.bond_length
    diagonals = [(0.7071, 0.7071), (-0.7071, 0.7071), (-0.7071, -0.7071), (0.7071, -0.7071)]  # y down
    out = []
    for p in lay.species.values():
        ox, oy = p.offset
        labeled = set(p.drawer.labels) if p.drawer else set()
        for a in p.mol.atoms:
            if not a.map:
                continue
            bonds = []
            for j, _ in p.mol.neighbors(a.idx):
                b = p.mol.atoms[j]
                n = math.hypot(b.x - a.x, b.y - a.y) or 1.0
                bonds.append(((b.x - a.x) / n, (b.y - a.y) / n))
            vx, vy = sum(u for u, _ in bonds), sum(v for _, v in bonds)
            n = math.hypot(vx, vy)
            if a.idx in labeled or n < 1e-6:
                ux, uy = max(diagonals, key=lambda d: min((-(d[0] * u + d[1] * v) for u, v in bonds), default=1))
                d = 0.6 * bl
            else:
                ux, uy, d = -vx / n, -vy / n, 0.45 * bl
            out.append((str(a.map), a.x + ox + dx + d * ux, a.y + oy + dy + d * uy))
    return out


def maps_svg(lay: SchemeLayout) -> str:
    """A copy of the drawing's SVG with the map-number overlay (and room for it at the edges)."""
    svg, offset = to_svg(lay.root, padding=PADDING + 0.6 * lay.preset.lengths.bond_length)
    size = 0.3 * lay.preset.lengths.bond_length
    texts = [f'    <text x="{num(x)}" y="{num(y + 0.35 * size)}" font-family="{FONT_FAMILY}" font-size="{num(size)}" '
             f'font-weight="bold" text-anchor="middle" fill="{MAP_COLOR}" stroke="#fff" stroke-width="{num(size / 4)}" '
             f'paint-order="stroke">{m}</text>' for m, x, y in map_labels(lay, offset)]
    return svg.replace("</svg>", "\n".join(['  <g id="atom-maps">', *texts, "  </g>", "</svg>"]))


def render(spec_path: Path, out_dir: Path | None = None, px_per_pt: float = 4.0,
           with_crops: bool = False, with_maps: bool = False, markers: bool = True) -> list[Path]:
    spec = load_spec(spec_path)
    lay = build_scheme(spec, spec_path.parent, spec_path.stem) if spec.steps else None
    root = lay.root if lay else render_spec(spec_path)
    svg, offset = to_svg(root, markers=markers)
    out_dir = out_dir or spec_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    svg_path = out_dir / f"{spec_path.stem}.svg"
    svg_path.write_text(svg)
    png_path = rasterize_file(svg_path, out_dir / f"{spec_path.stem}.png", px_per_pt)
    for f in faces_used(root):
        if "FALLBACK" in f:
            print(f"warning: labels measured with fallback font {f}", file=sys.stderr)
    b = root.bbox()
    print(f"{spec_path.stem}: {b[2] - b[0]:.1f} x {b[3] - b[1]:.1f} pt (ink)", file=sys.stderr)
    for a in lay.curved if lay else []:
        print(f"  curved arrow step {a.step} index {a.index} ({a.spec.from_} -> {a.spec.to}): bow {a.bow}",
              file=sys.stderr)
    paths = [svg_path, png_path]
    if with_crops and lay:
        paths += crops(svg, crop_boxes(lay, offset), out_dir, prefix=f"{spec_path.stem}_crop_")
    if with_maps and lay:
        maps_png = out_dir / f"{spec_path.stem}_maps.png"
        maps_png.write_bytes(rasterize(maps_svg(lay), px_per_pt))
        paths.append(maps_png)
    return paths


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("spec", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--px-per-pt", type=float, default=4.0)
    ap.add_argument("--crops", action="store_true",
                    help="also write zoomed crops of curved arrows, charged / dotted atoms, resonance arrows and bracket corners")
    ap.add_argument("--maps", action="store_true",
                    help="also write <stem>_maps.png with each mapped atom's map number (preview only)")
    ap.add_argument("--affinity", action="store_true",
                    help="write arrowheads as separate shapes, not SVG markers (Affinity drops markers on import)")
    args = ap.parse_args(argv)
    for p in render(args.spec, args.out, args.px_per_pt, args.crops, args.maps, not args.affinity):
        print(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
