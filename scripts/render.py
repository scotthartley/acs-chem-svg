"""Render a spec to SVG + PNG + CDXML: a scheme (steps), a catalytic cycle, or species side by side.

Layout is acs_chem_svg.figure.build_figure, shared with scripts/gallery.py. Writes
<stem>.svg and <stem>.png (4 px/pt) next to the spec, or in --out, and prints the drawing width.
--crops also writes zoomed crops (design.md §6.5, 16 px/pt) of every curved arrow, every
labeled atom with a charge or dots, every resonance arrow, both ends of every equilibrium arrow
and every bracket corner: <stem>_crop_arrow-<step>-<index>.png, <stem>_crop_<species>-atom<id>.png,
<stem>_crop_resonance-<step>-<k>.png, <stem>_crop_equilibrium-<step>-<tail|head>.png and
<stem>_crop_bracket-<step>-<left|right>-<top|bottom>.png. A catalytic cycle gets each arc arrow's
head (and an equilibrium arc's tail end), <stem>_crop_cycle-arc-<n>-<head|tail>.png, and each
tail, <stem>_crop_cycle-<in|out>-<n>-<i>.png; its radius and free arc, and any tail still
crowded, are printed.
--maps also writes <stem>_maps.png: the drawing with every mapped atom's map number in small red
type beside it, for choosing arrows and orientations. The overlay is preview only: it goes into
a copy of the SVG text, never into <stem>.svg.
--affinity writes for Affinity: each arrowhead as a separate filled shape instead of an SVG
<marker> (Affinity drops markers on import), and bold/italic text as a bare face name
("Helvetica Bold"), which Affinity needs and Inkscape and Chrome can't resolve. The default
output is for Inkscape and browsers; its markers follow edits.
<stem>.cdxml for ChemDraw is written too (acs_chem_svg.cdxml_writer, acs-style.md §13), in the
SVG's coordinates, from the same layout as the SVG: single species (Lewis structures included), or
a scheme (arrows, conditions, plus signs, numbers and <scheme><step>, §13.1; curved arrows and
dots, §13.2; abbreviations, resonance arrows and brackets, §13.3), or a catalytic cycle (arc arrows
and tails as ChemDraw arc arrows, §15.1). What it can't write is a "warning: no CDXML:" line on
stderr (the SVG and PNG are still written); what it writes as plain chemistry (atom maps) is
listed on stderr. --no-cdxml skips it (--cdxml is accepted and does nothing).

Usage: uv run python scripts/render.py spec.yaml [--out DIR] [--px-per-pt 4] [--crops] [--maps] [--affinity] [--no-cdxml]
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from acs_chem_svg.cdxml_writer import CdxmlUnsupported
from acs_chem_svg.cycle import CycleLayout
from acs_chem_svg.figure import build_figure
from acs_chem_svg.geometry import arc_point, pad_box, shift_box, union_box
from acs_chem_svg.layout import SchemeLayout
from acs_chem_svg.preview import crops, rasterize, rasterize_file
from acs_chem_svg.svg_writer import FONT_FAMILY, PADDING, faces_used, num, to_svg

MAP_COLOR = "#d01c5b"


def crop_boxes(lay: SchemeLayout, offset: tuple[float, float]) -> dict[str, tuple]:
    """name -> box in the SVG frame: each curved arrow, each labeled atom with a charge or dots
    (the union of its label runs, charge bars and dots), each resonance arrow, each end of an
    equilibrium arrow (12 pt long) and each bracket corner (6 pt around it)."""
    def moved(b):
        return shift_box(b, offset)
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
    for st in lay.steps:
        if st.kind != "forward":
            b = st.group.bbox()
            out[f"equilibrium-{st.n}-tail"] = moved((b[0] - 1, b[1] - 1, b[0] + 12, b[3] + 1))
            out[f"equilibrium-{st.n}-head"] = moved((b[2] - 12, b[1] - 1, b[2] + 1, b[3] + 1))
    for res in lay.resonance:
        for k, (_, _, g) in enumerate(res.arrows, 1):
            out[f"resonance-{res.n}-{k}"] = moved(g.bbox())
        for g in res.group.children:
            pts = g.children[0].points
            for corner, (x, y) in (("top", pts[1]), ("bottom", pts[2])):
                out[f"bracket-{res.n}-{g.id.rsplit('-', 1)[-1]}-{corner}"] = moved(pad_box((x, y, x, y), 6))
    return out


def cycle_crop_boxes(cyc: CycleLayout, offset: tuple[float, float]) -> dict[str, tuple]:
    """name -> box in the SVG frame: 12 pt round each arc arrow's tip (and an equilibrium arc's
    tail too, where its back shaft's head is), and each tail (its path, without its item)."""
    out = {}
    for a in cyc.arcs:
        ends = [("head", a.a1)] + ([("tail", a.a0)] if a.kind != "forward" else [])
        for end, ang in ends:
            x, y = arc_point(cyc.center, cyc.radius, ang)
            out[f"cycle-arc-{a.n}-{end}"] = shift_box(pad_box((x, y, x, y), 6), offset)
    for t in cyc.tails:
        out[f"cycle-{t.side}-{t.n}-{t.index}"] = shift_box(pad_box(t.group.children[0].bbox(), 2), offset)
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
           with_crops: bool = False, with_maps: bool = False, affinity: bool = False,
           with_cdxml: bool = True) -> list[Path]:
    fig = build_figure(spec_path)
    root, lay, cyc = fig.root, fig.scheme, fig.cycle
    svg, offset = to_svg(root, affinity=affinity)
    cdxml = None
    if with_cdxml:
        try:
            cdxml = fig.cdxml(offset)
        except CdxmlUnsupported as e:
            print(f"warning: no CDXML: {e}", file=sys.stderr)
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
    if cyc:
        print(f"  cycle: radius {cyc.radius:.2f} pt, free arc {cyc.free:.2f} pt", file=sys.stderr)
        for t in cyc.tails:
            if t.clash:
                print(f"  warning: cycle step {t.n} {t.side} tail {t.index}: {t.clash}", file=sys.stderr)
    paths = [svg_path, png_path]
    if cdxml is not None:
        cdxml_path = out_dir / f"{spec_path.stem}.cdxml"
        cdxml_path.write_text(cdxml)
        paths.append(cdxml_path)
    if with_crops and lay:
        paths += crops(svg, crop_boxes(lay, offset), out_dir, prefix=f"{spec_path.stem}_crop_")
    if with_crops and cyc:
        paths += crops(svg, cycle_crop_boxes(cyc, offset), out_dir, prefix=f"{spec_path.stem}_crop_")
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
                    help="also write zoomed crops of curved arrows, charged / dotted atoms, resonance and equilibrium arrows "
                         "and bracket corners")
    ap.add_argument("--maps", action="store_true",
                    help="also write <stem>_maps.png with each mapped atom's map number (preview only)")
    ap.add_argument("--affinity", action="store_true",
                    help="write for Affinity: arrowheads as separate shapes (it drops SVG markers), bold/italic as face names")
    ap.add_argument("--no-cdxml", dest="cdxml", action="store_false", help="don't write <stem>.cdxml for ChemDraw")
    ap.add_argument("--cdxml", dest="cdxml_legacy", action="store_true", help=argparse.SUPPRESS)  # the default; kept for old commands
    args = ap.parse_args(argv)
    paths = render(args.spec, args.out, args.px_per_pt, args.crops, args.maps, args.affinity, args.cdxml)
    for p in paths:
        print(p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
