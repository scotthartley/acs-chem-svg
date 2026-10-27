"""Compare our renderer with a ChemDraw reference on identical layouts (design.md §6.6).

Each fragment of the reference is rendered from its template (templates/m1/, written by
cdxml_to_molfile.py) and placed on ChemDraw's page frame. Both drawings are then
measured with the same code in measure_refs.py:
- Our atoms and bonds come from the Molecule, with ids equal to the CDXML ids.
- Our strokes become butt-capped rectangles.
- Our <text> goes through the same glyph loader.

Scheme references (scheme_*) also get our reaction arrows at ChemDraw's tail and head (resonance
arrows for ChemDraw's double-headed ones; the head variant from each arrow's CDXML HeadSize,
ARROW_HEADS), our brackets on ChemDraw's bracket lines, and
each caption line set by richtext with its first glyph at ChemDraw's (scheme_captions), so the
exact items (arrowheads, text advances, sub/superscripts) compare one to one. Mechanism
references (mech_*) get those too, plus our curved arrows drawn on ChemDraw's CurvePoints
(mechanism.curved_arrow) and our dots at ChemDraw's lone-pair midpoints and radical centers
(mech_dots): the exact items are the heads, strokes and dot sizes. Our own dot placement and
arrow layout (approximate) are what examples/m3_*.yaml show. Equilibrium references
(equilibrium_*, arrows_large_*) get our equilibrium arrows on ChemDraw's (ref_arrow) and the
captions only: their "B" is a lone boron atom standing in for a structure, with no template.

Prints the summarize() table (ChemDraw vs ours) and the per-bond end distances at labels.
With --spec, also renders an M2 spec with our own layout and prints its approximate scheme
items (arrow lengths, gaps, conditions offsets) next to the reference's; for scheme_wrap_* also
row widths, the row gap and alignment, and for resonance_* the arrow and bracket gaps.
Writes the side-by-side images to build/compare/:
    <ref>_ours.svg / .png     our drawing on the ChemDraw page
    <ref>_compare.png         ChemDraw | ours | overlay
    <ref>_zoom.png            per-fragment crops at 10 px/pt

Usage: uv run python scripts/compare.py m1_set_x1.0 [more refs] [--preset acs1996]
       uv run python scripts/compare.py scheme_multi_x1.0 --spec examples/m2_scheme_multi_x1.0.yaml
       uv run python scripts/compare.py scheme_wrap_x1.0 --spec examples/m4_wrap_x1.0.yaml
       uv run python scripts/compare.py resonance_x1.0 --spec examples/m4_resonance_x1.0.yaml
"""

from __future__ import annotations

import argparse
import copy
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import measure_refs as mr  # noqa: E402
from cdxml_to_molfile import fragments, template_path  # noqa: E402

from acs_chem_svg import chem, richtext  # noqa: E402
from acs_chem_svg.arrows import equilibrium_arrow, reaction_arrow, resonance_arrow  # noqa: E402
from acs_chem_svg.brackets import square_bracket  # noqa: E402
from acs_chem_svg.draw_mol import Drawer  # noqa: E402
from acs_chem_svg.preset import Preset, load_preset  # noqa: E402
from acs_chem_svg.preview import compare_png, rasterize_file, zoom_crops  # noqa: E402
from acs_chem_svg.geometry import add, mul, sub, unit  # noqa: E402
from acs_chem_svg.layout import build_scheme  # noqa: E402
from acs_chem_svg.mechanism import curved_arrow  # noqa: E402
from acs_chem_svg.scene import Circle, Group  # noqa: E402
from acs_chem_svg.spec import load_spec  # noqa: E402
from acs_chem_svg.svg_writer import to_svg  # noqa: E402

OUT = mr.REPO / "build" / "compare"
ARROW_HEADS = {"1000": "standard", "1500": "large"}  # CDXML HeadSize -> preset arrow.heads variant
ARROWS_ONLY = ("equilibrium_", "arrows_large_")  # references drawn without their structures


def ref_arrow(a: dict, preset: Preset, name: str):
    """Our arrow on a CDXML arrow (load_cdxml dict), with its head variant. ChemDraw's Head3D /
    Tail3D are the forward shaft's tip and butt; a balanced equilibrium arrow's ends (ours)
    reach equilibrium.stagger beyond them, and an unbalanced one's long shaft is the forward one."""
    variant = ARROW_HEADS.get(a["attrs"].get("HeadSize"))
    if a.get("equilibrium"):
        tail, head = a["tail"], a["head"]
        if a["attrs"].get("ArrowEquilibriumRatio"):
            return equilibrium_arrow(tail, head, preset, name, "products", variant)
        d = mul(unit(sub(head, tail)), preset.equilibrium.stagger)
        return equilibrium_arrow(sub(tail, d), add(head, d), preset, name, None, variant)
    kind = resonance_arrow if a["attrs"].get("ArrowheadTail") == "Full" else reaction_arrow
    return kind(a["tail"], a["head"], preset, name, variant)


def preset_for(name: str) -> str:
    for k in ("0.8", "0.6"):
        if name.endswith(f"_x{k}"):
            return f"acs1996_x{k}"
    return "acs1996"


@dataclass
class Rendered:
    name: str
    preset: Preset
    svg: str
    ref: mr.Doc
    ours: mr.Doc
    molecules: list[chem.Molecule]


def render_reference(name: str, preset_name: str | None = None) -> Rendered:
    ref = mr.load(name)
    preset = load_preset(preset_name or preset_for(name))
    groups, mols, drawers = [], [], []
    for n, _, _ in [] if name.startswith(ARROWS_ONLY) else fragments(ref):
        mol = chem.from_molfile(template_path(name, n), f"{name}__{n}", preset.lengths.bond_length)
        d = Drawer(mol, preset)
        g = d.draw(f"frag{n}")
        g.children = [c for c in g.children if not c.id.endswith(("-radicals", "-lone-pairs"))]  # ChemDraw's dots below
        groups.append(g)
        mols.append(mol)
        drawers.append(d)
    if name.startswith(("scheme_", "mech_", "resonance_") + ARROWS_ONLY):  # M1 references keep their approved renders (no captions)
        groups += [ref_arrow(a, preset, f"arrow-{k}") for k, a in enumerate(ref.arrows, 1)]
        groups += scheme_captions(ref, preset)
        groups += [square_bracket(b["x"], b["y_top"], b["y_bottom"], "left" if b["lips_point"] == "right" else "right",
                                  preset, f"bracket-{k}") for k, b in enumerate(mr._brackets(copy.deepcopy(ref)), 1)]
    if name.startswith("mech_"):
        groups += [curved_arrow(*curve_points(c), "two_electron" if c["head"] == "Full" else "fishhook", preset,
                                f"curved-arrow-{k}") for k, c in enumerate(ref.curves, 1)]
        groups.append(mech_dots(ref, preset))
    # templates hold ChemDraw's SVG-frame coordinates (measure_refs atoms), so no offset
    svg, (dx, dy) = to_svg(Group(f"compare-{name}", groups), size=ref.svg_size, offset=(0.0, 0.0),
                           affinity=True)  # load_rendered reads head paths, not <defs>
    pair_ids = {frozenset((b.a, b.b)): b.id for b in ref.bonds}
    atoms, bonds = {}, []
    for mol, d in zip(mols, drawers):
        for a in mol.atoms:
            label = None
            if a.idx in d.labels:
                label = a.abbrev or a.element + ("H" + (str(a.num_h) if a.num_h > 1 else "") if a.num_h else "")
            atoms[a.id] = mr.Atom(a.id, (a.x + dx, a.y + dy), a.element, label, a.charge, a.num_h, mol.name,
                                  nickname=a.abbrev is not None)
        for b in mol.bonds:
            ia, ib = mol.atoms[b.a].id, mol.atoms[b.b].id
            display = {"wedge": "WedgeBegin", "hash": "WedgedHashBegin"}.get(b.stereo, "")
            bonds.append(mr.Bond(pair_ids[frozenset((ia, ib))], ia, ib, b.order, display, frag=mol.name))
    ours = mr.load_rendered(svg, name + "_ours", atoms, bonds, copy.deepcopy(ref.settings), ref.arrows, ref.curves)
    return Rendered(name, preset, svg, ref, ours, mols)


def curve_points(c: dict) -> tuple:
    """The drawn cubic of a CDXML curve: anchor 0, its out-handle, anchor 1's in-handle, anchor 1.
    (ChemDraw's fishhooks add a third anchor on top of the second, a zero-length retrace.)"""
    p = c["pts"]
    return p[1], p[2], p[3], p[4]


def mech_dots(ref: mr.Doc, preset: Preset) -> Group:
    """Our dots (preset radius and spacing) at ChemDraw's lone-pair midpoints, along each pair's
    own direction, and at its radical centers."""
    pairs, singles = mr._dot_sets(ref)
    lp, out = preset.lone_pair, []
    for a, b in pairs:
        mid, d = mul(add(a.center, b.center), 0.5), unit(sub(b.center, a.center))
        out += [Circle(add(mid, mul(d, -lp.dot_spacing / 2)), lp.dot_radius),
                Circle(add(mid, mul(d, lp.dot_spacing / 2)), lp.dot_radius)]
    out += [Circle(c.center, lp.dot_radius) for c in singles]
    return Group("dots", out)


def richtext_source(runs) -> str:
    """CDXML <s> runs -> richtext input. Formula runs (face 96) are auto-subscripted like ours;
    italic runs (face bit 2) are wrapped in *...*; plain runs have digits and signs escaped."""
    out = []
    for text, _, face in runs:
        f = int(face or 0)
        if f == 96:
            out.append(text)
            continue
        t = "".join("\\" + ch if ch.isdigit() or ch in "+\u2212*\\" else ch for ch in text)
        out.append(f"*{t}*" if f & 2 else t)
    return "".join(out)


def scheme_captions(ref: mr.Doc, preset: Preset) -> list[Group]:
    """Every reference caption set by richtext, one line per CDXML line, each with its first
    glyph origin where ChemDraw put it. Bold captions (compound numbers) are set bold."""
    out = []
    for k, cap in enumerate(ref.captions, 1):
        glyphs = mr._caption_glyphs(ref, cap)
        if not glyphs:
            continue
        lines = mr.text_lines(glyphs)
        src = richtext_source(cap["runs"]).split("\n")
        if len(src) != len(lines):
            raise ValueError(f"{ref.name}: caption {cap['text']!r}: {len(src)} CDXML lines, {len(lines)} drawn")
        runs = []
        for text, ln in zip(src, lines):
            laid = richtext.layout_line(text, ln["x0"], ln["base"], ln["size"], preset, f"caption:{k}",
                                        style="bold" if cap["bold"] else "regular")
            runs += laid.runs
        out.append(Group(f"caption-{k}", runs))
    return out


def label_bond_ends(m: dict) -> dict[tuple[str, str, float], dict]:
    """(atom id, bond id, line v rounded) -> label-bond entry from measure()."""
    out = {}
    for e in m["labels"]:
        for b in e.get("bonds", []):
            out[(e["atom"], b["bond"], round(b["line_v"], 1))] = dict(b, label=e["label"], run=e["glyph_run"])
    return out


def report(r: Rendered) -> dict:
    m_ref, m_our = mr.measure(r.ref), mr.measure(r.ours)
    s_ref, s_our = mr.summarize(m_ref), mr.summarize(m_our)
    print(f"\n== {r.name} (preset {r.preset.name})")
    print(f"{'parameter':34}{'ChemDraw':>12}{'ours':>12}{'diff':>9}")
    for k in dict.fromkeys(list(s_ref) + list(s_our)):
        a, b = s_ref.get(k), s_our.get(k)
        diff = f"{b - a:+.3f}" if isinstance(a, (int, float)) and isinstance(b, (int, float)) else ""
        print(f"{k:34}{str(a):>12}{str(b):>12}{diff:>9}")
    ends_ref, ends_our = label_bond_ends(m_ref), label_bond_ends(m_our)
    print(f"\n{'label bond end':40}{'dir':>7}{'ChemDraw':>10}{'ours':>9}{'diff':>9}")
    for key in sorted(ends_ref, key=lambda k: (int(k[0]), int(k[1]), k[2])):
        e = ends_ref[key]
        o = ends_our.get(key)
        tag = f"{e['run']} atom {key[0]} bond {key[1]} v {key[2]:+.1f}"
        if o is None:
            print(f"{tag:40}{e['dir_deg']:>7}{e['end_dist_from_atom']:>10}{'—':>9}")
            continue
        d = o["end_dist_from_atom"] - e["end_dist_from_atom"]
        print(f"{tag:40}{e['dir_deg']:>7}{e['end_dist_from_atom']:>10}{o['end_dist_from_atom']:>9}{d:>+9.3f}"
              f"{'  !' if abs(d) > 0.1 else ''}")
    return {"ref": m_ref, "ours": m_our}


def images(r: Rendered) -> list[Path]:
    OUT.mkdir(parents=True, exist_ok=True)
    svg_path = OUT / f"{r.name}_ours.svg"
    svg_path.write_text(r.svg)
    rasterize_file(svg_path, OUT / f"{r.name}_ours.png", 4.0)
    ref_svg = (mr.REF_DIR / f"{r.name}.svg").read_text()
    w, h = r.ref.svg_size
    cmp = compare_png(ref_svg, r.svg, w, h, OUT / f"{r.name}_compare.png",
                      labels=(f"ChemDraw 26.1 ({r.name})", f"acs-chem-svg ({r.preset.name})"))
    boxes = [b for fid, b in r.ref.fragments.items()
             if any(a.frag == fid and not a.nickname for a in r.ref.atoms.values())]
    if r.name.startswith(ARROWS_ONLY):  # both ends of every arrow
        boxes = [(x - 12, a["tail"][1] - 6, x + 12, a["tail"][1] + 6) for a in r.ref.arrows
                 for x in sorted((a["tail"][0], a["head"][0]))]
    zoom = zoom_crops(ref_svg, r.svg, w, h, boxes, OUT / f"{r.name}_zoom.png")
    return [svg_path, cmp, zoom]


# --------------------------------------------------------------------------- our layout vs ChemDraw's (approx)

APPROX_KEYS = ["length", "gap_left_struct_to_tail", "gap_head_to_right_struct", "above_axis_minus_last_base",
               "below_first_base_minus_axis"]


def spec_doc(spec_path: Path, ref: mr.Doc) -> tuple[mr.Doc, str]:
    """Render an M2 spec with our layout and load it as a Doc (atoms from the placed molecules,
    arrows from the layout), so measure_refs._scheme measures it like the reference."""
    lay = build_scheme(load_spec(spec_path), spec_path.parent, spec_path.stem)
    svg, (dx, dy) = to_svg(lay.root, affinity=True)  # head paths, as above
    atoms, bonds = {}, []
    for name, pl in lay.species.items():
        d = Drawer(pl.mol, lay.preset)
        ox, oy = pl.offset[0] + dx, pl.offset[1] + dy
        for a in pl.mol.atoms:
            label = None
            if a.idx in d.labels:
                label = a.abbrev or a.element + ("H" + (str(a.num_h) if a.num_h > 1 else "") if a.num_h else "")
            atoms[f"{name}:{a.id}"] = mr.Atom(f"{name}:{a.id}", (a.x + ox, a.y + oy), a.element, label, a.charge,
                                             a.num_h, name, nickname=a.abbrev is not None)
        for b in pl.mol.bonds:
            display = {"wedge": "WedgeBegin", "hash": "WedgedHashBegin"}.get(b.stereo, "")
            bonds.append(mr.Bond(f"{name}:{b.id}", f"{name}:{pl.mol.atoms[b.a].id}", f"{name}:{pl.mol.atoms[b.b].id}",
                                 b.order, display, frag=name))
    arrows = [{"tail": (st.tail[0] + dx, st.tail[1] + dy), "head": (st.head[0] + dx, st.head[1] + dy), "attrs": {},
               "equilibrium": st.kind != "forward"} for st in lay.steps]
    arrows += [{"tail": (t[0] + dx, t[1] + dy), "head": (h[0] + dx, h[1] + dy), "attrs": {"ArrowheadTail": "Full"}}
               for res in lay.resonance for t, h, _ in res.arrows]
    return mr.load_rendered(svg, spec_path.stem, atoms, bonds, copy.deepcopy(ref.settings), arrows), svg


def _row(label, a, b):
    diff = f"{b - a:+.3f}" if isinstance(a, (int, float)) and isinstance(b, (int, float)) else ""
    print(f"{label:40}{'' if a is None else round(a, 3):>10}{'' if b is None else round(b, 3):>10}{diff:>9}")


def resonance_report(ours: mr.Doc, ref_name: str) -> None:
    """Resonance placement (acs-style.md §4.3), ours vs ChemDraw: arrow lengths and gaps,
    bracket gaps."""
    r, o = mr.resonance_summary(ref_name), mr.resonance_summary(ours)
    print(f"\n== our layout vs {ref_name}: resonance placement (approx; no tolerance)")
    print(f"{'item':40}{'ChemDraw':>10}{'ours':>10}{'diff':>9}")
    for k, (ra, oa) in enumerate(zip(r["arrows"], o["arrows"]), 1):
        for key in ("length", "gap_left_ink_to_end", "gap_end_to_right_ink"):
            _row(f"arrow {k} {key}", ra[key], oa[key])
    for key in ("bracket gap_x", "bracket gap_above", "bracket gap_below"):
        _row(key, r[key]["median"], o[key]["median"])


def wrap_report(ours: mr.Doc, ref_name: str) -> None:
    """Multi-row items (acs-style.md §3.4), ours vs ChemDraw: row widths, the row gap and the
    left alignment."""
    r, o = mr.wrap_summary(ref_name), mr.wrap_summary(ours)
    print(f"\n== our layout vs {ref_name}: rows (approx; no tolerance)")
    print(f"{'item':40}{'ChemDraw':>10}{'ours':>10}{'diff':>9}")
    for k, (rr, orow) in enumerate(zip(r["rows"], o["rows"]), 1):
        _row(f"row {k} width", rr["width"], orow["width"])
        print(f"{'row ' + str(k) + ' sequence':40}{len(rr['sequence']):>10}{len(orow['sequence']):>10}")
    key = "ink_gap (row bottom incl. numbers -> next row top incl. conditions)"
    for k, (pr, po) in enumerate(zip(r["between_rows"], o["between_rows"]), 1):
        _row(f"rows {k}-{k + 1} ink gap", pr[key], po[key])
        _row(f"rows {k}-{k + 1} left edge diff", pr["left_edge_diff (row2 - row1)"], po["left_edge_diff (row2 - row1)"])


def spec_report(spec_path: Path, ref_name: str, ref_arrows: list[int] | None = None) -> None:
    """Approximate scheme items, ours vs ChemDraw, per arrow. `ref_arrows` are 1-based indices
    of the reference arrows (top to bottom, left to right) the spec's steps stand for.
    Resonance and multi-row references add their own items."""
    ref = mr.load(ref_name)
    ours, _ = spec_doc(spec_path, ref)
    if ref_name.startswith("resonance_"):
        resonance_report(ours, ref_name)
        return
    if ref_name.startswith("scheme_wrap_"):
        wrap_report(ours, ref_name)
    m_ref, m_our = mr.measure(ref)["scheme"], mr.measure(ours)["scheme"]
    ra = sorted(m_ref["arrows"], key=lambda a: (round(a["axis_y"]), a["tail"][0]))
    oa = m_our["arrows"]
    pick = [ra[i - 1] for i in ref_arrows] if ref_arrows else ra[:len(oa)]
    print(f"\n== {spec_path.name} (our layout) vs {ref_name} (approx items; no tolerance)")
    print(f"{'item':40}{'ChemDraw':>10}{'ours':>10}{'diff':>9}")
    row = _row
    for k, (r, o) in enumerate(zip(pick, oa), 1):
        for key in APPROX_KEYS:
            row(f"step {k} {key}", r.get(key), o.get(key))
        for side in ("above", "below"):
            for j, (lr, lo) in enumerate(zip(r[side], o[side]), 1):
                row(f"step {k} {side} {j} center - tail/notch mid", lr["center_minus_tail_notch_mid"],
                    lo["center_minus_tail_notch_mid"])
        for j, (pr, po) in enumerate(zip(r["pitches"], o["pitches"]), 1):
            row(f"step {k} line pitch {j}", pr["pitch"], po["pitch"])
    for k, (pr, po) in enumerate(zip(m_ref["pluses"], m_our["pluses"]), 1):
        for key in ("gap_left_struct", "gap_right_struct", "ink_center_minus_axis"):
            row(f"plus {k} {key}", pr[key], po[key])
    for rr, ro in zip(m_ref["number_rows"], m_our["number_rows"]):
        for t, gr, go in zip(rr["texts"], rr["gap_row_bottom_to_ink_top"], ro["gap_row_bottom_to_ink_top"]):
            row(f"number {t}: row bottom -> ink top", gr, go)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="+")
    ap.add_argument("--preset", help="override the preset (default: from the reference name)")
    ap.add_argument("--no-images", action="store_true")
    ap.add_argument("--spec", type=Path, help="also compare this M2 spec's own layout with the (first) reference")
    ap.add_argument("--ref-arrows", help="1-based reference arrows the spec's steps stand for, e.g. 4 or 1,2,3")
    args = ap.parse_args(argv)
    if args.spec:
        spec_report(args.spec, args.names[0], [int(x) for x in args.ref_arrows.split(",")] if args.ref_arrows else None)
    for name in args.names:
        r = render_reference(name, args.preset)
        report(r)
        if not args.no_images:
            for p in images(r):
                print(p.relative_to(mr.REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
