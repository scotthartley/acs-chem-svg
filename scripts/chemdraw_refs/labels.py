"""Abbreviation and NH label parameters (acs-style.md §2, §3.5)."""

from __future__ import annotations

import math
import statistics
from pathlib import Path

from lxml import etree

from acs_chem_svg.geometry import unit

from .core import (
    ACS_PRESETS, Doc, Glyph, REF_DIR, _glyph_chain, _kern, _line_from_shape, assign, frame, label_glyphs, load, stats,
)
from .measurements import _r


ABBREV_REFS = ["abbrev_x1.0", "abbrev_x0.8", "abbrev_x0.6", "abbrev_2_x1.0", "abbrev_2_x0.8", "abbrev_2_x0.6",
               "ligands_x1.0", "ligands_x0.8", "ligands_x0.6",
               "metals_x1.0", "metals_x0.8", "metals_x0.6"]  # metals_x*: PPh3


def abbrev_summary(name: str, ref_dir: Path = REF_DIR, preset_name: str | None = None) -> dict:
    """§3.5: per abbreviation (abbrev_x*; PPh3 from metals_x*), each drawn form: displayed text, ChemDraw's
    LabelJustification, which glyph's advance center sits on the atom (its offset), baseline −
    atom y, and each bond end with its residual against the §2.1 model (preset values) on the
    anchor glyphs label_glyphs picks. `preset_name`: the preset to model with (default: the ACS
    preset of the file name's _x<k>)."""
    from acs_chem_svg import draw_mol
    from acs_chem_svg.preset import load_preset

    doc = load(name, ref_dir)
    asg = assign(doc)
    if preset_name is None:
        k = name.rsplit("_x", 1)[1]
        preset_name = next(p for p, s in ACS_PRESETS.items() if s == k)
    preset = load_preset(preset_name)
    L, r = preset.lengths, preset.ratios
    margin = float(doc.settings["MarginWidth"])
    forms, offsets, baselines, residuals, pitches = [], [], [], [], {}
    for a in doc.atoms.values():
        if not a.nickname:
            continue
        elem, whole = label_glyphs(doc, a)
        text = "".join(g.ch for g in whole)
        right = elem[-1] is whole[-1] and elem[0] is not whole[0]
        c = whole[-1] if right else whole[0]
        # a subscript at the centered end (N₃, P(o-tol)₃, Lₙ) clips as if its advance center were on the atom,
        # wherever it is drawn (N₃: drawn 0.70 right of that, as if full size). Not modeled: dppf above Pd,
        # whose bond ChemDraw stops below the p descenders (1.9 pt farther; x0.8 / x0.6: no bond drawn)
        # A first glyph kerned with the next (T|s in Ts) clips where its kerned advance would be
        # centered: -kern / 2 right of where the SVG draws it (draw_mol.Label.internal_dx)
        clip_glyphs = elem
        dx0 = (c.x + c.adv / 2 - a.p[0]) if c.size < float(doc.settings["LabelSize"]) else 0.0
        if not right and len(whole) > len(elem):
            dx0 += _kern(elem[-1], whole[len(elem)]) / 2
        outline = [[(x - a.p[0] - dx0, y - a.p[1]) for x, y in cont] for g in clip_glyphs for cont in g.outline()]
        ends = []
        for b in doc.bonds:
            if a.id not in (b.a, b.b) or not asg.get(b.id, {}).get("lines"):
                continue
            fr = frame(doc, b.b if b.a == a.id else b.a, a.id)
            end = statistics.median(fr.L - _line_from_shape(ln.shape, fr).u1 for ln in asg[b.id]["lines"])
            model = draw_mol.clip_distance_band(outline, (-fr.u[0], -fr.u[1]), r.label_band_per_size * preset.text.label_size, margin,
                                                r.label_clearance_per_margin,
                                                (L.label_outline_offset_x, L.label_outline_offset_y))
            ends.append({"dir_deg": round(math.degrees(math.atan2(fr.u[1], -fr.u[0])), 1), "end": _r(end),
                         "model_minus_end": _r(model - end)})
            residuals.append(model - end)
        off = c.x + c.adv / 2 - a.p[0]
        offsets.append((right, off))
        baselines.append(c.base - a.p[1])
        form = {"atom": a.id, "key": a.label, "text": text, "justification": a.justification, "node_type": a.node_type,
                "centered": c.ch, "anchor": "".join(g.ch for g in elem), "centered_dx": _r(off),
                "baseline_minus_atom_y": _r(c.base - a.p[1]), "bonds": ends}
        if a.alignment in ("Above", "Below"):
            form["alignment"] = a.alignment
            form["tail_line"] = tail = _tail_line(doc, elem[0], a.alignment)
            pitches.setdefault(a.alignment, []).append(tail["pitch"])
        forms.append(form)
    return {"forms": forms,
            "tail_line_pitch": {k: stats(v) for k, v in pitches.items()},
            "centered_dx_left_justified": stats(o for rt, o in offsets if not rt),
            "centered_dx_right_justified": stats(o for rt, o in offsets if rt),
            "baseline_minus_atom_y": stats(baselines),
            "bond_end_model_minus_end": stats(residuals),
            "bond_ends_over_0.1": sum(abs(x) > 0.1 for x in residuals)}


def _tail_line(doc: Doc, anchor: Glyph, alignment: str) -> dict:
    """A label set over two lines (ring NTs with LabelAlignment Above / Below, abbrev_2_x*): the
    tail's glyph chain on the nearest baseline above / below the anchor glyph, its text, its first
    glyph origin x minus the anchor's, and the baseline pitch."""
    sign = -1 if alignment == "Above" else 1
    cands = [g for g in doc.glyphs if g.box and abs(g.x - anchor.x) < anchor.size
             and 0.5 * anchor.size < sign * (g.base - anchor.base) < 1.3 * anchor.size]
    seed = min(cands, key=lambda g: (abs(g.x - anchor.x), abs(g.base - anchor.base)))
    whole = _glyph_chain(doc, [seed])
    return {"text": "".join(g.ch for g in whole), "origin_dx": _r(whole[0].x - anchor.x),
            "pitch": _r(abs(whole[0].base - anchor.base))}


NH_REFS = ["nh_x1.0", "nh_x0.8", "nh_x0.6"]


def nh_summary(name: str, ref_dir: Path = REF_DIR) -> dict:
    """§2: every NH of nh_x* (pyrrolidine turned in 15° steps, an acyclic R-NH-R at a peak, a
    valley and turned 15° / 30°, aziridine up and down). Per atom: the open direction (opposite
    the sum of its bond unit vectors; deg, counterclockwise from +x, y up) and its angle from
    vertical, the bond angle, ChemDraw's LabelAlignment (Left / Right: H beside; Above / Below),
    the H glyph origin and baseline relative to N's, and the CDXML <t> p (y relative to the
    usual first-line y, 3.59 at 10 pt) and LineStarts. Summary: the largest off-vertical angle
    with H above / below, the smallest with H beside, and the line pitch."""
    doc = load(name, ref_dir)
    root = etree.parse(str(ref_dir / f"{name}.cdxml")).getroot()
    nodes = {n.get("id"): n for n in root.iter("n")}
    atoms, pitch, stored = [], [], {}
    for a in doc.atoms.values():
        if a.element != "N" or a.label != "NH":
            continue
        nbs = [doc.atoms[b.b if b.a == a.id else b.a] for b in doc.bonds if a.id in (b.a, b.b)]
        us = [unit((o.p[0] - a.p[0], -(o.p[1] - a.p[1]))) for o in nbs]
        open_deg = math.degrees(math.atan2(-sum(u[1] for u in us), -sum(u[0] for u in us))) % 360
        off_vertical = min(abs((open_deg - 90 + 180) % 360 - 180), abs((open_deg - 270 + 180) % 360 - 180))
        bond_angle = math.degrees(math.acos(max(-1.0, min(1.0, us[0][0] * us[1][0] + us[0][1] * us[1][1]))))

        def nearest(ch):
            return min((g for g in doc.glyphs if g.ch == ch and g.box),
                       key=lambda g: math.dist((g.x + g.adv / 2, (g.box[1] + g.box[3]) / 2), a.p))
        n_g, h_g = nearest("N"), nearest("H")
        t = nodes[a.id].find("t")
        align = t.get("LabelAlignment")
        ty = float(t.get("p").split()[1]) - float(nodes[a.id].get("p").split()[1])
        if align in ("Above", "Below"):
            pitch.append(abs(h_g.base - n_g.base))
            stored.setdefault(align, []).append(ty)
        atoms.append({"atom": a.id, "open_deg": round(open_deg, 1), "off_vertical": round(off_vertical, 1),
                      "bond_angle": round(bond_angle, 1), "alignment": align, "line_starts": t.get("LineStarts"),
                      "h_dx": _r(h_g.x - n_g.x), "h_dbase": _r(h_g.base - n_g.base), "t_dy": _r(ty)})
    vertical = [r["off_vertical"] for r in atoms if r["alignment"] in ("Above", "Below")]
    beside = [r["off_vertical"] for r in atoms if r["alignment"] in ("Left", "Right")]
    ok = bool(vertical) and all((r["alignment"] == "Above") == (r["off_vertical"] <= max(vertical) and 0 < r["open_deg"] < 180)
                                for r in atoms if r["alignment"] != "Below")
    return {"atoms": sorted(atoms, key=lambda r: (r["bond_angle"], r["open_deg"])),
            "vertical_max_off_vertical": max(vertical, default=None), "beside_min_off_vertical": min(beside, default=None),
            "above_iff_open_upward": ok,
            "line_pitch": stats(pitch), "t_dy": {k: stats(v) for k, v in stored.items()}}
