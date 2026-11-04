"""Scheme, multi-row, resonance and equilibrium parameters per preset (acs-style.md §3, §4.3)."""

from __future__ import annotations

import statistics

from acs_chem_svg.geometry import add, length, mul, sub, unit

from .core import Doc, load, stats, union
from .measurements import _caption_glyphs, _r, _structures, measure, text_lines


# --------------------------------------------------------------------------- scheme summary (§3)

SCHEME_REFS = {"1.0": ["scheme_basic_x1.0", "scheme_multi_x1.0"], "0.8": ["scheme_basic_x0.8", "scheme_multi_x0.8"],
               "0.6": ["scheme_basic_x0.6", "scheme_multi_x0.6"]}


def scheme_summary(results: dict) -> dict:
    """Scheme parameters per preset (acs-style.md §3) from measure() results of
    scheme_basic_x<k> + scheme_multi_x<k>. Medians unless noted; `n` is the sample count."""
    def med(vals):
        vals = [v for v in vals if v is not None]
        return {"median": round(statistics.median(vals), 3), "min": round(min(vals), 3),
                "max": round(max(vals), 3), "n": len(vals)} if vals else None

    out = {}
    for k, names in SCHEME_REFS.items():
        ms = [results[n]["scheme"] for n in names if n in results]
        if not ms:
            continue
        arrows = [a for m in ms for a in m["arrows"]]
        lines = [(a, ln) for a in arrows for side in ("above", "below") for ln in a[side]]
        centered = [ln for _, ln in lines if not ln["left_justified"]]
        stepped = [ln for _, ln in lines if ln["left_justified"]]
        widest = [a["widest"] for a in arrows if "widest" in a and not a["widest"]["left_justified"]]
        pitches = [p for a in arrows for p in a["pitches"]]
        size = max(ln["size"] for _, ln in lines)
        plain = [p["pitch"] - p["size"] for p in pitches if not p["sub_above"] and not p["sub_below"]]
        rows = [r for m in ms for r in m["number_rows"]]
        ring_rows = [a for a in arrows if a["structures"] and all(s["ring_centers_minus_axis"] for s in a["structures"])]
        ink_rows = [a for a in arrows if a not in ring_rows]
        out[k] = {
            "arrow_length": med([a["length"] for a in arrows]),
            "conditions_pad (widest line: tail -> origin, end of advance -> notch)":
                med([w["pad_tail"] for w in widest] + [w["pad_notch"] for w in widest]),
            "conditions_center_minus_tail_notch_mid (advance)": med([ln["center_minus_tail_notch_mid"] for ln in centered]),
            "conditions_center_minus_tail_notch_mid (ink)": med([ln["ink_center_minus_tail_notch_mid"] for ln in centered]),
            "multistep x0 - tail (glyph origin)": med([ln["x0_minus_tail"] for ln in stepped]),
            "multistep ink left - tail": med([ln["ink_left_minus_tail"] for ln in stepped]),
            "above: axis - last baseline": med([a.get("above_axis_minus_last_base") for a in arrows]),
            "above: axis - ink bottom": med([a.get("above_axis_minus_ink_bottom") for a in arrows]),
            "below: first baseline - axis": med([a.get("below_first_base_minus_axis") for a in arrows]),
            "below: ink top - axis": med([a.get("below_ink_top_minus_axis") for a in arrows]),
            "line pitch - size, no subscripts": med(plain),
            "line pitch - size - plain, sub on upper line only": med(
                [p["pitch"] - p["size"] - plain[0] for p in pitches if p["sub_above"] and not p["sub_below"]] if plain else []),
            "line pitch - size - plain, sub on both lines": med(
                [p["pitch"] - p["size"] - plain[0] for p in pitches if p["sub_above"] and p["sub_below"]] if plain else []),
            "subscript size / size": med([ln["sub"]["size_ratio"] for _, ln in lines if "sub" in ln]),
            "subscript drop / size": med([ln["sub"]["drop_over_size"] for _, ln in lines if "sub" in ln]),
            "superscript size / size": med([ln["sup"]["size_ratio"] for _, ln in lines if "sup" in ln]),
            "superscript raise / size": med([ln["sup"]["raise_over_size"] for _, ln in lines if "sup" in ln]),
            "structure ring center - axis (rows where every structure has a ring)":
                med([s["ring_centers_minus_axis"][0] for a in ring_rows for s in a["structures"]]),
            "structure ink center - axis (other rows)": med([s["ink_center_minus_axis"] for a in ink_rows for s in a["structures"]]),
            "arrow gap (structure ink -> tail, head -> structure ink)":
                med([a["gap_left_struct_to_tail"] for a in arrows] + [a["gap_head_to_right_struct"] for a in arrows]),
            "plus gap (structure ink <-> plus ink)": med([p[g] for m in ms for p in m["pluses"]
                                                          for g in ("gap_left_struct", "gap_right_struct")]),
            "plus ink center - axis": med([p["ink_center_minus_axis"] for m in ms for p in m["pluses"]]),
            "number: row's lowest structure ink -> number ink top": med([g for r in rows for g in r["gap_row_bottom_to_ink_top"]]),
            "number: baseline spread within a row (max - min)": med([max(r["bases"]) - min(r["bases"]) for r in rows]),
            "number: center - structure ink center x": med([abs(c) for r in rows for c in r["center_minus_struct_center_x"]]),
            "caption size": size,
        }
    return out


WRAP_REFS = {"1.0": "scheme_wrap_x1.0", "0.8": "scheme_wrap_x0.8", "0.6": "scheme_wrap_x0.6"}


def wrap_summary(name: str | Doc) -> dict:
    """Multi-row scheme parameters (acs-style.md §3.4) from one scheme_wrap_x<k> file (or a
    loaded Doc: compare.py measures our layout the same way). Every structure, arrow and caption
    joins the row of the nearest arrow axis; a row's ink is the union of its structure ink,
    arrow boxes (axis ± head half-width) and caption ink."""
    doc = load(name) if isinstance(name, str) else name
    m = measure(doc)
    structs = _structures(doc)
    arrows = sorted(m["arrows"], key=lambda a: min(a["tail"][0], a["head"][0]))
    axes = sorted({round(a["tail"][1], 2) for a in arrows})

    def nearest(y):
        return min(axes, key=lambda ax: abs(y - ax))

    items = {ax: [] for ax in axes}  # (kind, label, ink)
    for fid, s in structs.items():
        items[nearest((s["ink"][1] + s["ink"][3]) / 2)].append(("structure", fid, s["ink"]))
    for a in arrows:
        hw = max(h["half_width_left"] + h["half_width_right"] for h in a["heads"].values() if h) / 2
        x0, x1 = sorted((a["tail"][0], a["head"][0]))
        items[nearest(a["tail"][1])].append(("arrow", f"{a['length']:.2f}", (x0, a["tail"][1] - hw, x1, a["tail"][1] + hw)))
    for cap in doc.captions:
        glyphs = _caption_glyphs(doc, cap)
        if glyphs:
            ink = union(ln["ink"] for ln in text_lines(glyphs))
            # numbers sit below their row: assign by the structure above, i.e. the nearest axis above the ink top
            above = [ax for ax in axes if ax <= ink[1]]
            bold = all(ln["bold"] for ln in text_lines(glyphs))
            ax = max(above) if bold and above else nearest((ink[1] + ink[3]) / 2)
            items[ax].append(("number" if bold else "caption", cap["text"], ink))
    rows = []
    for ax in axes:
        its = sorted(items[ax], key=lambda it: it[2][0])
        ink = union(it[2] for it in its)
        sink = union(it[2] for it in its if it[0] == "structure")
        row_arrows = [it for it in its if it[0] == "arrow"]
        rows.append({"axis_y": ax, "ink": ink, "structure_ink": sink, "width": ink[2] - ink[0],
                     "sequence": [f"{k}:{lbl}" for k, lbl, _ in its if k in ("structure", "arrow")],
                     "first": its[0][0], "last": [it for it in its if it[0] in ("structure", "arrow")][-1][0],
                     "arrow_lengths": [float(lbl) for _, lbl, _ in row_arrows],
                     "first_arrow_tail_minus_row_left": row_arrows[0][2][0] - ink[0] if row_arrows else None})
    pairs = []
    for r1, r2 in zip(rows, rows[1:]):
        pairs.append({"axis_pitch": r2["axis_y"] - r1["axis_y"],
                      "ink_gap (row bottom incl. numbers -> next row top incl. conditions)": r2["ink"][1] - r1["ink"][3],
                      "structure_ink_gap": r2["structure_ink"][1] - r1["structure_ink"][3],
                      "left_edge_diff (row2 - row1)": r2["ink"][0] - r1["ink"][0],
                      "right_edge_diff (row2 - row1)": r2["ink"][2] - r1["ink"][2],
                      "center_diff (row2 - row1)": (r2["ink"][0] + r2["ink"][2] - r1["ink"][0] - r1["ink"][2]) / 2})
    return _r({"file": doc.name, "overall_ink": union(r["ink"] for r in rows), "rows": rows, "between_rows": pairs})


RESONANCE_REFS = {"1.0": "resonance_x1.0", "0.8": "resonance_x0.8", "0.6": "resonance_x0.6"}


def _charge_box(c):
    """Ink box of a filled charge (bars of bar_length x bar_thickness; '+' adds the vertical bar)."""
    (x, y), h, t = c["center"], c["bar_length"] / 2, c["bar_thickness"] / 2
    return (x - h, y - (h if c["sign"] == "+" else t), x + h, y + (h if c["sign"] == "+" else t))


def resonance_summary(name: str | Doc) -> dict:
    """Resonance parameters (acs-style.md §4.3) from one resonance_x<k> file: the double-headed
    arrows (exact heads; hand-drawn length), structure ink (bonds, labels and charges: each
    charge joins the structure of its atom) -> arrow end gaps, the bracket gaps to the content
    ink, and the contributors' ring centers against the arrow axis. `name` may be a loaded Doc."""
    doc = load(name) if isinstance(name, str) else name
    m = measure(doc)
    structs = _structures(doc)
    frag_of = {a.id: a.frag for a in doc.atoms.values()}
    for c in m["charges"]:
        atom = min(doc.atoms.values(), key=lambda a: length(sub(a.p, c["center"])))
        s = structs[frag_of[atom.id]]
        s["ink_without_charges"] = s.get("ink_without_charges", s["ink"])
        s["ink"] = union([s["ink"], _charge_box(c)])
    order = sorted(structs.values(), key=lambda s: s["ink"][0])
    arrows = sorted(m["arrows"], key=lambda a: min(a["tail"][0], a["head"][0]))
    recs = []
    for a in arrows:
        x0, x1 = sorted((a["tail"][0], a["head"][0]))
        left = max((s for s in order if s["ink"][2] <= x0), key=lambda s: s["ink"][2])
        right = min((s for s in order if s["ink"][0] >= x1), key=lambda s: s["ink"][0])
        heads = [h for h in a["heads"].values() if h]
        recs.append({"length": a["length"], "axis_y": a["tail"][1],
                     "heads": {"length": stats([h["length"] for h in heads]),
                               "half_width": stats([(h["half_width_left"] + h["half_width_right"]) / 2 for h in heads]),
                               "inset": stats([h["inset_from_tip"] for h in heads])},
                     "shaft": a.get("shaft"),
                     "gap_left_ink_to_end": x0 - left["ink"][2], "gap_end_to_right_ink": right["ink"][0] - x1,
                     "gap_end_to_right_ink_without_charges": right.get("ink_without_charges", right["ink"])[0] - x1})
    brackets = m["brackets"]
    out = {"arrows": recs,
           "arrow_length": stats([r["length"] for r in recs]),
           "arrow_gap (structure ink incl. charges <-> arrow end)":
               stats([r["gap_left_ink_to_end"] for r in recs] + [r["gap_end_to_right_ink"] for r in recs]),
           "ring_center_minus_axis": [round(c[1] - recs[0]["axis_y"], 3) for s in order for c in s["ring_centers"]],
           "brackets": [{k: b[k] for k in ("lips_point", "height", "stroke_width", "lip_length", "gap_to_content_x",
                                            "extends_above_content", "extends_below_content")} for b in brackets],
           "bracket gap_x": stats([b["gap_to_content_x"] for b in brackets]),
           "bracket gap_above": stats([b["extends_above_content"] for b in brackets]),
           "bracket gap_below": stats([b["extends_below_content"] for b in brackets])}
    return _r(out)


EQUILIBRIUM_REFS = {"1.0": "equilibrium_x1.0", "0.8": "equilibrium_x0.8", "0.6": "equilibrium_x0.6",
                    "1.0 large": "arrows_large_x1.0"}


def equilibrium_summary(name: str | Doc) -> dict:
    """Equilibrium arrows and head variants (acs-style.md §3.1, §3.6) from one reference: per
    arrow its CDXML head attributes (HeadSize / ArrowheadCenterSize / ArrowheadWidth say which
    variant), then for equilibrium arrows the kind (balanced, or unbalanced favoring products or
    reactants: the long shaft's tip on the right or the left), shaft spacing and stagger, the
    short shaft's ratio and centering, and the half heads (forward shaft: the one pointing at
    Head3D; back: the other). Full heads (forward and ↔ arrows) as _arrow_entry measures them.
    Conditions (hand-placed) are provenance only: from the outer shafts and the arrow's center."""
    doc = load(name) if isinstance(name, str) else name
    m = measure(doc)
    by_axis = {round(a["axis_y"], 1): a for a in m["scheme"]["arrows"]}
    recs = []
    for a in sorted(m["arrows"], key=lambda a: a["tail"][1]):
        c = a["cdxml"]
        rec = {"head_variant": f"HeadSize {c['HeadSize']} / CenterSize {c['ArrowheadCenterSize']} / "
                               f"Width {c['ArrowheadWidth']}"}
        e = a.get("equilibrium")
        if not e:
            hs = [h for h in a["heads"].values() if h]
            rec.update({"kind": "resonance" if len(hs) == 2 else "forward", "length": a["length"],
                        "head_length": stats([h["length"] for h in hs]),
                        "head_half_width": stats([(h["half_width_left"] + h["half_width_right"]) / 2 for h in hs]),
                        "head_inset": stats([h["inset_from_tip"] for h in hs]), "shaft": a.get("shaft")})
            recs.append(rec)
            continue
        sh = e["shafts"]
        fwd = next(r for r in sh if r["direction"] == "forward")
        back = next(r for r in sh if r["direction"] == "back")
        rec["kind"] = e["kind"]
        if e["kind"] == "unbalanced":
            long = max(sh, key=lambda r: r["length"])
            tip_x = add(a["tail"], mul(unit(sub(a["head"], a["tail"])), long["tip_u"]))[0]
            rec["kind"] += " favoring " + ("products" if tip_x > (a["tail"][0] + a["head"][0]) / 2 else "reactants")
            rec.update({k: e[k] for k in ("long_length", "short_length", "short_center_minus_long_center",
                                          "long_over_short_butt_to_notch")})
        else:
            rec.update({k: e[k] for k in ("stagger_tail", "stagger_head")})
        rec.update({
            "spacing (center to center)": e["spacing"], "shafts_center_minus_axis": e["center_v"],
            "shaft_width": [r["width"] for r in sh], "barb": [r["barb"] for r in sh],
            "half_head_length": [r["head"]["length"] for r in sh],
            "half_head_barb_from_center": [r["head"]["barb_from_center"] for r in sh],
            "half_head_tip_from_center": [r["head"]["tip_from_center"] for r in sh],
            "half_head_inset": [r["head"]["inset_from_tip"] for r in sh],
            "shaft_end_from_tip": [r["u_end_from_tip"] for r in sh],
            "forward_back_control_per_barb": fwd["head"]["back_control_per_barb"],
            "back_back_control_from_center": back["head"]["back_control_from_center"],
            "back_control_shift (forward - back)": round(fwd["head"]["back_control_from_center"]
                                                         - back["head"]["back_control_from_center"], 3),
        })
        s = by_axis.get(round(min(a["tail"][1], a["head"][1]), 1)) or min(
            m["scheme"]["arrows"], key=lambda r: abs(r["axis_y"] - a["tail"][1]))
        rec["conditions (provenance)"] = {
            "above: top shaft - last baseline": s.get("above_axis_minus_last_base"),
            "below: first baseline - bottom shaft": s.get("below_first_base_minus_axis"),
            "line centers - arrow center": [ln["center_minus_tail_notch_mid"] for ln in s["above"] + s["below"]],
            "widest pad (each end)": s.get("widest") and round((s["length"] - s["widest"]["advance"]) / 2, 3)}
        recs.append(rec)
    return _r({"arrows": recs})
