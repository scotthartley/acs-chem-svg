"""The bond-to-label clip-model fit (style-spec.md §2.1) and the key-parameter summary table."""

from __future__ import annotations

import math
import statistics
from pathlib import Path

import pypdfium2

from .core import REF_DIR, _line_from_shape, all_names, assign, frame, label_glyphs, load


# --------------------------------------------------------------------------- clip-model fit


def clip_rows(names=None) -> list[dict]:
    """One row per drawn bond end at a label: the measured distance from the label atom to
    where the bond stops (plain lines: line end; wedges and hashes: the wide end), with the
    direction from the label atom toward the bond's other atom (unit, y down) and the element
    glyph outline as ChemDraw placed it, relative to the atom (y down)."""
    rows = []
    for n in names or all_names():
        doc = load(n)
        asg = assign(doc)
        margin = float(doc.settings["MarginWidth"])
        for b in doc.bonds:
            r = asg.get(b.id)
            if not r:
                continue
            for atom_id in (b.a, b.b):
                atom = doc.atoms[atom_id]
                if not atom.label:
                    continue
                elem, _ = label_glyphs(doc, atom)
                if not elem:
                    continue
                other = b.b if atom_id == b.a else b.a
                fr = frame(doc, other, atom_id)  # u runs toward the label
                if "stereo" in r:
                    narrow = b.a if b.display.endswith("Begin") else b.b
                    if narrow == atom_id or not r["shapes"]:
                        continue
                    end = fr.L - max(u for _, loc in r["shapes"] for u, _ in loc)
                    kind = "hash" if "Hash" in b.display else "wedge"
                elif r.get("lines"):
                    end = statistics.median(fr.L - _line_from_shape(ln.shape, fr).u1 for ln in r["lines"])
                    kind = "line"
                else:
                    continue
                outline = [[(x - atom.p[0], y - atom.p[1]) for x, y in c] for g in elem for c in g.outline()]
                rows.append({"file": n, "atom": atom_id, "bond": b.id, "kind": kind, "outline": outline,
                             "element": "".join(g.ch for g in elem), "size": elem[0].size, "margin": margin,
                             "dir": (-fr.u[0], -fr.u[1]), "end": end,
                             "dir_deg": round(math.degrees(math.atan2(fr.u[1], -fr.u[0])), 1)})
    return rows


def fit_label_clip(rows=None) -> dict:
    """Fit the bond-to-label clip model used by draw_mol.clip_distance (style-spec.md §2.1).

    A bond stops c x MarginWidth beyond the furthest point, along the bond, of the element
    glyph outline (as ChemDraw placed it) within ±band of the bond axis, with
    band = beta x label size and the outline shifted by (dx, dy) pt. Grid search: fewest
    rows off by more than 0.1 pt, then the worst residual among the rest, then rms.
    """
    from chem_svg_cdxml import draw_mol

    rows = rows if rows is not None else clip_rows()
    best = None
    for beta in [x / 1000 for x in range(170, 211, 5)]:
        for dx in [x / 100 for x in range(-5, 6)]:
            for dy in [x / 100 for x in range(-5, 6)]:
                ext = [draw_mol.clip_distance_band(r["outline"], r["dir"], beta * r["size"], r["margin"], 0.0,
                                                   (dx, dy)) for r in rows]
                for c in [x / 10000 for x in range(9900, 10301, 25)]:
                    res = [e + c * r["margin"] - r["end"] for e, r in zip(ext, rows)]
                    inside = [abs(x) for x in res if abs(x) <= 0.2]  # worst of the rows the model can fit
                    key = (sum(abs(x) > 0.1 for x in res), round(max(inside), 3),
                           math.sqrt(statistics.fmean(x * x for x in res)))
                    if best is None or key < best[0]:
                        best = (key, beta, dx, dy, c, res)
    (n_bad, worst, rms), beta, dx, dy, c, res = best
    for r, x in zip(rows, res):
        r["residual"] = round(x, 3)
    return {"band_per_size": beta, "offset": (dx, dy), "per_margin": c, "n": len(rows),
            "n_over_0.1": n_bad, "worst": worst, "rms": round(rms, 3),
            "rows": [{k: v for k, v in r.items() if k not in ("outline", "dir")} for r in rows]}


# --------------------------------------------------------------------------- summary


def summarize(m: dict) -> dict:
    """Flat key-parameter table from one measure() result (medians, pt)."""
    def med(vals):
        vals = [v for v in vals if v is not None]
        return round(statistics.median(vals), 3) if vals else None

    labels = [e for e in m["labels"] if "bonds" in e]
    lb = [b for e in labels for b in e["bonds"]]
    out = {
        "bond_length": (m["bond_length"] or {}).get("median"),
        "line_width": (m["line_width"] or {}).get("median"),
        "bond_spacing_ring": (m["bond_spacing_ring"] or {}).get("median"),
        "ring_inner_inset": (m["ring_inner_line"]["inset_unlabeled_ends"] or {}).get("median"),
        "acyclic_double_spacing": med([e["spacing"] for e in m["acyclic_double"]]),
        "triple_spacing": med([(e["v"][2] - e["v"][0]) / 2 for e in m["triple"]]),
        "label_baseline_below_atom": med([e["baseline_minus_atom_y"] for e in labels]),
        "label_anchor_dx": med([abs(e["anchor_advance_center_x_minus_atom"]) for e in labels
                                if len(e["element"]) == 1]),
        "end_C=O_vertical_from_O": med([b["end_dist_from_atom"] for e in labels for b in e["bonds"]
                                        if e["label"] == "O" and b["order"] == 2 and abs(b["dir_deg"]) == 90]),
        "end_OH_30deg_from_O": med([b["end_dist_from_atom"] for e in labels for b in e["bonds"]
                                    if e["label"] == "OH" and e["glyph_run"] == "OH" and b["dir_deg"] == -30]),
        "end_horizontal_into_N_from_left": med([b["end_dist_from_atom"] for e in labels for b in e["bonds"]
                                                if e["element"] == "N" and b["dir_deg"] == 0]),
        "gap_to_glyph_outline": med([b["gap_to_elem_outline"] for b in lb]),
        "axial_gap_hull_over_margin": med([b["axial_gap_to_hull_over_margin"] for e in labels
                                           for b in e["bonds"] if e["label"] not in ("H", "I")]),
        "subscript_size_ratio": med([o["size"] / e["size"] for e in labels for o in e["other_glyphs"]
                                     if o["ch"].isdigit()]),
        "subscript_drop": med([o["dbaseline"] for e in labels for o in e["other_glyphs"] if o["ch"].isdigit()]),
    }
    for st in m["stereo"]:
        if "Hash" in st["display"]:
            out["hash_count"] = st["count"]
            out["hash_pitch"] = (st["pitch"] or {}).get("median")
            out["hash_line_thickness"] = (st["line_thickness"] or {}).get("median")
            out["hash_last_len"] = st["hashes"][-1]["len_v"] if st["hashes"] else None
        else:
            out["wedge_wide_end_width"] = st["wide_end_width"]
            out["wedge_end_u"] = st["u_range"][1]
    if m["charges"]:
        out["charge_bar_length"] = med([c["bar_length"] for c in m["charges"]])
        out["charge_bar_thickness"] = med([c["bar_thickness"] for c in m["charges"]])
    heads = [h for a in m["arrows"] for h in a["heads"].values()]
    if heads:
        out["arrowhead_length"] = med([h["length"] for h in heads])
        out["arrowhead_half_width"] = med([(h["half_width_left"] + h["half_width_right"]) / 2 for h in heads])
        out["arrowhead_inset"] = med([h["inset_from_tip"] for h in heads])
        out["arrow_shaft_width"] = med([a["shaft"]["width"] for a in m["arrows"] if "shaft" in a])
    return out



def pdf_char_origins(name: str, ref_dir: Path = REF_DIR) -> tuple[str, list[float]]:
    """(page text, x origin of each character) from the reference PDF. These are ChemDraw's own
    glyph positions, independent of our kerning (the SVG gives one origin per run only)."""
    import ctypes

    import pypdfium2.raw as raw

    tp = pypdfium2.PdfDocument(str(ref_dir / f"{name}.pdf"))[0].get_textpage()
    xs = []
    for i in range(tp.count_chars()):
        x, y = ctypes.c_double(), ctypes.c_double()
        raw.FPDFText_GetCharOrigin(tp.raw, i, x, y)
        xs.append(x.value)
    return tp.get_text_range(), xs
