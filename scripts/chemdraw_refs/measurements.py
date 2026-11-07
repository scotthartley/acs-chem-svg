"""Per-reference measurements (measure()): bonds, labels, charges, dots, arrows, captions, schemes."""

from __future__ import annotations

import math
import re
import statistics

from chem_svg_cdxml.fontmetrics import load_face
from chem_svg_cdxml.geometry import (
    add, convex_hull, dot, length, mul, ray_box_gap, ray_polygon, seg_seg_dist, sub, unit,
)

from .core import Doc, _box_center, _line_from_shape, assign, frame, label_glyphs, neighbors, stats, union


# --------------------------------------------------------------------------- measurements


def measure(doc: Doc) -> dict:
    out: dict = {"file": doc.name, "settings": doc.settings}
    out["units"] = {"svg_viewbox": doc.svg_size, "pdf_page_pt": doc.page_pt,
                    "svg_unit_is_pt": bool(doc.page_pt) and all(
                        abs(a - b) < 0.01 for a, b in zip(doc.svg_size, doc.page_pt))}
    asg = assign(doc)
    bl = float(doc.settings["BondLength"])
    sp_nominal = bl * float(doc.settings["BondSpacing"]) / 100

    # bond lengths (CDXML atom coordinates)
    out["bond_length"] = stats(frame(doc, b.a, b.b).L for b in doc.bonds)

    # line widths of all plain lines
    widths = [ln.width for r in asg.values() for ln in r.get("lines", [])]
    out["line_width"] = stats(widths)

    # double / triple bonds
    ring_inner, ring_spacing, acyclic, triples = [], [], [], []
    for b in doc.bonds:
        r = asg.get(b.id)
        if not r or "lines" not in r:
            continue
        lines, fr = r["lines"], r["frame"]
        A, B = doc.atoms[b.a], doc.atoms[b.b]
        if b.order == 2 and len(lines) == 2:
            v = [ln.v for ln in lines]
            spacing = v[1] - v[0]
            if b.ring:
                main = min(lines, key=lambda ln: abs(ln.v))
                inner = max(lines, key=lambda ln: abs(ln.v))
                ring_spacing.append(spacing)
                ring_inner.append({
                    "bond": b.id, "spacing": round(spacing, 3),
                    "inset_a": round(inner.u0, 3), "inset_b": round(fr.L - inner.u1, 3),
                    "main_u": (round(main.u0, 3), round(main.u1, 3)), "L": round(fr.L, 3),
                    "a_labeled": bool(A.label), "b_labeled": bool(B.label),
                })
            else:
                acyclic.append(_acyclic_entry(doc, b, fr, lines, spacing))
        elif b.order == 3 and len(lines) == 3:
            triples.append({"bond": b.id, "v": [round(ln.v, 3) for ln in lines],
                            "u": [(round(ln.u0, 3), round(fr.L - ln.u1, 3)) for ln in lines]})
    out["bond_spacing_ring"] = stats(ring_spacing)
    out["bond_spacing_nominal"] = round(sp_nominal, 4)
    out["ring_inner_line"] = {
        "inset_unlabeled_ends": stats([e["inset_a"] for e in ring_inner if not e["a_labeled"]]
                                      + [e["inset_b"] for e in ring_inner if not e["b_labeled"]]),
        "inset_formula_s_tan30": round(sp_nominal * math.tan(math.radians(30)), 3),
        "entries": ring_inner,
    }
    out["acyclic_double"] = acyclic
    out["triple"] = triples

    # stereo
    out["stereo"] = [_stereo_entry(doc, b, asg[b.id]) for b in doc.bonds if b.id in asg and "stereo" in asg[b.id]]

    # labels
    out["labels"] = [e for a in doc.atoms.values() if a.label for e in [_label_entry(doc, a, asg)] if e]
    out["charges"] = _charges(doc)
    out["dots"] = _dots(doc)
    out["arrows"] = [_arrow_entry(doc, a) for a in doc.arrows]
    out["curved_arrows"] = [_curve_entry(doc, c) for c in doc.curves]
    out["scheme"] = _r(_scheme(doc, out["arrows"]))
    out["brackets"] = _brackets(doc)
    unclaimed = [s for s in doc.shapes if s.kind == "poly" and not s.owner]
    out["unclaimed_polys"] = [{"bbox": tuple(round(v, 2) for v in s.bbox), "n": len(s.pts),
                               "curved": s.curved} for s in unclaimed]
    return out


def _side(doc, fr, atom_id, exclude):
    """Signed perpendicular sides of an atom's other neighbours in frame fr."""
    return [fr.local(doc.atoms[n].p)[1] for n in neighbors(doc, atom_id) if n != exclude]


def _acyclic_entry(doc, b, fr, lines, spacing):
    A, B = doc.atoms[b.a], doc.atoms[b.b]
    sides_a, sides_b = _side(doc, fr, b.a, b.b), _side(doc, fr, b.b, b.a)
    v = [ln.v for ln in lines]
    if abs(v[0] + v[1]) < 0.15 * spacing:
        placement = "centered"
    else:
        off = max(lines, key=lambda ln: abs(ln.v)).v
        subs = sides_a + sides_b
        toward = sum(1 for s in subs if s * off > 0)
        away = sum(1 for s in subs if s * off < 0)
        placement = f"offset (second line toward {toward} / away from {away} substituents)"
    return {
        "bond": b.id,
        "a": f"{A.label or 'C'}(deg {len(neighbors(doc, b.a))})",
        "b": f"{B.label or 'C'}(deg {len(neighbors(doc, b.b))})",
        "v": [round(x, 3) for x in v], "spacing": round(spacing, 3), "placement": placement,
        "subst_sides_a": [round(s, 2) for s in sides_a], "subst_sides_b": [round(s, 2) for s in sides_b],
        "ends": [(round(ln.u0, 3), round(fr.L - ln.u1, 3)) for ln in lines], "L": round(fr.L, 3),
    }


def _stereo_entry(doc, b, r):
    fr = r["frame"]
    shapes = r["shapes"]
    entry = {"bond": b.id, "display": b.display, "L": round(fr.L, 3),
             "narrow_atom": doc.atoms[fr_id(doc, fr)].label or "C"}
    end_atom = doc.atoms[b.b if b.display.endswith("Begin") else b.a]
    entry["wide_atom"] = end_atom.label or "C"
    if "Hash" in b.display:
        hashes = []
        for sh, loc in shapes:
            us, vs = [p[0] for p in loc], [p[1] for p in loc]
            hashes.append({"u": round(statistics.fmean(us), 3),
                           "thickness_u": round(max(us) - min(us), 3),
                           "len_v": round(max(vs) - min(vs), 3)})
        hashes.sort(key=lambda h: h["u"])
        gaps = [h2["u"] - h1["u"] for h1, h2 in zip(hashes, hashes[1:])]
        # each hash line is a parallelogram: its true thickness is measured
        # perpendicular to its long edge, not along u
        thick = []
        for sh, _ in shapes:
            e = max(((sh.pts[i], sh.pts[(i + 1) % len(sh.pts)]) for i in range(len(sh.pts))),
                    key=lambda ab: length(sub(ab[1], ab[0])))
            d = unit(sub(e[1], e[0]))
            nn = (-d[1], d[0])
            proj = [dot(p, nn) for p in sh.pts]
            thick.append(max(proj) - min(proj))
        # extrapolate hash length to u = L (linear in u)
        if len(hashes) >= 2:
            k = (hashes[-1]["len_v"] - hashes[0]["len_v"]) / (hashes[-1]["u"] - hashes[0]["u"])
            at_L = hashes[-1]["len_v"] + k * (fr.L - hashes[-1]["u"])
            at_0 = hashes[0]["len_v"] - k * hashes[0]["u"]
        else:
            at_L = at_0 = None
        entry.update({"count": len(hashes), "hashes": hashes, "pitch": stats(gaps),
                      "line_thickness": stats(thick),
                      "len_extrapolated_at_narrow_atom": round(at_0, 3) if at_0 is not None else None,
                      "len_extrapolated_at_wide_atom": round(at_L, 3) if at_L is not None else None})
    else:
        pts = [p for sh, loc in shapes for p in loc]
        umax, umin = max(p[0] for p in pts), min(p[0] for p in pts)
        wide = [p[1] for p in pts if p[0] > umax - 0.3]
        narrow = [p[1] for p in pts if p[0] < umin + 0.3]
        w_end = max(wide) - min(wide)
        entry.update({"u_range": (round(umin, 3), round(umax, 3)),
                      "wide_end_width": round(w_end, 3),
                      "narrow_end_width": round(max(narrow) - min(narrow), 3),
                      "width_extrapolated_at_wide_atom": round(w_end * fr.L / umax, 3) if umax else None,
                      "n_vertices": len(pts)})
    return entry


def fr_id(doc, fr):
    for aid, a in doc.atoms.items():
        if length(sub(a.p, fr.o)) < 1e-6:
            return aid
    return None


def _label_entry(doc, atom, asg):
    elem, whole = label_glyphs(doc, atom)
    if not elem:
        return {"atom": atom.id, "label": atom.label, "error": "element glyph not found"}
    face = load_face("Regular")
    ebox = union(g.box for g in elem)
    eadv_c = elem[0].x + sum(g.adv for g in elem) / 2
    entry = {
        "atom": atom.id, "label": atom.label, "element": atom.element, "nickname": atom.nickname,
        "glyph_run": "".join(g.ch for g in whole),
        "anchor_ink_center_minus_atom": tuple(round(v, 3) for v in sub(_box_center(ebox), atom.p)),
        "anchor_advance_center_x_minus_atom": round(eadv_c - atom.p[0], 3),
        "baseline_minus_atom_y": round(elem[0].base - atom.p[1], 3),
        "cap_height": round(face.cap_height(elem[0].size), 3),
        "size": elem[0].size,
    }
    # H, subscripts, charges relative to the element
    others = []
    for g in whole:
        if g in elem:
            continue
        others.append({"ch": g.ch, "size": round(g.size, 3),
                       "dx_origin_from_elem_origin": round(g.x - elem[0].x, 3),
                       "dbaseline": round(g.base - elem[0].base, 3)})
    entry["other_glyphs"] = others
    # bonds into this label
    bonds = []
    outlines = [c for g in elem for c in g.outline()]
    hull = convex_hull([p for c in outlines for p in c])
    margin = float(doc.settings["MarginWidth"])
    for b in doc.bonds:
        if atom.id not in (b.a, b.b) or b.id not in asg:
            continue
        r = asg[b.id]
        other = b.b if b.a == atom.id else b.a
        fr = frame(doc, other, atom.id)  # u runs toward the label
        lines = r.get("lines") or []
        if "stereo" in r:
            continue
        for ln in lines:
            ln2 = _line_from_shape(ln.shape, fr)
            end = fr.world(ln2.u1, ln2.v)
            c1, c2 = fr.world(ln2.u1, ln2.v - ln2.width / 2), fr.world(ln2.u1, ln2.v + ln2.width / 2)
            gap_box = ray_box_gap(end, fr.u, ebox)
            gap_outline = min(seg_seg_dist(c1, c2, p, q) for c in outlines
                              for p, q in zip(c, c[1:] + c[:1]))
            gap_hull = ray_polygon(fr.world(ln2.u1, 0.0), fr.u, hull)
            bonds.append({"bond": b.id, "order": b.order, "line_v": round(ln2.v, 3),
                          "dir_deg": round(math.degrees(math.atan2(-fr.u[1], fr.u[0])), 1),
                          "end_dist_from_atom": round(fr.L - ln2.u1, 3),
                          "gap_to_elem_inkbox_along_bond": round(gap_box, 3) if gap_box is not None else None,
                          "gap_to_elem_outline": round(gap_outline, 3),
                          "axial_gap_to_elem_hull": round(gap_hull, 3) if gap_hull is not None else None,
                          "axial_gap_to_hull_over_margin": round(gap_hull / margin, 3) if gap_hull is not None else None,
                          "other_end_labeled": bool(doc.atoms[other].label)})
    entry["bonds"] = bonds
    return entry


def _line_scale(doc) -> float:
    """How far a search window fitted to LineWidth 0.6 (every ACS scale) must grow on this file's thicker lines."""
    return max(1.0, float(doc.settings.get("LineWidth") or 0) / 0.6)


def _small_bars(doc):
    # a sign's bars. 6 x 1 pt holds to LabelSize 10 and LineWidth 0.6 (every ACS scale). Past that the limits grow
    # with the style; 0.9 x LineWidth keeps out hash lines and the like, which are as thick as a line.
    line_width = float(doc.settings.get("LineWidth") or 0)
    label_size = float(doc.settings.get("LabelSize") or 0)
    max_length, max_thickness = max(6.0, 0.6 * label_size), max(1.0, 0.9 * line_width)
    bars = []
    for s in doc.shapes:
        if s.kind != "poly" or s.owner or s.curved or len(s.pts) != 4 or s.group.endswith("-bonds"):
            continue  # ours: a crossing's oblique cut end (draw_mol.cut_caps) is a short polygon among the bonds
        x0, y0, x1, y1 = s.bbox
        w, h = x1 - x0, y1 - y0
        if 1.0 < max(w, h) < max_length and min(w, h) < max_thickness:
            bars.append((s, "h" if w > h else "v"))
    return bars


def _charges(doc):
    bars = _small_bars(doc)
    used, out = set(), []
    for i, (s, o) in enumerate(bars):
        if i in used:
            continue
        mate = None
        if o == "h":
            for j, (t, o2) in enumerate(bars):
                if o2 == "v" and j not in used and length(sub(s.centroid, t.centroid)) < 0.3:
                    mate = j
        if o == "v":
            continue
        used.add(i)
        sign = "+" if mate is not None else "-"
        s.owner = "charge"
        if mate is not None:
            used.add(mate)
            bars[mate][0].owner = "charge"
        c = s.centroid
        x0, y0, x1, y1 = s.bbox
        atom = min(doc.atoms.values(), key=lambda a: length(sub(a.p, c)))
        rec = {"sign": sign, "center": (round(c[0], 3), round(c[1], 3)),
               "bar_length": round(x1 - x0, 3), "bar_thickness": round(y1 - y0, 3),
               "atom_label": atom.label or "C", "center_minus_atom": tuple(round(v, 3) for v in sub(c, atom.p))}
        if atom.label:
            elem, whole = label_glyphs(doc, atom)
            if elem:
                eb = union(g.box for g in elem)
                lb = union(g.box for g in whole)
                rec["center_minus_elem_inkbox_topright"] = (round(c[0] - eb[2], 3), round(c[1] - eb[1], 3))
                rec["center_minus_elem_inkbox_center"] = tuple(round(v, 3) for v in sub(c, _box_center(eb)))
                rec["label_run"] = "".join(g.ch for g in whole)
                rec["label_inkbox"] = tuple(round(v, 2) for v in lb)
        else:
            nb = neighbors(doc, atom.id)
            vecs = [unit(sub(doc.atoms[n].p, atom.p)) for n in nb]
            rec["bond_dirs_deg"] = [round(math.degrees(math.atan2(-v[1], v[0])), 1) for v in vecs]
            d = sub(c, atom.p)
            rec["charge_dir_deg"] = round(math.degrees(math.atan2(-d[1], d[0])), 1)
            rec["charge_dist"] = round(length(d), 3)
        out.append(rec)
    return out


def _dots(doc):
    circles = [s for s in doc.shapes if s.kind == "circle"]
    if not circles:
        return {}
    pairs, singles, used = [], [], set()
    for i, a in enumerate(circles):
        if i in used:
            continue
        mate = None
        for j, b in enumerate(circles):
            if j != i and j not in used and length(sub(a.center, b.center)) < 4.5:
                mate = j
                break
        atom = min(doc.atoms.values(), key=lambda at: length(sub(at.p, a.center)))
        if mate is None:
            used.add(i)
            singles.append({"center_minus_atom": tuple(round(v, 3) for v in sub(a.center, atom.p)),
                            "atom": atom.label or "C", "dist": round(length(sub(a.center, atom.p)), 3)})
            continue
        used |= {i, mate}
        b = circles[mate]
        mid = mul(add(a.center, b.center), 0.5)
        atom = min(doc.atoms.values(), key=lambda at: length(sub(at.p, mid)))
        rec = {"atom": atom.label or "C", "dot_spacing": round(length(sub(a.center, b.center)), 3),
               "mid_minus_atom": tuple(round(v, 3) for v in sub(mid, atom.p)),
               "mid_dist_from_atom": round(length(sub(mid, atom.p)), 3)}
        elem, _ = label_glyphs(doc, atom) if atom.label else ([], [])
        if elem:
            eb = union(g.box for g in elem)
            # gap from the nearer dot edge to the element ink box, measured on the side the pair sits
            d = sub(mid, _box_center(eb))
            if abs(d[0]) * (eb[3] - eb[1]) > abs(d[1]) * (eb[2] - eb[0]):  # left/right side
                gap = (min(a.center[0], b.center[0]) - eb[2]) if d[0] > 0 else (eb[0] - max(a.center[0], b.center[0]))
                side = "right" if d[0] > 0 else "left"
            else:
                gap = (min(a.center[1], b.center[1]) - eb[3]) if d[1] > 0 else (eb[1] - max(a.center[1], b.center[1]))
                side = "below" if d[1] > 0 else "above"
            rec.update({"side": side, "gap_dot_edge_to_elem_inkbox": round(gap - a.r, 3),
                        "mid_minus_elem_inkbox_center": tuple(round(v, 3) for v in d)})
        pairs.append(rec)
    return {"radius": stats([c.r for c in circles]), "pairs": pairs, "singles": singles}


def _head_geometry(sh, tip, axis):
    """Arrowhead polygon geometry relative to its tip and axis direction (unit,
    pointing toward the tip)."""
    n = (-axis[1], axis[0])
    loc = [(dot(sub(p, tip), axis), dot(sub(p, tip), n)) for p in sh.pts]
    back = min(u for u, _ in loc)
    vs = [v for _, v in loc]
    # inset point: control/vertex on the axis (|v| small) farthest from the tip, excluding the tip
    on_axis = [u for u, v in loc if abs(v) < 0.05 and u < -0.5]
    return {"length": round(-back, 3), "half_width_left": round(max(vs), 3),
            "half_width_right": round(-min(vs), 3),
            "inset_from_tip": round(-min(on_axis), 3) if on_axis else None}


def _arrow_entry(doc, arr):
    if arr.get("equilibrium"):
        return _equilibrium_entry(doc, arr)
    head, tail = arr["head"], arr["tail"]
    axis = unit(sub(head, tail))
    n = (-axis[1], axis[0])
    rec = {"tail": tuple(round(v, 2) for v in tail), "head": tuple(round(v, 2) for v in head),
           "length": round(length(sub(head, tail)), 3), "cdxml": arr["attrs"]}
    heads = []
    for sh in doc.shapes:
        if sh.kind == "poly" and sh.curved and not sh.owner:
            for tip, ax, nm in ((head, axis, "head"), (tail, mul(axis, -1), "tail")):
                if any(length(sub(p, tip)) < 0.05 for p in sh.pts):
                    sh.owner = f"arrow-{nm}"
                    heads.append((nm, _head_geometry(sh, tip, ax)))
    rec["heads"] = dict(heads)
    outlines = [s.width for s in doc.shapes if s.kind == "outline"]
    rec["head_outline_stroke"] = stats(outlines)
    # shaft: an unclaimed straight polygon along the axis
    scale = _line_scale(doc)
    for sh in doc.shapes:
        if sh.kind != "poly" or sh.curved or sh.owner:
            continue
        loc = [(dot(sub(p, tail), axis), dot(sub(p, tail), n)) for p in sh.pts]
        us, vs = [p[0] for p in loc], [p[1] for p in loc]
        if max(abs(v) for v in vs) < 0.6 * scale and min(us) > -1 and max(us) < rec["length"] + 1 and max(us) - min(us) > 3:
            sh.owner = "arrow-shaft"
            rec["shaft"] = {"width": round(max(vs) - min(vs), 3), "u_start": round(min(us), 3),
                            "u_end_from_tip": round(rec["length"] - max(us), 3)}
    return rec


def _half_head_geometry(lp, tip, d, vc, o):
    """Half arrowhead (one barb) in its shaft's frame: `lp` the polygon's (u, v) in the arrow frame,
    `tip` its (u, v), d = +1 / -1 the shaft's direction along u, vc the shaft center's v and
    o = +1 / -1 the side away from the other shaft. back = distance behind the tip, side =
    offset from the shaft's center line, positive away from the other shaft (where the barb is)."""
    loc = [(d * (tip[0] - u), o * (v - vc)) for u, v in lp]
    tip_side = o * (tip[1] - vc)
    flat = [b for b, s in loc if abs(s - tip_side) < 0.05]
    barb = max(s for _, s in loc)
    inset = max(flat)
    ctl = [s for b, s in loc if abs(b - inset) < 0.05 and tip_side + 0.05 < s < barb - 0.05]
    return {"length": round(max(b for b, _ in loc), 3), "barb_from_center": round(barb, 3),
            "tip_from_center": round(tip_side, 3), "inset_from_tip": round(inset, 3),
            "back_control_from_center": round(ctl[0], 3) if ctl else None,
            "back_control_per_barb": round(ctl[0] / barb, 4) if ctl else None}


def _equilibrium_entry(doc, arr):
    """Equilibrium arrow (style-spec.md §3.6): two straight shafts parallel to Tail3D -> Head3D,
    each with a filled half head (curved polygon) at one end. Per shaft: center offset v from the
    CDXML axis, width, tip and butt along the axis (u from Tail3D), direction (forward: tip
    toward Head3D) and its half head (_half_head_geometry). Per arrow: shaft spacing (center to
    center), balanced (equal shafts: the stagger of each butt behind the other shaft's tip) or
    unbalanced (the short shaft's center against the long one's, and long / short butt-to-notch
    ratio), and the ink extent along the axis."""
    head, tail = arr["head"], arr["tail"]
    axis = unit(sub(head, tail))
    n = (-axis[1], axis[0])
    L = length(sub(head, tail))
    scale = _line_scale(doc)
    shafts, heads = [], []
    for sh in doc.shapes:
        if sh.kind != "poly" or sh.owner:
            continue
        lp = [(dot(sub(p, tail), axis), dot(sub(p, tail), n)) for p in sh.pts]
        us, vs = [q[0] for q in lp], [q[1] for q in lp]
        if max(abs(v) for v in vs) > 4 * scale or min(us) < -3 * scale or max(us) > L + 3 * scale:
            continue
        if sh.curved:
            heads.append((sh, lp))
        elif max(vs) - min(vs) < scale and max(us) - min(us) > 3:
            shafts.append((sh, lp))
    recs = []
    for sh, lp in shafts:
        us, vs = [q[0] for q in lp], [q[1] for q in lp]
        vc, u0, u1 = (max(vs) + min(vs)) / 2, min(us), max(us)
        o = 1 if vc > 0 else -1
        for hs, hp in heads:
            hu = [q[0] for q in hp]
            if hs.owner or o * statistics.fmean(q[1] for q in hp) < 0:
                continue
            if min(hu) < u1 < max(hu) + 0.5 * scale:
                d, tip = 1, max(hp, key=lambda q: q[0])
            elif min(hu) - 0.5 * scale < u0 < max(hu):
                d, tip = -1, min(hp, key=lambda q: q[0])
            else:
                continue
            sh.owner, hs.owner = "arrow-shaft", "arrow-head"
            butt = u0 if d == 1 else u1
            hg = _half_head_geometry(hp, tip, d, vc, o)
            recs.append({"v": round(vc, 3), "width": round(max(vs) - min(vs), 3),
                         "direction": "forward" if d == 1 else "back",
                         "tip_u": round(tip[0], 3), "butt_u": round(butt, 3),
                         "length": round(abs(tip[0] - butt), 3),
                         "u_end_from_tip": round(d * (tip[0] - (u1 if d == 1 else u0)), 3),
                         "head": hg, "barb": "outer" if hg["barb_from_center"] > 0 else "inner",
                         "y": add(tail, mul(n, vc))[1]})
            break
    recs.sort(key=lambda r: r["v"])
    rec = {"tail": tuple(round(v, 2) for v in tail), "head": tuple(round(v, 2) for v in head),
           "length": round(L, 3), "cdxml": arr["attrs"], "equilibrium": {"shafts": recs}, "heads": {}}
    if len(recs) != 2:
        rec["equilibrium"]["error"] = f"{len(recs)} shafts with heads"
        return rec
    e = rec["equilibrium"]
    a, b = recs
    lo = min(min(r["tip_u"], r["butt_u"]) for r in recs)
    hi = max(max(r["tip_u"], r["butt_u"]) for r in recs)
    rec["ends"] = (add(tail, mul(axis, lo)), add(tail, mul(axis, hi)))
    rec["outer_y"] = tuple(sorted(r["y"] for r in recs))
    e.update({"spacing": round(b["v"] - a["v"], 3), "center_v": round((a["v"] + b["v"]) / 2, 3),
              "extent_u": (round(lo, 3), round(hi, 3))})
    fwd = next((r for r in recs if r["direction"] == "forward"), None)
    back = next((r for r in recs if r["direction"] == "back"), None)
    if abs(a["length"] - b["length"]) < 1.0:
        e["kind"] = "balanced"
        if fwd and back:
            e["stagger_tail"] = round(fwd["butt_u"] - back["tip_u"], 3)  # forward butt behind the back tip
            e["stagger_head"] = round(fwd["tip_u"] - back["butt_u"], 3)
    else:
        long, short = (a, b) if a["length"] > b["length"] else (b, a)
        e["kind"] = "unbalanced"
        e["long"] = long["direction"]
        e["short_length"], e["long_length"] = short["length"], long["length"]
        e["short_center_minus_long_center"] = round((short["tip_u"] + short["butt_u"]) / 2
                                                    - (long["tip_u"] + long["butt_u"]) / 2, 3)
        e["long_over_short_butt_to_notch"] = round((long["length"] - long["head"]["inset_from_tip"])
                                                   / (short["length"] - short["head"]["inset_from_tip"]), 4)
    return rec


def _curve_entry(doc, c):
    pts = c["pts"]
    stroke = None
    end = pts[-2] if len(pts) % 3 == 0 else pts[-1]
    # ChemDraw CurvePoints: (in-handle, anchor, out-handle) per anchor
    anchors = pts[1::3]
    best = None
    for sh in doc.shapes:
        if sh.kind == "stroke" and sh.curved:
            d = min(length(sub(sh.pts[0], anchors[0])), length(sub(sh.pts[-1], anchors[0])))
            if best is None or d < best[0]:
                best = (d, sh)
    rec = {"head_type": c["head"], "curve_type": c["type"],
           "anchors": [tuple(round(v, 2) for v in a) for a in anchors],
           "handles": [tuple(round(v, 2) for v in p) for i, p in enumerate(pts) if i % 3 != 1]}
    if best:
        stroke = best[1]
        rec["stroke_width"] = round(stroke.width, 3)
        rec["svg_path_pts"] = [tuple(round(v, 2) for v in p) for p in stroke.pts]
        p_end = stroke.pts[3] if len(stroke.pts) >= 4 else stroke.pts[-1]
        # end tangent from the last control point that isn't the end (ChemDraw's C–Br arrow has c2 = p3)
        c_last = next(p for p in reversed(stroke.pts[:3]) if length(sub(p_end, p)) > 1e-6)
        tangent = unit(sub(p_end, c_last))
        # the curve's last anchor is the stroke end (the head's notch), not the tip (§13.2)
        rec["last_anchor_minus_stroke_end"] = tuple(round(v, 3) for v in sub(anchors[-1], p_end))
        # head polygon: unclaimed straight poly with a vertex near the stroke end
        heads = [s for s in doc.shapes if s.kind == "poly" and not s.owner and not s.curved
                 and min(length(sub(p, p_end)) for p in s.pts) < 3.0 and len(s.pts) in (3, 4)]
        if heads:
            h = min(heads, key=lambda s: length(sub(s.centroid, p_end)))
            h.owner = "curve-head"
            tip = max(h.pts, key=lambda p: dot(sub(p, p_end), tangent))
            ax = tangent
            n = (-ax[1], ax[0])
            loc = [(dot(sub(p, tip), ax), dot(sub(p, tip), n)) for p in h.pts]
            pe = (dot(sub(p_end, tip), ax), dot(sub(p_end, tip), n))
            rec["head"] = {"n_vertices": len(h.pts), "length": round(-min(u for u, _ in loc), 3),
                           # stroke end in the head frame: fishhooks end beside the head's axis, toward the barb
                           "stroke_end_local": (round(pe[0], 3), round(pe[1], 3)),
                           "extent_left": round(max(v for _, v in loc), 3),
                           "extent_right": round(-min(v for _, v in loc), 3),
                           "tip_minus_stroke_end": tuple(round(v, 3) for v in sub(tip, p_end)),
                           "pts_local": [(round(u, 3), round(v, 3)) for u, v in loc]}
        # chord / bow (approximate: hand-shaped)
        p0, p3 = stroke.pts[0], p_end
        ch = sub(p3, p0)
        L = length(ch)
        nn = (-ch[1] / L, ch[0] / L)
        bow = max(abs(dot(sub(p, p0), nn)) for p in stroke.pts[:4])
        rec["chord"] = round(L, 3)
        rec["control_offset_max_over_chord"] = round(bow / L, 3)
    return rec


MULTISTEP = re.compile(r"^\d+\.")  # "1. n-BuLi, ..." : a numbered stage, set left-justified


def _caption_glyphs(doc, cap):
    """Glyphs of one caption: given directly (our renders) or found inside the CDXML text box."""
    if "glyphs" in cap:
        return [g for g in cap["glyphs"] if g.box]
    bb = cap["bbox"]
    return [g for g in doc.glyphs if g.box and bb[0] - 1 <= _box_center(g.box)[0] <= bb[2] + 1
            and bb[1] - 1 <= _box_center(g.box)[1] <= bb[3] + 1]


def text_lines(glyphs) -> list[dict]:
    """Group a caption's glyphs into lines, top to bottom. Full-size glyphs set the baselines;
    smaller ones (sub- and superscripts) join the nearest line. Advances: x0 is the first
    glyph origin, x1 the last glyph's origin + advance."""
    main = max(g.size for g in glyphs)
    bases = []
    for b in sorted(g.base for g in glyphs if abs(g.size - main) < 0.01):
        if not bases or b - bases[-1] > 0.05:
            bases.append(b)
    lines = [{"base": b, "glyphs": []} for b in bases]
    for g in glyphs:
        min(lines, key=lambda ln: abs(ln["base"] - g.base))["glyphs"].append(g)
    out = []
    for ln in lines:
        gl = sorted(ln["glyphs"], key=lambda g: g.x)
        base = ln["base"]
        small = [g for g in gl if g.size < main - 0.01]
        subs = [g for g in small if g.base > base + 0.01]
        sups = [g for g in small if g.base < base - 0.01]
        rec = {"text": "".join(g.ch for g in gl), "base": base, "size": main,
               "x0": gl[0].x, "x1": gl[-1].x + gl[-1].adv, "ink": union(g.box for g in gl),
               "has_sub": bool(subs), "has_sup": bool(sups), "bold": any(g.bold for g in gl),
               "italic": "".join(g.ch for g in gl if g.italic)}
        if subs:
            rec["sub"] = {"size_ratio": subs[0].size / main, "drop_over_size": (subs[0].base - base) / main}
        if sups:
            rec["sup"] = {"size_ratio": sups[0].size / main, "raise_over_size": (base - sups[0].base) / main}
        out.append(rec)
    return out


def _rings(doc, atom_ids) -> list[list[str]]:
    """Smallest cycle through each ring bond, as atom-id lists (no duplicates)."""
    adj = {}
    for b in doc.bonds:
        if b.a in atom_ids:
            adj.setdefault(b.a, set()).add(b.b)
            adj.setdefault(b.b, set()).add(b.a)
    rings = {}
    for b in doc.bonds:
        if not b.ring or b.a not in atom_ids:
            continue
        prev, queue = {b.a: None}, [b.a]
        while queue and b.b not in prev:
            u = queue.pop(0)
            for w in adj[u]:
                if w not in prev and {u, w} != {b.a, b.b}:
                    prev[w] = u
                    queue.append(w)
        if b.b in prev:
            path, u = [], b.b
            while u is not None:
                path.append(u)
                u = prev[u]
            rings.setdefault(frozenset(path), path)
    return list(rings.values())


def _structures(doc) -> dict:
    """fragment id -> {ink, rings} for fragments drawn with atoms (formula text excluded)."""
    out = {}
    for fid in dict.fromkeys(a.frag for a in doc.atoms.values()):
        if fid in doc.text_frags:
            continue
        ids = {a.id for a in doc.atoms.values() if a.frag == fid}
        bond_ids = {b.id for b in doc.bonds if b.frag == fid}
        boxes = [s.bbox for s in doc.shapes if s.owner in bond_ids]
        for aid in ids:
            if doc.atoms[aid].label:
                boxes += [g.box for g in label_glyphs(doc, doc.atoms[aid])[1]]
        if not boxes:
            continue
        centers = [(statistics.fmean(doc.atoms[i].p[0] for i in r), statistics.fmean(doc.atoms[i].p[1] for i in r))
                   for r in _rings(doc, ids)]
        out[fid] = {"ink": union(boxes), "ring_centers": centers}
    return out


def _rows(structs: dict) -> list[set]:
    """Structures grouped into rows by overlapping vertical ink extents."""
    rows: list[set] = []
    for fid, s in sorted(structs.items(), key=lambda kv: kv[1]["ink"][1]):
        y0, y1 = s["ink"][1], s["ink"][3]
        hit = [r for r in rows if any(structs[o]["ink"][1] < y1 and y0 < structs[o]["ink"][3] for o in r)]
        merged = set().union(*hit) | {fid} if hit else {fid}
        rows = [r for r in rows if r not in hit] + [merged]
    return rows


def _scheme(doc, arrows: list[dict]):
    """Scheme layout (style-spec.md §3): each caption is assigned to an arrow (horizontal overlap,
    then nearest in y), a plus sign or a compound number; structures are grouped in rows.
    `arrows` are the _arrow_entry records (they carry the measured head inset)."""
    structs = _structures(doc)
    rows = _rows(structs)

    def row_of(y):
        return next((r for r in rows if any(structs[f]["ink"][1] <= y <= structs[f]["ink"][3] for f in r)), set())

    arrow_recs = []
    for a in arrows:
        tail, head = a["tail"], a["head"]
        inset = (a["heads"].get("head") or {}).get("inset_from_tip") or 0.0
        if "ends" in a:  # equilibrium: its ink extent, left to right, centered (no notch)
            tail, head = sorted(a["ends"])
        sgn = 1 if head[0] >= tail[0] else -1
        row = row_of(tail[1])
        left = [structs[f]["ink"] for f in row if structs[f]["ink"][2] <= min(tail[0], head[0])]
        right = [structs[f]["ink"] for f in row if structs[f]["ink"][0] >= max(tail[0], head[0])]
        arrow_recs.append({
            "tail": tail, "head": head, "axis_y": tail[1], "length": a["length"], "notch": head[0] - sgn * inset,
            "mid": (tail[0] + head[0]) / 2 if "ends" in a else None,
            # conditions are measured from the outer shafts of an equilibrium arrow
            "axis_above": a["outer_y"][0] if "outer_y" in a else tail[1],
            "axis_below": a["outer_y"][1] if "outer_y" in a else tail[1],
            "row": row, "above": [], "below": [],
            "gap_left_struct_to_tail": min(tail[0], head[0]) - max(b[2] for b in left) if left else None,
            "gap_head_to_right_struct": min(b[0] for b in right) - max(tail[0], head[0]) if right else None,
            "structures": [{"ink_center_minus_axis": (structs[f]["ink"][1] + structs[f]["ink"][3]) / 2 - tail[1],
                            "ring_centers_minus_axis": sorted((c[1] - tail[1] for c in structs[f]["ring_centers"]),
                                                              key=abs)}
                           for f in sorted(row, key=lambda f: structs[f]["ink"][0])],
        })
    pluses, numbers, other = [], [], []
    for cap in doc.captions:
        glyphs = _caption_glyphs(doc, cap)
        if not glyphs:
            continue
        lines = text_lines(glyphs)
        ink = union(ln["ink"] for ln in lines)
        cx, cy = _box_center(ink)
        if cap["text"].strip() == "+":
            row = row_of(cy)
            left = [structs[f]["ink"] for f in row if structs[f]["ink"][2] < ink[0]]
            right = [structs[f]["ink"] for f in row if structs[f]["ink"][0] > ink[2]]
            axis = next((r["axis_y"] for r in arrow_recs if r["row"] == row), None)
            pluses.append({"gap_left_struct": ink[0] - max(b[2] for b in left) if left else None,
                           "gap_right_struct": min(b[0] for b in right) - ink[2] if right else None,
                           "ink_center_minus_axis": cy - axis if axis is not None else None})
        elif all(ln["bold"] for ln in lines):
            above = [f for f in structs if structs[f]["ink"][3] < ink[1] and structs[f]["ink"][0] < cx < structs[f]["ink"][2]]
            if not above:
                continue
            fid = max(above, key=lambda f: structs[f]["ink"][3])
            sb = structs[fid]["ink"]
            row = next(r for r in rows if fid in r)
            numbers.append({"text": cap["text"], "row": row, "base": lines[0]["base"], "ink": ink,
                            "advance": lines[0]["x1"] - lines[0]["x0"],
                            "gap_own_struct_to_ink_top": ink[1] - sb[3],
                            "gap_row_bottom_to_ink_top": ink[1] - max(structs[f]["ink"][3] for f in row),
                            "center_minus_struct_center_x": cx - (sb[0] + sb[2]) / 2})
        else:
            cands = [r for r in arrow_recs if ink[0] < max(r["tail"][0], r["head"][0])
                     and ink[2] > min(r["tail"][0], r["head"][0])]
            if not cands:
                other.append(cap["text"])
                continue
            r = min(cands, key=lambda r: abs(cy - r["axis_y"]))
            for ln in lines:
                ln["caption"] = cap["text"]
                ln["formula"] = cap.get("formula", False)
                ln["left_justified"] = bool(MULTISTEP.match(lines[0]["text"]))
            (r["above"] if cy < r["axis_y"] else r["below"]).append(lines)
    for r in arrow_recs:
        r["above"] = [ln for lines in sorted(r["above"], key=lambda ls: ls[0]["base"]) for ln in lines]
        r["below"] = [ln for lines in sorted(r["below"], key=lambda ls: ls[0]["base"]) for ln in lines]
        tail_x, notch = r["tail"][0], r["notch"]
        mid = r.pop("mid") or (tail_x + notch) / 2
        for side in ("above", "below"):
            for ln in r[side]:
                ln["center_minus_tail_notch_mid"] = (ln["x0"] + ln["x1"]) / 2 - mid
                ln["ink_center_minus_tail_notch_mid"] = (ln["ink"][0] + ln["ink"][2]) / 2 - mid
                ln["x0_minus_tail"] = ln["x0"] - tail_x
                ln["ink_left_minus_tail"] = ln["ink"][0] - tail_x
                ln["notch_minus_x1"] = notch - ln["x1"]
        if r["above"]:
            last = r["above"][-1]
            r["above_axis_minus_last_base"] = r["axis_above"] - last["base"]
            r["above_axis_minus_ink_bottom"] = r["axis_above"] - max(ln["ink"][3] for ln in r["above"])
        if r["below"]:
            first = r["below"][0]
            r["below_first_base_minus_axis"] = first["base"] - r["axis_below"]
            r["below_ink_top_minus_axis"] = min(ln["ink"][1] for ln in r["below"]) - r["axis_below"]
        lines = r["above"] + r["below"]
        if lines:
            w = max(lines, key=lambda ln: ln["x1"] - ln["x0"])
            r["widest"] = {"text": w["text"], "advance": w["x1"] - w["x0"], "ink_width": w["ink"][2] - w["ink"][0],
                           "pad_tail": w["x0"] - tail_x, "pad_notch": notch - w["x1"],
                           "left_justified": w["left_justified"]}
        r["pitches"] = [{"pitch": b["base"] - a["base"], "size": a["size"], "sub_above": a["has_sub"],
                         "sub_below": b["has_sub"], "texts": (a["text"], b["text"])}
                        for side in ("above", "below") for a, b in zip(r[side], r[side][1:])
                        if a["caption"] == b["caption"]]
    by_row = {}
    for n in numbers:
        by_row.setdefault(frozenset(n["row"]), []).append(n)
    number_rows = [{"texts": [n["text"] for n in ns], "bases": [n["base"] for n in ns],
                    "gap_row_bottom_to_ink_top": [n["gap_row_bottom_to_ink_top"] for n in ns],
                    "center_minus_struct_center_x": [n["center_minus_struct_center_x"] for n in ns],
                    "numbers": ns}
                   for ns in by_row.values()]
    for r in arrow_recs:
        r["row"] = sorted(r["row"])
    for n in numbers:
        n["row"] = sorted(n["row"])
    return {"arrows": arrow_recs, "pluses": pluses, "number_rows": number_rows, "unassigned": other}


def _r(v, nd=3):
    """Round floats (recursively) for JSON output."""
    if isinstance(v, float):
        return round(v, nd)
    if isinstance(v, (list, tuple)):
        return type(v)(_r(x, nd) for x in v)
    if isinstance(v, dict):
        return {k: _r(x, nd) for k, x in v.items()}
    return v


def _brackets(doc):
    """Square brackets: ChemDraw draws each as three stroked segments."""
    segs = [s for s in doc.shapes if s.kind == "stroke" and len(s.pts) == 2]
    out = []
    for v in segs:
        (x0, y0), (x1, y1) = v.pts
        if abs(x1 - x0) > 0.01 or abs(y1 - y0) < 20:
            continue
        lips = [h for h in segs if h is not v and abs(h.pts[0][1] - h.pts[1][1]) < 0.01
                and min(length(sub(p, q)) for p in h.pts for q in v.pts) < 0.01]
        for sh in [v, *lips]:
            sh.owner = "bracket"
        lip_len = [abs(h.pts[0][0] - h.pts[1][0]) for h in lips]
        lip_dir = "right" if lips and max(p[0] for h in lips for p in h.pts) > x0 + 0.01 else "left"
        out.append({"x": round(x0, 3), "y_top": round(min(y0, y1), 3), "y_bottom": round(max(y0, y1), 3),
                    "height": round(abs(y1 - y0), 3), "stroke_width": round(v.width, 3),
                    "lip_length": stats(lip_len), "lips_point": lip_dir})
    if out:
        content = [s.bbox for s in doc.shapes if s.owner and s.owner != "bracket"]
        content += [g.box for g in doc.glyphs if g.box]
        cb = union(content)
        for b in out:
            b["content_ink"] = tuple(round(v, 2) for v in cb)
            b["extends_above_content"] = round(cb[1] - b["y_top"], 3)
            b["extends_below_content"] = round(b["y_bottom"] - cb[3], 3)
            b["gap_to_content_x"] = round(cb[0] - b["x"], 3) if b["lips_point"] == "right" else round(b["x"] - cb[2], 3)
    return out
