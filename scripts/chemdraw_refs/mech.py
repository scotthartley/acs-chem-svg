"""Mechanism parameters per preset: curved arrows, lone pairs, radicals (acs-style.md §4)."""

from __future__ import annotations

import math

from acs_chem_svg.geometry import add, dot, length, mul, sub, unit

from .core import _box_center, frame, label_glyphs, load, stats, union
from .measurements import measure


MECH_REFS = {"1.0": "mech_basic_x1.0", "0.8": "mech_basic_x0.8", "0.6": "mech_basic_x0.6"}


def _box_gap(p, box):
    """Distance from p to an axis-aligned box (0 inside)."""
    dx = max(box[0] - p[0], 0.0, p[0] - box[2])
    dy = max(box[1] - p[1], 0.0, p[1] - box[3])
    return math.hypot(dx, dy)


def _dot_sets(doc):
    """Lone pairs (mutually nearest circles) and single dots (radicals), each with its atom:
    the nearest atom center to the pair midpoint / dot."""
    circles = [s for s in doc.shapes if s.kind == "circle"]

    def nearest(i):
        return min((j for j in range(len(circles)) if j != i),
                   key=lambda j: length(sub(circles[i].center, circles[j].center)), default=None)
    mate = {i: nearest(i) for i in range(len(circles))}
    spacing = min(length(sub(circles[i].center, circles[j].center)) for i, j in mate.items() if j is not None)
    pairs, singles, seen = [], [], set()
    for i, j in mate.items():
        if i in seen:
            continue
        a = circles[i]
        if j is not None and mate[j] == i and length(sub(a.center, circles[j].center)) < 1.2 * spacing:
            seen |= {i, j}
            pairs.append((a, circles[j]))
        else:
            seen.add(i)
            singles.append(a)
    return pairs, singles


def _elem_box(doc, atom):
    elem, _ = label_glyphs(doc, atom) if atom.label else ([], [])
    return union(g.box for g in elem) if elem else None


def mech_summary(name: str) -> dict:
    """Mechanism parameters from one mech_basic file (acs-style.md §4): exact dot and
    curved-arrow geometry, and the hand-placed positions (lone-pair gaps, radical distances,
    curved-arrow tails, heads and shape)."""
    doc = load(name)
    m = measure(doc)
    pairs, singles = _dot_sets(doc)
    out = {"dot_radius": stats([c.r for p in pairs for c in p] + [c.r for c in singles]),
           "dot_spacing": stats([length(sub(a.center, b.center)) for a, b in pairs])}
    lp = []
    for a, b in pairs:
        mid = mul(add(a.center, b.center), 0.5)
        atom = min(doc.atoms.values(), key=lambda at: length(sub(at.p, mid)))
        box = _elem_box(doc, atom)
        if box is None:
            continue
        gap = min(_box_gap(a.center, box), _box_gap(b.center, box)) - a.r
        d = sub(mid, _box_center(box))
        side = (("right" if d[0] > 0 else "left") if abs(d[0]) * (box[3] - box[1]) > abs(d[1]) * (box[2] - box[0])
                else ("below" if d[1] > 0 else "above"))
        lp.append({"atom": atom.label, "side": side, "gap_dot_edge_to_elem_inkbox": round(gap, 3),
                   "mid_offset_from_elem_center_along_side": round(abs(d[0]) if side in ("left", "right") else abs(d[1]), 3)})
    out["lone_pairs"] = lp
    out["lone_pair_gap"] = stats([e["gap_dot_edge_to_elem_inkbox"] for e in lp])
    rad = []
    for c in singles:
        atom = min(doc.atoms.values(), key=lambda at: length(sub(at.p, c.center)))
        box = _elem_box(doc, atom)
        d = sub(c.center, atom.p)
        rec = {"atom": atom.label or "C", "dir_deg": round(math.degrees(math.atan2(-d[1], d[0])), 1)}
        if box:
            rec["gap_dot_edge_to_elem_inkbox"] = round(_box_gap(c.center, box) - c.r, 3)
        else:
            rec["center_dist_from_atom"] = round(length(d), 3)
        rad.append(rec)
    out["radicals"] = rad
    lp_mids = [mul(add(a.center, b.center), 0.5) for a, b in pairs]
    bond_mids = [mul(add(doc.atoms[b.a].p, doc.atoms[b.b].p), 0.5) for b in doc.bonds]
    arrows = []
    for c in m["curved_arrows"]:
        pts = c["svg_path_pts"]
        p0, c1, c2, p3 = pts[:4]
        tip = add(p3, c["head"]["tip_minus_stroke_end"])
        chord = length(sub(p3, p0))
        u = unit(sub(p3, p0))
        nrm = (-u[1], u[0])
        rec = {"type": "fishhook" if c["head"]["n_vertices"] == 3 else "two_electron",
               "chord": round(chord, 3), "control_offset_over_chord": c["control_offset_max_over_chord"],
               "handle_lengths_over_chord": [round(length(sub(c1, p0)) / chord, 3), round(length(sub(c2, p3)) / chord, 3)],
               "handle_offsets_over_chord": [round(dot(sub(c1, p0), nrm) / chord, 3), round(dot(sub(c2, p0), nrm) / chord, 3)],
               "handle_along_over_chord": [round(dot(sub(c1, p0), u) / chord, 3), round(dot(sub(c2, p0), u) / chord, 3)]}
        near = min(range(len(lp_mids)), key=lambda i: length(sub(lp_mids[i], p0)))
        k = min(range(len(bond_mids)), key=lambda i: length(sub(bond_mids[i], p0)))
        if length(sub(lp_mids[near], p0)) < length(sub(bond_mids[k], p0)):  # from a lone pair, else a bond
            a, b = pairs[near]
            owner = min(doc.atoms.values(), key=lambda at: length(sub(at.p, lp_mids[near])))
            box = _elem_box(doc, owner)
            out_dir = unit(sub(lp_mids[near], _box_center(box) if box else owner.p))
            rec["tail_from_lone_pair_mid"] = round(length(sub(p0, lp_mids[near])), 3)
            rec["tail_gap_to_dot_edge"] = round(min(length(sub(p0, a.center)), length(sub(p0, b.center))) - a.r, 3)
            # along the pair's outward direction, beyond the dots' outer edge
            rec["tail_beyond_dot_edge_outward"] = round(dot(sub(p0, lp_mids[near]), out_dir) - a.r, 3)
        else:
            bd = doc.bonds[k]
            fr = frame(doc, bd.a, bd.b)
            along, off = fr.local(p0)
            rec["tail_from_bond_mid"] = round(length(sub(p0, bond_mids[k])), 3)
            rec["tail_offset_from_bond_axis"] = round(abs(off), 3)
            rec["tail_along_from_bond_mid"] = round(along - fr.L / 2, 3)
        atom = min(doc.atoms.values(), key=lambda at: length(sub(at.p, tip)))
        box = _elem_box(doc, atom)
        rec["head_target"] = atom.label or "C"
        if box:
            rec["head_tip_gap_to_elem_inkbox"] = round(_box_gap(tip, box), 3)
        else:
            rec["head_tip_from_atom"] = round(length(sub(tip, atom.p)), 3)
        arrows.append(rec)
    out["curved_arrows"] = arrows
    ca = m["curved_arrows"]
    out["curved_arrow_stroke"] = stats([c["stroke_width"] for c in ca])
    full = [c["head"] for c in ca if c["head"]["n_vertices"] == 4]
    half = [c["head"] for c in ca if c["head"]["n_vertices"] == 3]
    out["full_head"] = {"length": stats([h["length"] for h in full]),
                        "half_width": stats([(h["extent_left"] + h["extent_right"]) / 2 for h in full]),
                        "inset": stats([-max(u for u, v in h["pts_local"] if abs(v) < 0.05 and u < -0.5) for h in full])}
    out["full_head"]["stroke_end_off_axis"] = stats([h["stroke_end_local"][1] for h in full])
    out["fishhook_head"] = {"stroke_end_off_axis_toward_barb": stats(
                                [h["stroke_end_local"][1] * (1 if h["extent_left"] > h["extent_right"] else -1) for h in half]),
                            "length": stats([h["length"] for h in half]),
                            "barb": stats([max(h["extent_left"], h["extent_right"]) for h in half]),
                            "inset": stats([-max(u for u, v in h["pts_local"] if abs(v) < 0.05 and u < -0.5) for h in half])}
    return out
