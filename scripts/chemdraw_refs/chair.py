"""Chair templates and bond crossings (style-spec.md §14.1, §14.2)."""

from __future__ import annotations

import copy
import math
import statistics
from pathlib import Path

from lxml import etree

from chem_svg_cdxml.geometry import dot, length, seg_point_dist, sub

from .core import Doc, REF_DIR, angle_between, assign, frame, load, neighbors, stats


CHAIR_REFS = ["chair_x1.0", "chair_x0.8", "chair_x0.6"]


def _bond_z(name: str, ref_dir: Path = REF_DIR) -> tuple[dict, set]:
    """CDXML bond id -> Z, and the CrossingBonds pairs ChemDraw stores (frozensets of bond ids), for
    drawn bonds (a nickname's internal bonds cross each other: scheme_multi_x*)."""
    root = etree.parse(str(ref_dir / f"{name}.cdxml")).getroot()
    drawn = [b for b in root.iter("b") if b.getparent().getparent().tag != "n"]  # not nickname internals (OBn)
    z = {b.get("id"): int(b.get("Z")) for b in drawn}
    pairs = {frozenset((b.get("id"), c)) for b in drawn for c in (b.get("CrossingBonds") or "").split() if c in z}
    return z, pairs


def _poly_dist(p, pts) -> float:
    return min(seg_point_dist(p, a, b) for a, b in zip(pts, pts[1:] + pts[:1]))


def crossing_summary(name: str, doc: Doc | None = None, ref_dir: Path = REF_DIR) -> dict:
    """Bond crossings (§14.2). Every CrossingBonds pair: the bond with the higher Z is drawn whole
    (front), the other (rear) is cut. Per side of the front bond's axis, the rear polygon corners
    nearest the axis (the cut ends; unclaimed slivers in the rear's band count, so the stub left
    between a label and a gap shows): their distance to the front axis (`beyond_axis`, the same for
    both corners of an end when the cut is parallel to the front bond) and to the front bond's ink.
    `half_extent` is the front bond's drawn half-width at the crossing: lw/2 (single), the outer
    line's center offset (double), the wedge outline's half-width there (wedge). `doc`: measure
    this drawing instead (compare.py's render of the templates, same bond ids), with front and
    rear as in the reference."""
    if doc:
        doc = copy.deepcopy(doc)  # assign() claims shapes; measure() may have claimed them already
        for sh in doc.shapes:
            sh.owner = None
    else:
        doc = load(name, ref_dir)
    asg = assign(doc)
    z, pairs = _bond_z(name, ref_dir)
    bonds = {b.id: b for b in doc.bonds}
    lw, m = float(doc.settings["LineWidth"]), float(doc.settings["MarginWidth"])
    spacing = float(doc.settings["BondLength"]) * float(doc.settings["BondSpacing"]) / 100
    rows = []
    for pair in sorted(pairs, key=lambda p: sorted(p)):
        b1, b2 = (bonds[i] for i in sorted(pair))
        front, rear = (b1, b2) if z[b1.id] > z[b2.id] else (b2, b1)
        ff, rf = frame(doc, front.a, front.b), frame(doc, rear.a, rear.b)
        angle = math.degrees(math.acos(min(1.0, abs(dot(ff.u, rf.u)))))
        fa = asg.get(front.id, {})
        front_ink = [ln.shape.pts for ln in fa.get("lines", [])] + [sh.pts for sh, _ in fa.get("shapes", [])]
        # the crossing point, in the front frame
        den = ff.u[0] * rf.u[1] - ff.u[1] * rf.u[0]
        d = sub(rf.o, ff.o)
        t = (d[0] * rf.u[1] - d[1] * rf.u[0]) / den
        if fa.get("lines"):
            half = max(max(abs(ln.v) for ln in fa["lines"]), lw / 2)
        else:  # wedge: outline half-width at the crossing, between its two ends
            loc = [q for _, lc in fa["shapes"] for q in lc]
            near = [q for q in loc if q[0] < ff.L / 2]
            far = [q for q in loc if q[0] >= ff.L / 2]
            if near and far:
                (u0, h0), (u1, h1) = (max(near, key=lambda q: abs(q[1])), max(far, key=lambda q: abs(q[1])))
                half = abs(h0) + (abs(h1) - abs(h0)) * (t - u0) / (u1 - u0)
            else:
                half = None
        rear_polys = [ln.shape for ln in asg.get(rear.id, {}).get("lines", [])]
        for sh in doc.shapes:
            if sh.kind == "poly" and not sh.owner:
                loc = [rf.local(q) for q in sh.pts]
                if all(-2 <= u <= rf.L + 2 and abs(v) <= spacing + lw for u, v in loc):
                    rear_polys.append(sh)
        sides = {}
        for sh in rear_polys:
            for q in sh.pts:
                v = ff.local(q)[1]
                sides.setdefault(1 if v > 0 else -1, []).append((abs(v), q))
        ends = {}
        for sd, pts in sides.items():
            near = min(v for v, _ in pts)
            corners = [(v, q) for v, q in pts if v < near + 0.05]
            ends["+" if sd > 0 else "-"] = {
                "beyond_axis": stats([v for v, _ in corners]),
                "beyond_ink": stats([_poly_dist(q, ink) for _, q in corners for ink in [min(
                    front_ink, key=lambda k: _poly_dist(q, k))]]) if front_ink else None,
                "beyond_half_extent_per_margin": round((statistics.median(v for v, _ in corners) - half) / m, 3)
                if half is not None else None,
            }
        rows.append({"front": front.id, "rear": rear.id, "front_z": z[front.id], "rear_z": z[rear.id],
                     "front_order": front.order, "front_display": front.display or "plain",
                     "rear_order": rear.order, "angle_deg": round(angle, 2),
                     "front_half_extent": round(half, 3) if half is not None else None,
                     "rear_pieces": len(rear_polys), "cut": ends})
    return {"margin_width": m, "line_width": lw, "crossings": rows}


def _ring6(doc, frag: str) -> list[str] | None:
    """The six-ring atom ids of a fragment, in ring order starting at the first CDXML atom."""
    ring = [a for a, at in doc.atoms.items() if at.frag == frag
            and sum(1 for b in doc.bonds if b.ring and a in (b.a, b.b)) == 2]
    if len(ring) != 6:
        return None
    order = [ring[0]]
    while len(order) < 6:
        nxt = [n for n in neighbors(doc, order[-1]) if n in ring and n not in order]
        order.append(nxt[0] if len(order) == 1 else nxt[-1])
    return order


def chair_summary(name: str, ref_dir: Path = REF_DIR) -> dict:
    """Chairs (§14.1). Fragments of chair_x* in CDXML order: two bare chairs (ChemDraw's
    Conformers templates: `chair` and `flipped`), axial methyls on all six atoms, equatorial methyls
    on all six, β-D-glucopyranose, an axial OAc, then the crossing cases (crossing_summary).
    Per template: vertex offsets from vertex 1 (pt, y down) and bond lengths. Per vertex of a
    substituted chair (matched to a template by its offsets): the substituent's direction (deg,
    counterclockwise from +x, y up) and length, and for equatorials the angle to the ring bond one
    removed (the bonds parallel to it in a textbook chair). Crossing bonds in the chairs: which is
    in front (higher Z) and the mean screen y of each (larger y = lower on the page)."""
    doc = load(name, ref_dir)
    z, pairs = _bond_z(name, ref_dir)
    bl = float(doc.settings["BondLength"])
    chairs = [(fid, r) for fid in doc.fragments for r in [_ring6(doc, fid)] if r]
    templates = {}
    for key, (fid, ring) in zip(("chair", "flipped"), chairs[:2]):
        o = doc.atoms[ring[0]].p
        templates[key] = {"fragment": fid, "atoms": ring,
                          "vertices": [tuple(round(v, 3) for v in sub(doc.atoms[a].p, o)) for a in ring],
                          "bond_lengths": [round(frame(doc, a, b).L, 3) for a, b in zip(ring, ring[1:] + ring[:1])]}

    def match(ring):
        """(template, vertex offset of ring[0] in it, direction): the template whose vertices the
        ring's positions repeat (within 0.03 pt: two 0.01 pt CDXML roundings) under a translation."""
        pts = [doc.atoms[a].p for a in ring]
        for key, t in templates.items():
            tv = t["vertices"]
            for start in range(6):
                for dr in (1, -1):
                    o = sub(pts[0], tv[start])
                    if all(length(sub(sub(pts[k], o), tv[(start + dr * k) % 6])) < 0.03 for k in range(6)):
                        return key, start, dr
        return None

    def ang(v):
        return round(math.degrees(math.atan2(-v[1], v[0])), 2)

    subs = []
    for fid, ring in chairs[2:]:
        mt = match(ring)
        if not mt:
            subs.append({"fragment": fid, "error": "no template match"})
            continue
        key, start, dr = mt
        verts = []
        for k, a in enumerate(ring):
            vtx = (start + dr * k) % 6
            for s in neighbors(doc, a):
                if s in ring:
                    continue
                v = sub(doc.atoms[s].p, doc.atoms[a].p)
                far = [(ring[(k + 1) % 6], ring[(k + 2) % 6]), (ring[(k - 1) % 6], ring[(k - 2) % 6])]
                par = min(angle_between(v, sub(doc.atoms[q].p, doc.atoms[p].p)) for p, q in far)
                verts.append({"vertex": vtx + 1, "atom": a, "substituent": s, "label": doc.atoms[s].label,
                              "dir_deg": ang(v), "length": round(length(v), 3),
                              "vertical": abs(v[0]) < 0.02, "up": v[1] < 0,
                              "deg_from_ring_bond_one_removed": round(par, 2),
                              "element": doc.atoms[a].element})
        hetero = [{"vertex": (start + dr * k) % 6 + 1, "element": doc.atoms[a].element}
                  for k, a in enumerate(ring) if doc.atoms[a].element != "C"]
        cross = []
        for p in sorted(pairs, key=sorted):
            ids = sorted(p)
            bs = [next(b for b in doc.bonds if b.id == i) for i in ids]
            if not any(b.frag == fid for b in bs):
                continue
            front = max(bs, key=lambda b: z[b.id])
            cross.append({"front": front.id, "rear": next(b.id for b in bs if b is not front),
                          "mean_y": {b.id: round((doc.atoms[b.a].p[1] + doc.atoms[b.b].p[1]) / 2, 2) for b in bs},
                          "ring_bond": {b.id: b.ring for b in bs}})
        subs.append({"fragment": fid, "template": key, "ring_atom_1_at_vertex": start + 1, "direction": dr,
                     "substituents": verts, "heteroatoms": hetero, "crossings": cross})
    return {"bond_length": bl, "templates": templates,
            "templates_per_bond_length": {k: [tuple(round(c / bl, 4) for c in v) for v in t["vertices"]]
                                          for k, t in templates.items()},
            "substituted": subs}
