"""Scheme sanity checks on the scene (design.md §6.4).

M2:
- Ink boxes of scheme elements (species, plus signs, arrows, each conditions line, numbers)
  must not overlap. M4 adds resonance arrows and brackets (each bracket as its bar and two
  lips, so the contributors inside it don't count as overlaps).
- Conditions must not run past their arrow's tail or head.

M3 (check_mechanism). Each curved arrow is sampled into a polyline (mechanism.SAMPLES
segments), and these must hold:
- the arrow does not cross or touch a bond, a label, a charge, another curved arrow, a
  reaction or resonance arrow, a bracket, conditions or a number. Plus signs are allowed:
  mech_basic's 2e⁻ arrow passes over its '+'. The step right after a resonance bracket may
  cross that bracket (its arrows reach the last contributor from outside);
- its head and tail stay at least margin_width from label glyph boxes (each text run's ink box);
- lone-pair and radical dots do not overlap bonds, labels, charges or other dots, and no bond
  runs between the two dots of a pair;
- dots are not drawn on an unlabeled atom with at most one bond (a bare line end, like ⁻:C≡N
  drawn without its C): show_carbons labels it.
"""

from __future__ import annotations

import math

from .geometry import point_box_dist, seg_point_dist, segment_box_dist, segment_distance, union_box
from .layout import SchemeLayout
from .mechanism import clashes, species_parts
from .scene import Polyline


def _overlap(a, b, tol: float = 0.0) -> bool:
    return a[0] < b[2] - tol and b[0] < a[2] - tol and a[1] < b[3] - tol and b[1] < a[3] - tol


def elements(lay: SchemeLayout) -> list[tuple[str, tuple]]:
    out = [(f"species-{n}", p.group.ink_bbox()) for n, p in lay.species.items()]
    out += [(f"plus-{i}", r.bbox()) for i, r in enumerate(lay.pluses, 1)]
    for st in lay.steps:
        out.append((f"arrow-step-{st.n}", st.group.bbox()))
        for side in ("above", "below"):
            for k, ln in enumerate(getattr(st, side), 1):
                out.append((f"conditions-step-{st.n}-{side}-{k}", ln.ink()))
    out += [(f"number-{n}", union_box(r.bbox() for r in runs)) for n, runs in lay.numbers.items()]
    for res in lay.resonance:
        out += [(g.id, g.bbox()) for _, _, g in res.arrows]
        for g in res.group.children:
            line = g.children[0]
            out += [(f"{g.id}-{k}", Polyline([a, b], line.width).ink_bbox())
                    for k, (a, b) in enumerate(zip(line.points, line.points[1:]), 1)]
    return out


def check_scheme(lay: SchemeLayout) -> list[str]:
    """Problems as readable strings; empty when the scheme passes."""
    problems = []
    els = elements(lay)
    for i, (na, a) in enumerate(els):
        for nb, b in els[i + 1:]:
            if na.rsplit("-", 1)[0] == nb.rsplit("-", 1)[0] and na.startswith("bracket-"):
                continue  # one bracket's bar and lips meet at its corners
            if _overlap(a, b):
                problems.append(f"{na} overlaps {nb}")
    for st in lay.steps:
        for side in ("above", "below"):
            for k, ln in enumerate(getattr(st, side), 1):
                ink = ln.ink()
                if ink[0] < st.tail[0] or ink[2] > st.head[0]:
                    problems.append(f"conditions-step-{st.n}-{side}-{k} ({ln.text!r}) is wider than its arrow")
    return problems


def check_mechanism(lay: SchemeLayout) -> list[str]:
    """Problems with curved arrows and dots as readable strings (module docstring)."""
    parts = {n: species_parts(p.group) for n, p in lay.species.items()}
    problems = []
    sw = lay.preset.curved_arrow.stroke_width
    scheme_boxes = [(n, b) for n, b in elements(lay) if not n.startswith(("species-", "plus-"))]
    for k, a in enumerate(lay.curved):
        problems += [f"{a.group.id} {c}" for c in clashes(a, parts, [o for o in lay.curved[:k] if o.step == a.step], lay.preset)]
        for el, box in scheme_boxes:
            if el.startswith(f"bracket-step-{a.step - 1}-"):
                continue  # the step after a bracket reaches into it
            if any(segment_box_dist(p, q, box) < sw / 2 for p, q in zip(a.polyline(), a.polyline()[1:])):
                problems.append(f"{a.group.id} crosses {el}")
    for sp, pt in parts.items():
        for k, (c, kind) in enumerate(pt.dots):
            atom_id = c.source.split(':', 1)[-1]
            what = f"{kind} on atom {atom_id} of {sp}"
            hit = next((src for a, b, hw, src in pt.bonds if seg_point_dist(c.center, a, b) < c.r + hw), None)
            if hit:
                problems.append(f"{what} overlaps a bond ({hit}){_unlabeled_hint(lay, sp, atom_id)}")
            for box, text, _ in pt.labels:
                if point_box_dist(c.center, box) < c.r:
                    problems.append(f"{what} overlaps the label {text!r}")
            if any(seg_point_dist(c.center, a, b) < c.r + hw for a, b, hw, _ in pt.charges):
                problems.append(f"{what} overlaps the charge")
            mate = pt.dots[k + 1][0] if kind == "lone pair" and k % 2 == 0 and k + 1 < len(pt.dots) else None
            if mate is not None:  # a bond running between the two dots of a pair
                hit = next((src for a, b, hw, src in pt.bonds if segment_distance(c.center, mate.center, a, b) < hw), None)
                if hit:
                    problems.append(f"{what} overlaps a bond ({hit})")
            for o, okind in pt.dots[k + 1:]:
                if math.dist(c.center, o.center) < c.r + o.r:
                    problems.append(f"{what} overlaps the {okind} on atom {o.source.split(':', 1)[-1]}")
    for sp, p in lay.species.items():
        if p.drawer is None:
            continue
        for i in sorted(p.drawer.dot_sets()):
            a = p.mol.atoms[i]
            if i not in p.drawer.labels and len(p.mol.neighbors(i)) <= 1:  # a bare line end reads as no atom
                problems.append(f"dots on atom {a.id} of {sp} sit on a line end with no label: "
                                f"add show_carbons: [{a.id}] to {sp}")
    return list(dict.fromkeys(problems))


def _unlabeled_hint(lay, sp: str, atom_id: str) -> str:
    """For dots on a carbon drawn without a label (bonds meet at its point): how to label it."""
    p = lay.species[sp]
    a = next((a for a in p.mol.atoms if a.id == atom_id), None)
    if a is None or (p.drawer and a.idx in p.drawer.labels):
        return ""
    return f"; the carbon has no label, so its dots sit on the bonds: add show_carbons: [{atom_id}] to {sp}"
