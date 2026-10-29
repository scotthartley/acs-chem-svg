"""Scheme sanity checks on the scene (design.md §6.4).

M2:
- Ink boxes of scheme elements (species, plus signs, arrows, each conditions line, numbers)
  must not overlap. M4 adds resonance arrows and brackets (each bracket as its bar and two
  lips, so the contributors inside it don't count as overlaps).
- Conditions must not run past their arrow's tail or head.
- Inside each species (M1 rows too: check_labels), an atom label must not overlap another
  atom's label or a bond that isn't on that atom (draw_mol flips a crowded CH₃ to H₃C first),
  and drawn bond lines of bonds that share no atom must not touch: a crossing ChemDraw-style,
  with the rear bond cut (draw_mol, acs-style.md §14.2), leaves a gap and passes.

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

Lewis structures (check_lewis, M1):
- a shared pair's dots keep at least 2 × lone_pair.gap from both element ink boxes, so they
  read as farther out than lone pairs (a wide `labels:` symbol such as R' can break this);
- an atom's assigned bond electrons stay on its half of the gap between the ink boxes;
- dots of different sets (lone pairs, radicals, bond electrons) are at least dot_spacing apart.
"""

from __future__ import annotations

import math

from .chem import Molecule
from .draw_mol import Drawer
from .geometry import point_box_dist, seg_point_dist, segment_box_dist, segment_distance, union_box
from .layout import SchemeLayout
from .mechanism import clashes, species_parts
from .scene import Group, Polyline

LABEL_TOL = 0.05  # pt; boxes that only touch, or a bond end grazing a label, don't count


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
    for n, p in lay.species.items():
        problems += check_labels(n, p.group, p.mol)
    return problems


def check_labels(name: str, group: Group, mol: Molecule) -> list[str]:
    """Labels of one drawn species that overlap another atom's label, or a bond not on their atom."""
    pt = species_parts(group)
    problems = []
    for k, (box, _, src) in enumerate(pt.labels):
        atom_id = src.split(":", 1)[-1]
        i = next((a.idx for a in mol.atoms if a.id == atom_id), None)
        own = {b.id for _, b in mol.neighbors(i)} if i is not None else set()
        for obox, _, osrc in pt.labels[k + 1:]:
            if osrc != src and _overlap(box, obox, LABEL_TOL):
                problems.append(f"species-{name}: {src} label overlaps {osrc} label")
        for a, b, hw, bsrc in pt.bonds:
            if own.isdisjoint(bsrc.split(":", 1)[-1].split(",")) and segment_box_dist(a, b, box) < hw - LABEL_TOL:
                problems.append(f"species-{name}: {src} label overlaps {bsrc}")
    atoms_of = {b.id: {b.a, b.b} for b in mol.bonds}

    def atoms(src):
        return set().union(*(atoms_of.get(i, set()) for i in src.split(":", 1)[-1].split(",")))
    segs = [(a, b, hw, src, atoms(src)) for a, b, hw, src in pt.bonds]
    for k, (a, b, hw, src, at) in enumerate(segs):
        for c, d, hw2, src2, at2 in segs[k + 1:]:
            if src2 != src and at.isdisjoint(at2) and segment_distance(a, b, c, d) < hw + hw2 - LABEL_TOL:
                problems.append(f"species-{name}: {src} crosses {src2} without a gap")
    return list(dict.fromkeys(problems))


def check_lewis(name: str, d: Drawer) -> list[str]:
    """Problems with a Lewis structure's dots (module docstring); empty for other species."""
    if d.mol.lewis is None:
        return []
    lp = d.preset.lone_pair
    problems = []
    bond_sets = d.lewis_dot_sets()
    for ds in bond_sets:
        b = d.mol.bonds[ds.bond]
        fr, u0, u1 = d.bond_gap(b)
        ids = f"{d.mol.atoms[b.a].id}-{d.mol.atoms[b.b].id}"
        for p in ds.dots:
            u = (p[0] - fr.o[0]) * fr.u[0] + (p[1] - fr.o[1]) * fr.u[1]
            if ds.kind == "shared":
                near = min(u - u0, u1 - u) - lp.dot_radius
                if near < 2 * lp.gap - LABEL_TOL:
                    problems.append(f"species-{name}: a shared pair of bond {ids} is {near:.2f} pt from a label "
                                    f"(at least {2 * lp.gap:.2f}, 2 × lone_pair.gap, so it reads as a bond)")
                    break
            else:
                mid = (u0 + u1) / 2
                past = (u + lp.dot_radius - mid) if ds.atom == b.a else (mid - (u - lp.dot_radius))
                if past > LABEL_TOL:
                    problems.append(f"species-{name}: the electrons of bond {ids} assigned to atom "
                                    f"{d.mol.atoms[ds.atom].id} cross the middle of the bond by {past:.2f} pt")
                    break
    sets = [ds for i in sorted(d.dot_sets()) for ds in d.dot_sets()[i]] + bond_sets
    for k, s1 in enumerate(sets):
        for s2 in sets[k + 1:]:
            gap = min(math.dist(p, q) for p in s1.dots for q in s2.dots)
            if gap < lp.dot_spacing - LABEL_TOL:
                problems.append(f"species-{name}: dots of atom {d.mol.atoms[s1.atom].id} ({s1.kind}) and atom "
                                f"{d.mol.atoms[s2.atom].id} ({s2.kind}) are {gap:.2f} pt apart (at least {lp.dot_spacing:.2f})")
    return list(dict.fromkeys(problems))


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
        if p.drawer is not None:
            problems += check_line_end_dots(sp, p.drawer)
    return list(dict.fromkeys(problems))


def check_line_end_dots(name: str, d) -> list[str]:
    """Dots (a radical, a carbene's pair) on an unlabeled terminal atom: a bare line end reads as no atom."""
    return [f"dots on atom {d.mol.atoms[i].id} of {name} sit on a line end with no label: "
            f"add show_carbons: [{d.mol.atoms[i].id}] to {name}"
            for i in sorted(d.dot_sets()) if i not in d.labels and len(d.mol.neighbors(i)) <= 1]


def _unlabeled_hint(lay, sp: str, atom_id: str) -> str:
    """For dots on a carbon drawn without a label (bonds meet at its point): how to label it."""
    p = lay.species[sp]
    a = next((a for a in p.mol.atoms if a.id == atom_id), None)
    if a is None or (p.drawer and a.idx in p.drawer.labels):
        return ""
    return f"; the carbon has no label, so its dots sit on the bonds: add show_carbons: [{atom_id}] to {sp}"
