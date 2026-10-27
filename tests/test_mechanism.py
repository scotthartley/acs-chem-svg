"""Mechanism layout (draw_mol dots, mechanism.py arrows, acs-style.md §4.2) and checks."""

import math
from functools import lru_cache
from pathlib import Path

import pytest
from validate import validate

from acs_chem_svg import chem
from acs_chem_svg.draw_mol import Drawer
from acs_chem_svg.geometry import add, angle_deg, dot, length, mul, perp, point_box_dist, sub, unit
from acs_chem_svg.geometry_checks import check_mechanism
from acs_chem_svg.layout import build_scheme
from acs_chem_svg.mechanism import ARC_HANDLES, FARTHER, MISS_DEG, aims, miss, species_parts
from acs_chem_svg.preset import load_preset
from acs_chem_svg.scene import Circle
from acs_chem_svg.spec import Spec, load_spec
from acs_chem_svg.svg_writer import to_svg

EX = Path(__file__).resolve().parents[1] / "examples"
M3 = sorted(p.stem for p in EX.glob("m3_*.yaml"))
P = load_preset("acs1996")


@lru_cache(maxsize=None)
def scheme(stem):
    return build_scheme(load_spec(EX / f"{stem}.yaml"), EX, stem)


def sides(smiles, tails=None, flip_h=False):
    m = chem.from_smiles(smiles, "x", P.lengths.bond_length, flip_h=flip_h)
    i = next(a.idx for a in m.atoms if a.map == 1)
    d = Drawer(m, P, show_lone_pairs={i}, tails={i: tails} if tails else None)
    return [(s.kind, s.side) for s in d.dot_sets()[i]]


@pytest.mark.parametrize("stem", M3)
def test_examples_validate(stem):
    assert validate(EX / f"{stem}.yaml") == []


def test_lone_pair_sides_like_chemdraw():
    # mech_basic: H-O⁻ with its bond on the left: top, right, bottom; Cl·: pairs top, right,
    # bottom and the radical on the left; Br⁻: all four sides
    assert sides("[H][O-:1]") in ([("pair", "top"), ("pair", "right"), ("pair", "bottom")],
                                  [("pair", "top"), ("pair", "bottom"), ("pair", "left")])
    assert sides("[Cl:1]") == [("pair", "top"), ("pair", "right"), ("pair", "bottom"), ("radical", "left")]
    assert sides("[Br-:1]") == [("pair", s) for s in ("top", "right", "bottom", "left")]


def test_lone_pairs_spread_out_and_avoid_bonds():
    # acetone O (bond below): left and right, not top + one side
    got = sides("CC(C)=[O:1]")
    assert sorted(s for _, s in got) == ["left", "right"]


def test_tail_set_faces_its_target():
    assert ("radical", "right") in sides("[Br:1]", {"radical": (1.0, 0.0)})
    assert ("pair", "left") in sides("[Br-:1]", {"pair": (-1.0, 0.0)})


def test_dots_are_exact_size_and_gap():
    lay = scheme("m3_sn2_x1.0")
    d = lay.species["hydroxide"].drawer
    i = next(a.idx for a in lay.species["hydroxide"].mol.atoms if a.map == 1)
    box = d.labels[i].elem_ink_box()
    for ds in d.dot_sets()[i]:
        assert length(sub(ds.dots[0], ds.dots[1])) == pytest.approx(P.lone_pair.dot_spacing)
        assert min(point_box_dist(p, box) for p in ds.dots) - P.lone_pair.dot_radius == pytest.approx(P.lone_pair.gap)
    circles = [leaf for leaf in lay.species["hydroxide"].group.leaves() if isinstance(leaf, Circle)]
    assert len(circles) == 6 and all(c.r == P.lone_pair.dot_radius for c in circles)


def test_radical_distance_on_unlabeled_atom():
    m = chem.from_smiles("C[C:1](C)C", "x", P.lengths.bond_length)
    i = next(a.idx for a in m.atoms if a.map == 1)
    (ds,) = Drawer(m, P).dot_sets()[i]  # radicals are drawn without being asked
    assert ds.kind == "radical" and length(sub(ds.center, m.atoms[i].p)) == pytest.approx(P.radical.distance)
    gaps = [math.degrees(math.acos(max(-1, min(1, dot(ds.direction, sub(m.atoms[j].p, m.atoms[i].p)) / P.lengths.bond_length))))
            for j, _ in m.neighbors(i)]
    assert min(gaps) == pytest.approx(60, abs=1)  # on a bisector of the 120° gaps


def test_sn2_arrow_bows_up_over_the_plus_and_keeps_its_gaps():
    lay = scheme("m3_sn2_x1.0")
    a = lay.curved[0]
    assert a.spec.from_ == "lp(1)" and a.c1[1] < a.p0[1] and a.c2[1] < a.p3[1]  # up (y down)
    ref = lay.species["hydroxide"]
    pair = min((ds for ds in ref.drawer.dot_sets()[next(x.idx for x in ref.mol.atoms if x.map == 1)]),
               key=lambda ds: length(sub((ds.center[0] + ref.offset[0], ds.center[1] + ref.offset[1]), a.p0)))
    mid = (pair.center[0] + ref.offset[0], pair.center[1] + ref.offset[1])
    assert length(sub(a.p0, mid)) == pytest.approx(P.lone_pair.dot_radius + P.curved_arrow.tail_gap)
    mebr = lay.species["mebr"]
    c = next(x.idx for x in mebr.mol.atoms if x.map == 2)
    b = mebr.drawer.labels[c].elem_ink_box()
    box = (b[0] + mebr.offset[0], b[1] + mebr.offset[1], b[2] + mebr.offset[0], b[3] + mebr.offset[1])
    assert point_box_dist(a.tip, box) == pytest.approx(P.curved_arrow.head_gap, abs=0.75)


def test_two_fishhooks_from_one_bond_take_opposite_sides():
    lay = scheme("m3_radical_initiation_x1.0")
    p = lay.species["br2"]
    a1, a2 = (x for x in p.mol.atoms if x.map in (1, 2))
    mid = ((a1.x + a2.x) / 2 + p.offset[0], (a1.y + a2.y) / 2 + p.offset[1])
    n = (-(a2.y - a1.y), a2.x - a1.x)
    s1, s2 = (dot(sub(arr.p0, mid), n) for arr in lay.curved)
    assert s1 * s2 < 0
    for arr in lay.curved:  # loops: chords shorter than a bond
        assert length(sub(arr.p3, arr.p0)) < P.lengths.bond_length
        assert len(arr.head) == 3


def test_overrides_flip_the_bow():
    import yaml
    spec = yaml.safe_load((EX / "m3_sn2_x1.0.yaml").read_text())
    spec["overrides"] = {"arrows": [{"step": 1, "index": 0, "bow": "right"}]}
    from acs_chem_svg.spec import Spec
    lay = build_scheme(Spec.model_validate(spec), EX, "sn2_flipped")
    a = lay.curved[0]
    assert a.c1[1] > a.p0[1]  # tail -> head runs rightward, so its right-hand side is down


def test_svg_groups_and_circles():
    svg, _ = to_svg(scheme("m3_sn2_x1.0").root)
    for gid in ("curved-arrows-step-1", "curved-arrow-step-1-0", "species-hydroxide-lone-pairs"):
        assert f'<g id="{gid}">' in svg
    assert svg.count("<circle ") == 6 + 6 + 8  # HO⁻ 3 pairs, the Br of MeBr (an arrow target) 3, Br⁻ 4
    assert 'fill="none" stroke="#000" stroke-width="0.6" stroke-linejoin="miter" stroke-linecap="butt"' in svg


def test_bond_to_atom_arrow_lands_on_the_atom():
    # C=O -> O (carbonyl addition): the head points at the O, not a curl closing below its tail
    lay = scheme("m3_carbonyl_addition_x1.0")
    a = next(x for x in lay.curved if x.step == 1 and x.spec.from_ == "bond(4=5)")
    assert a.target_atom == ("acetone", next(x.idx for x in lay.species["acetone"].mol.atoms if x.map == 5))
    assert aims(a, P.lone_pair.dot_radius)
    assert length(sub(a.p3, a.p0)) > P.curved_arrow.head_length


def test_short_bond_to_atom_arrow_points_at_the_o():
    # water's O-H -> O (inbox 2026-10-26, item 1): the head points at the O's center, not its lower pair
    lay = scheme("m3_carbonyl_addition_x1.0")
    a = next(x for x in lay.curved if x.step == 2 and x.spec.from_ == "bond(7-8)")
    assert miss(a) <= MISS_DEG and aims(a, P.lone_pair.dot_radius)


def _map_point(lay, sp, m):
    p = lay.species[sp]
    return add(next(x.p for x in p.mol.atoms if x.map == m), p.offset)


@pytest.mark.parametrize("stem,step,index,sp", [
    ("m4_nitration_x1.0", 4, 0, "sigma_a"), ("m4_nitration_x1.0", 4, 1, "sigma_b"),  # sigma complex (item 4)
])
def test_hooks_leave_and_land_square(stem, step, index, sp):
    # a bond -> adjacent bond arrow leaves its bond along the normal and lands across the target bond
    lay = scheme(stem)
    a = next(x for x in lay.curved if x.step == step and x.index == index)
    s1, s2 = (_map_point(lay, sp, m) for m in a.spec.source.atoms)
    _, t1, t2 = a.target
    for tangent, (p, q) in ((sub(a.c1, a.p0), (s1, s2)), (sub(a.p3, a.c2), (t1, t2))):
        assert abs(dot(unit(tangent), perp(unit(sub(q, p))))) >= math.cos(math.radians(20)), a.group.id


@pytest.mark.parametrize("stem,step,index", [("m3_e2_x1.0", 1, 1), ("m4_nitration_x1.0", 5, 1)])
def test_ch_hook_goes_through_the_forming_pi_bond_and_back(stem, step, index):
    # C-H -> C-C (user choice 2026-10-26, option D): no room inside the H-C-C angle, so the
    # arrow passes through the target bond, curls back, and its head points at the bond from the
    # far side, turned toward the shared atom; that crossing isn't a clash. Ring C=C -> C-C
    # arrows keep the near side
    lay = scheme(stem)
    a = next(x for x in lay.curved if x.step == step and x.index == index)
    _, t1, t2 = a.target
    n = perp(unit(sub(t2, t1)))
    n = n if dot(sub(a.p0, t1), n) > 0 else mul(n, -1)  # toward the tail's side
    assert dot(sub(a.tip, t1), n) < 0  # the tip is on the far side
    head = unit(sub(a.p3, a.c2))
    assert dot(head, n) > 0  # pointing back at the bond
    (shared,) = set(a.spec.source.atoms) & set(a.spec.target.atoms)
    c = t1 if a.spec.target.atoms[0] == shared else t2  # a.target holds the target atoms in spec order
    mid = mul(add(t1, t2), 0.5)
    assert dot(head, unit(sub(c, mid))) > 0  # turned toward the shared atom
    assert 10 <= angle_deg(head, n) <= 50
    assert a.crosses and check_mechanism(lay) == []
    ring = next(x for x in scheme("m4_nitration_x1.0").curved if x.step == 4 and x.index == 0)
    assert not ring.crosses


@pytest.mark.parametrize("stem", ["m4_nitration_x1.0", "m4_nitration_x0.6"])
@pytest.mark.parametrize("index", [0, 1])
def test_ring_hook_inside_the_ring(stem, index):
    # ring C=C -> adjacent ring bond (user decision 2026-10-27): inside the ring as a circular
    # 120° arc, both ends slid toward the bonds' far atoms; outside it wrapped round the vertex
    lay = scheme(stem)
    a = next(x for x in lay.curved if x.step == 4 and x.index == index)
    sp = lay.species[("sigma_a", "sigma_b")[index]]
    (r,) = sp.mol.rings
    ring = [add(sp.mol.atoms[i].p, sp.offset) for i in r]
    c = mul((sum(p[0] for p in ring), sum(p[1] for p in ring)), 1 / len(ring))
    _, t1, t2 = a.target
    s1, s2 = (add(next(x.p for x in sp.mol.atoms if x.map == m), sp.offset) for m in a.spec.source.atoms)
    for end, (p, q) in ((a.p0, (s1, s2)), (a.tip, (t1, t2))):
        assert length(sub(end, c)) < length(sub(mul(add(p, q), 0.5), c))
    t0, t3 = unit(sub(a.c1, a.p0)), unit(sub(a.p3, a.c2))
    assert angle_deg(t0, t3) == pytest.approx(120, abs=25)
    th, L = math.radians(angle_deg(t0, t3)), length(sub(a.p3, a.p0))
    circ = L * 2 / 3 * math.tan(th / 4) / math.sin(th / 2)
    for h in (length(sub(a.c1, a.p0)), length(sub(a.p3, a.c2))):
        assert h == pytest.approx(circ, rel=0.1)
    assert not a.crosses and check_mechanism(lay) == []


def test_bond_tail_at_the_drawn_bond_midpoint():
    # H-C -> C (radical propagation, item 3): the H label shortens the drawn bond; the tail sits
    # beside the middle of what is drawn, and the loop stays compact
    lay = scheme("m3_radical_propagation_x1.0")
    a = next(x for x in lay.curved if x.spec.to == "atom(3)")
    h, c = _map_point(lay, "isobutane", 6), _map_point(lay, "isobutane", 3)
    u = unit(sub(c, h))
    drawn = [dot(sub(q, h), u) for p, q2, _, _ in species_parts(lay.species["isobutane"].group).bonds
             if abs(dot(unit(sub(q2, p)), u)) > 0.999 and abs(dot(sub(p, h), perp(u))) < 0.1 for q in (p, q2)]
    assert dot(sub(a.p0, h), u) == pytest.approx((min(drawn) + max(drawn)) / 2, abs=1.0)
    b = a.group.bbox()
    assert max(b[2] - b[0], b[3] - b[1]) < 0.75 * P.lengths.bond_length


@pytest.mark.parametrize("h", [0.5, 0.531])
def test_arc_curvature_rises_evenly(h):
    # item 5: ChemDraw's arc dips mid-way and spikes at the head; ours never drops more than 5%
    # below its running maximum, and the head end is no sharper than 1.5x the middle
    c1 = (ARC_HANDLES[0][0], ARC_HANDLES[0][1] * h)
    c2 = (ARC_HANDLES[1][0], ARC_HANDLES[1][1] * h)

    def kappa(t):
        s = 1 - t
        d1 = [3 * (s * s * c1[i] + 2 * s * t * (c2[i] - c1[i]) + t * t * ((1, 0)[i] - c2[i])) for i in (0, 1)]
        d2 = [6 * (s * (c2[i] - 2 * c1[i]) + t * ((1, 0)[i] - 2 * c2[i] + c1[i])) for i in (0, 1)]
        return abs(d1[0] * d2[1] - d1[1] * d2[0]) / math.hypot(*d1) ** 3
    ks = [kappa(k / 200) for k in range(201)]
    top = 0.0
    for k in ks:
        top = max(top, k)
        assert k >= 0.95 * top
    assert ks[-1] <= 1.5 * ks[100]
    if h == 0.531:  # ChemDraw's tangent angles, kept
        assert angle_deg(c1, (1, 0)) == pytest.approx(63.8, abs=0.5)
        assert angle_deg(sub((1, 0), c2), (1, 0)) == pytest.approx(72.8, abs=0.5)


@pytest.mark.parametrize("stem", M3)
def test_heads_stop_clear_of_the_label_and_its_dots(stem):
    # the tip keeps head_gap (or FARTHER) from the nearest of the element ink box and the
    # atom's dots, so it ends above a lone pair, never between its dots, and points into the box
    lay = scheme(stem)
    ca, r = lay.preset.curved_arrow, lay.preset.lone_pair.dot_radius
    for a in lay.curved:
        if a.target[0] != "box":
            continue
        tip = add(a.p3, mul(unit(sub(a.p3, a.c2)), ca.head_inset))  # on the axis: a fishhook's sits fishhook_offset aside
        d = min([point_box_dist(tip, a.target[1])] + [length(sub(tip, q)) - r for q in a.target_dots])
        # roots hold a FARTHER stop exactly; an aiming end (a short loop whose end tangent can't
        # meet the stop, water's O-H -> O above an arrow) stays between the nearest and farthest
        assert ca.head_gap - 0.01 <= d <= ca.head_gap + max(FARTHER) + 0.01, (a.group.id, d)
        assert aims(a, r), a.group.id


def test_target_atoms_keep_their_lone_pairs_in_place():
    # the arrows are decorations: target pairs sit where they would with no arrow, and heads clear them
    for stem, sp, m in (("m3_carbonyl_addition_x1.0", "acetone", 5), ("m3_e2_x1.0", "etbr", 5),
                        ("m3_sn2_x0.6", "mebr", 3), ("m3_radical_initiation_x1.0", "br2", 1)):
        lay = scheme(stem)
        p = lay.species[sp]
        i = next(x.idx for x in p.mol.atoms if x.map == m)
        plain = Drawer(p.mol, lay.preset, show_lone_pairs={i}).dot_sets()[i]
        assert [(d.side, d.dots) for d in p.drawer.dot_sets()[i] if d.kind == "pair"] == \
               [(d.side, d.dots) for d in plain if d.kind == "pair"]
    assert validate(EX / "m3_carbonyl_addition_x1.0.yaml") == []  # no head touches them


def test_fishhooks_meet_where_the_bond_forms():
    # tBu· + Br2: electron(3) and bond(7-8) both go to between(3,7), past the plus sign, on one side
    lay = scheme("m3_radical_propagation_2_x1.0")
    a, b = (x for x in lay.curved if x.spec.to == "between(3,7)")
    plus = lay.pluses[0].bbox()
    for x in (a, b):
        assert plus[0] - 3 < x.tip[0] < plus[2] + 3 and x.tip[1] < plus[1]  # above the plus, not on an atom
    assert length(sub(a.tip, b.tip)) < 3 * P.curved_arrow.head_gap


def test_charge_moves_left_when_a_bond_takes_the_upper_right():
    m = chem.from_smiles("[O-:1][H]", "x", P.lengths.bond_length)
    i = next(a.idx for a in m.atoms if a.map == 1)
    d = Drawer(m, P)
    x0, _, x1, _ = d.labels[i].elem_ink_box()
    assert d.charge_center(i)[0] < x0
    m = chem.from_smiles("[H][O-:1]", "x", P.lengths.bond_length)
    i = next(a.idx for a in m.atoms if a.map == 1)
    d = Drawer(m, P)
    assert d.charge_center(i)[0] > d.labels[i].elem_ink_box()[2]


@pytest.mark.parametrize("rotate", range(0, 360, 30))
def test_ring_carbanion_pair_on_bisector_charge_beside(rotate):
    # inbox 2026-10-26: [CH-] in a ring with an arrow from its pair. The pair takes the outward
    # bisector; the ⊖ goes beside it, outside the ring, clear of the dots and the ring bonds
    spec = Spec.model_validate({
        "preset": "acs1996",
        "species": {"ion_adj": {"smiles": "[CH:1]1[CH:2]=[CH:3][CH+:4][CH-:5][CH:6]=1",
                                "orientation": {"rotate": rotate}},
                    "ion_para": {"smiles": "[CH-:1]1[CH:2]=[CH:3][CH+:4][CH:5]=[CH:6]1"}},
        "steps": [{"resonance": ["ion_adj", "ion_para"],
                   "arrows": [[{"from": "lp(5)", "to": "between(5,6)"}, {"from": "bond(6=1)", "to": "atom(1)"}]]}]})
    lay = build_scheme(spec, EX, "ring_carbanion")
    problems = check_mechanism(lay)
    assert not [p for p in problems if "overlaps the charge" in p or "overlaps a bond" in p], problems
    for name, mapno in (("ion_adj", 5), ("ion_para", 1)):
        d = lay.species[name].drawer
        i = next(a.idx for a in d.mol.atoms if a.map == mapno)
        gap, out = d._largest_gap(i)
        (pair,) = [s for s in d.dot_sets()[i] if s.kind == "pair"]
        assert angle_deg(pair.direction, out) < 1.0
        q = sub(d.charge_center(i), d.mol.atoms[i].p)
        assert angle_deg(unit(q), out) <= math.degrees(gap) / 2 + 1e-6  # in the exterior gap: outside the ring
        assert length(q) >= P.charge.unlabeled_distance - 1e-9


# ---------------------------------------------------------------- heads as markers


def _old_head(c1, c2, p3, kind, p):
    """The curved-arrow head as mechanism.head_polygon built it before markers."""
    from acs_chem_svg.geometry import perp
    ca = p.curved_arrow
    c_last = c2 if length(sub(p3, c2)) > 1e-6 else c1
    t = unit(sub(p3, c_last))
    nv = perp(t)
    tip = add(p3, mul(t, ca.head_inset))

    def at(back, side):
        return add(tip, add(mul(t, -back), mul(nv, side)))
    if kind == "two_electron":
        return [at(ca.head_length, ca.head_half_width), p3, at(ca.head_length, -ca.head_half_width), tip]
    bend = add(sub(p3, mul(c2, 2)), c1)
    side = -1.0 if dot(nv, bend) > 0 else 1.0
    notch = add(p3, mul(nv, -side * ca.fishhook_offset))
    tip = add(notch, mul(t, ca.head_inset))
    return [notch, add(tip, add(mul(t, -ca.head_length), mul(nv, side * ca.fishhook_barb))), tip]


CURVES = {"up": ((0.0, 0.0), (6.0, -22.0), (30.0, -22.0), (38.0, -2.0)),
          "down": ((0.0, 0.0), (6.0, 22.0), (30.0, 22.0), (38.0, 2.0)),
          "c2_on_end": ((0.0, 0.0), (10.0, -15.0), (30.0, 4.0), (30.0, 4.0))}  # ChemDraw's C–Br arrow


@pytest.mark.parametrize("kind", ["two_electron", "fishhook"])
@pytest.mark.parametrize("curve", sorted(CURVES))
def test_curved_head_marker_places_the_old_head(curve, kind):
    from acs_chem_svg.mechanism import curved_arrow, head_polygon
    p0, c1, c2, p3 = CURVES[curve]
    (stroke,) = curved_arrow(p0, c1, c2, p3, kind, P).children
    (head,) = stroke.head_leaves()
    old = _old_head(c1, c2, p3, kind, P)
    assert len(head.points) == len(old)
    assert all(length(sub(a, b)) < 1e-6 for a, b in zip(head.points, old))
    assert head_polygon(c1, c2, p3, kind, P) == head.points


def test_fishhook_markers_take_both_sides():
    from acs_chem_svg.mechanism import head_marker
    ids = {head_marker(*CURVES[c][1:], "fishhook", P).id for c in ("up", "down")}
    assert ids == {"acs1996-head-fishhook-left", "acs1996-head-fishhook-right"}


PHENOL = Path(__file__).resolve().parent / "specs" / "phenol_x_abstraction.yaml"


@pytest.mark.parametrize("bow", [None, "left", "right"])
def test_two_fishhooks_from_a_ring_double_bond(bow, tmp_path):
    """inbox 2026-10-26: phenoxy -> o-radical. The fishhook to atom(3) inside the ring starts
    past the double bond's inner line; on the same side as its partner (bow: left), its tail
    slides along the bond toward atom 3. Either way nothing crosses."""
    import yaml
    data = yaml.safe_load(PHENOL.read_text())
    if bow:
        data["overrides"] = {"arrows": [{"step": 2, "index": 2, "bow": bow}]}
    path = tmp_path / "phenol.yaml"
    path.write_text(yaml.safe_dump(data))
    assert validate(path) == []
    lay = build_scheme(load_spec(path), PHENOL.parent, "phenol")
    one, two = (a for a in lay.curved if a.step == 2 and a.index in (1, 2))
    sp = lay.species["phenoxy"]
    c3, c4 = (next(add(a.p, sp.offset) for a in sp.mol.atoms if a.map == m) for m in (3, 4))
    mid, u = mul(add(c3, c4), 0.5), unit(sub(c3, c4))
    n = (-u[1], u[0])
    off = [dot(sub(a.p0, mid), n) for a in (one, two)]
    gap = P.curved_arrow.bond_tail_gap
    if (off[0] > 0) == (off[1] > 0):  # same side: the second tail moved toward atom 3
        assert dot(sub(two.p0, mid), u) == pytest.approx(gap)
    else:  # the inner tail is past the inner line, not at bond_tail_gap
        assert max(abs(v) for v in off) > gap + 1.0
