"""Mechanism layout (draw_mol dots, mechanism.py arrows, acs-style.md §4.2) and checks."""

import math
from functools import lru_cache
from pathlib import Path

import pytest
from validate import validate

from acs_chem_svg import chem
from acs_chem_svg.draw_mol import Drawer
from acs_chem_svg.geometry import add, dot, length, mul, point_box_dist, sub, unit
from acs_chem_svg.layout import build_scheme
from acs_chem_svg.mechanism import FARTHER, aims
from acs_chem_svg.preset import load_preset
from acs_chem_svg.scene import Circle
from acs_chem_svg.spec import load_spec
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
