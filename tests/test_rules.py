"""Unit tests of the pure drawing rules in draw_mol (references/acs-style.md §1, §2)."""

import math

import pytest
from cdxml_to_molfile import template_path

from acs_chem_svg import chem
from acs_chem_svg.draw_mol import Drawer, element_x, h_side
from acs_chem_svg.preset import load_preset
from acs_chem_svg.scene import Polyline

P = load_preset("acs1996")
TOL = 0.01


def drawer(ref, n, preset=P):
    return Drawer(chem.from_molfile(template_path(ref, n), f"{ref}__{n}", preset.lengths.bond_length), preset)


def doubles(d):
    return [d.bond_geometry(b) for b in d.mol.bonds if b.order == 2]


def second_line_insets(g):
    """(inset at a, inset at b) of an offset/ring second line."""
    (v, u0, u1), = g.lines
    return round(u0, 3), round(g.frame.L - u1, 3)


# --------------------------------------------------------------------------- §1.1 double bonds

INSET = P.lengths.ring_inner_inset  # 1.497 = s tan 30


@pytest.mark.parametrize("n, name", [(1, "propene"), (2, "2-methylpropene"), (3, "(E)-2-butene"),
                                     (4, "(Z)-2-butene"), (5, "acetaldehyde, H drawn"), (6, "acetone"),
                                     (7, "cyclohexanone"), (8, "methylenecyclohexane")])
def test_double_bond_placement(n, name):
    """The eight double bonds of double_bonds_x1.0 (acs-style.md §1.1 table)."""
    g, = doubles(drawer("double_bonds_x1.0", n))
    expected = {1: "offset", 2: "centered", 3: "offset", 4: "offset"}.get(n, "centered")
    assert g.placement == expected, name
    if expected == "offset":
        insets = sorted(second_line_insets(g))
        want = {1: [0, INSET], 3: [0, INSET], 4: [INSET, INSET]}[n]
        assert insets == pytest.approx(want, abs=TOL), name
        assert abs(g.lines[0][0]) == pytest.approx(P.lengths.bond_spacing)
    else:
        assert sorted(v for v, _, _ in g.lines) == pytest.approx([-P.lengths.bond_spacing / 2,
                                                                  P.lengths.bond_spacing / 2])


def test_offset_line_goes_toward_substituents():
    """Propene: the second line sits on the side of the methyl."""
    d = drawer("double_bonds_x1.0", 1)
    g, = doubles(d)
    b = next(b for b in d.mol.bonds if b.order == 2)
    me = next(j for i in (b.a, b.b) for j, _ in d.mol.neighbors(i) if j not in (b.a, b.b))
    side = (d.mol.atoms[me].x - g.frame.o[0]) * g.frame.n[0] + (d.mol.atoms[me].y - g.frame.o[1]) * g.frame.n[1]
    assert math.copysign(1, side) == math.copysign(1, g.lines[0][0])


def test_centered_lines_run_into_degree3_join():
    """C(3)=CH2: both lines run past the C(3) atom to the neighbor centerlines (u = -h/tan 60)."""
    d = drawer("double_bonds_x1.0", 2)
    g, = doubles(d)
    b = next(b for b in d.mol.bonds if b.order == 2)
    deg3_is_a = d.mol.degree(b.a) == 3
    want = -(P.lengths.bond_spacing / 2) / math.tan(math.radians(60))  # -0.748
    for v, u0, u1 in g.lines:
        run_on = u0 if deg3_is_a else g.frame.L - u1
        square = g.frame.L - u1 if deg3_is_a else u0
        assert run_on == pytest.approx(want, abs=TOL)
        assert square == pytest.approx(0, abs=1e-9)


@pytest.mark.parametrize("ref, n", [("m1_set_x1.0", 5), ("m1_set_x1.0", 4)])
def test_carbonyl_and_nitronium_centered(ref, n):
    for g in doubles(drawer(ref, n)):
        assert g.placement == "centered"


def test_ring_inner_line_inset():
    """Benzene: inner lines inside the ring, inset s tan 30 at both ends."""
    d = drawer("m1_set_x1.0", 1)
    for g in doubles(d):
        assert g.placement == "ring"
        assert second_line_insets(g) == pytest.approx((INSET, INSET), abs=TOL)
        assert P.lengths.bond_spacing * math.tan(math.radians(30)) == pytest.approx(INSET, abs=0.002)


# --------------------------------------------------------------------------- chains / joins


def test_chain_is_one_mitered_polyline():
    d = drawer("double_bonds_x1.0", 9)  # hexane
    group = d.draw("hexane")
    lines = [x for x in group.leaves() if isinstance(x, Polyline)]
    assert len(lines) == 1 and len(lines[0].points) == 6 and lines[0].join == "miter"


def test_ring_main_lines_close():
    d = drawer("m1_set_x1.0", 1)
    closed = [x for x in d.draw("benzene").leaves() if isinstance(x, Polyline) and x.closed]
    assert len(closed) == 1 and len(closed[0].points) == 6


# --------------------------------------------------------------------------- §2 labels


@pytest.mark.parametrize("n, side, runs", [(2, "left", ["H", "2", "N"]), (1, "right", ["NH", "2"]),
                                           (4, "left", ["HO"]), (3, "right", ["OH"])])
def test_h_side_and_runs(n, side, runs):
    d = drawer("labels_x1.0", n)
    i, = d.labels
    assert h_side(d.mol, i) == side
    assert [r.text for r in d.labels[i].runs] == runs


def test_left_label_set_right_to_left():
    """H2N: H at -11.393 and 2 at -4.171 from the N origin (labels_x1.0)."""
    d = drawer("labels_x1.0", 2)
    lab, = d.labels.values()
    h, two, n = lab.runs
    assert (h.x - n.x, two.x - n.x) == pytest.approx((-11.393, -4.171), abs=0.005)
    assert two.size == pytest.approx(0.75 * n.size) and two.y - n.y == pytest.approx(2.30)


@pytest.mark.parametrize("preset", ["acs1996", "acs1996_x0.8", "acs1996_x0.6"])
def test_label_anchor_and_baseline(preset):
    p = load_preset(preset)
    d = drawer("labels_x1.0", 4, p)
    i, lab = next(iter(d.labels.items()))
    a = d.mol.atoms[i]
    ox, base = lab.elem_origin  # HO: the O is centered on the atom, the H set to its left
    assert ox + d.face.advance("O", lab.size) / 2 == pytest.approx(a.x)
    assert base - a.y == pytest.approx(0.3585 * p.text.label_size)
    assert lab.runs[0].text == "HO" and lab.runs[0].x + d.face.advance("H", lab.size) == pytest.approx(ox)


# --------------------------------------------------------------------------- §2.2 charges


def test_charge_on_unlabeled_carbon():
    """Isopropyl cation: the + sits on the bisector of the largest gap, unlabeled_distance away."""
    mol = chem.from_smiles("C[CH+]C", "ipr", P.lengths.bond_length)
    d = Drawer(mol, P)
    i = next(a.idx for a in mol.atoms if a.charge)
    assert i not in d.labels
    cx, cy = d.charge_center(i)
    a = mol.atoms[i]
    assert math.dist((cx, cy), a.p) == pytest.approx(P.charge.unlabeled_distance)
    dirs = [math.atan2(mol.atoms[j].y - a.y, mol.atoms[j].x - a.x) for j, _ in mol.neighbors(i)]
    bis = math.atan2(-(math.sin(dirs[0]) + math.sin(dirs[1])), -(math.cos(dirs[0]) + math.cos(dirs[1])))
    assert math.atan2(cy - a.y, cx - a.x) == pytest.approx(bis, abs=1e-6)
    bars = d.charge_bars(i)
    assert len(bars) == 2 and all(b.width == P.charge.bar_thickness for b in bars)
    assert all(math.dist(*b.points) == pytest.approx(P.charge.bar_length) for b in bars)


def test_charge_above_nh3():
    d = drawer("labels_x1.0", 5)
    i, = d.labels
    a = d.mol.atoms[i]
    assert d.charge_center(i) == pytest.approx((a.x, a.y - P.charge.above_label))


# --------------------------------------------------------------------------- §1.2 stereo


@pytest.mark.parametrize("preset, count", [("acs1996", 3), ("acs1996_x0.8", 2), ("acs1996_x0.6", 2)])
def test_hash_count_and_last_length(preset, count):
    p = load_preset(preset)
    ref = {"acs1996": "m1_set_x1.0", "acs1996_x0.8": "m1_set_x0.8", "acs1996_x0.6": "m1_set_x0.6"}[preset]
    d = drawer(ref, 3, p)
    g = next(d.bond_geometry(b) for b in d.mol.bonds if b.stereo == "hash")
    assert len(g.hashes) == count
    assert 2 * g.hashes[-1][1] == pytest.approx(p.lengths.wedge_wide_end)
    assert g.hashes[0][0] - p.lengths.hash_line_width / 2 == pytest.approx(p.lengths.hash_spacing)


def test_wedge_outline_starts_at_half_line_width():
    d = drawer("m1_set_x1.0", 3)
    g = next(d.bond_geometry(b) for b in d.mol.bonds if b.stereo == "wedge")
    apex, s1, w1, w2, s2 = g.polygon
    assert math.dist(w1, w2) == pytest.approx(P.lengths.wedge_wide_end)
    fr = g.frame
    for s, w in ((s1, w1), (s2, w2)):  # side line extrapolated to the stereocenter: |v| = lw/2
        su = (s[0] - fr.o[0]) * fr.u[0] + (s[1] - fr.o[1]) * fr.u[1]
        sv = (s[0] - fr.o[0]) * fr.n[0] + (s[1] - fr.o[1]) * fr.n[1]
        wu = (w[0] - fr.o[0]) * fr.u[0] + (w[1] - fr.o[1]) * fr.u[1]
        wv = (w[0] - fr.o[0]) * fr.n[0] + (w[1] - fr.o[1]) * fr.n[1]
        v0 = sv - (wv - sv) / (wu - su) * su
        assert abs(v0) == pytest.approx(P.lengths.line_width / 2, abs=1e-6)


# --------------------------------------------------------------------------- chem input


def test_kekule_override():
    m = chem.from_smiles("[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1", "bz", 14.4, kekule=["2=3", "4=5", "6=1"])
    dbl = {frozenset((m.atoms[b.a].map, m.atoms[b.b].map)) for b in m.bonds if b.order == 2}
    assert dbl == {frozenset((2, 3)), frozenset((4, 5)), frozenset((6, 1))}


def test_smiles_bond_length_and_ez():
    m = chem.from_smiles("C/C=C\\C", "z", 11.52)
    for b in m.bonds:
        assert math.dist(m.atoms[b.a].p, m.atoms[b.b].p) == pytest.approx(11.52, abs=0.01)
    c1, c2, c3, c4 = (a.p for a in m.atoms)
    # Z: both methyls on the same side of the C=C axis
    side = lambda p: (c3[0] - c2[0]) * (p[1] - c2[1]) - (c3[1] - c2[1]) * (p[0] - c2[0])  # noqa: E731
    assert side(c1) * side(c4) > 0


def test_spec_validation(tmp_path):
    from pydantic import ValidationError

    from acs_chem_svg.spec import Spec
    with pytest.raises(ValidationError):
        Spec.model_validate({"species": {"x": {"smiles": "C", "molfile": "a.mol"}}})
    with pytest.raises(ValidationError):
        Spec.model_validate({"preset": "acs2000", "species": {"x": {"smiles": "C"}}})
    with pytest.raises(ValidationError):
        Spec.model_validate({"species": {"x": {"smiles": "c1ccccc1", "kekule": ["1-2"]}}})


@pytest.mark.parametrize("which, labeled", [("none", set()), ("terminal", {"1", "4"}), ([4], {"4"})])
def test_show_carbons(which, labeled):
    """Explicit carbon labels: CC(=O)C has carbons 1, 2, 4; O (3) is always labeled."""
    m = chem.show_carbons(chem.from_smiles("CC(=O)C", "acetone", P.lengths.bond_length), which)
    d = Drawer(m, P)
    assert {m.atoms[i].id for i in d.labels} == labeled | {"3"}
    if "1" in labeled:
        assert "".join(r.text for r in d.labels[0].runs) in ("CH3", "H3C")


def test_show_carbons_rejects_unknown_ids():
    m = chem.from_smiles("CCO", "ethanol", P.lengths.bond_length)
    with pytest.raises(ValueError, match="no carbon with id 3"):
        chem.show_carbons(m, [3])
    with pytest.raises(NotImplementedError, match="not terminal"):
        chem.show_carbons(m, [2])


def test_molfile_carbon_alias():
    """stereo_x1.0 fragment 18: ChemDraw's CH3 label arrives as a molfile alias."""
    d = drawer("stereo_x1.0", 18)
    c = next(i for i, a in enumerate(d.mol.atoms) if a.show_label)
    assert c in d.labels and d.mol.atoms[c].num_h == 3
    with pytest.raises(ValueError, match="condensed label"):
        chem.from_molfile(template_path("stereo_x1.0", 18).read_text().replace("CH3\n", "CH2Cl\n"),
                          "x", P.lengths.bond_length)
    with pytest.raises(ValueError, match="attaches through O, not C"):  # an abbreviation on the wrong atom
        chem.from_molfile(template_path("stereo_x1.0", 18).read_text().replace("CH3\n", "OMe\n"),
                          "x", P.lengths.bond_length)
    with pytest.raises(ValueError, match="has 3 H"):
        chem.from_molfile(template_path("stereo_x1.0", 18).read_text().replace("CH3\n", "CH2\n"),
                          "x", P.lengths.bond_length)


def test_spec_show_carbons():
    from pydantic import ValidationError

    from acs_chem_svg.spec import Spec
    assert Spec.model_validate({"species": {"x": {"smiles": "CC", "show_carbons": [1, "2"]}}})
    with pytest.raises(ValidationError):
        Spec.model_validate({"species": {"x": {"smiles": "CC", "show_carbons": "all"}}})


# --------------------------------------------------------------------------- new references (stereo_x*, clearance_x*, double_bonds_2)


def test_two_letter_label_right_justified_when_bond_points_right():
    """clearance_x1.0: -Cl puts the C on the atom; Cl- puts the l on the atom."""
    seen = set()
    for n in range(1, 85):
        d = drawer("clearance_x1.0", n)
        cl = [i for i, lab in d.labels.items() if lab.element == "Cl"]
        if not cl:
            continue
        i = cl[0]
        lab, a = d.labels[i], d.mol.atoms[i]
        (j, _), = d.mol.neighbors(i)
        dx = d.mol.atoms[j].x - a.x
        anchor = "l" if dx > 1e-6 else "C"
        x0 = lab.elem_origin[0] + (d.face.advance("C", lab.size) if anchor == "l" else 0)
        assert x0 + d.face.advance(anchor, lab.size) / 2 == pytest.approx(a.x, abs=1e-6)
        seen.add(anchor)
    assert seen == {"C", "l"}


def test_e_alkene_tie_goes_to_negative_v():
    """double_bonds_2 fragments 1-4: (E)-2-butene at four rotations, second line at v < 0."""
    for n in (1, 2, 3, 4):
        g, = doubles(drawer("double_bonds_2", n))
        assert g.placement == "offset" and g.lines[0][0] < 0


@pytest.mark.parametrize("n, size", [(10, 6), (11, 6)])
def test_fused_ring_bond_goes_into_ring_with_more_double_bonds(n, size):
    """Naphthalene (10) and indene (11): the shared double bond's second line is in a
    six-membered ring (indene: not the smaller five-membered ring)."""
    d = drawer("double_bonds_2", n)
    for b in d.mol.bonds:
        rings = [r for r in d.mol.rings if b.a in r and b.b in r]
        if b.order != 2 or len(rings) < 2:
            continue
        g = d.bond_geometry(b)
        chosen = next(r for r in rings if d.ring_side(b, g.frame) == (1 if g.lines[0][0] > 0 else -1)
                      and len(r) == size)
        assert len(chosen) == size


def test_five_ring_inset_uses_bisector():
    """Indene's five-membered ring double bond: inset 1.88 at 108° (s cot 54°)."""
    d = drawer("double_bonds_2", 11)
    five = next(r for r in d.mol.rings if len(r) == 5)
    b = next(b for b in d.mol.bonds if b.order == 2 and b.a in five and b.b in five
             and not any(b.a in r and b.b in r for r in d.mol.rings if len(r) == 6))
    assert second_line_insets(d.bond_geometry(b)) == pytest.approx((1.884, 1.884), abs=0.01)


@pytest.mark.parametrize("n", [13, 14])
def test_wedge_wide_end_notched_at_atom_with_bonds(n):
    """stereo_x1.0: a wedge into a ring carbon (13) or chain carbon (14) runs 0.55 pt past
    the atom and is cut back to the atom center (in 14 the cut is straight: ChemDraw's
    corners and the center are collinear within 0.01 pt)."""
    d = drawer("stereo_x1.0", n)
    b = next(b for b in d.mol.bonds if b.stereo == "wedge")
    g = d.bond_geometry(b)
    loc = [((p[0] - g.frame.o[0]) * g.frame.u[0] + (p[1] - g.frame.o[1]) * g.frame.u[1]) for p in g.polygon]
    assert max(loc) == pytest.approx(g.frame.L + 0.55, abs=0.03)
    assert any(math.dist(p, d.mol.atoms[b.b].p) < 1e-9 for p in g.polygon)


def test_chain_bond_fills_wedge_notch():
    """stereo_x1.0 fragment 14 (butane): the one plain bond at the wedge's wide end starts
    behind the atom, at the wedge's far corner (ChemDraw: 1.51 pt), so the join is solid."""
    d = drawer("stereo_x1.0", 14)
    wide = next(b for b in d.mol.bonds if b.stereo == "wedge").b
    b = next(b for b in d.mol.bonds if not b.stereo and wide in (b.a, b.b))
    u0, u1 = d.bond_geometry(b).main
    assert (-u0 if b.a == wide else u1 - d.frame(b).L) == pytest.approx(1.511, abs=0.02)


@pytest.mark.parametrize("n", [5, 6])
def test_terminal_label_centered_between_double_lines(n):
    """double_bonds_2 acetaldehyde: the O label moves bond_spacing/2 (1.30 pt) toward the
    second line, so it sits centered between the two lines; the lines don't move."""
    d = drawer("double_bonds_2", n)
    b = next(b for b in d.mol.bonds if b.order == 2)
    o = b.a if b.a in d.labels else b.b
    lab, fr = d.labels[o], d.frame(b)
    g = d.bond_geometry(b)
    (v, _, _), = g.lines
    a = d.mol.atoms[o]
    shift = (lab.elem_origin[0] - a.x - element_x(lab.face, "O", lab.size, lab.h_side),
             lab.elem_origin[1] - a.y - P.ratios.label_baseline_per_size * lab.size)
    assert (shift[0] * fr.n[0] + shift[1] * fr.n[1]) == pytest.approx(v / 2, abs=TOL)
    assert (shift[0] * fr.u[0] + shift[1] * fr.u[1]) == pytest.approx(0, abs=TOL)
    assert abs(v) == pytest.approx(P.lengths.bond_spacing)


@pytest.mark.parametrize("n, count", [(1, 5), (2, 4), (3, 2), (4, 6), (5, 8), (6, 12)])
def test_hash_count_at_unlabeled_ends(n, count):
    """stereo_x1.0 fragments 1-6: hashes of 1, 0.75, 0.56, 1.2, 1.5, 2.2 x bond length."""
    d = drawer("stereo_x1.0", n)
    g = next(d.bond_geometry(b) for b in d.mol.bonds if b.stereo == "hash")
    assert len(g.hashes) == count
