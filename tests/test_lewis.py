"""Lewis structures (`lewis:` on a species): lattice layout, labels, lone pairs, bond electrons."""

import math

import pytest
import yaml

from acs_chem_svg import chem
from acs_chem_svg.draw_mol import Drawer, element_x
from acs_chem_svg.fontmetrics import load_face
from acs_chem_svg.geometry_checks import check_lewis
from acs_chem_svg.preset import load_preset
from acs_chem_svg.scene import Circle, Group, Polyline
from acs_chem_svg.spec import Spec

PRESET = load_preset("acs1996")
L = PRESET.lewis.bond_length
SHARED = ["ClC(Cl)(Cl)Cl", "CO", "O", "O=C=O", "C#N", "N#N", "[NH4+]", "[C-]#N", "[N]=O", "CC(=O)O",
          "[OH3+]", "[O-][N+](=O)[O-]"]


def mol(smiles, mode="shared"):
    return chem.from_smiles(smiles, "m", L, lewis=mode)


def drawer(smiles, mode="shared", preset=PRESET):
    return Drawer(chem.from_smiles(smiles, "m", preset.lewis.bond_length, lewis=mode), preset)


# --------------------------------------------------------------------------- layout


@pytest.mark.parametrize("smiles", SHARED)
def test_bonds_on_the_lattice(smiles):
    m = mol(smiles)
    for b in m.bonds:
        dx, dy = (abs(m.atoms[b.b].x - m.atoms[b.a].x), abs(m.atoms[b.b].y - m.atoms[b.a].y))
        assert min(dx, dy) == pytest.approx(0, abs=1e-6)
        assert max(dx, dy) == pytest.approx(L, abs=1e-6)
    pts = [a.p for a in m.atoms]
    assert len({(round(x, 3), round(y, 3)) for x, y in pts}) == len(pts)


@pytest.mark.parametrize("smiles", ["ClC(Cl)(Cl)Cl", "CO", "O", "O=C=O", "C#N", "[NH4+]", "[O-][N+](=O)[O-]"])
def test_no_unbonded_neighbors_where_avoidable(smiles):
    m = mol(smiles)
    bonded = {frozenset((b.a, b.b)) for b in m.bonds}
    for a in m.atoms:
        for o in m.atoms[a.idx + 1:]:
            if abs(math.dist(a.p, o.p) - L) < 1e-6:
                assert frozenset((a.idx, o.idx)) in bonded, smiles


def test_acetic_acid_keeps_one_unavoidable_contact_and_its_lone_pairs_avoid_it():
    m = mol("CC(=O)O")
    bonded = {frozenset((b.a, b.b)) for b in m.bonds}
    close = [(a.idx, o.idx) for a in m.atoms for o in m.atoms[a.idx + 1:]
             if abs(math.dist(a.p, o.p) - L) < 1e-6 and frozenset((a.idx, o.idx)) not in bonded]
    assert len(close) == 1
    d = Drawer(m, PRESET)
    for i, j in close:
        for ds in d.dot_sets().get(i, []) + d.dot_sets().get(j, []):
            other = j if ds.atom == i else i
            toward = ((m.atoms[other].x - m.atoms[ds.atom].x) / L, (m.atoms[other].y - m.atoms[ds.atom].y) / L)
            assert ds.direction[0] * toward[0] + ds.direction[1] * toward[1] < 0.5


def test_methanol_is_straight():
    m = mol("CO")
    c, o = m.atoms[0], m.atoms[1]
    h_o = next(a for a in m.atoms if a.element == "H" and any(b.other(a.idx) == 1 for _, b in m.neighbors(a.idx)))
    assert c.y == o.y == h_o.y
    h_c = [a for a in m.atoms if a.element == "H" and a is not h_o]
    assert sorted((round(a.x - c.x, 3), round(a.y - c.y, 3)) for a in h_c) == sorted([(-L, 0), (0, -L), (0, L)])


@pytest.mark.parametrize("smiles", ["O", "O=C=O"])
def test_linear_triatomics(smiles):
    assert len({round(a.y, 6) for a in mol(smiles).atoms}) == 1


def test_ccl4_is_a_cross():
    m = mol("ClC(Cl)(Cl)Cl")
    c = m.atoms[1]
    assert sorted((round(a.x - c.x), round(a.y - c.y)) for a in m.atoms if a.element == "Cl") == \
        sorted([(round(L), 0), (-round(L), 0), (0, round(L)), (0, -round(L))])


@pytest.mark.parametrize("smiles,msg", [("C1CC1", "rings are not supported"), ("ClP(Cl)(Cl)(Cl)Cl", "has 5 neighbors"),
                                        ("CC(C)(C)C", "no square-lattice layout")])
def test_layout_errors(smiles, msg):
    with pytest.raises(ValueError, match=msg):
        mol(smiles)


# --------------------------------------------------------------------------- labels and lone pairs


@pytest.mark.parametrize("smiles", SHARED)
def test_every_atom_labeled_without_h_runs(smiles):
    d = drawer(smiles)
    assert set(d.labels) == set(range(len(d.mol.atoms)))
    assert all(a.num_h == 0 for a in d.mol.atoms)
    assert all(r.text in ("C", "H", "N", "O", "Cl") for lab in d.labels.values() for r in lab.runs)


@pytest.mark.parametrize("smiles", SHARED)
def test_every_lone_pair_drawn(smiles):
    d = drawer(smiles)
    for a in d.mol.atoms:
        sets = d.dot_sets().get(a.idx, [])
        assert sum(ds.kind == "pair" for ds in sets) == a.lone_pairs
        assert sum(ds.kind == "radical" for ds in sets) == a.radicals


def test_lone_pair_counts():
    assert mol("O").atoms[0].lone_pairs == 2
    assert [a.lone_pairs for a in mol("ClC(Cl)(Cl)Cl").atoms if a.element == "Cl"] == [3, 3, 3, 3]


# --------------------------------------------------------------------------- shared pairs


def _u(fr, p):
    return (p[0] - fr.o[0]) * fr.u[0] + (p[1] - fr.o[1]) * fr.u[1]


def _v(fr, p):
    return (p[0] - fr.o[0]) * fr.n[0] + (p[1] - fr.o[1]) * fr.n[1]


@pytest.mark.parametrize("smiles", SHARED)
def test_shared_columns(smiles):
    d = drawer(smiles)
    lp = PRESET.lone_pair
    for ds in d.lewis_dot_sets():
        b = d.mol.bonds[ds.bond]
        fr, u0, u1 = d.bond_gap(b)
        us = sorted({round(_u(fr, p), 6) for p in ds.dots})
        assert len(us) == b.order and len(ds.dots) == 2 * b.order
        assert all(abs(abs(_v(fr, p)) - lp.dot_spacing / 2) < 1e-9 for p in ds.dots)  # across the bond
        assert all(math.isclose(y - x, lp.dot_spacing) for x, y in zip(us, us[1:]))
        assert (us[0] + us[-1]) / 2 == pytest.approx((u0 + u1) / 2)  # centered between the ink boxes


@pytest.mark.parametrize("name", ["acs1996", "acs1996_x0.8", "acs1996_x0.6"])
@pytest.mark.parametrize("smiles", SHARED)
def test_shared_pairs_farther_out_than_lone_pairs(name, smiles):
    """The distinctness rule: a shared dot is farther from each label than that atom's lone pairs."""
    preset = load_preset(name)
    d = drawer(smiles, preset=preset)
    lp = preset.lone_pair
    assert check_lewis("m", d) == []
    for ds in d.lewis_dot_sets():
        b = d.mol.bonds[ds.bond]
        fr, u0, u1 = d.bond_gap(b)
        near = min(min(_u(fr, p) - u0, u1 - _u(fr, p)) for p in ds.dots) - lp.dot_radius
        assert near >= 2 * lp.gap - 1e-9
        assert near > lp.gap


@pytest.mark.parametrize("name", ["acs1996", "acs1996_x0.8", "acs1996_x0.6"])
def test_bond_length_is_the_smallest_that_fits(name):
    """lewis.bond_length (preset comment): the smallest 0.5 pt value where three columns fit between
    two Cl or Br labels (first letter centered on both atoms) with 2 x lone_pair.gap to each box."""
    p = load_preset(name)
    face, size, lp = load_face("Regular"), p.text.label_size, p.lone_pair
    need = 0.0
    for el in ("Cl", "Br"):
        b = face.text_ink_bounds(el, size)
        ex = element_x(face, el, size, "right")
        need = max(need, (ex + b[2]) - (ex + b[0]) + 2 * lp.dot_spacing + 2 * lp.dot_radius + 4 * lp.gap)
    assert p.lewis.bond_length == math.ceil(2 * need) / 2


# --------------------------------------------------------------------------- assigned electrons


def owned(d):
    """atom index -> electrons it owns: lone pairs, radicals and assigned bond electrons."""
    out = {a.idx: 2 * a.lone_pairs + a.radicals for a in d.mol.atoms}
    for ds in d.lewis_dot_sets():
        out[ds.atom] += len(ds.dots)
    return out


def valence(el):
    from rdkit import Chem
    return Chem.GetPeriodicTable().GetNOuterElecs(el)


@pytest.mark.parametrize("mode", ["formal_charge", "oxidation_state"])
def test_assigned_keeps_label_positions(mode):
    for smiles in ("O=C=O", "C#N", "CO"):
        a, b = drawer(smiles), drawer(smiles, mode)
        assert [x.elem_origin for x in a.labels.values()] == [x.elem_origin for x in b.labels.values()]


@pytest.mark.parametrize("smiles", ["O=C=O", "C#N", "[NH4+]", "[C-]#N", "[N]=O"])
def test_formal_charge_counts(smiles):
    d = drawer(smiles, "formal_charge")
    e = owned(d)
    assert check_lewis("m", d) == []
    for a in d.mol.atoms:
        assert valence(a.element) - e[a.idx] == a.charge, (smiles, a.element)


def test_oxidation_states():
    d = drawer("CO", "oxidation_state")
    e = owned(d)
    assert check_lewis("m", d) == []
    states = {(a.element, valence(a.element) - e[a.idx]) for a in d.mol.atoms}
    assert states == {("C", -2), ("O", -2), ("H", 1)}
    d = drawer("O=C=O", "oxidation_state")
    assert valence("C") - owned(d)[1] == 4


def test_assigned_columns_pack_from_the_owner():
    d = drawer("C#N", "formal_charge")
    lp = PRESET.lone_pair
    cn = next(b for b in d.mol.bonds if b.order == 3)
    fr, u0, u1 = d.bond_gap(cn)
    for ds in (s for s in d.lewis_dot_sets() if s.bond == cn.idx):
        edge, sgn = (u0, 1) if ds.atom == cn.a else (u1, -1)
        offs = sorted(round(sgn * (_u(fr, p) - edge), 6) for p in ds.dots)
        first = lp.gap + lp.dot_radius
        assert offs == pytest.approx([first, first, first + lp.dot_spacing])  # a pair, then one dot on the axis


# --------------------------------------------------------------------------- spec and SVG


def test_svg_groups_and_no_bond_lines():
    g = drawer("O=C=O").draw("co2")
    ids = {c.id for c in g.children}
    assert "species-co2-bonding-pairs" in ids and "species-co2-bonds" not in ids
    assert not [x for x in g.leaves() if isinstance(x, Polyline)]
    g = drawer("O=C=O", "formal_charge").draw("co2")
    grp = next(c for c in g.children if c.id == "species-co2-bond-electrons")
    assert isinstance(grp, Group) and all(isinstance(c, Circle) for c in grp.children)


def _spec(**sp):
    return Spec.model_validate({"species": {"x": {"smiles": "O=C=O", **sp}}})


def test_spec_lewis_true_means_shared():
    assert _spec(lewis=True).species["x"].lewis_mode == "shared"
    assert _spec(lewis={"electrons": "oxidation_state"}).species["x"].lewis_mode == "oxidation_state"
    assert _spec(lewis=False).species["x"].lewis_mode is None


@pytest.mark.parametrize("extra,msg", [
    ({"show_carbons": "terminal"}, "lewis can't be combined with show_carbons"),
    ({"show_lone_pairs": [1]}, "lewis can't be combined with show_lone_pairs"),
    ({"abbreviate": ["OMe"]}, "lewis can't be combined with abbreviate"),
    ({"orientation": {"rotate": 45}}, "multiple of 90"),
])
def test_spec_lewis_conflicts(extra, msg):
    with pytest.raises(ValueError, match=msg):
        _spec(lewis=True, **extra)


def test_spec_lewis_rotate_and_flip_ok():
    _spec(lewis=True, orientation={"rotate": 90, "flip_h": True})


def test_spec_lewis_not_in_schemes():
    with pytest.raises(ValueError, match="lewis species are M1-only for now"):
        Spec.model_validate(yaml.safe_load("""
species:
  a: {smiles: "O", lewis: true}
  b: {smiles: "[OH-]"}
steps: [{from: a, to: b}]
"""))


# --------------------------------------------------------------------------- bonds: lines


def lines_drawer(smiles, preset=PRESET):
    return Drawer(chem.from_smiles(smiles, "m", preset.lewis.bond_length, lewis="shared", lewis_lines=True), preset)


@pytest.mark.parametrize("smiles", SHARED)
def test_lines_draw_bonds_not_shared_pairs(smiles):
    """bonds: lines keeps the dots layout and lone pairs, draws each bond as lines and no shared pair."""
    d, dots = lines_drawer(smiles), drawer(smiles)
    assert d.lewis_dot_sets() == []
    assert [a.p for a in d.mol.atoms] == [a.p for a in dots.mol.atoms]
    assert d.dot_sets() == dots.dot_sets()
    g = {c.id.rsplit("-", 1)[-1]: c for c in d.draw("m").children}
    assert "bonds" in g and "bonding" not in "".join(g)
    assert check_lewis("m", d) == []


def test_lines_need_shared_electrons():
    with pytest.raises(ValueError, match="bonds: lines"):
        Spec.model_validate(yaml.safe_load(
            'preset: acs1996\nspecies: {co2: {smiles: "O=C=O", lewis: {electrons: formal_charge, bonds: lines}}}'))
