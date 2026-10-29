"""Chairs (chair.py) and bond crossings (draw_mol, acs-style.md §14)."""

import math

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from acs_chem_svg import chair, chem
from acs_chem_svg.draw_mol import Drawer, cut_pieces
from acs_chem_svg.preset import load_preset
from acs_chem_svg.scene import Polygon, Polyline

PRESET = load_preset("acs1996")
BL = PRESET.lengths.bond_length


def _faces_3d(smi):
    """Substituent faces from an ETKDG embedding: the side of the ring's normal (oriented so the
    ring order runs counterclockwise seen from it) relative to the geminal partner."""
    m = Chem.MolFromSmiles(smi)
    ring = chair.chair_rings(m)[0]
    mh = Chem.AddHs(m)
    AllChem.EmbedMolecule(mh, randomSeed=7)
    conf = mh.GetConformer()

    def p(i):
        return list(conf.GetAtomPosition(i))
    c = [sum(p(i)[k] for i in ring) / 6 for k in range(3)]
    n = [0.0, 0.0, 0.0]
    for a, b in zip(ring, ring[1:] + ring[:1]):
        pa, pb = [p(a)[k] - c[k] for k in range(3)], [p(b)[k] - c[k] for k in range(3)]
        n = [n[0] + pa[1] * pb[2] - pa[2] * pb[1], n[1] + pa[2] * pb[0] - pa[0] * pb[2], n[2] + pa[0] * pb[1] - pa[1] * pb[0]]
    out = {}
    for i in ring:
        nb = [x.GetIdx() for x in mh.GetAtomWithIdx(i).GetNeighbors() if x.GetIdx() not in ring]
        h = [sum((p(s)[k] - p(i)[k]) * n[k] for k in range(3)) for s in nb]
        for s, hv, other in zip(nb, h, h[::-1]):
            if s < m.GetNumAtoms():
                out[s] = 1 if hv > other else -1
    return m, ring, out


@pytest.mark.parametrize("smi", [
    "C[C@H]1CC[C@@H](C)CC1", "C[C@H]1CC[C@H](C)CC1", "C[C@H]1C[C@@H](C)CCC1", "C[C@H]1C[C@H](C)CCC1",
    "C[C@H]1CCCC[C@@H]1C", "C[C@H]1CCCC[C@H]1C", "CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O",
    "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"])
def test_faces_match_a_3d_embedding(smi):
    m, ring, want = _faces_3d(smi)
    got = chair.faces(m, ring)
    assert {s: f for s, f in got.items() if f is not None} == {s: want[s] for s, f in got.items() if f is not None}


def _axial(smi, **kw):
    """Map numbers of the substituents drawn axial."""
    p = chem.prepare(chem.parse_smiles(smi), "t", BL, chair=chem.ChairSpec(PRESET.chair, **kw))
    mol = p.molecule()
    ring = set(p.chair_ring)
    out = set()
    for a in mol.atoms:
        if a.idx in ring or not a.map:
            continue
        r = next(j for j, _ in mol.neighbors(a.idx) if j in ring)
        d = (a.x - mol.atoms[r].x, a.y - mol.atoms[r].y)
        if abs(d[0]) < 1e-6:
            out.add(a.map)
    return out


@pytest.mark.parametrize("smi, axial", [
    ("[CH3:1][C@H]1CC[C@H]([CH3:2])CC1", set()),  # trans-1,4: diequatorial
    ("[CH3:1][C@H]1CC[C@@H]([CH3:2])CC1", "one"),  # cis-1,4: one axial
    ("[CH3:1][C@H]1C[C@@H]([CH3:2])CCC1", set()),  # cis-1,3: diequatorial
    ("[CH3:1][C@H]1C[C@H]([CH3:2])CCC1", "one"),  # trans-1,3
    ("CC(C)(C)[C@H]1CC[C@@H]([OH:1])CC1", {1}),  # cis-4-tBu: tBu (4.9) equatorial, OH axial
    ("[CH3:1][C@H]1CC[C@@H]([F:2])CC1", {2}),  # cis: F (0.25) axial rather than CH3 (1.74)
])
def test_flip_puts_the_bulkier_group_equatorial(smi, axial):
    got = _axial(smi)
    assert (len(got) == 1) if axial == "one" else got == axial


def test_1_2_pair_cis_one_axial_trans_diequatorial():
    got = [_axial(s) for s in ("[CH3:1][C@H]1CCCC[C@H]1[CH3:2]", "[CH3:1][C@H]1CCCC[C@@H]1[CH3:2]")]
    assert sorted(len(g) for g in got) == [0, 1]


def test_flip_override_and_axial_request():
    assert _axial("[CH3:1]C1CCCCC1") == set()
    assert _axial("[CH3:1]C1CCCCC1", flip=True) == {1}
    assert _axial("[CH3:1][C@H]1CC[C@@H]([F:2])CC1", flip=True) == {1}
    assert _axial("[CH3:1][C@H]1CC[C@@H]([F:2])CC1", axial=[1]) == {1}
    with pytest.raises(ValueError, match="no ring flip"):
        _axial("[CH3:1][C@H]1CC[C@H]([CH3:2])CC1", axial=[1], equatorial=[2])


def test_geminal_pair_one_axial_one_equatorial():
    got = _axial("[CH3:1][C:2]1([OH:3])CCCCC1")
    assert got == {3}  # OH (0.87) axial, CH3 (1.74) equatorial


def test_chair_geometry_is_the_template():
    p = chem.prepare(chem.parse_smiles("C1CCCCC1"), "t", BL, chair=chem.ChairSpec(PRESET.chair))
    mol = p.molecule()
    pts = [mol.atoms[i].p for i in p.chair_ring]
    lengths = sorted(round(math.dist(a, b), 2) for a, b in zip(pts, pts[1:] + pts[:1]))
    ref = [tuple(v) for v in PRESET.chair.vertices]
    want = sorted(round(math.dist(a, b), 2) for a, b in zip(ref, ref[1:] + ref[:1]))
    assert lengths == want


@pytest.mark.parametrize("name", ["acs1996", "acs1996_x0.8", "acs1996_x0.6"])
def test_flipped_chair_is_the_exact_mirror(name):
    c = load_preset(name).chair
    (v, up, eq), (fv, fup, feq) = chair.template(c, False), chair.template(c, True)
    assert fv[0] == (0.0, 0.0)
    for k in range(6):
        j = (4 - k) % 6
        assert fv[k][1] == v[j][1] and fup[k] == up[j] and feq[k] == (-eq[j][0], eq[j][1])
    # mirror about one vertical line: x + x' is the same for every vertex pair
    assert max(fv[k][0] + v[(4 - k) % 6][0] for k in range(6)) - min(fv[k][0] + v[(4 - k) % 6][0] for k in range(6)) < 1e-9


def test_depth_follows_the_ring_atom():
    p = chem.prepare(chem.parse_smiles("[CH3:1][CH2:2]C1CCCCC1"), "t", BL, chair=chem.ChairSpec(PRESET.chair))
    mol = p.molecule()
    ring_atom = next(j for j, _ in mol.neighbors(mol.by_id("2").idx) if j in p.chair_ring)
    assert mol.by_id("1").z == mol.by_id("2").z == mol.atoms[ring_atom].z == mol.atoms[ring_atom].y
    assert not any(b.stereo for b in mol.bonds)


def test_wrong_rings_are_refused():
    fused = "[CH2:1]1[CH2:2][CH2:3][CH:4]2CCCC[CH:5]2[CH2:6]1"
    for smi, ring, msg in (("C1=CCCCC1", None, "drawn flat"), ("c1ccccc1", None, "drawn flat"),
                           (fused, None, "more than one six-membered ring"), (fused, [1, 2, 3, 4, 5, 6], "fused"),
                           ("C1CCCC1", None, "no six-membered ring")):
        with pytest.raises(ValueError, match=msg):
            chem.prepare(chem.parse_smiles(smi), "t", BL, chair=chem.ChairSpec(PRESET.chair, ring=ring))


# --------------------------------------------------------------------------- crossings


def _cross(order_rear=1, angle=90.0, front=None):
    """A horizontal bond (rear by default) crossed at its middle by a bond at `angle`."""
    rear = "C=C" if order_rear == 2 else "CC"
    m = chem.from_smiles(f"{rear}.CC", "x", BL, rotate_deg=0)
    a = math.radians(angle)
    pts = [(0.0, 0.0), (BL, 0.0), (BL / 2 - BL / 2 * math.cos(a), BL / 2 * math.sin(a)),
           (BL / 2 + BL / 2 * math.cos(a), -BL / 2 * math.sin(a))]
    for atom, (x, y) in zip(m.atoms, pts):
        atom.x, atom.y = x, y
    if front is not None:
        m.front = [front]
    return m


@pytest.mark.parametrize("angle", [30.0, 45.0, 60.0, 74.0, 90.0])
def test_rear_bond_cut_parallel_to_the_front_bond(angle):
    m = _cross(angle=angle)
    d = Drawer(m, PRESET)
    d.draw("x")
    geoms = {b.idx: d.bond_geometry(b) for b in m.bonds}
    d.apply_crossings(geoms)
    assert geoms[0].cuts and not geoms[1].cuts  # the later bond is in front
    lw, c = PRESET.lengths.line_width, PRESET.lengths.margin_width
    pieces = cut_pieces(geoms[0].frame, 0.0, 0.0, BL, geoms[0].cuts, lw)
    assert len(pieces) == 2
    s = math.sin(math.radians(angle))
    k = 1 / math.tan(math.radians(angle))
    for pc, cap in ((pieces[0], pieces[0].cap1), (pieces[1], pieces[1].cap0)):
        mid = sum(cap) / 2
        assert abs(abs(mid - BL / 2) - (lw / 2 + c) / s) < 1e-6  # centerline cut on the strip edge
        assert abs(abs(cap[0] - cap[1]) - lw * k) < 1e-6  # corners on a line parallel to the front bond


def test_front_list_and_double_rear_bond():
    m = _cross(order_rear=2, front=1)
    g = Drawer(m, PRESET).draw("x")
    bonds = next(c for c in g.children if c.id.endswith("-bonds"))
    rear_src = f"bond:{m.bonds[0].id}"
    # the rear double bond's two lines are each cut in two (perpendicular: no oblique caps)
    lines = [c for c in bonds.children if c.source == rear_src]
    assert len(lines) == 4 and all(isinstance(c, (Polyline, Polygon)) for c in lines)
    m2 = _cross(order_rear=2, front=1)
    d = Drawer(m2, PRESET)
    assert d.in_front(m2.bonds[1], m2.bonds[0]) and not d.in_front(m2.bonds[0], m2.bonds[1])


def _axial_template():
    """chair_x1.0's six axial methyls (templates/chair/, committed), with ChemDraw's front bonds:
    two crossings, axial bonds over and under ring bonds."""
    import json
    from pathlib import Path
    t = Path(__file__).resolve().parents[1] / "templates" / "chair" / "chair_x1.0__3.mol"
    mol = chem.from_molfile(t, "t", BL)
    ids = {a.id: a.idx for a in mol.atoms}
    pairs = {frozenset(ids[i] for i in p) for p in json.loads(t.with_suffix(".front.json").read_text())}
    mol.front = [b.idx for b in mol.bonds if frozenset((b.a, b.b)) in pairs]
    return mol


def test_cut_ring_bond_keeps_the_miters():
    """Axial bonds crossing ring bonds: each cut ring bond ends one chain and starts another, the
    ring's joins stay at its atoms, and no line runs into a front bond's strip."""
    mol = _axial_template()
    d = Drawer(mol, PRESET)
    assert len(d.crossings()) == 2
    g = d.draw("t")
    bonds = next(c for c in g.children if c.id.endswith("-bonds"))
    chains = [c for c in bonds.children if isinstance(c, Polyline)]
    ring = [c for c in chains if len(c.points) > 2]
    assert ring and not any(c.closed for c in ring)  # the ring is cut, so it is no longer one closed line
    assert sum(len(c.points) - 2 for c in ring) == 6  # a miter join at every ring atom
    from acs_chem_svg.geometry_checks import check_labels
    assert not [p for p in check_labels("t", g, mol) if "crosses" in p]


def test_cdxml_puts_rear_bonds_below_front_bonds():
    from acs_chem_svg.cdxml_writer import z_order
    mol = _axial_template()
    d = Drawer(mol, PRESET)
    pairs = d.crossings()
    assert len(pairs) == 2
    rank = z_order(mol, d)
    assert all(rank[r] < rank[f] for f, r in pairs)
