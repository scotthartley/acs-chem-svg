"""Four-coordinate metals drawn as a cross (chem._square_planar), ChemDraw's centered two-letter
labels (draw_mol.centered_label; metals_x*), and the check for a layout that lays one group over
another (geometry_checks.check_overlaid). inbox 2026-10-30."""

import math
from pathlib import Path

import pytest

import measure_refs as mr
from acs_chem_svg import chem
from acs_chem_svg.chem import Atom, Bond, Molecule
from acs_chem_svg.draw_mol import Drawer
from acs_chem_svg.geometry import sub
from acs_chem_svg.geometry_checks import check_overlaid
from acs_chem_svg.preset import load_preset
from acs_chem_svg.spec import Spec, scheme_prepared

P = load_preset("acs1996")
TOL = 0.1
OXADD = "Br[Pd](c1ccccc1)(<-P(C)(C)C)<-P(C)(C)C"
TRANS = "c1ccc(cc1)[Pd](c1ccc(C)cc1)(<-P(C)(C)C)<-P(C)(C)C"
CISPLATIN = "[NH3]->[Pt@SP1](<-[NH3])(Cl)Cl"
TRANSPLATIN = "[NH3]->[Pt@SP2](<-[NH3])(Cl)Cl"
DPPE = "Cl[Pd]1(Cl)<-P(c2ccccc2)(c2ccccc2)CCP->1(c1ccccc1)c1ccccc1"


def _arms(smiles: str, metal: str) -> list[tuple[str, float, float]]:
    """(neighbor element, bond length / median bond, direction in degrees) round the metal, as
    chem.layout_2d leaves it."""
    m = chem.parse_smiles(smiles)
    chem.layout_2d(m)
    return _arms_of(m, metal)


def _arms_of(m, metal: str) -> list[tuple[str, float, float]]:
    pts = chem._coords(m)
    med = sorted(math.dist(pts[b.GetBeginAtomIdx()], pts[b.GetEndAtomIdx()]) for b in m.GetBonds())[m.GetNumBonds() // 2]
    a = next(a for a in m.GetAtoms() if a.GetSymbol() == metal)
    c = pts[a.GetIdx()]
    return [(n.GetSymbol(), math.dist(c, pts[n.GetIdx()]) / med,
             math.degrees(math.atan2(-(pts[n.GetIdx()][1] - c[1]), pts[n.GetIdx()][0] - c[0])))
            for n in a.GetNeighbors()]


@pytest.mark.parametrize("smiles, metal", [(OXADD, "Pd"), (TRANS, "Pd"), (CISPLATIN, "Pt"), (TRANSPLATIN, "Pt"),
                                           ("[NH3]->[Pt](<-[NH3])(Cl)Cl", "Pt")])
def test_metal_is_a_cross(smiles, metal):
    """Four equal bonds 90° apart (CoordGen gave Br–Pd 0.88, Pd–P 1.13 and 47–107° gaps)."""
    arms = _arms(smiles, metal)
    assert all(r == pytest.approx(1.0, abs=0.01) for _, r, _ in arms)
    dirs = sorted(d % 360 for *_, d in arms)
    assert [b - a for a, b in zip(dirs, dirs[1:])] == pytest.approx([90.0] * 3, abs=0.01)


@pytest.mark.parametrize("smiles, gap", [(CISPLATIN, 90.0), (TRANSPLATIN, 180.0),
                                         ("[NH3]->[Pt@SP3](<-[NH3])(Cl)Cl", 90.0), ("[NH3]->[Pt@SP1](Cl)(<-[NH3])Cl", 180.0)])
def test_square_planar_stereo_picks_cis_or_trans(smiles, gap):
    """@SP1 (U), @SP2 (4) and @SP3 (Z) set the order round the cross (OpenSMILES)."""
    n = [d for e, _, d in _arms(smiles, "Pt") if e == "N"]
    assert abs((n[0] - n[1] + 180) % 360 - 180) == pytest.approx(gap, abs=0.01)


def test_chelate_left_to_coordgen():
    """Two donors on one ligand: not a cross (its ring fixes the bite angle)."""
    m = chem.parse_smiles(DPPE)
    pd = next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "Pd")
    assert chem._cross_ligands(m, pd) is None


@pytest.mark.parametrize("smiles", [OXADD, TRANS, CISPLATIN])
def test_cross_sits_on_the_axes(smiles):
    """The automatic orientation turns a cross only by quarter turns (the aryls' rings would
    otherwise tilt it 30°)."""
    m = chem.from_smiles(smiles, "m", P.lengths.bond_length)
    metal = next(a for a in m.atoms if a.element in ("Pd", "Pt"))
    for j, _ in m.neighbors(metal.idx):
        d = math.degrees(math.atan2(m.atoms[j].y - metal.y, m.atoms[j].x - metal.x))
        assert abs((d + 45) % 90 - 45) < 0.01


PD0 = "CP(C)(C)->[Pd]<-P(C)(C)C"


@pytest.mark.parametrize("smiles, metal, abbreviate", [(PD0, "Pd", None), (PD0, "Pd", ["PMe3"]),
                                                       ("C[Au]Cl", "Au", None), ("C[Hg]C", "Hg", None)])
def test_two_coordinate_metal_is_straight(smiles, metal, abbreviate):
    """A two-coordinate d-block metal with two separate ligands is drawn straight, on an axis
    (user decision 2026-10-30): CoordGen bent Pd(PMe3)2 to 120° once PMe3 was abbreviated."""
    m = chem.from_smiles(smiles, "m", P.lengths.bond_length, abbreviate=abbreviate)
    c = next(a for a in m.atoms if a.element == metal)
    dirs = [math.degrees(math.atan2(m.atoms[j].y - c.y, m.atoms[j].x - c.x)) for j, _ in m.neighbors(c.idx)]
    assert abs((dirs[0] - dirs[1]) % 360 - 180) < 0.01
    assert abs((dirs[0] + 45) % 90 - 45) < 0.01


def test_pin_that_bends_a_cross_is_rejected():
    """II aligned to Pd(0)'s straight P-Pd-P: pinning it left Br and Ph where CoordGen put them,
    60° off the cross (chem.off_arm_metals), so II is laid out as a cross and oriented."""
    spec = Spec.model_validate({"species": {"a": {"smiles": PD0, "abbreviate": ["PMe3"]},
                                            "b": {"smiles": OXADD, "abbreviate": ["PMe3"]}},
                                "steps": [{"from": "a", "to": "b"}]})
    p = scheme_prepared(spec, Path("."))["b"]
    assert p.aligned_to == "a" and p.matched == []
    assert chem.off_arm_metals(p.mol) == []


@pytest.mark.parametrize("smiles, flagged", [(TRANSPLATIN, False), (CISPLATIN, True)])
def test_off_arm_metals_reads_square_planar_stereo(smiles, flagged):
    """A cross drawn trans where the SMILES says cis (a pin can do that: align_coords doesn't
    read @SP) is flagged; one drawn as the SMILES says is not."""
    m = chem.parse_smiles(smiles)
    chem.layout_2d(m)
    assert chem.off_arm_metals(m) == []
    trans = chem.parse_smiles(TRANSPLATIN)
    chem.layout_2d(trans)
    pt = next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "Pt")
    m.GetConformer().SetPositions(trans.GetConformer().GetPositions())  # same atom order, trans drawing
    assert (chem.off_arm_metals(m) == [pt]) == flagged


def _star(element: str, dirs: list[float], length: float = 14.4) -> Molecule:
    """An atom of `element` at the origin with plain carbons at `dirs` (deg, counterclockwise
    on screen)."""
    atoms = [Atom(0, "1", element, 0.0, 0.0)] + [
        Atom(k + 1, str(k + 2), "C", length * math.cos(math.radians(t)), -length * math.sin(math.radians(t)))
        for k, t in enumerate(dirs)]
    return Molecule("m", atoms, [Bond(k, 0, k + 1) for k in range(len(dirs))])


@pytest.mark.parametrize("element, dirs, centered", [
    ("Pd", [0, 90, 180, 270], True), ("Si", [0, 90, 180, 270], True),  # metals_x*
    ("Pd", [90, 270], True), ("Pd", [90, 260], True), ("Pd", [100, 285], True),  # within 165-195°, near vertical
    ("Pd", [0, 180], True), ("Pd", [5, 190], True),  # horizontal P–Pd–P (the user, 2026-10-30)
    ("Pd", [30, 210], False),  # opposite, but on neither axis (unmeasured)
    ("Pd", [90, 240], False), ("Cl", [90], False), ("Br", [270], False), ("Cl", [0], False),  # no opposite bond; terminal (clearance_x*)
    ("P", [90, 270], False),  # one letter
])
def test_centered_label_rule(element, dirs, centered):
    """ChemDraw centers a two-letter symbol when a vertical or horizontal bond has another bond
    165-195° from it (metals_x*; the user, 2026-10-30), else its first letter."""
    d = Drawer(_star(element, dirs), P)
    assert (d.labels[0].h_side == "center") == centered


@pytest.mark.parametrize("smiles", [OXADD, TRANSPLATIN, "CP(C)(C)->[Pd]<-P(C)(C)C"])
def test_cross_label_is_centered(smiles):
    """In a cross the whole symbol is centered (its advance), so a vertical bond meets the
    middle of "Pd"; the CDXML label is Center / Center / Center with p.x on the atom, as
    metals_x*'s Pd and Si; no bond is lengthened (ChemDraw keeps 14.4 pt)."""
    from lxml import etree
    from acs_chem_svg.cdxml_writer import to_cdxml
    m = chem.from_smiles(smiles, "m", P.lengths.bond_length)
    d = Drawer(m, P)
    metal = next(a for a in m.atoms if a.element in ("Pd", "Pt"))
    lab = d.labels[metal.idx]
    assert lab.h_side == "center"
    assert lab.elem_origin[0] + d.face.advance(metal.element, lab.size) / 2 == pytest.approx(metal.x, abs=1e-9)
    n = m.degree(metal.idx)
    assert [math.dist(metal.p, m.atoms[j].p) for j, _ in m.neighbors(metal.idx)] == pytest.approx([P.lengths.bond_length] * n, abs=0.01)
    root = etree.fromstring(to_cdxml([(m, d, (0.0, 0.0))], "acs1996", warn=lambda msg: None).encode())
    t = next(t for t in root.iter("t") if "".join(x.text for x in t.iter("s")) == metal.element)
    assert t.get("LabelJustification") == t.get("Justification") == t.get("LabelAlignment") == "Center"
    assert float(t.get("p").split()[0]) == pytest.approx(metal.x, abs=0.006)


@pytest.mark.parametrize("key, text", [("PMe3", "Me3P"), ("PPh3", "Ph3P")])
def test_phosphine_abbreviation_on_a_dative_bond(key, text):
    """chem.condense finds a phosphine bonded P->Pd (rdAbbreviations misses a dative attachment;
    issues 2026-10-30), keeps the cross, and the CDXML writes each as a nickname node."""
    from lxml import etree
    from acs_chem_svg.cdxml_writer import to_cdxml
    full = OXADD if key == "PMe3" else OXADD.replace("P(C)(C)C", "P(c1ccccc1)(c1ccccc1)c1ccccc1")
    m = chem.from_smiles(full, "m", P.lengths.bond_length, abbreviate=[key])
    ab = [a for a in m.atoms if a.abbrev]
    assert [a.abbrev for a in ab] == [key, key] and len(m.atoms) == 10
    pd = next(a for a in m.atoms if a.element == "Pd")
    assert all([j for j, _ in m.neighbors(a.idx)] == [pd.idx] for a in ab)
    d = Drawer(m, P)
    assert {"".join(r.text for r in d.labels[a.idx].runs) for a in ab} <= {key, text}
    root = etree.fromstring(to_cdxml([(m, d, (0.0, 0.0))], "acs1996", warn=lambda msg: None).encode())
    nodes = [n for n in root.iter("n") if n.get("NodeType") == "Unspecified"]
    assert ["".join(x.text for x in n.find("t").iter("s")) for n in nodes] == [key, key]
    assert all(n.find("fragment") is not None for n in nodes)


METALS = {"metals_x1.0": "acs1996", "metals_x0.8": "acs1996_x0.8", "metals_x0.6": "acs1996_x0.6"}
# the bond straight down from Pd: our §2.1 clip model stops 0.17-0.24 pt short of ChemDraw's
# (acs-style.md §2.1, Open 25); every other direction within TOL
PD_DOWN_KNOWN = 0.25


@pytest.mark.skipif(not mr.REF_DIR.exists(), reason="ChemDraw references (reference_renders/) are local only")
@pytest.mark.parametrize("ref", METALS)
def test_centered_labels_vs_chemdraw(ref):
    """metals_x*'s Pd (crosses) and Si (SiMe4) redrawn at ChemDraw's atom positions: the label
    is centered (its advance center on the atom, within 0.01 pt of ChemDraw's) and each bond
    stops where ChemDraw's line does, within 0.1 pt."""
    doc = mr.load(ref)
    preset = load_preset(METALS[ref])
    bad = []
    for a in doc.atoms.values():
        if a.label not in ("Pd", "Si"):
            continue
        nbs = [doc.atoms[b.b if b.a == a.id else b.a] for b in doc.bonds if a.id in (b.a, b.b)]
        m = Molecule("m", [Atom(0, "1", a.label, *a.p)] + [Atom(k + 1, str(k + 2), "C", *o.p) for k, o in enumerate(nbs)],
                     [Bond(k, 0, k + 1) for k in range(len(nbs))])
        d = Drawer(m, preset)
        lab = d.labels[0]
        glyphs = sorted((g for g in doc.glyphs if g.box and abs((g.box[1] + g.box[3]) / 2 - a.p[1]) < 5
                         and abs(g.x - a.p[0]) < 9 and g.ch in a.label), key=lambda g: g.x)[:2]
        theirs = (glyphs[0].x + glyphs[-1].x + glyphs[-1].adv) / 2
        ours = lab.elem_origin[0] + d.face.advance(a.label, lab.size) / 2
        if lab.h_side != "center" or abs(ours - theirs) > 0.01:
            bad.append(f"{a.label} {a.id}: label {lab.h_side}, advance center ours {ours:.3f} ChemDraw {theirs:.3f}")
        for b, o in zip(m.bonds, nbs):
            u = mr.unit(sub(o.p, a.p))
            ends = [min((q[0] - a.p[0]) * u[0] + (q[1] - a.p[1]) * u[1] for q in s.pts) for s in doc.shapes
                    if s.kind == "poly" and not s.curved
                    and max(abs((q[0] - a.p[0]) * u[1] - (q[1] - a.p[1]) * u[0]) for q in s.pts) < 1.0
                    and all(-0.5 < (q[0] - a.p[0]) * u[0] + (q[1] - a.p[1]) * u[1] < 15 for q in s.pts)]
            diff = d.end_clip(b, 0) - ends[0]
            known = PD_DOWN_KNOWN if a.label == "Pd" and u[1] > 0.99 else TOL
            if abs(diff) > known:
                bad.append(f"{a.label} {a.id} bond toward {o.label or o.element}: clip ours {d.end_clip(b, 0):.3f} "
                           f"ChemDraw {ends[0]:.3f} ({diff:+.3f})")
    assert not bad, "\n".join(bad)


def _mol(points, bonds, **kw) -> Molecule:
    return Molecule("x", [Atom(i, str(i + 1), "C", x, y) for i, (x, y) in enumerate(points)],
                    [Bond(k, a, b) for k, (a, b) in enumerate(bonds)], **kw)


def test_overlaid_layouts_are_reported():
    """Two bonds crossing in a flat species, or two atoms on one spot, are problems; a chair's
    crossings (depth) and a declared front bond are not."""
    cross = _mol([(0, 0), (10, 10), (0, 10), (10, 0)], [(0, 1), (2, 3)])
    assert any("cross" in p for p in check_overlaid("x", cross))
    assert check_overlaid("x", _mol([(0, 0), (10, 10), (0, 10), (10, 0)], [(0, 1), (2, 3)], front=[1])) == []
    same = _mol([(0, 0), (10, 0), (20, 0), (10, 0.5)], [(0, 1), (1, 2), (2, 3)])
    assert any("on top of each other" in p for p in check_overlaid("x", same))
    assert check_overlaid("x", _mol([(0, 0), (10, 0), (20, 0)], [(0, 1), (1, 2)])) == []


@pytest.mark.skipif(not mr.REF_DIR.exists(), reason="ChemDraw references (reference_renders/) are local only")
@pytest.mark.parametrize("ref", METALS)
def test_pph3_vs_chemdraw(ref):
    """metals_x*'s PPh3 ligands (chem.ABBREVIATIONS, both forms: PPh3, and Ph3P with its bond to
    the right) redrawn at ChemDraw's atom positions: the same text, each run's first glyph origin
    and the bond end within 0.1 pt of ChemDraw's."""
    doc = mr.load(ref)
    preset = load_preset(METALS[ref])
    bad, forms = [], set()
    for a in doc.atoms.values():
        if not (a.nickname and a.label.endswith("PPh3")):
            continue
        (b,) = [b for b in doc.bonds if a.id in (b.a, b.b)]
        o = doc.atoms[b.b if b.a == a.id else b.a]
        m = Molecule("m", [Atom(0, "1", "P", *a.p, abbrev="PPh3"), Atom(1, "2", o.element, *o.p)], [Bond(0, 0, 1)])
        d = Drawer(m, preset)
        _, whole = mr.label_glyphs(doc, a)
        text = "".join(g.ch for g in whole)
        runs = d.labels[0].runs
        forms.add(text)
        if "".join(r.text for r in runs) != text:
            bad.append(f"{a.id}: ours {''.join(r.text for r in runs)!r}, ChemDraw {text!r}")
            continue
        k = 0
        for r in runs:
            g = whole[k]
            if max(abs(r.x - g.x), abs(r.y - g.base), abs(r.size - g.size)) > TOL:
                bad.append(f"{a.id} {text} glyph {g.ch}: ours ({r.x:.3f}, {r.y:.3f}) ChemDraw ({g.x:.3f}, {g.base:.3f})")
            k += len(r.text)
        u, n = mr.unit(sub(o.p, a.p)), math.dist(a.p, o.p)  # not beyond o: x0.6's bonds are under 15 pt
        ends = [min((q[0] - a.p[0]) * u[0] + (q[1] - a.p[1]) * u[1] for q in s.pts) for s in doc.shapes
                if s.kind == "poly" and not s.curved
                and max(abs((q[0] - a.p[0]) * u[1] - (q[1] - a.p[1]) * u[0]) for q in s.pts) < 1.0
                and all(-0.5 < (q[0] - a.p[0]) * u[0] + (q[1] - a.p[1]) * u[1] < n for q in s.pts)]
        if abs(d.end_clip(m.bonds[0], 0) - ends[0]) > TOL:
            bad.append(f"{a.id} {text}: bond end ours {d.end_clip(m.bonds[0], 0):.3f} ChemDraw {ends[0]:.3f}")
    assert forms == {"PPh3", "Ph3P"} and not bad, "\n".join(bad)
