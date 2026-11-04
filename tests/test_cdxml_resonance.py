"""CDXML export, phase 4a (acs-style.md §13.3): resonance arrows and brackets (resonance_x*),
wrapped schemes (scheme_wrap_x*), abbreviation nodes; and every example spec (ALL_SPECS).
"""

import pytest
from lxml import etree
from rdkit import Chem

import measure_refs as mr
from acs_chem_svg.draw_mol import Drawer
from acs_chem_svg.preset import load_preset
from acs_chem_svg.spec import load_spec

from cdxml_helpers import (
    PARSER, TOL, _check_structure, _floats, _spec_smiles, _text, check_scheme_round_trip, exported,
    needs_dtd, needs_refs, reference_pairs, reference_scheme_export, scheme_pairs,
)


RESONANCE_REFS = [f"resonance_x{s}" for s in ("1.0", "0.8", "0.6")]
RESONANCE_KNOWN = {
    # the C+ signs on unlabeled carbons are hand-placed (charge.unlabeled_distance, a round
    # default): 1.52 below ChemDraw's on the two outer contributors at x0.8 and x0.6
    "charge center": 1.6, "fragment box (with a charge)": 1.6,
    # ChemDraw writes one of resonance_x1.0's three sign boxes mirrored (dx -7.5; preset comment)
    ("resonance_x1.0", "charge box dx", "node 31"): 15.01,
}
ALL_SPECS = sorted(str(p.relative_to(mr.REPO)) for p in (mr.REPO / "examples").glob("m*/*.yaml"))


@needs_refs
@pytest.mark.parametrize("ref", RESONANCE_REFS)
def test_resonance_reference_diff(ref):
    """resonance_x* rebuilt (compare.reference_cdx_scheme: templates with their NO2 nicknames,
    ↔ arrows on ChemDraw's Head3D / Tail3D, brackets on its bars, its steps) and compared with
    its file: arrows (ArrowType Resonance, Full at both ends, box), brackets (BoundingBox,
    LipSize), one step per ↔ in one <scheme>, and the nodes, charges and NO2 nicknames (charged
    inner N+ and O-) as for the M1 references."""
    bad, nums = scheme_pairs(ref)
    nbad, nnums = reference_pairs(ref, reference_scheme_export(ref))
    bad += nbad
    for cat, key, a, b in [(c, k, a, b) for c, (_, k), a, b in nums] + nnums:
        if abs(b - a) > RESONANCE_KNOWN.get((ref, cat, key), RESONANCE_KNOWN.get(cat, TOL)):
            bad.append(f"{cat} {key}: ChemDraw {a} ours {b} (diff {b - a:+.3f})")
    assert not bad, "\n".join(bad)
    cats = {c for c, *_ in nums + nnums}
    assert {"arrow points", "arrow box", "bracket box", "nickname atoms"} <= cats
    root = etree.fromstring(reference_scheme_export(ref).encode(), PARSER)
    assert [len(s.findall("step")) for s in root.iter("scheme")] == [2]
    assert {g.get("ArrowType") for g in root.iter("graphic") if g.get("SupersededBy")} == {"Resonance"}


@needs_dtd
@pytest.mark.parametrize("ref", RESONANCE_REFS)
def test_resonance_reference_exports_validate(ref):
    _check_structure(reference_scheme_export(ref))


def test_resonance_objects():
    """m4_resonance_x1.0: two ↔ arrows (Full at both ends), one <step> each from contributor k to
    k + 1, and two free bracket <graphic>s: no <bracketedgroup>, no objecttag (user decision
    2026-10-28); "[" lists its bar's bottom first, "]" its top."""
    root = etree.fromstring(exported("examples/m4/m4_resonance_x1.0.yaml").encode(), PARSER)
    page = root.find("page")
    arrows = list(page.iter("arrow"))
    assert len(arrows) == 2 and all(a.get("ArrowheadHead") == a.get("ArrowheadTail") == "Full" for a in arrows)
    frags = [f.get("id") for f in page.findall("fragment")]
    steps = [(st.get("ReactionStepReactants"), st.get("ReactionStepProducts")) for st in root.iter("step")]
    assert steps == [(frags[0], frags[1]), (frags[1], frags[2])]
    brackets = [g for g in page.findall("graphic") if g.get("GraphicType") == "Bracket"]
    assert [g.get("LipSize") for g in brackets] == ["60", "60"]
    (x0, b0, _, t0), (x1, t1, _, b1) = (_floats(g.get("BoundingBox")) for g in brackets)
    assert x0 < x1 and b0 > t0 and b1 > t1
    assert root.find(".//bracketedgroup") is None and root.find(".//objecttag") is None


def test_wrapped_sequence_is_one_scheme():
    """m4_wrap_x1.0: its two rows are one sequence, so one <scheme>; the step that starts row 2
    cites its arrow, which lies below row 1's."""
    root = etree.fromstring(exported("examples/m4/m4_wrap_x1.0.yaml").encode(), PARSER)
    page = root.find("page")
    assert [len(s.findall("step")) for s in root.iter("scheme")] == [len(load_spec(mr.REPO / "examples/m4/m4_wrap_x1.0.yaml").steps)]
    graphics = {g.get("id"): g for g in page.findall("graphic")}
    ys = [_floats(graphics[st.get("ReactionStepArrows")].get("BoundingBox"))[1] for st in root.iter("step")]
    assert len(set(round(y) for y in ys)) == 2 and ys == sorted(ys)


def test_abbreviation_nodes():
    """m1_abbrev_x1.0: each abbreviation is a nickname node (no Element) holding ChemDraw's
    group, its <t> the key text in both forms: MeO2C right-justified is CO2Me; Ph and TMS are
    NodeType Nickname with IsNickname, the formula-like ones NodeType Fragment."""
    root = etree.fromstring(exported("examples/m1/m1_abbrev_x1.0.yaml").encode(), PARSER)
    nodes = {_text(n.find("t")): n for n in root.iter("n") if n.get("NodeType") in ("Fragment", "Nickname")}
    assert set(nodes) == {"CO2Me", "OMe", "CN", "NO2", "CF3", "TMS", "Ph", "NHAc", "SMe"}
    for key, n in nodes.items():
        assert n.get("Element") is None and n.get("NeedsClean") == "yes"
        assert (n.get("NodeType") == "Nickname") == (key in ("TMS", "Ph")) == (n.find("t").get("IsNickname") == "yes")
        inner = n.find("fragment")
        ecp = [m for m in inner.findall("n") if m.get("NodeType") == "ExternalConnectionPoint"]
        assert len(ecp) == 1 and any(_floats(m.get("p")) == _floats(n.get("p")) for m in inner.findall("n"))
    assert nodes["CO2Me"].find("t").get("Justification") == "Right"  # drawn MeO2C


@needs_dtd
@pytest.mark.parametrize("spec", ALL_SPECS)
def test_all_spec_exports_validate(spec):
    """Every example spec exports and passes the DTD, id and citation checks."""
    _check_structure(exported(spec))


@pytest.mark.parametrize("spec", ALL_SPECS)
def test_all_spec_round_trips(spec):
    """Every example spec reads back with RDKit as its species (and steps) and raises no
    CdxmlUnsupported."""
    if "_cycle_" in spec:
        pytest.skip("catalytic cycles: test_cdxml_cycle.py::test_cycle_round_trip")
    s = load_spec(mr.REPO / spec)
    if any(s.registry.get(t) is not None and s.registry.get(t).kind == "text" for t in s.labels.values()):
        pytest.skip("text labels have no chemistry: test_abbrev2_example_round_trip")
    if s.steps:
        check_scheme_round_trip(spec)
    else:
        got = sorted(Chem.MolToSmiles(m) for m in Chem.MolsFromCDXML(exported(spec)))
        assert got == sorted(_spec_smiles(mr.REPO / spec).values())


@pytest.mark.parametrize("spec", ALL_SPECS)
def test_signs_sit_nearest_their_atom(spec):
    """ChemDraw binds each sign symbol to the nearest atom on open and resets its charge from it
    (phase 4b gate), so every drawn charge is nearer its own atom than any other."""
    from acs_chem_svg.spec import molecules
    sp = load_spec(mr.REPO / spec)
    for mol in molecules(sp, (mr.REPO / spec).parent):
        d = Drawer(mol, load_preset(sp.preset))
        for a in mol.atoms:
            if d.charge_sign(a.idx) is not None:
                assert d._nearest_atom(d.charge_center(a.idx)) == a.idx, (mol.name, a.id)
