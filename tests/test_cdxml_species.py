"""CDXML export, phase 1 (acs-style.md §13): M1 species.

- DTD: the references are validated first, to show that a ChemDraw file passes; then our exports.
- RDKit round trip: every exported species reads back (Chem.MolsFromCDXML) as its spec's
  canonical SMILES, so stereo, charges and H counts survive.
- Reference diff: each M1 reference's templates at ChemDraw's coordinates, exported and compared
  node by node with the reference CDXML: elements, H counts, charges, bond orders and wedges,
  label text and justification, and within 0.1 pt label p and BoundingBox and the fragment box.
  The charge sign's center is hand-placed (CLAUDE.md), so it and a box it extends only have to
  be within CHARGE_TOL of ChemDraw's.
- Header templates (presets/cdxml/, scripts/cdxml_headers.py) against the presets and the other
  references at each scale.
"""

import math

import pytest
import yaml
from lxml import etree
from rdkit import Chem
from test_m1_vs_chemdraw import REFS

import measure_refs as mr
from acs_chem_svg.cdxml_writer import HEADER_DIR
from acs_chem_svg.figure import spec_cdxml
from acs_chem_svg.geometry import sub
from acs_chem_svg.preset import load_preset
from acs_chem_svg.spec import load_spec

from cdxml_helpers import (
    ABBREV_REFS, CHARGE_TOL, PARSER, SPECS, STYLE_ATTRS, TOL, _canonical, _check_structure, _text, _xy, _xy4,
    dtd_errors, exported, needs_dtd, needs_refs, reference_export, reference_pairs,
)


@needs_dtd
@pytest.mark.parametrize("ref", mr.all_names() if mr.REF_DIR.exists() else [])
def test_chemdraw_references_validate(ref):
    """Calibration: ChemDraw's own files pass, up to the DTD's known gaps."""
    assert dtd_errors(etree.parse(str(mr.REF_DIR / f"{ref}.cdxml"), PARSER)) == []


@needs_dtd
@pytest.mark.parametrize("spec", SPECS)
def test_spec_exports_validate(spec):
    _check_structure(exported(spec))


@needs_dtd
@pytest.mark.parametrize("ref", REFS + ABBREV_REFS)
def test_reference_exports_validate(ref):
    _check_structure(reference_export(ref))


@pytest.mark.parametrize("spec", SPECS)
def test_rdkit_round_trip(spec):
    species = yaml.safe_load((mr.REPO / spec).read_text())["species"]
    got = [Chem.MolToSmiles(m) for m in Chem.MolsFromCDXML(exported(spec))]
    assert got == [_canonical(sp["smiles"]) for sp in species.values()]


def test_round_trip_sees_wedges():
    """The round trip is sensitive to stereo: alanine with its hash turned to a wedge reads
    back as the other enantiomer."""
    text = exported("examples/m1/m1_smiles_x1.0.yaml")
    assert text.count('Display="WedgedHashBegin"') == 1
    flipped = Chem.MolToSmiles(Chem.MolsFromCDXML(text.replace("WedgedHashBegin", "WedgeBegin"))[2])
    assert flipped == _canonical("C[C@H](N)C(=O)O")


def test_only_atom_maps_warn():
    """Since phase 4b every drawn charge, dot and Lewis bond is written; only atom maps warn."""
    for spec in ("examples/m1/m1_charges_x1.0.yaml", "examples/m1/m1_lewis_assigned_x1.0.yaml"):
        msgs = []
        spec_cdxml(mr.REPO / spec, warn=msgs.append)
        assert all(m.endswith("atom maps not written") for m in msgs), msgs


def _species_frags(text: str) -> list:
    return etree.fromstring(text.encode()).find("page").findall("fragment")


def _symbols(frag, kind=None) -> list:
    return [g for g in frag.findall("graphic") if g.get("GraphicType") == "Symbol"
            and (kind is None or g.get("SymbolType") == kind)]


def test_charge_objects():
    """§13.4 on m1_charges_x1.0: circled signs are ChemDraw's CirclePlus / CircleMinus bound to
    the node, box dx cdxml.circled_box_dx; 2+ / 2- and the +3 label are one free <t> with
    InterpretChemically="no" and no symbol (ChemDraw rebinds any sign symbol), the node
    keeping its SMILES Charge; CH2's two
    unpaired electrons are two bound Electron symbols with Radical="Singlet", as ChemDraw wrote
    them."""
    text = exported("examples/m1/m1_charges_x1.0.yaml")
    root = etree.fromstring(text.encode())
    frags = dict(zip(load_spec(mr.REPO / "examples/m1/m1_charges_x1.0.yaml").species, _species_frags(text)))
    cx = load_preset("acs1996").cdxml
    nitro = _symbols(frags["nitromethane"])
    assert sorted(g.get("SymbolType") for g in nitro) == ["CircleMinus", "CirclePlus"]
    for g in nitro:
        x0, _, x1, _ = _xy4(g.get("BoundingBox"))
        assert round(x0 - x1, 3) == cx.circled_box_dx and g.find("represent").get("attribute") == "Charge"
    assert [g.get("SymbolType") for g in _symbols(frags["plain"])] == ["Plus", "Minus"] or \
        sorted(g.get("SymbolType") for g in _symbols(frags["plain"])) == ["Minus", "Plus"]
    texts = root.find("page").findall("t")
    for name, charge in (("iron", "2"), ("oxide", "-2"), ("iron_iii", "2")):
        (n,) = frags[name].findall("n")
        assert n.get("Charge") == charge and not _symbols(frags[name])
    assert sorted("".join(x.text for x in t.iter("s")) for t in texts) == ["+3", "2+", "2\u2212"]
    assert {t.get("InterpretChemically") for t in texts} == {"no"}
    (ch2,) = frags["methylene"].findall("n")
    assert ch2.get("Radical") == "Singlet"
    assert [g.find("represent").get("object") for g in _symbols(frags["methylene"], "Electron")] == [ch2.get("id")] * 2


def test_lewis_objects():
    """§13.4: every atom labeled; dots form: bonds Visible="no" and each shared column a free
    LonePair; lines form: visible bonds, no shared pairs; assigned pairs free LonePairs, an odd
    assigned electron a filled oval (ChemDraw makes any Electron near a node a radical);
    oxidation states free <t>s, Charge unchanged."""
    dots = _species_frags(exported("examples/m1/m1_lewis_x1.0.yaml"))
    lines = _species_frags(exported("examples/m1/m1_lewis_lines_x1.0.yaml"))
    for fd, fl in zip(dots, lines):
        assert all(n.find("t") is not None for n in fd.findall("n"))
        assert all(b.get("Visible") == "no" for b in fd.findall("b"))
        assert all(b.get("Visible") is None for b in fl.findall("b"))
        free = [g for g in _symbols(fd, "LonePair") if g.find("represent") is None]
        assert len(free) == sum(int(b.get("Order", "1")) for b in fd.findall("b"))
        assert all(g.find("represent") is not None for g in _symbols(fl))
    assigned = etree.fromstring(exported("examples/m1/m1_lewis_assigned_x1.0.yaml").encode())
    frags = dict(zip(load_spec(mr.REPO / "examples/m1/m1_lewis_assigned_x1.0.yaml").species,
                     assigned.find("page").findall("fragment")))
    assert all(n.get("Radical") is None for f in frags.values() for n in f.findall("n"))
    assert not _symbols(frags["hcn"], "Electron")
    ovals = [g for g in frags["hcn"].findall("graphic") if g.get("GraphicType") == "Oval"]
    assert len(ovals) == 4 and all(g.get("OvalType") == "Circle Filled" for g in ovals)
    assert all(n.get("Charge") is None for n in frags["methanol"].findall("n"))
    assert not [g for g in _symbols(frags["methanol"]) if g.get("SymbolType") in ("Plus", "Minus")]
    texts = sorted("".join(x.text for x in t.iter("s")) for t in assigned.find("page").findall("t"))
    assert texts == sorted(["2\u2212"] * 3 + ["+"] * 6)  # CH3OH: C, O 2-, four H +; H2O: O 2-, two H +


@needs_refs
@pytest.mark.parametrize("scale", ["x1.0", "x0.8", "x0.6"])
def test_charges_lewis_reference_keeps_our_lewis(scale, tmp_path):
    """charges_lewis_x* is the user's edit of cdxml_gate.py --lewis-draft (2026-10-28): its Lewis
    rows kept every atom and dot where the draft put them. Our draft today still puts them there
    (nodes and LonePair / Electron dot centers within 0.015 pt, the files' two decimals), so the reference approves them."""
    import cdxml_gate
    ours = etree.parse(str(cdxml_gate.lewis_draft(tmp_path)[["x1.0", "x0.8", "x0.6"].index(scale)])).getroot()
    ref = etree.parse(str(mr.REF_DIR / f"charges_lewis_{scale}.cdxml")).getroot()
    skip = len(cdxml_gate.DRAFT_CHARGES)  # row 1: the user's charges

    def items(root):  # per structure, nodes and dot centers relative to its node centroid (rows may shift)
        out = []
        for f in root.find("page").findall("fragment")[skip:]:
            pts_ = [_xy(n.get("p")) for n in f.findall("n")]
            o = (sum(q[0] for q in pts_) / len(pts_), sum(q[1] for q in pts_) / len(pts_))
            for g in f.findall("graphic"):
                if g.get("SymbolType") == "LonePair":
                    b = _xy4(g.get("BoundingBox"))
                    pts_ += [b[:2], b[2:]]
                elif g.get("SymbolType") == "Electron" or g.get("GraphicType") == "Oval":  # odd electrons: ovals now
                    pts_.append(_xy(g.get("BoundingBox")))
            out.append([sub(q, o) for q in pts_])
        return out
    got_all, want_all = items(ours), items(ref)
    assert len(got_all) == len(want_all)
    for got, want in zip(got_all, want_all):
        assert len(got) == len(want)
        assert all(min(math.dist(w, g) for g in got) < 0.015 for w in want)


@needs_refs
@pytest.mark.parametrize("ref", REFS + ABBREV_REFS)
def test_reference_diff(ref):
    bad, nums = reference_pairs(ref)
    tol = {"charge center": CHARGE_TOL, "fragment box (with a charge)": CHARGE_TOL}
    bad += [f"{c} [{k}]: ChemDraw {a} ours {b} (diff {b - a:+.3f})" for c, k, a, b in nums if abs(b - a) > tol.get(c, TOL)]
    assert not bad, "\n".join(bad)
    assert any(c == "label p" for c, *_ in nums) or ref == "double_bonds_x1.0"
    if ref.startswith("abbrev_"):  # all 21 keys, both forms
        assert sum(c == "nickname atoms" for c, *_ in nums) > 300


def _header(preset: str):
    return etree.parse(str(HEADER_DIR / f"{preset}.xml"), PARSER).getroot()


@pytest.mark.parametrize("preset", ["acs1996", "acs1996_x0.8", "acs1996_x0.6"])
def test_header_matches_preset(preset):
    root, p = _header(preset), load_preset(preset)
    assert float(root.get("BondLength")) == pytest.approx(p.lengths.bond_length, abs=0.005)
    assert float(root.get("LineWidth")) == p.lengths.line_width
    assert float(root.get("BoldWidth")) == pytest.approx(p.lengths.bold_width)
    assert float(root.get("HashSpacing")) == pytest.approx(p.lengths.hash_spacing)
    assert float(root.get("MarginWidth")) == pytest.approx(p.lengths.margin_width)
    assert float(root.get("BondSpacing")) / 100 == p.ratios.bond_spacing_fraction
    assert float(root.get("LabelSize")) == p.text.label_size
    assert float(root.get("CaptionSize")) == p.text.caption_size
    assert float(root.get("ChainAngle")) == p.angles.chain_angle
    fonts = {f.get("id"): f.get("name") for f in root.iter("font")}
    assert fonts[root.get("LabelFont")] == fonts[root.get("CaptionFont")] == p.font_family
    assert all(root.get(k) is None for k in ("CreationProgram", "Name", "BoundingBox"))
    assert root.find("page").get("id") is None
    assert len(root.find("page")) == 0


@needs_refs
@pytest.mark.parametrize("preset,scale", [("acs1996", "1.0"), ("acs1996_x0.8", "0.8"), ("acs1996_x0.6", "0.6")])
def test_header_matches_other_references(preset, scale):
    """Every reference at that scale has the same style settings as the template."""
    root = _header(preset)
    for ref in [n for n in mr.all_names() if n.endswith(f"_x{scale}")] + (["double_bonds_2"] if scale == "1.0" else []):
        r = etree.parse(str(mr.REF_DIR / f"{ref}.cdxml"), PARSER).getroot()
        assert {a: r.get(a) for a in STYLE_ATTRS} == {a: root.get(a) for a in STYLE_ATTRS}, ref


@needs_refs
def test_headers_are_current():
    """presets/cdxml/ is what scripts/cdxml_headers.py writes."""
    import cdxml_headers
    for preset, ref in cdxml_headers.SOURCES.items():
        assert (HEADER_DIR / f"{preset}.xml").read_bytes() == cdxml_headers.header(ref), preset


def test_charge_box_dx_follows_label_size():
    """cdxml.charge_box_dx = 0.075 x LabelSize^2 in 19 of 20 ChemDraw charge symbols
    (measure_refs.py --cdxml-symbols); the presets hold the measured values."""
    for name in ("acs1996", "acs1996_x0.8", "acs1996_x0.6"):
        p = load_preset(name)
        assert p.cdxml.charge_box_dx == pytest.approx(0.075 * p.text.label_size ** 2)


def test_creation_program_is_ours():
    root = etree.fromstring(exported(SPECS[0]).encode(), PARSER)
    assert root.get("CreationProgram").startswith("acs-chem-svg ")


def test_page_grows_to_fit():
    """m1_smiles_x1.0 (706 pt wide) gets two pages across, as ChemDraw expands it; m1_charges fits one."""
    wide = etree.fromstring(exported("examples/m1/m1_smiles_x1.0.yaml").encode(), PARSER).find("page")
    assert (wide.get("BoundingBox"), wide.get("WidthPages"), wide.get("HeightPages")) == ("0 0 1080 720", "2", "1")
    one = etree.fromstring(exported("examples/m1/m1_charges_x1.0.yaml").encode(), PARSER).find("page")
    assert (one.get("BoundingBox"), one.get("WidthPages")) == ("0 0 540 720", "1")


def test_vertical_bond_label_is_left_justified():
    """Alanine's OH sits on a vertical bond (to ~1e-4 in the RDKit layout): OH, not HO, as
    ChemDraw re-saved it."""
    root = etree.fromstring(exported("examples/m1/m1_smiles_x1.0.yaml").encode(), PARSER)
    alanine = root.find("page").findall("fragment")[2]
    oh = [t for t in alanine.iter("t") if _text(t) == "OH"]
    assert len(oh) == 1 and oh[0].get("LabelJustification") == "Left"
