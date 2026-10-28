"""CDXML export, phase 1 (acs-style.md §13): M1 species.

- DTD: reference_renders/cdxml/CDXML.dtd (Revvity's file, fetched once from the DOCTYPE URL;
  ChemDraw.app doesn't ship it). Local only, like the references: those tests skip without it. It dates from 2003, so libxml2 needs one enumeration loosened (AS lists +, - and ?, which
  aren't name tokens), and ChemDraw 26 writes two attributes it doesn't declare (DTD_GAPS). The
  references are validated first, to show that a ChemDraw file passes; then our exports.
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

import io
import re
from functools import lru_cache

import pytest
import yaml
from compare import render_reference
from lxml import etree
from rdkit import Chem
from test_m1_vs_chemdraw import REFS

import measure_refs as mr
from acs_chem_svg.cdxml_writer import HEADER_DIR, CdxmlUnsupported, spec_cdxml, to_cdxml
from acs_chem_svg.draw_mol import Drawer
from acs_chem_svg.preset import load_preset

DTD_PATH = mr.REPO / "reference_renders" / "cdxml" / "CDXML.dtd"  # Revvity's: local only, never committed
DTD_GAPS = {"No declaration for attribute ShowPrimeLabels of element CDXML",
            "No declaration for attribute IsNickname of element t"}
TOL = 0.1
CHARGE_TOL = 1.0  # hand-placed sign centers (preset charge.labeled_offset_*: round defaults)
PAIR_TOL = 0.01
SPECS = [f"examples/m1_{k}_x{s}.yaml" for k in ("smiles", "charges") for s in ("1.0", "0.8", "0.6")]
STYLE_ATTRS = ("BondLength", "LineWidth", "BoldWidth", "BondSpacing", "HashSpacing", "MarginWidth", "LabelSize",
               "CaptionSize", "LabelFont", "CaptionFont", "LabelFace", "ChainAngle")
PARSER = etree.XMLParser(load_dtd=False, no_network=True)
needs_dtd = pytest.mark.skipif(not DTD_PATH.exists(), reason="CDXML.dtd (reference_renders/cdxml/) is local only")
needs_refs = pytest.mark.skipif(not mr.REF_DIR.exists(), reason="ChemDraw references (reference_renders/) are local only")


@lru_cache(maxsize=None)
def dtd() -> etree.DTD:
    text = DTD_PATH.read_text()

    def loosen(m):  # an enumeration with a value that isn't a name token -> CDATA
        return m[0] if all(re.fullmatch(r"[\w.:-]+", v.strip()) for v in m[1].split("|")) else "CDATA"
    return etree.DTD(io.StringIO(re.sub(r"\(([^()]*\|[^()]*)\)(?=\s+(?:\"|#))", loosen, text)))


def dtd_errors(root) -> list[str]:
    d = dtd()
    d.validate(root)
    return [e.message for e in d.error_log if e.message not in DTD_GAPS]


@lru_cache(maxsize=None)
def exported(spec: str) -> str:
    return spec_cdxml(mr.REPO / spec, warn=lambda m: None)


@lru_cache(maxsize=None)
def reference_export(ref: str) -> str:
    """The reference's templates exported at ChemDraw's CDXML coordinates (the templates hold
    its SVG frame, svg = cdxml + Doc.offset)."""
    r = render_reference(ref)
    off = (-r.ref.offset[0], -r.ref.offset[1])
    return to_cdxml([(m, Drawer(m, r.preset), off) for m in r.molecules], r.preset.name, f"{ref}.cdxml",
                    warn=lambda m: None)


# --------------------------------------------------------------------------- DTD


@needs_dtd
@pytest.mark.parametrize("ref", mr.all_names() if mr.REF_DIR.exists() else [])
def test_chemdraw_references_validate(ref):
    """Calibration: ChemDraw's own files pass, up to the DTD's known gaps."""
    assert dtd_errors(etree.parse(str(mr.REF_DIR / f"{ref}.cdxml"), PARSER)) == []


def _check_structure(text: str) -> None:
    root = etree.fromstring(text.encode(), PARSER)
    assert dtd_errors(root) == []
    # object ids; font ids are their own namespace (clearance_x1.0 has font 21 and bond 21)
    ids = [el.get("id") for el in root.find("page").iter() if el.get("id")]
    assert len(ids) == len(set(ids)), "duplicate ids"
    nodes = {n.get("id") for n in root.iter("n")}
    for b in root.iter("b"):
        assert {b.get("B"), b.get("E")} <= nodes, f"bond {b.get('id')} points to a missing node"
    for r in root.iter("represent"):
        assert r.get("object") in nodes, f"represent points to missing node {r.get('object')}"


@needs_dtd
@pytest.mark.parametrize("spec", SPECS)
def test_spec_exports_validate(spec):
    _check_structure(exported(spec))


@needs_dtd
@pytest.mark.parametrize("ref", REFS)
def test_reference_exports_validate(ref):
    _check_structure(reference_export(ref))


# --------------------------------------------------------------------------- RDKit round trip


def _canonical(smiles: str) -> str:
    m = Chem.MolFromSmiles(smiles)
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(m)


@pytest.mark.parametrize("spec", SPECS)
def test_rdkit_round_trip(spec):
    species = yaml.safe_load((mr.REPO / spec).read_text())["species"]
    got = [Chem.MolToSmiles(m) for m in Chem.MolsFromCDXML(exported(spec))]
    assert got == [_canonical(sp["smiles"]) for sp in species.values()]


def test_round_trip_sees_wedges():
    """The round trip is sensitive to stereo: alanine with its hash turned to a wedge reads
    back as the other enantiomer."""
    text = exported("examples/m1_smiles_x1.0.yaml")
    assert text.count('Display="WedgedHashBegin"') == 1
    flipped = Chem.MolToSmiles(Chem.MolsFromCDXML(text.replace("WedgedHashBegin", "WedgeBegin"))[2])
    assert flipped == _canonical("C[C@H](N)C(=O)O")


@pytest.mark.parametrize("spec", ["examples/m1_abbrev_x1.0.yaml", "examples/m1_lewis_x1.0.yaml",
                                  "examples/m2_suzuki_x1.0.yaml"])
def test_unsupported_specs_are_refused(spec):
    with pytest.raises(CdxmlUnsupported):
        spec_cdxml(mr.REPO / spec, warn=lambda m: None)


def test_omissions_warn():
    msgs = []
    spec_cdxml(mr.REPO / "examples/m1_charges_x1.0.yaml", warn=msgs.append)
    assert "enolate: lone-pair dots omitted from CDXML" in msgs
    assert "nitromethane: circled charge written as a plain sign" in msgs
    assert "iron_iii: charge label '+3' written as the formal charge +2" in msgs
    assert len(msgs) == len(set(msgs))


# --------------------------------------------------------------------------- reference diff


def _xy(s):
    return tuple(float(v) for v in s.split()[:2])


def _xy4(s):
    return tuple(float(v) for v in s.split())


def _nodes(root):
    """Top-level fragment nodes (not nickname internals) by id."""
    return {n.get("id"): n for f in root.iter("fragment") if f.getparent().tag == "page" for n in f.findall("n")}


def _bonds(root, nodes):
    """frozenset(node ids) -> (Order, display, narrow-end node id); WedgeEnd -> Begin from E."""
    out = {}
    for b in root.iter("b"):
        if b.get("B") not in nodes:
            continue
        disp = b.get("Display", "")
        narrow = None
        if disp:
            narrow = b.get("B") if disp.endswith("Begin") else b.get("E")
            disp = disp.replace("End", "Begin")
        out[frozenset((b.get("B"), b.get("E")))] = (b.get("Order", "1"), disp, narrow)
    return out


def _charges(root):
    return {g.find("represent").get("object"): (g.get("SymbolType"), [float(v) for v in g.get("BoundingBox").split()])
            for g in root.iter("graphic") if g.get("GraphicType") == "Symbol" and g.find("represent") is not None}


def _text(t):
    return "".join(s.text or "" for s in t.iter("s"))


def reference_pairs(ref: str) -> tuple[list[str], list[tuple[str, str, float, float]]]:
    """(mismatches, [(category, key, ChemDraw, ours)] for numbers)."""
    theirs = etree.parse(str(mr.REF_DIR / f"{ref}.cdxml"), PARSER).getroot()
    ours = etree.fromstring(reference_export(ref).encode(), PARSER)
    tn, on = _nodes(theirs), _nodes(ours)
    to_ref = {}
    bad, nums = [], []
    for oid, o in on.items():
        p = _xy(o.get("p"))
        match = [tid for tid, t in tn.items() if max(abs(a - b) for a, b in zip(_xy(t.get("p")), p)) < PAIR_TOL]
        if len(match) != 1:
            bad.append(f"node at {p}: {len(match)} ChemDraw nodes")
            continue
        to_ref[oid] = match[0]
    for oid, tid in to_ref.items():
        o, t = on[oid], tn[tid]
        key = f"node {tid}"
        for attr, default in (("Element", "6"), ("Charge", "0")):
            if o.get(attr, default) != t.get(attr, default):
                bad.append(f"{key} {attr}: ChemDraw {t.get(attr, default)} ours {o.get(attr, default)}")
        ot, tt = o.find("t"), t.find("t")
        if (ot is None) != (tt is None):
            bad.append(f"{key}: label {'missing' if ot is None else 'extra'}")
            continue
        if tt is None:
            continue
        if o.get("NumHydrogens") != t.get("NumHydrogens"):
            bad.append(f"{key} NumHydrogens: ChemDraw {t.get('NumHydrogens')} ours {o.get('NumHydrogens')}")
        for attr in ("LabelJustification", "Justification"):
            if ot.get(attr) != tt.get(attr):
                bad.append(f"{key} {attr}: ChemDraw {tt.get(attr)} ours {ot.get(attr)}")
        if _text(ot) != _text(tt):
            bad.append(f"{key} text: ChemDraw {_text(tt)} ours {_text(ot)}")
        for k, (a, b) in enumerate(zip(_xy(tt.get("p")), _xy(ot.get("p")))):
            nums.append(("label p", f"{key} {'xy'[k]} {_text(tt)}", a, b))
        for k, (a, b) in enumerate(zip(_xy4(tt.get("BoundingBox")), _xy4(ot.get("BoundingBox")))):
            nums.append(("label box", f"{key} [{k}] {_text(tt)}", a, b))
    # every node of a fragment we drew is paired (templates skip only nickname fragments)
    drawn = {tn[t].getparent() for t in to_ref.values()}
    missing = [tid for tid, t in tn.items() if t.getparent() in drawn and tid not in to_ref.values()]
    if missing:
        bad.append(f"ChemDraw nodes not exported: {missing}")
    tb = _bonds(theirs, tn)
    for pair, (order, disp, narrow) in _bonds(ours, on).items():
        tpair = frozenset(to_ref.get(i) for i in pair)
        if tpair not in tb:
            bad.append(f"bond {sorted(tpair)}: not in ChemDraw's file")
            continue
        t_order, t_disp, t_narrow = tb[tpair]
        if (order, disp) != (t_order, t_disp) or (narrow and to_ref[narrow] != t_narrow):
            bad.append(f"bond {sorted(tpair)}: ChemDraw {t_order} {t_disp} from {t_narrow}, "
                       f"ours {order} {disp} from {to_ref.get(narrow)}")
    for f in ours.iter("fragment"):  # the drawn ink box; a hand-placed charge sign moves it
        t = tn[to_ref[f.find("n").get("id")]].getparent()
        cat = "fragment box (with a charge)" if f.find("graphic") is not None else "fragment box"
        for k, (a, b) in enumerate(zip(_xy4(t.get("BoundingBox")), _xy4(f.get("BoundingBox")))):
            nums.append((cat, f"fragment {t.get('id')} [{k}]", a, b))
    tc, oc = _charges(theirs), _charges(ours)
    if sorted(tc) != sorted(to_ref[i] for i in oc):
        bad.append(f"charge symbols on ChemDraw nodes {sorted(tc)}, ours {sorted(to_ref[i] for i in oc)}")
    for oid, (sym, box) in oc.items():
        tsym, tbox = tc.get(to_ref[oid], (None, None))
        if tsym != sym:
            bad.append(f"charge on node {to_ref[oid]}: ChemDraw {tsym} ours {sym}")
            continue
        for k in range(2):
            nums.append(("charge center", f"node {to_ref[oid]} {'xy'[k]}", tbox[k], box[k]))
        nums.append(("charge box dx", f"node {to_ref[oid]}", tbox[0] - tbox[2], box[0] - box[2]))
    return bad, nums


@needs_refs
@pytest.mark.parametrize("ref", REFS)
def test_reference_diff(ref):
    bad, nums = reference_pairs(ref)
    tol = {"charge center": CHARGE_TOL, "fragment box (with a charge)": CHARGE_TOL}
    bad += [f"{c} [{k}]: ChemDraw {a} ours {b} (diff {b - a:+.3f})" for c, k, a, b in nums if abs(b - a) > tol.get(c, TOL)]
    assert not bad, "\n".join(bad)
    assert any(c == "label p" for c, *_ in nums) or ref == "double_bonds_x1.0"


# --------------------------------------------------------------------------- header templates


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
    wide = etree.fromstring(exported("examples/m1_smiles_x1.0.yaml").encode(), PARSER).find("page")
    assert (wide.get("BoundingBox"), wide.get("WidthPages"), wide.get("HeightPages")) == ("0 0 1080 720", "2", "1")
    one = etree.fromstring(exported("examples/m1_charges_x1.0.yaml").encode(), PARSER).find("page")
    assert (one.get("BoundingBox"), one.get("WidthPages")) == ("0 0 540 720", "1")


def test_vertical_bond_label_is_left_justified():
    """Alanine's OH sits on a vertical bond (to ~1e-4 in the RDKit layout): OH, not HO, as
    ChemDraw re-saved it."""
    root = etree.fromstring(exported("examples/m1_smiles_x1.0.yaml").encode(), PARSER)
    alanine = root.find("page").findall("fragment")[2]
    oh = [t for t in alanine.iter("t") if _text(t) == "OH"]
    assert len(oh) == 1 and oh[0].get("LabelJustification") == "Left"
