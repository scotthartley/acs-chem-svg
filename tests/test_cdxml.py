"""CDXML export (acs-style.md §13): phase 1, M1 species; phase 2 (§13.1), M2 schemes; phase 3
(§13.2), M3 mechanisms (mech_basic_x* round trip: curves, dot symbols; M3 specs).

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
- Phase 2, scheme round trip: each scheme reference rebuilt as a CdxScheme at ChemDraw's
  coordinates (templates; arrows from its Head3D / Tail3D; captions set by richtext at its glyph
  origins; formula fragments such as K2CO3 as captions; its steps) and compared with its file:
  arrow points, boxes and head attributes, caption p, box, LineStarts and text, and step
  membership by paired ids, all within 0.1 pt except KNOWN. Abbreviated species (scheme_multi's
  OMe / OBn) aren't written in phase 2: they're left out of the export and of the step check.
  equilibrium_x*'s A and B (no templates) are one-atom stand-ins at their node positions.
- Phase 2 specs: structure (DTD, ids, step citations), RDKit (species, each step's reactants and
  products by fragment id, ReactionsFromCDXMLBlock on one-step specs) and refusals.
"""

import io
import math
import re
from functools import lru_cache

import pytest
import yaml
from compare import preset_for, reference_cdx_scheme, render_reference
from lxml import etree
from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions
from test_m1_vs_chemdraw import REFS

import measure_refs as mr
from acs_chem_svg.cdxml_writer import HEADER_DIR, CdxmlUnsupported, spec_cdxml, to_cdxml
from acs_chem_svg.draw_mol import Drawer
from acs_chem_svg.geometry import sub
from acs_chem_svg.preset import load_preset
from acs_chem_svg.spec import load_spec, resolve_molfile

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


STEP_KINDS = {"ReactionStepReactants": ("fragment",), "ReactionStepProducts": ("fragment",),
              "ReactionStepArrows": ("graphic",),  # above / below: captions, and species a mechanism step adds
              "ReactionStepObjectsAboveArrow": ("t", "fragment"), "ReactionStepObjectsBelowArrow": ("t", "fragment")}


def _check_structure(text: str) -> None:
    root = etree.fromstring(text.encode(), PARSER)
    assert dtd_errors(root) == []
    # object ids; font ids are their own namespace (clearance_x1.0 has font 21 and bond 21)
    els = {el.get("id"): el for el in root.find("page").iter() if el.get("id")}
    ids = [el.get("id") for el in root.find("page").iter() if el.get("id")]
    assert len(ids) == len(set(ids)), "duplicate ids"
    nodes = {n.get("id") for n in root.iter("n")}
    for b in root.iter("b"):
        assert {b.get("B"), b.get("E")} <= nodes, f"bond {b.get('id')} points to a missing node"
    for r in root.iter("represent"):
        assert r.get("object") in nodes, f"represent points to missing node {r.get('object')}"
    for st in root.iter("step"):  # every id a step cites exists and is of the right kind
        for attr, tags in STEP_KINDS.items():
            for i in (st.get(attr) or "").split():
                assert i in els and els[i].tag in tags, f"step {st.get('id')} {attr}: {i} is not a <{tags}>"
    for g in root.iter("graphic"):
        if g.get("SupersededBy"):
            a = els.get(g.get("SupersededBy"))
            assert a is not None and a.tag == "arrow" and a.get("Z") == g.get("Z")
    page = root.find("page")  # the references' page order
    order = ["fragment", "t", "graphic", "curve", "scheme", "arrow"]
    tags = [order.index(el.tag) for el in page]
    assert tags == sorted(tags)


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
                                  "examples/m2_scheme_multi_x1.0.yaml", "examples/m2_sequence_x1.0.yaml",
                                  "examples/m4_wrap_x1.0.yaml", "examples/m4_resonance_x1.0.yaml"])
def test_unsupported_specs_are_refused(spec):
    with pytest.raises(CdxmlUnsupported):
        spec_cdxml(mr.REPO / spec, warn=lambda m: None)


def test_omissions_warn():
    msgs = []
    spec_cdxml(mr.REPO / "examples/m1_charges_x1.0.yaml", warn=msgs.append)
    assert not any("dots" in m for m in msgs)  # written since phase 3
    assert "nitromethane: circled charge written as a plain sign" in msgs
    assert "iron_iii: charge label '+3' written as the formal charge +2" in msgs
    assert len(msgs) == len(set(msgs))


# --------------------------------------------------------------------------- reference diff


def _xy(s):
    return tuple(float(v) for v in s.split()[:2])


def _xy4(s):
    return tuple(float(v) for v in s.split())


def _nodes(root):
    """Fragment nodes (not nickname internals) by id; mech_basic's tBu-Cl is in a <group>."""
    return {n.get("id"): n for f in root.iter("fragment") if f.getparent().tag != "n" for n in f.findall("n")}


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
            for g in root.iter("graphic") if g.get("SymbolType") in ("Plus", "Minus")}


def _text(t):
    return "".join(s.text or "" for s in t.iter("s"))


def reference_pairs(ref: str, text: str | None = None) -> tuple[list[str], list[tuple[str, str, float, float]]]:
    """(mismatches, [(category, key, ChemDraw, ours)] for numbers); `text` defaults to
    reference_export(ref)."""
    theirs = etree.parse(str(mr.REF_DIR / f"{ref}.cdxml"), PARSER).getroot()
    ours = etree.fromstring((text or reference_export(ref)).encode(), PARSER)
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
    for f in ours.iter("fragment"):  # the drawn ink box; a hand-placed charge sign or dot moves it
        t = tn[to_ref[f.find("n").get("id")]].getparent()
        syms = {g.get("SymbolType") for g in f.findall("graphic")}
        cat = ("fragment box (with dots)" if syms & {"LonePair", "Electron"}
               else "fragment box (with a charge)" if syms else "fragment box")
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


# --------------------------------------------------------------------------- phase 2: schemes

SCHEME_REFS = [f"{k}_x{s}" for k in ("scheme_basic", "scheme_multi", "equilibrium") for s in ("1.0", "0.8", "0.6")]
SCHEME_SPECS = ([f"examples/m2_{k}_x{s}.yaml" for k in ("scheme_basic", "equilibrium") for s in ("1.0", "0.8", "0.6")]
                + ["examples/m2_stages_x1.0.yaml", "examples/m2_stages3_x1.0.yaml", "examples/m2_suzuki_x1.0.yaml"])
# ChemDraw values the rules (acs-style.md §13.1) don't reproduce, with their measured size
KNOWN = {
    ("scheme_basic_x0.6", "arrow box", 1): 0.35,  # pad above the axis 0.18, not 0.5 (scheme_wrap_x0.6 too)
    ("*", "stacked caption box", 3): 0.35,  # bottom of a caption of 2+ lines: 0.04-0.32 below the ink, unexplained
    ("*", "formula superscript box", 1): 3.6,  # H3O+ (face 96): box top 0.354 x size above the drawn "+"
}


def _floats(s):
    return [float(v) for v in s.split()]


@lru_cache(maxsize=None)
def reference_scheme_export(ref: str) -> str:
    """compare.reference_cdx_scheme written out."""
    return to_cdxml(reference_cdx_scheme(ref), preset_for(ref), f"{ref}.cdxml", warn=lambda m: None)


def _near(p, q, tol=PAIR_TOL):
    return max(abs(a - b) for a, b in zip(p, q)) < tol


def scheme_pairs(ref: str) -> tuple[list[str], list[tuple[str, str, float, float]]]:
    """(mismatches, [(category, key, ChemDraw, ours)]) for the scheme objects."""
    theirs = etree.parse(str(mr.REF_DIR / f"{ref}.cdxml"), PARSER).getroot().find("page")
    ours = etree.fromstring(reference_scheme_export(ref).encode(), PARSER).find("page")
    bad, nums = [], []
    # arrows, paired by Head3D
    o2t = {}
    t_arrows = list(theirs.iter("arrow"))
    for oa in ours.iter("arrow"):
        m = [ta for ta in t_arrows if _near(_floats(ta.get("Head3D")), _floats(oa.get("Head3D")))]
        if len(m) != 1:
            bad.append(f"arrow at {oa.get('Head3D')}: {len(m)} ChemDraw arrows")
            continue
        ta = m[0]
        og = next(g for g in ours.iter("graphic") if g.get("SupersededBy") == oa.get("id"))
        tg = next(g for g in theirs.iter("graphic") if g.get("SupersededBy") == ta.get("id"))
        o2t[og.get("id")] = tg.get("id")
        key = f"arrow {ta.get('id')}"
        for attr in ("FillType", "ArrowheadHead", "ArrowheadTail", "ArrowheadType", "HeadSize", "ArrowheadCenterSize",
                     "ArrowheadWidth", "ArrowShaftSpacing", "ArrowEquilibriumRatio"):
            if oa.get(attr) != ta.get(attr):
                bad.append(f"{key} {attr}: ChemDraw {ta.get(attr)} ours {oa.get(attr)}")
        for attr in ("ArrowType", "HeadSize"):
            if og.get(attr) != tg.get(attr):
                bad.append(f"{key} graphic {attr}: ChemDraw {tg.get(attr)} ours {og.get(attr)}")
        for attr in ("Head3D", "Tail3D", "Center3D", "MajorAxisEnd3D", "MinorAxisEnd3D", "BoundingBox"):
            cat = "arrow box" if attr == "BoundingBox" else "arrow points"
            for k, (a, b) in enumerate(zip(_floats(ta.get(attr)), _floats(oa.get(attr)))):
                nums.append((cat, (k, f"{key} {attr}[{k}]"), a, b))
        for k, (a, b) in enumerate(zip(_floats(tg.get("BoundingBox")), _floats(og.get("BoundingBox")))):
            nums.append(("arrow points", (k, f"{key} graphic box[{k}]"), a, b))
    # captions, paired by p: their <t> or formula fragment (text frags by drawn origin: the label's p)
    doc = mr.load(ref)
    tobj = {t.get("id"): t for t in theirs.iter("t") if t.getparent().tag != "n"}  # mech_basic: in <group>s
    for f in theirs.findall("fragment"):
        if f.get("id") in doc.text_frags:
            tobj[f.get("id")] = f.find("n").find("t")
    # the lines of their stacked captions as drawn: we write one <t> per line where a line has a
    # sub- or superscript (acs-style.md §13.1), each paired with its line: p on its drawn
    # baseline, box its ink
    their_lines = {}
    for cap in doc.captions:
        if "\n" not in cap["text"]:
            continue
        p = sub(cap["p"], doc.offset)
        tid = next(i for i, t in tobj.items() if _near(_floats(t.get("p")), p))
        for k, ln in enumerate(mr.text_lines(mr._caption_glyphs(doc, cap))):
            ink = (*sub(ln["ink"][:2], doc.offset), *sub(ln["ink"][2:], doc.offset))
            their_lines[(tid, k)] = ((p[0], ln["base"] - doc.offset[1]), ink, cap["text"].split("\n")[k])
    split = {}
    for ot in ours.findall("t"):
        m = [i for i, t in tobj.items() if _near(_floats(t.get("p")), _floats(ot.get("p")))]
        if len(m) == 1 and "\n" in _text(tobj[m[0]]) and "\n" not in _text(ot):
            m = []  # the first line of a split caption: paired with its line below
        if not m:
            lm = [k for k, (q, _, _) in their_lines.items() if _near(q, _floats(ot.get("p")))]
            if len(lm) == 1:
                (tid, k), (q, ink, text) = lm[0], their_lines[lm[0]]
                o2t[ot.get("id")] = tid
                split.setdefault(tid, []).append(k)
                key = f"caption {tid} line {k + 1} {text!r}"
                if _text(ot) != text:
                    bad.append(f"{key} text: ours {_text(ot)!r}")
                if ot.get("Justification") != tobj[tid].get("Justification"):
                    bad.append(f"{key} Justification: ChemDraw {tobj[tid].get('Justification')} ours {ot.get('Justification')}")
                for j, (a, b) in enumerate(zip(q, _floats(ot.get("p")))):
                    nums.append(("caption p", (j, f"{key} p[{j}]"), a, b))
                for j, (a, b) in enumerate(zip(ink, _floats(ot.get("BoundingBox")))):
                    nums.append(("caption box", (j, f"{key} box[{j}]"), a, b))
                continue
        if len(m) != 1:
            bad.append(f"caption {_text(ot)!r} at {ot.get('p')}: {len(m)} ChemDraw captions")
            continue
        tt = tobj[m[0]]
        o2t[ot.get("id")] = m[0]
        key = f"caption {m[0]} {_text(tt)!r}"
        if _text(ot) != _text(tt):
            bad.append(f"{key} text: ours {_text(ot)!r}")
        for attr in ("Justification", "CaptionJustification", "LineStarts"):
            if ot.get(attr) != tt.get(attr):
                bad.append(f"{key} {attr}: ChemDraw {tt.get(attr)} ours {ot.get(attr)}")
        for k, (a, b) in enumerate(zip(_floats(tt.get("p")), _floats(ot.get("p")))):
            nums.append(("caption p", (k, f"{key} p[{k}]"), a, b))
        stacked = "\n" in _text(tt)
        formula_sup = any(s.get("face") == "96" and re.search(r"[A-Za-z0-9)][+\u2212-]$", s.text or "")
                          for s in tt.iter("s"))
        cat = "stacked caption box" if stacked else "formula superscript box" if formula_sup else "caption box"
        for k, (a, b) in enumerate(zip(_floats(tt.get("BoundingBox")), _floats(ot.get("BoundingBox")))):
            nums.append((cat, (k, f"{key} box[{k}]"), a, b))
    for tid, ks in split.items():
        if sorted(ks) != list(range(len(_text(tobj[tid]).split("\n")))):
            bad.append(f"caption {tid}: lines {sorted(ks)} written, ChemDraw has {_text(tobj[tid]).count(chr(10)) + 1}")
    # fragments, paired by their first node
    tnodes = {tuple(_floats(n.get("p"))[:2]): f.get("id") for f in theirs.iter("fragment") if f.getparent().tag != "n"
              for n in f.findall("n")}
    for f in ours.findall("fragment"):
        p = _floats(f.find("n").get("p"))
        m = [i for q, i in tnodes.items() if _near(q, p)]
        if len(m) == 1:
            o2t[f.get("id")] = m[0]
        else:
            bad.append(f"fragment {f.get('id')}: {len(m)} ChemDraw fragments")
    written = set(o2t.values())
    # a step citing a <group> (mech_basic's tBu-Cl with its fishhooks) cites the fragment in it
    in_group = {g.get("id"): [f.get("id") for f in g.iter("fragment")] for g in theirs.iter("group")}

    def their_ids(st, attr):
        return [j for i in (st.get(attr) or "").split() for j in in_group.get(i, [i])]
    # steps, paired by arrow; their ids we didn't write (abbreviated species) are left out
    tsteps = {st.get("ReactionStepArrows"): st for st in theirs.iter("step")}
    osteps = list(ours.iter("step"))
    if len(osteps) != len(tsteps):
        bad.append(f"{len(osteps)} steps, ChemDraw {len(tsteps)}")
    for st in osteps:
        tst = tsteps.get(o2t.get(st.get("ReactionStepArrows")))
        if tst is None:
            bad.append(f"step {st.get('id')}: no ChemDraw step with its arrow")
            continue
        for attr in STEP_KINDS:
            got = sorted({o2t.get(i, f"?{i}") for i in (st.get(attr) or "").split()})  # split lines: one id
            want = sorted(i for i in their_ids(tst, attr) if i in written)
            if got != want:
                bad.append(f"step {tst.get('id')} {attr}: ChemDraw {want} ours {got}")
    # the steps sit in the same <scheme>s
    group = [sorted(o2t.get(st.get("ReactionStepArrows")) for st in sch.findall("step")) for sch in ours.iter("scheme")]
    tgroup = [sorted(st.get("ReactionStepArrows") for st in sch.findall("step")) for sch in theirs.iter("scheme")]
    if sorted(group) != sorted(tgroup):
        bad.append(f"schemes: ChemDraw {tgroup} ours {group}")
    return bad, nums


@needs_refs
@pytest.mark.parametrize("ref", SCHEME_REFS)
def test_scheme_reference_diff(ref):
    bad, nums = scheme_pairs(ref)
    for cat, (k, key), a, b in nums:
        tol = KNOWN.get((ref, cat, k), KNOWN.get(("*", cat, k), TOL))
        if abs(b - a) > tol:
            bad.append(f"{cat} {key}: ChemDraw {a} ours {b} (diff {b - a:+.3f})")
    assert not bad, "\n".join(bad)
    assert any(c == "arrow points" for c, *_ in nums) and any(c == "caption p" for c, *_ in nums)


@needs_refs
@pytest.mark.parametrize("ref", ["scheme_basic_x1.0", "scheme_multi_x1.0", "equilibrium_x1.0"])
def test_scheme_reference_nodes(ref):
    """Phase 1's node and label checks on the scheme references' written species."""
    if ref.startswith("equilibrium_"):
        pytest.skip("stand-in species")
    bad, nums = reference_pairs(ref, reference_scheme_export(ref))
    bad += [f"{c} [{k}]: ChemDraw {a} ours {b}" for c, k, a, b in nums if abs(b - a) > TOL]
    assert not bad, "\n".join(bad)


@needs_dtd
@pytest.mark.parametrize("ref", SCHEME_REFS)
def test_scheme_reference_exports_validate(ref):
    _check_structure(reference_scheme_export(ref))


@needs_dtd
@pytest.mark.parametrize("spec", SCHEME_SPECS)
def test_scheme_spec_exports_validate(spec):
    _check_structure(exported(spec))


def _spec_smiles(spec_path) -> dict[str, str]:
    """name -> canonical SMILES (no maps) of each species as the spec gives it."""
    spec = load_spec(spec_path)
    out = {}
    for name, sp in spec.species.items():
        m = (Chem.MolFromSmiles(sp.smiles) if sp.smiles
             else Chem.MolFromMolFile(str(resolve_molfile(sp.molfile, spec_path.parent))))
        for a in m.GetAtoms():
            a.SetAtomMapNum(0)
        out[name] = Chem.MolToSmiles(m)
    return out


@pytest.mark.parametrize("spec", SCHEME_SPECS)
def test_scheme_rdkit_round_trip(spec):
    """Every species reads back as its spec's (charges and radicals too); each step's reactants
    and products, by fragment id, are the spec step's (reversed for an arrow favoring reactants,
    as ChemDraw stores it), and the species drawn above its arrow are cited there."""
    path = mr.REPO / spec
    text = exported(spec)
    smiles = _spec_smiles(path)
    root = etree.fromstring(text.encode(), PARSER)
    frags = [f.get("id") for f in root.find("page").findall("fragment")]
    got = dict(zip(frags, (Chem.MolToSmiles(m) for m in Chem.MolsFromCDXML(text))))
    s = load_spec(path)
    drawn = [n for st in s.steps for n in st.sources + st.drawn_products]
    assert sorted(got.values()) == sorted(smiles[n] for n in dict.fromkeys(drawn))
    steps = list(root.iter("step"))
    assert len(steps) == len(s.steps)
    for i, (st, el) in enumerate(zip(s.steps, steps), 1):
        over = s.over_names(i)
        want = ([smiles[n] for n in st.sources if n not in over], [smiles[n] for n in st.drawn_products])
        if st.arrow_type == "equilibrium_reactants":
            want = want[::-1]
        have = tuple([got[i] for i in el.get(k).split()] for k in ("ReactionStepReactants", "ReactionStepProducts"))
        assert have == want
        above = [got[j] for j in (el.get("ReactionStepObjectsAboveArrow") or "").split() if j in got]
        assert above == [smiles[n] for n in over]
    if len(s.steps) == 1:
        RDLogger.DisableLog("rdApp.*")  # RDKit looks for the <t> conditions among its fragments
        try:
            rxns = rdChemReactions.ReactionsFromCDXMLBlock(text)
        finally:
            RDLogger.EnableLog("rdApp.*")
        assert len(rxns) == 1
        # a reaction's SMILES keeps drawn H atoms ([H]C(C)(C)C): read back as plain molecules
        r, p = ([Chem.MolToSmiles(Chem.MolFromSmiles(x)) for x in side.split(".")]
                for side in rdChemReactions.ReactionToSmiles(rxns[0]).split(">>"))
        assert sorted(r) == sorted(want[0]) and sorted(p) == sorted(want[1])


def test_scheme_objects():
    """m2_scheme_basic_x1.0: the conditions are one centered <t> per side (two lines below, with
    LineStarts), the plus and numbers Left; the step cites them and the arrow's <graphic>."""
    root = etree.fromstring(exported("examples/m2_scheme_basic_x1.0.yaml").encode(), PARSER)
    texts = {_text(t): t for t in root.find("page").findall("t")}
    assert set(texts) == {"+", "K2CO3", "Acetone, reflux\n65%", "1", "2", "3"}
    below = texts["Acetone, reflux\n65%"]
    assert below.get("Justification") == "Center" and below.get("LineStarts") == "16 19"
    assert texts["+"].get("Justification") is None
    assert [s.get("face") for s in texts["K2CO3"].iter("s")] == [None, "32", None, "32"]
    assert all(s.get("font") == "310" and s.get("face") == "1" for n in "123" for s in texts[n].iter("s"))
    step = root.find(".//step")
    assert step.get("ReactionStepObjectsAboveArrow") == texts["K2CO3"].get("id")
    assert step.get("ReactionStepObjectsBelowArrow") == below.get("id")
    assert root.find(".//graphic").get("id") == step.get("ReactionStepArrows")


def test_stage_lines_are_left_captions():
    """m2_stages3_x1.0: one Left <t> per stage, two above and one below; *n* italic."""
    root = etree.fromstring(exported("examples/m2_stages3_x1.0.yaml").encode(), PARSER)
    page = root.find("page")
    step = root.find(".//step")
    above = [page.find(f"t[@id='{i}']") for i in step.get("ReactionStepObjectsAboveArrow").split()]
    below = [page.find(f"t[@id='{i}']") for i in step.get("ReactionStepObjectsBelowArrow").split()]
    assert [_text(t) for t in above] == ["1. n-BuLi, THF, \u221278 \u00b0C", "2. B(OMe)3"]
    assert [_text(t) for t in below] == ["3. H3O+"]
    assert all(t.get("Justification") is None for t in above + below)
    assert [(s.text, s.get("font"), s.get("face")) for s in above[0].iter("s")][1] == ("n", "307", "2")
    assert [s.get("face") for s in below[0].iter("s")][-1] == "64"


def test_new_row_starts_a_scheme():
    """m2_equilibrium_x1.0: three sequences, three <scheme>s; the arrow favoring reactants has
    Head3D on the left and its step reads right to left, as in equilibrium_x1.0."""
    root = etree.fromstring(exported("examples/m2_equilibrium_x1.0.yaml").encode(), PARSER)
    assert [len(s.findall("step")) for s in root.iter("scheme")] == [1, 1, 1]
    arrows = list(root.iter("arrow"))
    assert [a.get("ArrowEquilibriumRatio") for a in arrows] == [None, "300", "300"]
    h, t = (_floats(arrows[2].get(k))[0] for k in ("Head3D", "Tail3D"))
    assert h < t


def test_cdxml_head_attrs_match_preset_heads():
    """cdxml.arrow_heads against the drawn heads: HeadSize and ArrowheadCenterSize x LineWidth /
    100 are head_length and head_inset; ArrowShaftSpacing and ArrowEquilibriumRatio follow the
    preset as the writer derives them (300 and 300 in every reference)."""
    for name in ("acs1996", "acs1996_x0.8", "acs1996_x0.6"):
        p = load_preset(name)
        lw = p.lengths.line_width
        for variant, (size, center, _) in p.cdxml.arrow_heads.items():
            h = p.head(variant)
            assert size * lw / 100 == pytest.approx(h.head_length, abs=0.03)
            assert center * lw / 100 == pytest.approx(h.head_inset, abs=0.03)
        assert round(p.equilibrium.shaft_spacing / lw * 100) == 300
        assert round(p.ratios.equilibrium_long_per_short * 100) == 300


@needs_refs
def test_cdxml_arrow_conventions_hold_in_references():
    """The rules the writer applies (acs-style.md §13.1), on every scheme reference's arrows:
    Center3D the midpoint, MajorAxisEnd3D = center + (head - tail), MinorAxisEnd3D turned 90°,
    the box's half width from ArrowheadWidth / ArrowShaftSpacing, and the pad, all to 0.02 pt,
    except arrows_large_x1.0's balanced arrow (minor axis 28.8 long, a copied arrow's) and the
    x0.6 forward arrows padded 0.18 (scheme_basic_x0.6, scheme_wrap_x0.6)."""
    pad = 0.5
    for ref in mr.CDXML_SCHEME_REFS:
        for k, a in enumerate(mr.cdxml_scheme(ref)["arrows"]):
            assert max(map(abs, a["center_minus_mid"] + a["major_minus_(center + head - tail)"])) < 0.02, ref
            if not (ref == "arrows_large_x1.0" and k == 0):
                assert max(map(abs, a["minor_minus_(center + rot90(head - tail))"])) < 0.02, ref
            top, bottom = a["box_pad (top, bottom)"]
            eq = "ArrowShaftSpacing" in a["attrs"]
            both = a["attrs"].get("ArrowheadTail") == "Full"
            want = ((pad, pad) if a["head3d_right_of_tail3d"] else (0, 0)) if eq else (pad, pad if both else 0)
            if ref in ("scheme_basic_x0.6", "scheme_wrap_x0.6"):
                want = (0.175, 0)
            assert top == pytest.approx(want[0], abs=0.02) and bottom == pytest.approx(want[1], abs=0.02), (ref, k)


def test_stacked_scripts_are_split():
    """m2_suzuki_x1.0: its stacked conditions have subscripts, so each line is its own centered
    <t> (ChemDraw sets a stack of explicit script runs 11.775 apart: gate 2026-10-28), on our
    baselines; the step cites every line. A plain stack (m2_scheme_basic) stays one <t>."""
    root = etree.fromstring(exported("examples/m2_suzuki_x1.0.yaml").encode(), PARSER)
    page = root.find("page")
    step = root.find(".//step")
    for attr, want in (("ReactionStepObjectsAboveArrow", ["PhB(OH)2", "Pd(PPh3)4"]),
                       ("ReactionStepObjectsBelowArrow", ["1,4-dioxane/H2O (4:1)", "90 °C, 12 h"])):
        ts = [page.find(f"t[@id='{i}']") for i in step.get(attr).split()]
        assert [_text(t) for t in ts] == want
        assert all(t.get("Justification") == "Center" and t.get("LineStarts") is None for t in ts)
        (x0, y0), (x1, y1) = (_floats(t.get("p")) for t in ts)
        assert x0 == pytest.approx(x1, abs=0.01) and y1 > y0


# --------------------------------------------------------------------------- phase 3: mechanisms

MECH_REFS = [f"mech_basic_x{s}" for s in ("1.0", "0.8", "0.6")]
M3_SPECS = sorted(str(p.relative_to(mr.REPO)) for p in (mr.REPO / "examples").glob("m3_*.yaml"))
CURVE_TOL = 0.01  # anchors and inner handles: ChemDraw's numbers copied
MIRROR_TOL = 0.025  # end handles: we mirror its rounded points (2p0 - c1), it mirrors unrounded ones


def _dot_symbols(root, nodes: dict) -> dict:
    """node id (mapped by `nodes`) -> [(SymbolType, BoundingBox)] of its LonePair / Electron symbols."""
    out = {}
    for g in root.iter("graphic"):
        if g.get("SymbolType") in ("LonePair", "Electron"):
            r = g.find("represent")
            assert r.get("attribute") == "Radical"
            out.setdefault(nodes.get(r.get("object"), r.get("object")), []).append((g.get("SymbolType"), _floats(g.get("BoundingBox"))))
    return out


def _ccw(atom, box) -> bool:
    u, v = sub(box[:2], atom), sub(box[2:], atom)
    return u[0] * v[1] - u[1] * v[0] < 0


def mech_pairs(ref: str) -> tuple[list[str], list[tuple[str, str, float, float]]]:
    """(mismatches, [(category, key, ChemDraw, ours)]) for the curves and dot symbols."""
    theirs = etree.parse(str(mr.REF_DIR / f"{ref}.cdxml"), PARSER).getroot()
    ours = etree.fromstring(reference_scheme_export(ref).encode(), PARSER)
    bad, nums = [], []
    # curves, paired by first anchor
    tcurves = list(theirs.iter("curve"))
    ocurves = list(ours.iter("curve"))
    if len(ocurves) != len(tcurves):
        bad.append(f"{len(ocurves)} curves, ChemDraw {len(tcurves)}")
    for oc in ocurves:
        op = _floats(oc.get("CurvePoints"))
        m = [c for c in tcurves if _near(_floats(c.get("CurvePoints"))[2:4], op[2:4])]
        if len(m) != 1:
            bad.append(f"curve at {op[2:4]}: {len(m)} ChemDraw curves")
            continue
        tc = m[0]
        key = f"curve {tc.get('id')}"
        for attr in ("CurveType", "ArrowheadHead", "ArrowheadType"):
            if oc.get(attr) != tc.get(attr):
                bad.append(f"{key} {attr}: ChemDraw {tc.get(attr)} ours {oc.get(attr)}")
        tp = _floats(tc.get("CurvePoints"))
        if tc.get("CurveType") == "32":  # its fishhooks repeat the last anchor; we don't (§13.2)
            assert len(tp) == 18 and tp[12:] == tp[6:12], key
            tp = tp[:12]
        if len(tp) != len(op):
            bad.append(f"{key}: {len(op) // 2} points, ChemDraw {len(tp) // 2}")
            continue
        for k, (a, b) in enumerate(zip(tp, op)):
            cat = "curve end handle" if k // 2 in (0, 5) else "curve points"
            nums.append((cat, f"{key} [{k}]", a, b))
    # dot symbols, per atom (nodes paired by p)
    tn, on = _nodes(theirs), _nodes(ours)
    to_ref = {oid: tid for oid, o in on.items() for tid, t in tn.items() if _near(_floats(t.get("p")), _floats(o.get("p")))}
    tsym, osym = _dot_symbols(theirs, {}), _dot_symbols(ours, to_ref)
    if sorted(tsym) != sorted(osym):
        bad.append(f"dots on ChemDraw nodes {sorted(tsym)}, ours {sorted(osym)}")
    for nid, syms in osym.items():
        theirs_here = tsym.get(nid, [])
        atom = _floats(tn[nid].get("p"))
        for kind in ("LonePair", "Electron"):
            o = [b for k, b in syms if k == kind]
            t = [b for k, b in theirs_here if k == kind]
            if len(o) != len(t):
                bad.append(f"node {nid}: {len(o)} {kind}, ChemDraw {len(t)}")
                continue
            for ob in o:  # paired with ChemDraw's nearest of the same kind
                tb = min(t, key=lambda b: abs(b[0] - ob[0]) + abs(b[1] - ob[1]))
                key = f"node {nid} {kind}"
                for k in range(2):
                    nums.append(("dot center", f"{key} {'xy'[k]}", tb[k], ob[k]))
                if kind == "LonePair":
                    nums.append(("lone pair separation", key, math.dist(tb[:2], tb[2:]), math.dist(ob[:2], ob[2:])))
                    if not (_ccw(atom, ob) and _ccw(atom, tb)):
                        bad.append(f"{key}: not counterclockwise (ChemDraw {_ccw(atom, tb)}, ours {_ccw(atom, ob)})")
                else:
                    nums.append(("electron box", f"{key} dx", tb[0] - tb[2], ob[0] - ob[2]))
                    nums.append(("electron box", f"{key} dy", tb[3] - tb[1], ob[3] - ob[1]))
    return bad, nums


# Hand-placed: dot centers (lone_pair.gap is a round default), the charge signs the dots push
# out (Drawer.charge_center), and the fragment boxes they extend. Largest seen: 1.14 (O- sign, x1.0).
PLACE_TOL = 1.2
PLACED = ("dot center", "charge center", "fragment box (with dots)", "fragment box (with a charge)")
# what the rules don't reproduce: (ref or "*", category, key prefix) -> its measured size
MECH_KNOWN = {
    # tBu radical on an unlabeled C: ChemDraw's dot is 2.8-5.0 from ours (our largest-gap rule)
    ("*", "dot center", "node 70 Electron"): 5.1,
    # ChemDraw pads the arrow box 0.18 above the axis, not 0.5 (as scheme_basic_x0.6)
    ("mech_basic_x0.8", "arrow box", "arrow"): 0.35, ("mech_basic_x0.6", "arrow box", "arrow"): 0.35,
    # Cl. typed as text (no face 96): ChemDraw stores its box 0.1 higher than a drawn label's
    ("mech_basic_x1.0", "label box", "node 79"): 0.11,
}


def _mech_tol(ref, cat, key, default):
    return next((t for (r, c, k), t in MECH_KNOWN.items() if r in (ref, "*") and c == cat and key.startswith(k)),
                default)


@needs_refs
@pytest.mark.parametrize("ref", MECH_REFS)
def test_mech_reference_diff(ref):
    """mech_basic_x* rebuilt (compare.reference_cdx_scheme: templates, ChemDraw's CurvePoints,
    lone pairs on ChemDraw's atoms, its steps) and compared with its file: curves to 0.01 (end
    handles to MIRROR_TOL), CurveType and ArrowheadHead; per atom, dot symbol kinds and counts,
    LonePair separation and Electron dx exact, both counterclockwise; dot centers, the signs
    they push out and the boxes they extend within PLACE_TOL (hand-placed); the scheme objects
    and nodes as for the scheme references (steps citing a <group> cite its fragment)."""
    bad, nums = mech_pairs(ref)
    sbad, snums = scheme_pairs(ref)
    nbad, nnums = reference_pairs(ref, reference_scheme_export(ref))
    bad += sbad + nbad
    tols = {"curve points": CURVE_TOL, "curve end handle": MIRROR_TOL, "lone pair separation": 0.01,
            "electron box": 0.01, **dict.fromkeys(PLACED, PLACE_TOL)}
    for cat, key, a, b in nums + [(c, k, a, b) for c, (_, k), a, b in snums] + nnums:
        tol = _mech_tol(ref, cat, key, tols.get(cat, TOL))
        if abs(b - a) > tol:
            bad.append(f"{cat} {key}: ChemDraw {a} ours {b} (diff {b - a:+.3f})")
    assert not bad, "\n".join(bad)
    cats = {c for c, *_ in nums}
    assert {"curve points", "dot center", "lone pair separation", "electron box"} <= cats


@needs_refs
def test_curve_end_is_stroke_end():
    """In all 12 reference curves the last anchor is the SVG stroke end (the head's notch, our
    MechArrow.p3), not the head tip, so the writer writes p3 as it is."""
    n = 0
    for ref in MECH_REFS:
        for c in mr.measure(mr.load(ref))["curved_arrows"]:
            assert max(map(abs, c["last_anchor_minus_stroke_end"])) < 0.01, ref
            assert max(map(abs, c["head"]["tip_minus_stroke_end"])) > 0.5  # the tip is elsewhere
            n += 1
    assert n == 12


@needs_refs
def test_fishhook_side_matches_reference():
    """mechanism.fishhook_side gives +1 (cdxml_writer: HalfRight) on every reference fishhook,
    and the writer's ArrowheadHead equals ChemDraw's on every curve."""
    from compare import curve_points
    from acs_chem_svg.mechanism import fishhook_side
    for ref in MECH_REFS:
        hooks = [c for c in mr.load(ref).curves if c["type"] == "32"]
        assert len(hooks) == 2 and all(c["head"] == "HalfRight" for c in hooks)
        p = [curve_points(c) for c in hooks]
        assert all(fishhook_side(c1, c2, p3) == 1.0 for _, c1, c2, p3 in p)


def test_fishhook_left_is_half_left():
    """The references' fishhooks turn left on screen (counterclockwise, y down) and are
    HalfRight; one turning right is written HalfLeft, by symmetry (no reference has one; checked
    at the phase 3 gate)."""
    from acs_chem_svg.cdxml_writer import CdxCurve, CdxScheme
    turns_left = CdxCurve((0.0, 0.0), (10.0, 10.0), (20.0, 10.0), (30.0, 0.0), "fishhook")  # down, then up
    turns_right = CdxCurve((0.0, 0.0), (10.0, -10.0), (20.0, -10.0), (30.0, 0.0), "fishhook")
    heads = []
    for c in (turns_left, turns_right):
        root = etree.fromstring(to_cdxml(CdxScheme([], curves=[c]), "acs1996", warn=lambda m: None).encode(), PARSER)
        heads.append(root.find(".//curve").get("ArrowheadHead"))
    assert heads == ["HalfRight", "HalfLeft"]


@needs_refs
def test_dot_box_keys_match_references():
    """cdxml.electron_box_dx, lone_pair_box_reach and electron_box_reach against
    measure_refs.py --cdxml-symbols: every Electron box's dx; per fragment side bounded by a dot,
    the reach (more where other ink bounds it), each reach seen at least once."""
    for ref, preset in zip(MECH_REFS, ("acs1996", "acs1996_x0.8", "acs1996_x0.6")):
        cx = load_preset(preset).cdxml
        m = mr.cdxml_symbols(ref)
        assert all(e["dx"] == pytest.approx(cx.electron_box_dx, abs=0.01) and e["dy"] == 0 for e in m["electrons"])
        assert all(p["separation"] == pytest.approx(load_preset(preset).lone_pair.dot_spacing, abs=0.01)
                   for p in m["lone_pairs"])
        reach = {"pair": cx.lone_pair_box_reach, "electron": cx.electron_box_reach}
        seen = set()
        for r in m["fragment_reach_beyond_dot_centers"]:
            for side in ("left", "top", "right", "bottom"):
                v, kind = r[side]
                assert v >= reach[kind] - 0.01, (ref, r)
                if abs(v - reach[kind]) <= 0.01:
                    seen.add(kind)
        assert seen == {"pair", "electron"}, ref


@needs_dtd
@pytest.mark.parametrize("ref", MECH_REFS)
def test_mech_reference_exports_validate(ref):
    _check_structure(reference_scheme_export(ref))


@needs_dtd
@pytest.mark.parametrize("spec", M3_SPECS)
def test_mech_spec_exports_validate(spec):
    _check_structure(exported(spec))


@pytest.mark.parametrize("spec", M3_SPECS)
def test_mech_rdkit_round_trip(spec):
    test_scheme_rdkit_round_trip(spec)


def test_mech_objects():
    """m3_sn2_x1.0: one <curve> per curved arrow, free on the page after the graphics; LonePair
    symbols on the nucleophile and leaving group, bound by represent attribute="Radical"."""
    root = etree.fromstring(exported("examples/m3_sn2_x1.0.yaml").encode(), PARSER)
    page = root.find("page")
    curves = page.findall("curve")
    assert len(curves) == 2 and all(c.get("CurveType") == "8" and c.get("ArrowheadHead") == "Full" for c in curves)
    assert all(len(_floats(c.get("CurvePoints"))) == 12 for c in curves)
    assert {g.get("SymbolType") for g in root.iter("graphic") if g.get("GraphicType") == "Symbol"} >= {"LonePair"}
    root = etree.fromstring(exported("examples/m3_radical_initiation_x1.0.yaml").encode(), PARSER)
    hooks = root.findall(".//curve")
    # six points, as a 2e- arrow: no repeated last anchor (a stray node in ChemDraw; gate 2026-10-28)
    assert hooks and all(c.get("CurveType") == "32" and len(_floats(c.get("CurvePoints"))) == 12 for c in hooks)
    assert all(n.get("Radical") == "Doublet" for g in root.iter("graphic") if g.get("SymbolType") == "Electron"
               for n in root.iter("n") if n.get("id") == g.find("represent").get("object"))


def test_mechanism_maps_do_not_warn():
    """A mechanism's species are all mapped (its arrows need the maps), so the map warning is
    silenced there (user decision 2026-10-28); M1 species still warn."""
    msgs = []
    spec_cdxml(mr.REPO / "examples/m3_carbonyl_addition_x1.0.yaml", warn=msgs.append)
    assert msgs == []
    spec_cdxml(mr.REPO / "examples/m1_charges_x1.0.yaml", warn=msgs.append)
    assert any("atom maps not written" in m for m in msgs)


def test_added_reactant_is_above_arrow():
    """m3_carbonyl_addition_x1.0 step 2: water, added by the mechanism step and drawn above its
    arrow, is cited in ReactionStepObjectsAboveArrow, not as a reactant (user decision 2026-10-28)."""
    text = exported("examples/m3_carbonyl_addition_x1.0.yaml")
    root = etree.fromstring(text.encode(), PARSER)
    frags = [f.get("id") for f in root.find("page").findall("fragment")]
    smiles = dict(zip(frags, (Chem.MolToSmiles(m) for m in Chem.MolsFromCDXML(text))))
    st = list(root.iter("step"))[1]
    assert [smiles[i] for i in st.get("ReactionStepObjectsAboveArrow").split()] == ["O"]
    assert [smiles[i] for i in st.get("ReactionStepReactants").split()] == [_canonical("N#CC(C)([O-])C")]
