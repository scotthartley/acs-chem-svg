"""CDXML export, phase 4c (acs-style.md §13.3): general abbreviations (abbrev_2_x*, ligands_x*):
generated fragments, ring NTs, text labels.
"""

import math
import re

import pytest
import yaml
from lxml import etree
from rdkit import Chem

import measure_refs as mr
from acs_chem_svg.cdxml_writer import spec_cdxml
from acs_chem_svg.geometry import sub

from cdxml_helpers import (
    PARSER, TOL, _canonical, _check_structure, _nodes, _text, _xy, exported, needs_dtd, needs_refs,
    reference_export, reference_pairs,
)


# abbrev_2_x0.8 / x0.6 and ligands_x0.8 / x0.6: several fragments sit 1.3-2.3 pt apart in their CDXML
# and their SVG (moved between the two saves), so the CDXML diff runs on the x1.0 files
ABBREV2_REFS = ["abbrev_2_x1.0", "ligands_x1.0"]
# deliberate differences: ChemDraw's left-justified P(o-tol)3 and p-Tol are Unspecified with no
# fragment; ours keep their chemistry as a Fragment (only right-justified, where ChemDraw would
# scramble the text, ours are its Unspecified form too). IPr is text (ChemDraw read it as I + Pr),
# and IPr / P(o-tol)3 right-justified keep our reversal (user decision 2026-11-03): P(o-tol)3
# stored as (o-tol)3P, IPr stored as PrI, which ChemDraw shows as IPr (phase 4c label variants,
# rows B and C)
DECIDED_NODE = re.compile(r"NodeType: ChemDraw Unspecified ours Fragment|IPr: no inner fragment"
                          r"|text: ChemDraw P\(o-tol\)3 ours \(o-tol\)3P|text: ChemDraw IPr ours PrI")
DECIDED_LABELS = ("P(o-tol)3", "IPr")


@needs_refs
@pytest.mark.parametrize("ref", ABBREV2_REFS)
def test_abbrev2_reference_diff(ref):
    """Our CDXML of abbrev_2_x1.0 / ligands_x1.0 vs ChemDraw's: node types, label texts,
    justification and alignment, label p and boxes, inner fragments of every reference template;
    ring NTs's generated inner fragment by its structure (test_ring_nts_cdxml)."""
    bad, nums = reference_pairs(ref)
    bad = [b for b in bad if not DECIDED_NODE.search(b) and " NTs" not in b]
    bad += [f"{c} [{k}]: ChemDraw {a} ours {b} (diff {b - a:+.3f})" for c, k, a, b in nums
            if abs(b - a) > TOL and not any(lab in k for lab in DECIDED_LABELS)
            and not (c == "fragment box" and ref == "ligands_x1.0")]  # its P(o-tol)3 / IPr fragments
    assert not bad, "\n".join(bad)
    assert sum(c == "nickname atoms" for c, *_ in nums) > (300 if ref.startswith("abbrev") else 100)


def _abbrev_export(spec: dict) -> etree._Element:
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "s.yaml"
        p.write_text(yaml.safe_dump(spec))
        return etree.fromstring(spec_cdxml(p, warn=lambda m: None).encode(), PARSER)


def _label_node(root, text):
    (n,) = [n for n in _nodes(root).values() if n.find("t") is not None and _text(n.find("t")) == text]
    return n


def test_ring_nts_cdxml():
    """A composed label with two bonds (ring NTs, abbrev_2_x*): NodeType Fragment, one
    ExternalConnectionPoint per bond (ExternalConnectionNum 1, 2) on its bond, BondLength - 0.05
    out; the fragment's ConnectionOrder and the node's BondOrdering name them in the same order;
    RDKit reads the tosyl pyrrolidine back."""
    smiles = "Cc1ccc(cc1)S(=O)(=O)N1CCCC1"
    root = _abbrev_export({"species": {"a": {"smiles": smiles, "abbreviate": ["Ts"]}}})
    n = _label_node(root, "NTs")
    assert n.get("NodeType") == "Fragment"
    inner = n.find("fragment")
    ecps = [m for m in inner.findall("n") if m.get("NodeType") == "ExternalConnectionPoint"]
    assert [m.get("ExternalConnectionNum") for m in ecps] == ["1", "2"]
    assert inner.get("ConnectionOrder") == " ".join(m.get("id") for m in ecps)
    bonds = {b.get("id"): b for b in root.iter("b")}
    outer = [bonds[i] for i in n.get("BondOrdering").split()]
    nodes = _nodes(root)
    p = _xy(n.get("p"))
    for e, b in zip(ecps, outer):
        other = _xy(nodes[b.get("E") if b.get("B") == n.get("id") else b.get("B")].get("p"))
        u = sub(other, p)
        d = math.hypot(*u)
        want = (p[0] + u[0] / d * 14.35, p[1] + u[1] / d * 14.35)
        assert math.dist(_xy(e.get("p")), want) < 0.01
    text = etree.tostring(root, encoding="unicode")
    (m,) = Chem.MolsFromCDXML(text)
    assert Chem.MolToSmiles(m) == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


@pytest.mark.parametrize("smiles, abbreviate, key", [
    ("COCOc1ccccc1", ["OMOM"], "OMOM"),          # a reference template (abbrev_2_x*)
    ("CCOCOc1ccccc1", ["Xyz"], "OXyz"),          # a spec's own group, composed: generated
    ("COCCOCOc1ccccc1", ["MEM"], "OMEM"),        # a registry group without a reference: generated
])
def test_abbreviation_round_trip(smiles, abbreviate, key):
    """Template or generated, the nickname node reads back as the structure (RDKit). The label
    is left-justified here (orientation), so even OXyz, which ChemDraw can't read, keeps its
    fragment (right-justified: test_unreadable_right_label)."""
    spec = {"species": {"a": {"smiles": smiles, "abbreviate": abbreviate, "orientation": {"flip_h": True}}}}
    if "Xyz" in abbreviate:
        spec["abbreviations"] = {"Xyz": {"smiles": "*CC", "kind": "group"}}
    root = _abbrev_export(spec)
    n = _label_node(root, key)
    assert n.get("NodeType") == "Fragment" and n.find("fragment") is not None
    (m,) = Chem.MolsFromCDXML(etree.tostring(root, encoding="unicode"))
    assert Chem.MolToSmiles(m) == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_unreadable_right_label():
    """Right-justified, a label ChemDraw can't read (a spec's own group, XyzO) is ChemDraw's own
    form: Unspecified with its Warning, no inner fragment, the displayed text (label variants row C)."""
    root = _abbrev_export({"abbreviations": {"Xyz": {"smiles": "*CC", "kind": "group"}},
                           "species": {"a": {"smiles": "CCOCOc1ccccc1", "abbreviate": ["Xyz"]}}})
    n = _label_node(root, "XyzO")
    assert n.get("NodeType") == "Unspecified" and n.get("Warning") and n.find("fragment") is None
    assert n.find("t").get("LabelJustification") == "Right"


def test_text_label_cdxml():
    """A text label (ligands_x*: L, Lₙ, dppf) is Unspecified with no inner fragment; Lₙ's n is a
    letter subscript, italic (face 34), the L in the label face (96)."""
    root = _abbrev_export({"species": {"pd": {"smiles": "[*:1][Pd]([*:2])Cl"}}, "labels": {1: "Ln", 2: "dppf"}})
    for key in ("Ln", "dppf"):
        n = _label_node(root, key)
        assert n.get("NodeType") == "Unspecified" and n.find("fragment") is None
    runs = [(s.text, s.get("face")) for s in _label_node(root, "Ln").find("t").iter("s")]
    assert runs == [("L", "96"), ("n", "34")]


def test_italic_label_runs():
    """t-Bu's t is italic in the label face (98), the rest plain (96), as abbrev_2_x1.0 writes it."""
    root = _abbrev_export({"species": {"a": {"smiles": "CC(C)(C)c1ccccc1", "abbreviate": ["t-Bu"]}}})
    runs = [(s.text, s.get("face")) for s in _label_node(root, "t-Bu").find("t").iter("s")]
    assert runs == [("t", "98"), ("-Bu", "96")]


@pytest.mark.parametrize("preset", ["acs1996", "acs1996_x0.8", "acs1996_x0.6"])
def test_template_node_types_follow_the_registry(preset):
    """Every reference template is a registry group or formula, Nickname exactly when the
    registry says ChemDraw writes it as one."""
    from acs_chem_svg import abbrev
    from acs_chem_svg.cdxml_writer import load_nicknames
    reg = abbrev.default()
    for key, tpl in load_nicknames(preset).items():
        e = reg.get(key)
        assert e is not None and e.smiles, key
        assert (tpl.get("NodeType") == "Nickname") == e.nickname, key


def test_abbrev2_example_round_trip():
    """m1_abbrev2_x1.0 (one of each label kind): every species with chemistry reads back as its
    SMILES (RDKit expands templates and generated fragments alike); the text labels (Ln, dppf)
    have none, and RDKit reads their Unspecified nodes as carbons."""
    spec = "examples/m1/m1_abbrev2_x1.0.yaml"
    species = yaml.safe_load((mr.REPO / spec).read_text())["species"]
    got = dict(zip(species, (Chem.MolToSmiles(m) for m in Chem.MolsFromCDXML(exported(spec)))))
    for name, sp in species.items():
        if name != "ligands":
            assert got[name] == _canonical(sp["smiles"]), name
    assert got["ligands"] == "[CH3][Pd]([CH3])([Cl])[Cl]"


@needs_dtd
def test_abbrev2_example_validates():
    _check_structure(exported("examples/m1/m1_abbrev2_x1.0.yaml"))


def test_text_label_chemdraw_reads():
    """A text label ChemDraw reads as a formula (IPr: I + Pr) and reverses when right-justified is
    stored as the text whose reversal is ours: PrI, which ChemDraw shows as IPr (label variants,
    row B). Left-justified it stays IPr."""
    root = _abbrev_export({"species": {"a": {"smiles": "[*:1][Pd][*:2]"}}, "labels": {1: "IPr", 2: "IPr"}})
    got = sorted((_text(t), t.get("LabelJustification")) for n in _nodes(root).values()
                 if (t := n.find("t")) is not None and _text(t) in ("IPr", "PrI"))
    assert got == [("IPr", "Left"), ("PrI", "Right")]


@needs_refs
def test_typed_reversed_label_round_trip():
    """abbrev_2_x*'s Bpin typed right-justified as pinB (its `reversed:` text) is a template too;
    ChemDraw can't read Bpin, so ours right-justified is its own form: Unspecified, Warning, pinB."""
    root = etree.fromstring(reference_export("abbrev_2_x1.0").encode(), PARSER)
    got = sorted((_text(t), t.get("LabelJustification"), n.get("NodeType"), n.find("fragment") is None)
                 for n in _nodes(root).values() if (t := n.find("t")) is not None and "pin" in _text(t))
    assert got == [("Bpin", "Left", "Fragment", False), ("pinB", "Right", "Unspecified", True)]


def test_label_atom_has_no_lone_pair():
    """A label stands for its group: SnBu3's Sn (one drawn bond, RDKit radicals in the template)
    gets no carbene pair, in the SVG or the CDXML."""
    root = _abbrev_export({"species": {"a": {"smiles": "CCCC[Sn](CCCC)(CCCC)c1ccccc1", "abbreviate": ["Bu"]}}})
    assert not [g for g in root.iter("graphic") if g.get("SymbolType") == "LonePair"]
    import compare
    r = compare.render_reference("abbrev_2_x1.0") if (mr.REF_DIR / "abbrev_2_x1.0.cdxml").exists() else None
    if r is not None:
        assert not [a for m in r.molecules for a in m.atoms if a.abbrev and (a.singlet or a.radicals)]
