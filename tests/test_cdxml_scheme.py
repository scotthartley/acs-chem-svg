"""CDXML export, phase 2 (acs-style.md §13.1): M2 schemes.

- Scheme round trip: each scheme reference rebuilt as a CdxScheme at ChemDraw's coordinates
  (templates; arrows from its Head3D / Tail3D; captions set by richtext at its glyph origins;
  formula fragments such as K2CO3 as captions; its steps) and compared with its file: arrow
  points, boxes and head attributes, caption p, box, LineStarts and text, and step membership by
  paired ids, all within 0.1 pt except KNOWN; brackets by x (BoundingBox, LipSize).
  equilibrium_x*'s A and B (no templates) are one-atom stand-ins at their node positions.
- Specs: structure (DTD, ids, step citations), RDKit (species, each step's reactants and
  products by fragment id, ReactionsFromCDXMLBlock on one-step specs) and refusals.
"""

import pytest
from compare import preset_for
from lxml import etree

import measure_refs as mr
from acs_chem_svg.preset import load_preset

from cdxml_helpers import (
    PARSER, TOL, _check_structure, _floats, _text, check_scheme_round_trip, exported, needs_dtd, needs_refs,
    reference_pairs, reference_scheme_export, scheme_pairs,
)


SCHEME_REFS = [f"{k}_x{s}" for k in ("scheme_basic", "scheme_multi", "scheme_wrap", "equilibrium")
               for s in ("1.0", "0.8", "0.6")]
SCHEME_SPECS = ([f"examples/m2/m2_{k}_x{s}.yaml" for k in ("scheme_basic", "equilibrium") for s in ("1.0", "0.8", "0.6")]
                + ["examples/m2/m2_stages_x1.0.yaml", "examples/m2/m2_stages3_x1.0.yaml", "examples/m2/m2_suzuki_x1.0.yaml"])
# ChemDraw values the rules (acs-style.md §13.1) don't reproduce, with their measured size
KNOWN = {
    ("scheme_basic_x0.6", "arrow box", 1): 0.35,  # pad above the axis 0.18, not 0.5
    ("scheme_wrap_x0.6", "arrow box", 1): 0.35,  # the same
    ("*", "stacked caption box", 3): 0.35,  # bottom of a caption of 2+ lines: 0.04-0.32 below the ink, unexplained
    ("*", "formula superscript box", 1): 3.6,  # H3O+ (face 96): box top 0.354 x size above the drawn "+"
}


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
@pytest.mark.parametrize("ref", ["scheme_basic_x1.0", "scheme_multi_x1.0", "scheme_wrap_x1.0", "equilibrium_x1.0"])
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


@pytest.mark.parametrize("spec", SCHEME_SPECS)
def test_scheme_rdkit_round_trip(spec):
    """Every species reads back as its spec's (charges and radicals too); each step's reactants
    and products, by fragment id, are the spec step's (reversed for an arrow favoring reactants,
    as ChemDraw stores it), and the species drawn above its arrow are cited there."""
    check_scheme_round_trip(spec)


def test_scheme_objects():
    """m2_scheme_basic_x1.0: the conditions are one centered <t> per side (two lines below, with
    LineStarts), the plus and numbers Left; the step cites them and the arrow's <graphic>."""
    root = etree.fromstring(exported("examples/m2/m2_scheme_basic_x1.0.yaml").encode(), PARSER)
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
    root = etree.fromstring(exported("examples/m2/m2_stages3_x1.0.yaml").encode(), PARSER)
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
    root = etree.fromstring(exported("examples/m2/m2_equilibrium_x1.0.yaml").encode(), PARSER)
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
    x0.6 forward arrows padded 0.18 (scheme_basic_x0.6, scheme_wrap_x0.6). resonance_x*'s ↔
    arrows are padded cdxml.resonance_box_pad on both sides (arrows_large_x1.0's large one, 0.5)."""
    pad = 0.5
    for ref in mr.CDXML_SCHEME_REFS:
        if ref.startswith("resonance_"):
            rp = load_preset(preset_for(ref)).cdxml.resonance_box_pad
            for a in mr.cdxml_scheme(ref)["arrows"]:
                assert a["attrs"]["ArrowheadTail"] == "Full" and a["head3d_right_of_tail3d"], ref
                assert all(v == pytest.approx(rp, abs=0.02) for v in a["box_pad (top, bottom)"]), ref
            continue
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
    root = etree.fromstring(exported("examples/m2/m2_suzuki_x1.0.yaml").encode(), PARSER)
    page = root.find("page")
    step = root.find(".//step")
    for attr, want in (("ReactionStepObjectsAboveArrow", ["PhB(OH)2", "Pd(PPh3)4"]),
                       ("ReactionStepObjectsBelowArrow", ["1,4-dioxane/H2O (4:1)", "90 °C, 12 h"])):
        ts = [page.find(f"t[@id='{i}']") for i in step.get(attr).split()]
        assert [_text(t) for t in ts] == want
        assert all(t.get("Justification") == "Center" and t.get("LineStarts") is None for t in ts)
        (x0, y0), (x1, y1) = (_floats(t.get("p")) for t in ts)
        assert x0 == pytest.approx(x1, abs=0.01) and y1 > y0
