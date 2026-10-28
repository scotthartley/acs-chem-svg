"""Scheme layout rules (layout.py, acs-style.md §3) and alignment (chem.prepare)."""

import math
from functools import lru_cache
from pathlib import Path

import pytest

from acs_chem_svg import chem
from acs_chem_svg.arrows import notch_x
from acs_chem_svg.draw_mol import Drawer, draw_molecule
from acs_chem_svg.geometry_checks import check_scheme
from acs_chem_svg.layout import build_scheme, render_spec
from acs_chem_svg.preset import load_preset
from acs_chem_svg.scene import Group
from acs_chem_svg.spec import Spec, load_spec, molecules, scheme_prepared
from acs_chem_svg.svg_writer import to_svg

REPO = Path(__file__).resolve().parents[1]
EX = REPO / "examples"
SEQ = ["m2_sequence_x1.0", "m2_sequence_x0.8", "m2_sequence_x0.6"]


@lru_cache(maxsize=None)
def scheme(stem):
    return build_scheme(load_spec(EX / f"{stem}.yaml"), EX, stem)


def union(boxes):
    boxes = list(boxes)
    return (min(b[0] for b in boxes), min(b[1] for b in boxes), max(b[2] for b in boxes), max(b[3] for b in boxes))


@pytest.mark.parametrize("stem", SEQ)
def test_numbers_share_a_baseline_and_center(stem):
    lay = scheme(stem)
    s = lay.preset.scheme
    bases = {round(r.y, 6) for runs in lay.numbers.values() for r in runs}
    assert len(bases) == 1 and len(lay.numbers) == 5
    bottom = max(p.group.ink_bbox()[3] for p in lay.species.values())
    tops = [union(r.bbox() for r in runs)[1] for runs in lay.numbers.values()]
    assert min(tops) - bottom == pytest.approx(s.number_gap, abs=1e-6)
    for name, runs in lay.numbers.items():
        b, ink = lay.species[name].group.ink_bbox(), union(r.bbox() for r in runs)
        assert (ink[0] + ink[2]) / 2 == pytest.approx((b[0] + b[2]) / 2, abs=1e-6)
        assert all(r.weight == "bold" for r in runs)


@pytest.mark.parametrize("stem", SEQ)
def test_conditions_centered_on_tail_to_notch(stem):
    lay = scheme(stem)
    s = lay.preset.scheme
    for st in lay.steps:
        mid = (st.tail[0] + notch_x(st.tail, st.head, lay.preset)) / 2
        for ln in st.above + st.below:
            assert ln.x0 + ln.advance / 2 == pytest.approx(mid, abs=1e-6)
        assert st.above[-1].baseline == pytest.approx(st.tail[1] - s.conditions_above_baseline)
        assert st.below[0].baseline == pytest.approx(st.tail[1] + s.conditions_below_baseline)


def test_multistep_left_justified():
    lay = scheme("m2_stages_x1.0")
    st = lay.steps[0]
    s = lay.preset.scheme
    assert st.staged and [ln.text[:2] for ln in st.above + st.below] == ["1.", "2."]
    for ln in st.above + st.below:
        assert ln.ink()[0] - st.tail[0] == pytest.approx(s.multistep_indent, abs=1e-6)
    assert any(r.style == "italic" and r.text == "n" for r in st.above[0].runs)


def test_three_stages():
    lay = scheme("m2_stages3_x1.0")
    st = lay.steps[0]
    s = lay.preset.scheme
    assert [ln.text[:2] for ln in st.above] == ["1.", "2."] and [ln.text[:2] for ln in st.below] == ["3."]
    assert st.above[0].baseline < st.above[1].baseline < st.tail[1] < st.below[0].baseline
    assert st.above[-1].baseline == pytest.approx(st.tail[1] - s.conditions_above_baseline)
    for ln in st.above + st.below:
        assert ln.ink()[0] - st.tail[0] == pytest.approx(s.multistep_indent, abs=1e-6)
        assert ln.x0 + ln.advance + s.conditions_pad <= notch_x(st.tail, st.head, lay.preset) + 1e-6


def test_two_line_pitch():
    lay = scheme("m2_suzuki_x1.0")
    st = lay.steps[0]
    assert st.above[1].baseline - st.above[0].baseline == pytest.approx(10.975)  # both lines subscripted
    assert st.below[1].baseline - st.below[0].baseline == pytest.approx(10.825)  # upper line only


def _one_step(above, preset="acs1996"):
    s = Spec.model_validate({"preset": preset, "species": {"a": {"smiles": "CCO"}, "b": {"smiles": "CC=O"}},
                             "steps": [{"from": "a", "to": "b", "above": above}]})
    return build_scheme(s, REPO, "t")


def test_arrow_grows_with_text():
    p = load_preset("acs1996")
    short = _one_step("PCC").steps[0]
    assert short.head[0] - short.tail[0] == pytest.approx(p.scheme.arrow_min_length)
    long = _one_step("a much longer reagent description").steps[0]
    adv = long.above[0].advance
    assert long.head[0] - long.tail[0] == pytest.approx(adv + 2 * p.scheme.conditions_pad + p.head().head_inset)


@pytest.mark.parametrize("stem", SEQ + ["m2_suzuki_x1.0", "m2_stages_x1.0", "m2_scheme_basic_x1.0",
                                        "m2_scheme_multi_x1.0"])
def test_no_overlaps(stem):
    assert check_scheme(scheme(stem)) == []


def test_overlap_check_catches_problems():
    lay = _one_step("PCC")
    lay.steps[0].above[0].shift(-40.0)  # conditions pushed off the arrow, onto the reactant
    problems = check_scheme(lay)
    assert any("wider than its arrow" in p for p in problems)
    assert any(p.startswith("species-a overlaps conditions-step-1") for p in problems)


def test_ring_centers_on_the_axis():
    """scheme_multi_x*: every structure has a ring, so ring centers sit center_above_axis up."""
    lay = scheme("m2_sequence_x1.0")
    for p in lay.species.values():
        assert any(c[1] == pytest.approx(-lay.preset.scheme.center_above_axis, abs=1e-6) for c in p.ring_centers)


def test_alignment_holds_matched_atoms():
    """Each product's atoms matched to its reference land on the same pt coordinates (0.01)."""
    spec = load_spec(EX / "m2_sequence_x1.0.yaml")
    prepared = scheme_prepared(spec, EX)
    for name, p in prepared.items():
        if not p.aligned_to:
            continue
        assert len(p.matched) >= 7, name  # the ring and OBn; a carbon whose neighbors change is laid out afresh (chem.stable_core)
        ref = prepared[p.aligned_to].molecule()
        m = p.molecule()
        for ri, mi in p.matched:
            assert math.dist(ref.atoms[ri].p, m.atoms[mi].p) < 0.01, (name, ri, mi)
    assert {n: p.aligned_to for n, p in prepared.items()} == {
        "ester": None, "bnbr": None, "ether": "ester", "alcohol": "ether", "aldehyde": "alcohol"}


def test_alignment_opt_out():
    spec = Spec.model_validate({"species": {"a": {"smiles": "OCc1ccccc1"},
                                            "b": {"smiles": "O=Cc1ccccc1", "align": False}},
                                "steps": [{"from": "a", "to": "b"}]})
    assert scheme_prepared(spec, REPO)["b"].aligned_to is None


def test_abbreviation_label():
    m = chem.from_smiles("COc1ccccc1", "anisole", 14.4, abbreviate=["OMe"])
    ab = [a for a in m.atoms if a.abbrev]
    assert len(ab) == 1 and ab[0].element == "O" and len(m.atoms) == 7
    with pytest.raises(ValueError, match="no such terminal group"):
        chem.from_smiles("CCO", "ethanol", 14.4, abbreviate=["OMe"])


@pytest.mark.parametrize("smiles, key, text, anchor", [("COc1ccccc1", "OMe", "MeO", "O"),
                                                       ("COC(=O)c1ccccc1", "CO2Me", "MeO2C", "C"),
                                                       ("C[Si](C)(C)c1ccccc1", "TMS", "TMS", "S")])
def test_right_justified_abbreviation(smiles, key, text, anchor):
    """A group with its bond leaving to the right is set in ChemDraw's reversed text, the last
    letter centered on the atom (abbrev_x*; formerly NotImplementedError)."""
    preset = load_preset("acs1996")
    m = chem.from_smiles(smiles, "a", 14.4, abbreviate=[key])
    i = next(a.idx for a in m.atoms if a.abbrev)
    ring = [a for a in m.atoms if not a.abbrev]
    if m.atoms[i].x > min(a.x for a in ring):  # put the group on the left of the ring
        for a in m.atoms:
            a.x = -a.x
    lab = Drawer(m, preset).labels[i]
    assert lab.h_side == "left" and "".join(r.text for r in lab.runs) == text and lab.element == anchor
    last = lab.runs[-1]
    end = last.x + last.advance()
    adv = lab.face.advance(text[-1], lab.size)
    assert end - adv / 2 == pytest.approx(m.atoms[i].x, abs=1e-9)
    draw_molecule(m, preset)


@pytest.mark.parametrize("stem", ["m1_smiles_x1.0", "m1_smiles_x0.8", "m1_smiles_x0.6"])
def test_m1_specs_unchanged(stem):
    """M1 specs (no steps) render byte-identical to the M1 layout code."""
    path = EX / f"{stem}.yaml"
    spec = load_spec(path)
    preset = load_preset(spec.preset)
    groups, x = [], 0.0
    for mol in molecules(spec, path.parent):  # scripts/render.py as committed in M1
        g = draw_molecule(mol, preset)
        x0, y0, x1, y1 = g.bbox()
        g.translate(x - x0, -(y0 + y1) / 2)
        x += (x1 - x0) + preset.scheme.plus_gap
        groups.append(g)
    assert to_svg(render_spec(path))[0] == to_svg(Group(path.stem, groups))[0]


# ---------------------------------------------------------------- arrow heads and equilibrium arrows


def _step(arrow_type="forward", above="PCC", below=None, preset="acs1996", arrow_head=None):
    d = {"preset": preset, "species": {"a": {"smiles": "CCO"}, "b": {"smiles": "CC=O"}},
         "steps": [{"from": "a", "to": "b", "above": above, "below": below, "arrow_type": arrow_type}]}
    if arrow_head:
        d["arrow_head"] = arrow_head
    return build_scheme(Spec.model_validate(d), REPO, "t")


def _marker_ids(lay):
    return {m.id for leaf in lay.root.leaves() for m in (getattr(leaf, "marker_start", None),
                                                         getattr(leaf, "marker_end", None)) if m}


@pytest.mark.parametrize("preset, head", [("acs1996", "large"), ("acs1996_x0.8", "standard"),
                                          ("acs1996_x0.6", "standard")])
def test_default_head_per_preset(preset, head):
    p = load_preset(preset)
    assert p.arrow.default_head == head and p.head() == getattr(p.arrow.heads, head)
    assert _marker_ids(_step(preset=preset)) == {f"{preset}-head-{head}"}


def test_arrow_head_override():
    lay = _step(arrow_head="standard")
    assert _marker_ids(lay) == {"acs1996-head-standard"}
    p = lay.preset
    lay = _step("equilibrium", above="a much longer reagent description", arrow_head="standard")
    st = lay.steps[0]
    assert st.head[0] - st.tail[0] == pytest.approx(st.above[0].advance + 2 * p.scheme.conditions_pad
                                                     + 2 * p.head("standard").head_inset)
    with pytest.raises(ValueError, match="arrow_head"):
        Spec.model_validate({"species": {"a": {"smiles": "C"}}, "arrow_head": "huge"})


@pytest.mark.parametrize("arrow_type, lengths", [("equilibrium", "equal"), ("equilibrium_products", "forward"),
                                                 ("equilibrium_reactants", "back")])
def test_equilibrium_step(arrow_type, lengths):
    lay = _step(arrow_type, above="H2SO4", below="heat")
    p, st = lay.preset, lay.steps[0]
    fwd, back = st.group.children
    assert [m.marker_end.id for m in (fwd, back)] == ["acs1996-head-large-half-forward", "acs1996-head-large-half-back"]
    off = p.equilibrium.shaft_spacing / 2
    assert [q[1] for q in fwd.points + back.points] == pytest.approx([-off, -off, off, off])
    assert fwd.points[1][0] > fwd.points[0][0] and back.points[1][0] < back.points[0][0]
    # each shaft ends shaft_end_from_tip before its tip; the tips are the arrow's ends (or centered)
    tips = [h.points[0][0] for h in fwd.head_leaves() + back.head_leaves()]
    L = st.head[0] - st.tail[0]
    if lengths == "equal":
        assert tips == [pytest.approx(st.head[0]), pytest.approx(st.tail[0])]
        assert fwd.points[0][0] - st.tail[0] == pytest.approx(p.equilibrium.stagger)
    else:
        long, short = (fwd, back) if lengths == "forward" else (back, fwd)
        h = p.head()
        short_len = abs(short.head_leaves()[0].points[0][0] - short.points[0][0])
        assert short_len == pytest.approx(h.head_inset + (L - h.head_inset) / p.ratios.equilibrium_long_per_short)
        assert (short.head_leaves()[0].points[0][0] + short.points[0][0]) / 2 == pytest.approx((st.tail[0] + st.head[0]) / 2)
        assert long.head_leaves()[0].points[0][0] in (pytest.approx(st.tail[0]), pytest.approx(st.head[0]))
    # conditions: centered on the arrow, baselines from the outer shafts
    s = p.scheme
    for ln in st.above + st.below:
        assert ln.x0 + ln.advance / 2 == pytest.approx((st.tail[0] + st.head[0]) / 2)
    assert st.above[-1].baseline == pytest.approx(-off - s.conditions_above_baseline)
    assert st.below[0].baseline == pytest.approx(off + s.conditions_below_baseline)
    assert not check_scheme(lay)


def test_equilibrium_half_heads_barb_outward():
    """Half heads (acs-style.md §3.6): tip and notch on the shaft's inner edge, barb outward; the
    back shaft's back-edge control back_control_shift nearer its center line."""
    from acs_chem_svg.arrows import head_marker
    p = load_preset("acs1996_x0.8")
    h, w = p.head(), p.arrow.shaft_width
    for half, shift in (("forward", 0.0), ("back", p.equilibrium.back_control_shift)):
        m = head_marker(p, half=half)
        (_, tip), (_, notch), (_, ctl, corner, _) = m.segments
        e = h.shaft_end_from_tip
        assert tip == (e, w / 2) and notch == pytest.approx((e - h.head_inset, w / 2))
        assert corner == pytest.approx((e - h.head_length, -h.head_half_width))
        assert ctl[1] == pytest.approx(-(p.ratios.arrow_back_control_per_half_width * h.head_half_width - shift))
