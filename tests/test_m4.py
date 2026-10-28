"""M4: resonance steps, several sequences and wrapped rows (spec.py, electrons.py, layout.py,
arrows.resonance_arrow, brackets; acs-style.md §3.4, §4.3)."""

from functools import lru_cache
from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError
from rdkit import Chem
from validate import validate

from acs_chem_svg.arrows import resonance_arrow
from acs_chem_svg.brackets import square_bracket
from acs_chem_svg.electrons import check_resonance
from acs_chem_svg.layout import build_scheme
from acs_chem_svg.preset import load_preset
from acs_chem_svg.scene import TextRun
from acs_chem_svg.spec import ResonanceStep, Spec, load_spec

REPO = Path(__file__).resolve().parents[1]
EX = REPO / "examples"
SPECIES = {n: {"smiles": s} for n, s in {"a": "CCO", "b": "CC=O", "w": "O", "c": "CC(=O)O", "d": "CC", "e": "C=C"}.items()}
ALLYL = {"p": {"smiles": "[CH3:1][CH:2]=[CH2:3]"}, "x": {"smiles": "[CH2+:1][CH:2]=[CH2:3]"},
         "y": {"smiles": "[CH2:1]=[CH:2][CH2+:3]"}}
SHIFT = {"from": "bond(2=3)", "to": "between(1,2)"}


def spec(steps, species=None, **kw):
    return Spec.model_validate({"species": species or SPECIES, "steps": steps, **kw})


@lru_cache(maxsize=None)
def scheme(stem):
    return build_scheme(load_spec(EX / f"{stem}.yaml"), EX, stem)


# ---------------------------------------------------------------- schema


def test_by_products_are_left_behind():
    s = spec([{"from": "a", "to": ["b", "w"]}, {"from": "b", "to": "c"}],
             species={k: SPECIES[k] for k in "abwc"})
    assert s.order() == ["a", "b", "w", "c"] and s.kept(2) == 1


def test_extra_reactants_after_by_products():
    s = spec([{"from": "a", "to": ["b", "w"]}, {"from": ["b", "d"], "to": "c"}],
             species={k: SPECIES[k] for k in "abwdc"})
    assert s.order() == ["a", "b", "w", "d", "c"]


def test_new_row_starts_a_sequence():
    s = spec([{"from": "a", "to": "b"}, {"new_row": True, "from": "d", "to": "e"}],
             species={k: SPECIES[k] for k in "abde"})
    assert s.order() == ["a", "b", "d", "e"] and s.kept(2) == 0
    with pytest.raises(ValidationError, match="new_row needs `from`"):
        spec([{"from": "a", "to": "b"}, {"new_row": True, "to": "e"}], species={k: SPECIES[k] for k in "abe"})


def test_resonance_step_chains_on_its_last_contributor():
    s = spec([{"from": "p", "to": "x"}, {"resonance": ["x", "y"], "arrows": [[SHIFT]]}], species=ALLYL)
    assert isinstance(s.steps[1], ResonanceStep)
    assert s.order() == ["p", "x", "y"] and s.steps[1].products == ["y"]
    (a,) = s.arrow_sets()
    assert (a.step, a.sources, a.products, a.contributor) == (2, ("x",), ("y",), 0)


def test_resonance_rules():
    with pytest.raises(ValidationError, match="previous step's first product"):
        spec([{"from": "p", "to": "x"}, {"resonance": ["y", "x"]}], species=ALLYL)
    with pytest.raises(ValidationError, match="one list per contributor"):
        spec([{"resonance": ["x", "y"], "arrows": [[SHIFT], [SHIFT]]}], species={k: ALLYL[k] for k in "xy"})
    with pytest.raises(ValidationError):
        spec([{"resonance": ["x"]}], species={"x": ALLYL["x"]})
    with pytest.raises(ValidationError):  # no conditions on a resonance step
        spec([{"resonance": ["x", "y"], "above": "hv"}], species={k: ALLYL[k] for k in "xy"})


def test_resonance_after_a_step_with_by_products(tmp_path):
    sp = {**ALLYL, "w": {"smiles": "O"}}
    s = spec([{"from": "p", "to": ["x", "w"]}, {"resonance": ["x", "y"], "arrows": [[SHIFT]]}], species=sp)
    assert s.order() == ["p", "x", "w", "y"] and s.steps[1].sources == ["x"]
    lay = build_scheme(s, EX, "byproduct")
    (res,) = lay.resonance
    w = lay.species["w"].group.ink_bbox()
    assert w[0] > res.bar_x[1] and len(lay.pluses) == 1 and lay.pluses[0].x > res.bar_x[1]
    path = tmp_path / "byproduct.yaml"
    path.write_text(yaml.safe_dump({"species": sp, "steps": [{"from": "p", "to": ["x", "w"]},
                                                             {"resonance": ["x", "y"], "arrows": [[SHIFT]]}]}))
    assert validate(path) == []


def test_phenol_h_abstraction_by_x():
    path = REPO / "tests" / "specs" / "phenol_x_abstraction.yaml"
    assert validate(path) == []
    lay = build_scheme(load_spec(path), path.parent, "phenol")
    texts = {n: "".join(r.text for r in lay.species[n].group.leaves() if isinstance(r, TextRun)) for n in ("xrad", "hx")}
    assert texts == {"xrad": "X", "hx": "HX"}
    assert lay.species["hx"].group.ink_bbox()[0] > lay.resonance[0].bar_x[1]


def test_resonance_override_index_counts_all_lists():
    sp = {**{k: ALLYL[k] for k in "xy"}, "z": {"smiles": "[CH2+:1][CH:2]=[CH2:3]"}}
    back = {"from": "bond(1=2)", "to": "between(2,3)"}
    s = spec([{"resonance": ["x", "y", "z"], "arrows": [[SHIFT], [back]]}], species=sp,
             overrides={"arrows": [{"step": 1, "index": 1, "bow": "left"}]})
    assert [a.offset for a in s.arrow_sets()] == [0, 1]
    with pytest.raises(ValidationError, match="no arrow with index 2"):
        spec([{"resonance": ["x", "y", "z"], "arrows": [[SHIFT], [back]]}], species=sp,
             overrides={"arrows": [{"step": 1, "index": 2}]})


# ---------------------------------------------------------------- bookkeeping


def _mols(**smiles):
    return {k: Chem.MolFromSmiles(v) for k, v in smiles.items()}


def test_check_resonance_passes_and_names_the_contributor():
    ok = check_resonance(1, _mols(x="[CH2+:1][CH:2]=[CH2:3]", y="[CH2:1]=[CH:2][CH2+:3]"),
                         [[ResonanceStep.model_validate({"resonance": ["x", "y"], "arrows": [[SHIFT]]}).arrows[0][0]]])
    assert [c.problems for c in ok] == [[]]
    bad = check_resonance(2, _mols(x="[CH2+:1][CH:2]=[CH2:3]", y="[CH2-:1][CH:2]=[CH2:3]"), [])
    assert "contributor y (2 of 2) has total charge -1, but x has +1" in bad[0].problems[0]
    moved = check_resonance(2, _mols(x="[CH2+:1][CH:2]=[CH2:3]", y="[CH2:1]=[CH:2][CH3:3]"), [])
    assert "different connectivity" in moved[0].problems[0]


# ---------------------------------------------------------------- primitives


def test_resonance_arrow_has_a_head_at_each_end():
    p = load_preset()
    g = resonance_arrow((0.0, 0.0), (28.8, 0.0), p)
    (shaft,) = g.children
    h1, h2 = shaft.head_leaves()
    a = p.head()
    assert shaft.points == [(a.shaft_end_from_tip, 0.0), (pytest.approx(28.8 - a.shaft_end_from_tip), 0.0)]
    assert h1.points[0] == pytest.approx((0.0, 0.0)) and h2.points[0] == pytest.approx((28.8, 0.0))
    xs = [q[0] for q in h1.points]
    assert max(xs) == pytest.approx(a.head_length)


def test_square_bracket_lips_point_inward():
    p = load_preset()
    left = square_bracket(10.0, 0.0, 50.0, "left", p).children[0]
    right = square_bracket(90.0, 0.0, 50.0, "right", p).children[0]
    lip = p.bracket.lip_length
    assert left.points == [(10.0 + lip, 0.0), (10.0, 0.0), (10.0, 50.0), (10.0 + lip, 50.0)]
    assert right.points[0] == (90.0 - lip, 0.0) and right.width == p.bracket.stroke_width


# ---------------------------------------------------------------- layout


@pytest.mark.parametrize("stem", ["m4_resonance_x1.0", "m4_resonance_x0.8", "m4_resonance_x0.6"])
def test_resonance_gaps(stem):
    lay = scheme(stem)
    p = lay.preset
    (res,) = lay.resonance
    inks = [lay.species[n].group.ink_bbox() for n in res.contributors]
    for (tail, head, _), left, right in zip(res.arrows, inks, inks[1:]):
        assert head[0] - tail[0] == pytest.approx(p.resonance.arrow_length)
        assert tail[0] - left[2] == pytest.approx(p.resonance.arrow_gap)
        assert right[0] - head[0] == pytest.approx(p.resonance.arrow_gap)
    content = (inks[0][0], min(b[1] for b in inks), inks[-1][2], max(b[3] for b in inks))
    lb, rb = (g.children[0].points for g in res.group.children)
    assert content[0] - lb[1][0] == pytest.approx(p.bracket.gap_x)
    assert rb[1][0] - content[2] == pytest.approx(p.bracket.gap_x)
    assert content[1] - lb[1][1] == pytest.approx(p.bracket.gap_above)
    assert lb[2][1] - content[3] == pytest.approx(p.bracket.gap_below)


@pytest.mark.parametrize("stem", ["m4_wrap_x1.0", "m4_wrap_x0.8", "m4_wrap_x0.6"])
def test_wrapped_rows(stem):
    lay = scheme(stem)
    rows = [g for g in lay.root.children]
    assert [g.id for g in rows] == ["row-1", "row-2"]
    assert [st.row for st in lay.steps] == [0, 0, 1, 1, 1]
    r1, r2 = (g.ink_bbox() for g in rows)
    assert r2[1] - r1[3] == pytest.approx(lay.preset.scheme.row_gap)
    st = lay.steps[2]  # the continuation row starts with its arrow at the left margin
    assert st.tail[0] == pytest.approx(lay.species["ester"].group.ink_bbox()[0])
    assert max(r1[2], r2[2]) - r1[0] <= load_spec(EX / f"{stem}.yaml").layout.width
    assert {n for n in lay.numbers} == set("ester bnbr ether alcohol aldehyde enoate ester_h2".split())


def test_default_width_is_the_preset_wrap_width():
    s = load_spec(EX / "m4_wrap_x1.0.yaml")
    s.layout.width = None
    lay = build_scheme(s, EX, "wide")
    assert lay.preset.scheme.wrap_width == 504.0 and lay.width <= 504.0


def test_a_step_after_a_bracket_reaches_into_it():
    lay = scheme("m4_nitration_x1.0")
    (res,) = lay.resonance
    assert res.contributors == ["sigma_a", "sigma_b", "sigma_c"]
    deprot = [a for a in lay.curved if a.step == 5]
    assert len(deprot) == 2


@pytest.mark.parametrize("stem", ["m4_nitration_x1.0", "m4_nitration_x0.8", "m4_nitration_x0.6",
                                  "m4_resonance_x1.0", "m4_wrap_x1.0", "m4_perkin_x1.0"])
def test_m4_examples_validate(stem):
    assert validate(EX / f"{stem}.yaml") == []
