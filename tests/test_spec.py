"""Scheme spec validation (spec.py, references/scheme-schema.md)."""

import pytest
from pydantic import ValidationError

from acs_chem_svg.spec import Spec

SPECIES = {"a": {"smiles": "CCO"}, "b": {"smiles": "CC=O"}, "c": {"smiles": "CC(=O)O"}, "d": {"smiles": "CC"}}


def spec(steps, species=None, **kw):
    return Spec.model_validate({"species": species or SPECIES, "steps": steps, **kw})


def test_chain_and_defaults():
    s = spec([{"from": "a", "to": "b", "above": "PCC"}, {"to": "c"}],
             species={k: SPECIES[k] for k in "abc"}, numbers={"a": 1, "c": "4a"})
    assert s.steps[1].sources == ["b"]  # from: defaults to the previous `to`
    assert s.order() == ["a", "b", "c"]


def test_extra_reactant_joins_chain():
    s = spec([{"from": "a", "to": "b"}, {"from": ["b", "d"], "to": "c"}])
    assert s.order() == ["a", "b", "d", "c"]


def test_broken_chain():
    with pytest.raises(ValidationError, match="linear chain"):
        spec([{"from": "a", "to": "b"}, {"from": "c", "to": "d"}])


def test_first_step_needs_from():
    with pytest.raises(ValidationError, match="needs `from`"):
        spec([{"to": "b"}])


def test_unknown_names():
    with pytest.raises(ValidationError, match="unknown species 'x'"):
        spec([{"from": "a", "to": "x"}])
    with pytest.raises(ValidationError, match="numbers: unknown species 'x'"):
        spec([{"from": "a", "to": "b"}], species={k: SPECIES[k] for k in "ab"}, numbers={"x": 1})
    with pytest.raises(ValidationError, match="align_to: unknown species"):
        spec([{"from": "a", "to": "b"}], species={"a": SPECIES["a"], "b": {"smiles": "CC=O", "align_to": "q"}})


def test_unused_and_repeated_species():
    with pytest.raises(ValidationError, match="not used in any step: c, d"):
        spec([{"from": "a", "to": "b"}])
    # M5: a repeated name is drawn again, as a copy (test_m5.py)
    s = spec([{"from": "a", "to": "b"}, {"to": "a"}], species={k: SPECIES[k] for k in "ab"})
    assert s.order() == ["a", "b", "a~2"] and s.steps[1].products == ["a~2"]


def test_stages():
    s = spec([{"from": "a", "to": "b", "stages": ["*n*-BuLi", "H3O+"]}], species={k: SPECIES[k] for k in "ab"})
    assert s.steps[0].lines() == (["1. *n*-BuLi"], ["2. H3O+"])
    with pytest.raises(ValidationError, match="1 to 2 are supported"):
        spec([{"from": "a", "to": "b", "stages": ["x", "y", "z"]}], species={k: SPECIES[k] for k in "ab"})
    with pytest.raises(ValidationError, match="either stages or above/below"):
        spec([{"from": "a", "to": "b", "stages": ["x"], "above": "y"}], species={k: SPECIES[k] for k in "ab"})


def test_unknown_abbreviation():
    with pytest.raises(ValidationError, match="unknown abbreviation OTf"):
        Spec.model_validate({"species": {"a": {"smiles": "CCO", "abbreviate": ["OTf"]}}})


def test_m1_spec_has_no_steps():
    s = Spec.model_validate({"species": SPECIES})
    assert s.steps is None and s.order() == list(SPECIES)
    with pytest.raises(ValidationError, match="numbers need steps"):
        Spec.model_validate({"species": SPECIES, "numbers": {"a": 1}})


def test_curved_arrows_and_overrides():
    sp = {k: SPECIES[k] for k in "ab"}
    s = spec([{"from": "a", "to": "b", "arrows": [{"from": "lp(1)", "to": "between(1,2)"}]}], species=sp,
             overrides={"arrows": [{"step": 1, "index": 0, "bow": "left", "curvature": 0.6}]})
    assert s.steps[0].arrows[0].type == "two_electron" and s.arrow_override(1, 0).bow == "left"
    with pytest.raises(ValidationError, match="expected lp"):
        spec([{"from": "a", "to": "b", "arrows": [{"from": "lone(1)", "to": "atom(2)"}]}], species=sp)
    with pytest.raises(ValidationError, match="no arrow with index 1"):
        spec([{"from": "a", "to": "b", "arrows": [{"from": "lp(1)", "to": "atom(2)"}]}], species=sp,
             overrides={"arrows": [{"step": 1, "index": 1, "bow": "left"}]})
    with pytest.raises(ValidationError, match="raw|extra|Extra"):
        spec([{"from": "a", "to": "b", "arrows": [{"from": "lp(1)", "to": "atom(2)"}]}], species=sp,
             overrides={"arrows": [{"step": 1, "index": 0, "x": 3.0}]})  # no raw coordinates


def test_labels():
    sp = {"x": {"smiles": "[Cl:9]"}, "hx": {"smiles": "[ClH:9]"}}
    assert spec([{"from": "x", "to": "hx"}], species=sp, labels={9: "X"}).labels == {9: "X"}
    for bad in ("x", "X1", "Cl-", "XYZ", ""):
        with pytest.raises(ValidationError, match="plain symbol"):
            spec([{"from": "x", "to": "hx"}], species=sp, labels={9: bad})


def test_align_on_a_step_points_to_the_species():
    with pytest.raises(ValidationError, match=r"align belongs on the species \(species.<name>.align\), not on a step"):
        spec([{"from": "a", "to": "b", "align": False}], species={k: SPECIES[k] for k in "ab"})
    with pytest.raises(ValidationError, match="orientation belongs on the species"):
        spec([{"from": "a", "to": "b", "orientation": {"flip_h": True}}], species={k: SPECIES[k] for k in "ab"})
