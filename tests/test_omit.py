"""By-products and reagents in mechanisms (user feedback 2026-10-28): `omit` keeps a product for
the electron bookkeeping but leaves it out of the drawing; `over` puts reactants above a step's
arrow, the first step's too (spec.py, layout.py, scripts/validate.py notes)."""

from functools import lru_cache
from pathlib import Path

import yaml
from validate import notes

from acs_chem_svg.electrons import check_spec
from acs_chem_svg.layout import build_scheme
from acs_chem_svg.spec import Spec, load_spec

EX = Path(__file__).resolve().parents[1] / "examples"
ACID_BASE = {
    "acid": {"smiles": "[CH3:1][C:2](=[O:3])[O:4][H:5]"},
    "hydroxide": {"smiles": "[OH-:6]"},
    "acetate": {"smiles": "[CH3:1][C:2](=[O:3])[O-:4]"},
    "water": {"smiles": "[H:5][OH:6]"},
}
DEPROTONATE = [{"from": "lp(6)", "to": "between(6,5)"}, {"from": "bond(4-5)", "to": "atom(4)"}]
REPROTONATE = [{"from": "lp(4)", "to": "between(4,5)"}, {"from": "bond(5-6)", "to": "atom(6)"}]


@lru_cache(maxsize=None)
def nitration():
    return build_scheme(load_spec(EX / "m4_nitration_x1.0.yaml"), EX, "nitration")


def drawn(lay):
    def ids(g):
        yield getattr(g, "id", "")
        for c in getattr(g, "children", []):
            yield from ids(c)
    return {i.removeprefix("species-") for i in ids(lay.root) if i.startswith("species-") and i.count("-") == 1}


def axis_y(lay, step):
    return next(st for st in lay.steps if st.n == step).tail[1]


def test_omitted_products_are_not_drawn_and_the_bookkeeping_still_passes():
    lay = nitration()
    assert all(not c.problems for c in lay.checks)
    assert {c.step: c.products for c in lay.checks}[2] == ("nitronium", "water")
    assert "water" not in drawn(lay) and "sulfuric_2" not in drawn(lay)
    assert "nitronium" in drawn(lay)


def test_over_on_the_first_step():
    lay = nitration()
    sulfuric = lay.species["sulfuric"].group.ink_bbox()
    st = next(s for s in lay.steps if s.n == 1)
    assert sulfuric[3] < axis_y(lay, 1) and st.tail[0] < (sulfuric[0] + sulfuric[2]) / 2 < st.head[0]
    row1 = [p for p in lay.pluses if abs(p.bbox()[1] - axis_y(lay, 1)) < 20]
    assert not [p for p in row1 if p.bbox()[2] < st.tail[0]]  # no "+ H2SO4" inline
    assert [a for a in lay.curved if a.step == 1]  # its curved arrows reach it there


def test_an_omitted_product_added_later_is_drawn_over_that_arrow():
    # bisulfate: omitted after step 1, then the base of step 5
    lay = nitration()
    b = lay.species["bisulfate"].group.ink_bbox()
    assert b[3] < axis_y(lay, 5)


def test_an_omitted_product_taken_by_the_next_step_goes_over_its_arrow():
    spec = Spec.model_validate({"species": ACID_BASE, "steps": [
        {"from": ["acid", "hydroxide"], "over": ["hydroxide"], "to": ["acetate", "water"], "omit": ["water"],
         "arrows": DEPROTONATE},
        {"from": ["acetate", "water"], "to": ["acid", "hydroxide"], "omit": ["hydroxide"], "arrows": REPROTONATE}]})
    assert spec.order() == ["acid", "hydroxide", "acetate", "water", "acid~2"]
    assert spec.over_names(1) == ["hydroxide"] and spec.over_names(2) == ["water"]
    assert all(not c.problems for c in check_spec(spec, EX))
    lay = build_scheme(spec, EX, "acid_base")
    assert drawn(lay) == {"acid", "hydroxide", "acetate", "water", "acid~2"}
    assert lay.species["water"].group.ink_bbox()[3] < axis_y(lay, 2)
    assert lay.species["hydroxide"].group.ink_bbox()[3] < axis_y(lay, 1)
    assert not lay.pluses


def test_by_products_after_a_bracket_respect_omit():
    species = {"p": {"smiles": "[CH3:1][CH:2]=[CH2:3]"}, "x": {"smiles": "[CH2+:1][CH:2]=[CH2:3]"},
               "y": {"smiles": "[CH2:1]=[CH:2][CH2+:3]"}, "h": {"smiles": "[H-]"}}
    steps = [{"from": "p", "to": ["x", "h"], "above": "hν"}, {"resonance": ["x", "y"]}]
    shown = build_scheme(Spec.model_validate({"species": species, "steps": steps}), EX, "s")
    steps[0] = {**steps[0], "omit": ["h"]}
    hidden = build_scheme(Spec.model_validate({"species": species, "steps": steps}), EX, "s")
    assert "h" in drawn(shown) and len(shown.pluses) == 1
    assert "h" not in drawn(hidden) and not hidden.pluses


def test_validator_notes_products_that_are_drawn_and_never_used():
    data = yaml.safe_load((EX / "m4_nitration_x1.0.yaml").read_text())
    assert notes(Spec.model_validate(data)) == []
    for st in data["steps"]:
        for k in ("omit", "over", "above"):
            st.pop(k, None)
    data["species"]["bisulfate_2"] = data["species"]["bisulfate"]
    data["steps"][-1]["from"] = ["sigma_c", "bisulfate_2"]
    got = notes(Spec.model_validate(data))
    assert [n.split(" is drawn")[0] for n in got] == [
        "note: step 1: bisulfate", "note: step 2: water", "note: step 5: sulfuric_2"]
    assert notes(load_spec(EX / "m3_sn2_x1.0.yaml")) == []  # one mechanism step: Br- stays


def test_validator_notes_a_loss_below_an_equilibrium_arrow():
    # user review 2026-10-28: "−AcOH" under ⇌ reads as AcOH lost in the reverse direction
    species = {"a": {"smiles": "CC(=O)O"}, "b": {"smiles": "CC(=O)[O-]"}}
    below = Spec.model_validate({"species": species, "steps": [
        {"from": "a", "to": "b", "arrow_type": "equilibrium", "below": "−H2O"}]})
    (n,) = notes(below)
    assert "'−H2O' below an equilibrium arrow" in n
    above = Spec.model_validate({"species": species, "steps": [
        {"from": "a", "to": "b", "arrow_type": "equilibrium", "above": "−H2O", "below": "H2O"}]})
    assert notes(above) == []
