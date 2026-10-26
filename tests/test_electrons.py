"""Electron bookkeeping (electrons.py, design.md §6.3): arrow semantics, charges, radicals."""

from pathlib import Path

import pytest

from acs_chem_svg import chem, electrons
from acs_chem_svg.spec import CurvedArrow, load_spec, parsed_species

ROOT = Path(__file__).resolve().parents[1]


def step(reactants, products, arrows):
    return electrons.check_step(1, {f"r{i}": chem.parse_smiles(s) for i, s in enumerate(reactants)},
                                {f"p{i}": chem.parse_smiles(s) for i, s in enumerate(products)},
                                [CurvedArrow.model_validate(a) for a in arrows])


SN2 = (["[OH-:1]", "[CH3:2][Br:3]"], ["[OH:1][CH3:2]", "[Br-:3]"],
       [{"from": "lp(1)", "to": "between(1,2)"}, {"from": "bond(2-3)", "to": "atom(3)"}])


def test_sn2_moves_charge_and_lone_pairs():
    c = step(*SN2)
    assert c.problems == []
    assert c.result == "[Br-:3].[OH:1][CH3:2]"
    assert c.lp_tails == {1} and c.gained_lone_pairs == {3}


def test_formal_charges_follow_the_electrons():
    model = electrons.build_model({"r": chem.parse_smiles("[C-:1]#[N:2]")})
    k1, k2 = model.by_map[1], model.by_map[2]
    assert (model.charge(k1), model.charge(k2)) == (-1, 0)
    assert (model.atoms[k1].lone_pairs, model.atoms[k2].lone_pairs) == (1, 1)


def test_homolysis_gives_two_radicals():
    c = step(["[Br:1][Br:2]"], ["[Br:1]", "[Br:2]"],
             [{"type": "fishhook", "from": "bond(1-2)", "to": "atom(1)"},
              {"type": "fishhook", "from": "bond(1-2)", "to": "atom(2)"}])
    assert c.problems == [] and c.result == "[Br:1].[Br:2]"


def test_radical_abstraction_with_explicit_h():
    c = step(["[Br:1]", "[CH3:2][C:3]([CH3:4])([CH3:5])[H:6]"], ["[Br:1][H:6]", "[CH3:2][C:3]([CH3:4])[CH3:5]"],
             [{"type": "fishhook", "from": "electron(1)", "to": "between(1,6)"},
              {"type": "fishhook", "from": "bond(3-6)", "to": "between(1,6)"},
              {"type": "fishhook", "from": "bond(3-6)", "to": "atom(3)"}])
    assert c.problems == []


def test_explicit_and_implicit_h_compare_equal():
    # E2: the mapped H moves; the products write water with explicit H, ethene with implicit H
    c = step(["[H][O-:1]", "[H:4][CH2:3][CH2:2][Br:5]"], ["[H:4][O:1][H]", "[CH2:3]=[CH2:2]", "[Br-:5]"],
             [{"from": "lp(1)", "to": "between(1,4)"}, {"from": "bond(3-4)", "to": "between(2,3)"},
              {"from": "bond(2-5)", "to": "atom(5)"}])
    assert c.problems == []


@pytest.mark.parametrize("arrows, message", [
    ([{"from": "lp(2)", "to": "between(1,2)"}], "lp(2): atom 2 (C) has no lone pairs"),
    ([{"from": "bond(2=3)", "to": "atom(3)"}], "joined by a single bond, not double"),
    ([{"from": "bond(1-3)", "to": "atom(3)"}], "atoms 1 (O) and 3 (Br) are not bonded"),
    ([{"from": "electron(1)", "to": "atom(2)", "type": "fishhook"}], "electron(1): atom 1 (O) has no unpaired electron"),
    ([{"from": "lp(1)", "to": "between(1,2)", "type": "fishhook"}], "cannot start at a lone pair"),
    ([{"from": "lp(1)", "to": "between(1,2)"}], "atom 2 (C) ends with 10 electrons (octet exceeded)"),
    ([{"from": "lp(9)", "to": "atom(2)"}], "no atom with map number 9"),
    ([{"type": "fishhook", "from": "bond(2-3)", "to": "atom(3)"}], "ends with 1 electrons"),
])
def test_errors(arrows, message):
    c = step(SN2[0], SN2[1], arrows)
    assert any(message in p for p in c.problems), c.problems


def test_wrong_product():
    c = step(SN2[0], ["[OH:1][CH3:2]", "[Br:3]"], SN2[2])
    assert c.problems == ["step 1: the arrows give [Br-:3].[OH:1][CH3:2], but the products are [Br:3].[OH:1][CH3:2]"]


def test_hydrogen_duet():
    c = step(["[H-:1]", "[H:2][H:3]"], ["[H-:1]", "[H:2][H:3]"], [{"from": "lp(1)", "to": "between(1,2)"}])
    assert any("atom 2 (H) ends with 4 electrons" in p for p in c.problems)


def test_duplicate_map_numbers():
    c = step(["[OH-:1]", "[CH3:1][Br:3]"], SN2[1], SN2[2])
    assert c.problems == ["step 1: map number 1 (r0, r1) is used more than once in the step's reactants"]


def test_endpoint_vocabulary():
    with pytest.raises(ValueError, match="an arrow ends at atom"):
        CurvedArrow.model_validate({"from": "lp(1)", "to": "lp(2)"})
    with pytest.raises(ValueError, match="starts at lp"):
        CurvedArrow.model_validate({"from": "atom(1)", "to": "atom(2)"})
    assert CurvedArrow.model_validate({"from": "bond(1#2)", "to": "atom(2)"}).source.order == 3


SPECS = sorted(p for g in ("m3_*.yaml", "m4_*.yaml") for p in (ROOT / "examples").glob(g)
               if load_spec(p).arrow_sets())


@pytest.mark.parametrize("path", SPECS, ids=[p.stem for p in SPECS])
def test_suggested_products_round_trip(path):
    """Every suggestion (validate.py --suggest) parses and equals the declared products."""
    spec = load_spec(path)
    mols = parsed_species(spec, path.parent)
    checks = electrons.check_spec(spec, path.parent)
    assert len(checks) == len(spec.arrow_sets())
    for c in checks:
        assert c.problems == [] and c.suggested
        assert [n for n, _ in c.suggested] == list(c.products)
        pasted = [chem.parse_smiles(s) for _, s in c.suggested]
        assert electrons.canonical(electrons.combine(pasted)) == c.expected
        assert electrons.canonical(electrons.combine([mols[n] for n in c.products])) == c.expected


def test_suggestion_names_an_unmatched_fragment():
    c = step(["[OH-:1]", "[CH3:2][Br:3]"], ["[OH:1][CH3:2]"],
             [{"from": "lp(1)", "to": "between(1,2)"}, {"from": "bond(2-3)", "to": "atom(3)"}])
    assert c.suggested == [("p0", "[OH:1][CH3:2]"), ("?", "[Br-:3]")]


def test_unmapped_atoms_keep_their_implicit_h():
    # r0 trials (radical_bromination, carbonyl_addition): plain `C` methyls lost their H in the
    # comparison and printed as [C]
    c = step(["[Br:1]", "C[C:3](C)(C)[H:6]"], ["[Br:1][H:6]", "C[C:3](C)C"],
             [{"type": "fishhook", "from": "electron(1)", "to": "between(1,6)"},
              {"type": "fishhook", "from": "bond(3-6)", "to": "between(1,6)"},
              {"type": "fishhook", "from": "bond(3-6)", "to": "atom(3)"}])
    assert c.problems == [] and "[C]" not in c.result
    c = step(["[C-:1]#[N:2]", "C[C:4](=[O:5])C"], ["[N:2]#[C:1][C:4](C)(C)[O-:5]"],
             [{"from": "lp(1)", "to": "between(1,4)"}, {"from": "bond(4=5)", "to": "atom(5)"}])
    assert c.problems == []


def test_a_mismatch_names_missing_map_numbers():
    # r0 carbonyl_addition: atom 3 mapped in the reactant but not in the product
    c = step(["[C-:1]#[N:2]", "C[C:3](=[O:4])C"], ["CC(C)([O-:4])[C:1]#[N:2]"],
             [{"from": "lp(1)", "to": "between(1,3)"}, {"from": "bond(3=4)", "to": "atom(4)"}])
    assert "map number(s) 3 of the reactants are missing from the products" in c.problems[0]
    assert c.suggested[0][1] == "C[C:3](C)([C:1]#[N:2])[O-:4]"


def test_lp_to_another_atom_points_to_between():
    # r1 sn2: lp(1) -> atom(2) instead of between(1,2); the suggestion alone misled
    c = step(["[OH-:1]", "[CH3:2][Br:3]"], ["[OH:1][CH3:2]", "[Br-:3]"],
             [{"from": "lp(1)", "to": "atom(2)"}, {"from": "bond(2-3)", "to": "atom(3)"}])
    assert 'arrows[0] (lp(1) -> atom(2)) puts the electrons on atom 2 as a lone pair' in c.problems[0]
    assert 'write to: "between(1,2)"' in c.problems[0]


def test_suggested_water_draws_both_h():
    # r1 e2: [OH:1][H:2] drew as HO-H
    c = step(["[H][O-:1]", "[H:4][CH2:3][CH2:2][Br:5]"], ["[H:4][O:1][H]", "[CH2:3]=[CH2:2]", "[Br-:5]"],
             [{"from": "lp(1)", "to": "between(1,4)"}, {"from": "bond(3-4)", "to": "between(2,3)"},
              {"from": "bond(2-5)", "to": "atom(5)"}])
    assert c.suggested[0] == ("p0", "[O:1]([H:4])[H]")
