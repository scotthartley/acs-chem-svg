"""Deliberately broken specs (design.md §10 M3): each must fail validation with the message on
its first line (`# expect: ...`)."""

from pathlib import Path

import pytest
from validate import validate

BROKEN = sorted((Path(__file__).parent / "broken").glob("*.yaml"))


def expected(path: Path) -> str:
    first = path.read_text().splitlines()[0]
    assert first.startswith("# expect: "), path.name
    return first.removeprefix("# expect: ")


@pytest.mark.parametrize("path", BROKEN, ids=[p.stem for p in BROKEN])
def test_broken_spec_is_caught(path):
    problems = validate(path)
    assert any(expected(path) in p for p in problems), problems


def test_every_planned_failure_has_a_spec():
    assert {p.stem for p in BROKEN} >= {
        "wrong_product", "lp_without_lone_pair", "wrong_bond_order", "fishhook_from_lone_pair",
        "carbon_with_ten_electrons", "duplicate_map_number", "arrow_crosses_bond", "lone_pair_on_bond",
        "resonance_wrong_connectivity", "resonance_wrong_arrow", "from_skips_products", "resonance_first_not_product",
        "row_too_wide", "label_unknown_map"}


def test_wrong_product_prints_a_paste_ready_suggestion(capsys):
    import validate as v
    assert v.main([str(Path(__file__).parent / "broken" / "wrong_product.yaml")]) == 1
    out = capsys.readouterr().out
    assert 'suggested: bromine: {smiles: "[Br-:3]"}' in out
    assert 'suggested: methanol: {smiles: "[OH:1][CH3:2]"}' in out
