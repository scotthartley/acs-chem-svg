"""Abbreviation grammar (abbrev.py, acs-style.md §3.5): reversal, anchors and composition.

MEASURED holds ChemDraw's displayed right-justified text for every label in abbrev_x*,
abbrev_2_x*, ligands_x* and metals_x* (measure_refs.py --abbrev); test_measured_table re-reads it
from the references when they are present. reverse(key) must equal it, except the two user
decisions of 2026-11-03 (DECIDED)."""

import pytest

from acs_chem_svg import abbrev, chem
from acs_chem_svg.spec import Spec

MEASURED = {
    # abbrev_x*
    "OMe": "MeO", "OEt": "EtO", "OBn": "BnO", "OPh": "PhO", "OAc": "AcO", "SMe": "MeS", "NHAc": "AcHN",
    "CO2Me": "MeO2C", "CO2Et": "EtO2C", "CO2H": "HO2C", "CF3": "F3C", "CN": "NC", "NO2": "O2N",
    "Ph": "Ph", "Bn": "Bn", "Et": "Et", "Ac": "Ac", "TMS": "TMS", "TES": "TES", "TBS": "TBS", "TIPS": "TIPS",
    # abbrev_2_x*
    "OTf": "TfO", "OTs": "TsO", "OTBS": "TBSO", "NHBoc": "BocHN", "NHTs": "TsHN", "Boc": "Boc", "Ts": "Ts",
    "B(OH)2": "(HO)2B", "t-Bu": "t-Bu", "i-Pr": "i-Pr", "Me": "Me", "OBz": "BzO", "OPiv": "PivO",
    "OMOM": "MOMO", "OTHP": "THPO", "OPMB": "PMBO", "OTIPS": "TIPSO", "Cbz": "Cbz", "Fmoc": "Fmoc",
    "SPh": "PhS", "SO2Ph": "PhO2S", "NMe2": "Me2N", "N3": "N3", "CHO": "OHC", "CO2t-Bu": "t-BuO2C",
    "SnBu3": "Bu3Sn", "n-Bu": "n-Bu", "Cy": "Cy", "Mes": "Mes", "NTs": "TsN", "p-Tol": "p-Tol",
    "Bpin": "pinB",  # typed by the user as pinB (Bpin's reversed: override)
    # ligands_x*, metals_x*
    "PMe3": "Me3P", "PCy3": "Cy3P", "P(t-Bu)3": "(t-Bu)3P", "Ln": "Ln", "dppf": "dppf",
    "XPhos": "XPhos", "PPh3": "Ph3P", "IPr": "PrI", "P(o-tol)3": "P(o-tol)3",
}
# user decisions (2026-11-03): one rule, not ChemDraw's reading of labels it can't interpret
DECIDED = {"IPr": "IPr", "P(o-tol)3": "(o-tol)3P"}


@pytest.mark.parametrize("key", sorted(MEASURED))
def test_reverse_matches_chemdraw(key):
    assert abbrev.reverse(key) == DECIDED.get(key, MEASURED[key])


def test_measured_table():
    """MEASURED is what ChemDraw shows (re-read from the references when present)."""
    import measure_refs as mr
    if not (mr.REF_DIR / "abbrev_2_x1.0.cdxml").exists():
        pytest.skip("ChemDraw references (reference_renders/) are local only")
    shown = {}
    for ref in ("abbrev_x1.0", "abbrev_2_x1.0", "ligands_x1.0", "metals_x1.0"):
        for f in mr.abbrev_summary(ref)["forms"]:
            if f["justification"] == "Right" and not f.get("alignment"):
                shown[f["key"]] = f["text"]
    shown["Bpin"] = shown.pop("pinB")
    assert shown == MEASURED


def test_markup_keeps_italics_and_subscripts():
    assert abbrev.markup("CO2t-Bu", "left") == "*t*-BuO2C"
    assert abbrev.markup("P(t-Bu)3", "left") == "(*t*-Bu)3P"
    assert abbrev.markup("Ln", "right") == "L_*n*"


@pytest.mark.parametrize("key, element, left, right", [
    ("OMe", "O", "O", "O"), ("CO2Me", "C", "C", "C"), ("NO2", "N", "N", "N"), ("Ph", "C", "P", "h"),
    ("TMS", "Si", "T", "S"), ("N3", "N", "N", "3"), ("SnBu3", "Sn", "Sn", "Sn"), ("NTs", "N", "N", "N"),
    ("Bpin", "B", "B", "B"), ("dppf", "*", "d", "f"),
])
def test_anchor(key, element, left, right):
    """The glyphs on the atom: the attachment element at the centered end, else the end glyph."""
    assert abbrev.anchor(key, element, "right") == left
    assert abbrev.anchor(key, element, "left") == right


def test_unknown_token_names_the_label():
    with pytest.raises(ValueError, match=r"label 'OXyz': can't read 'Xyz'"):
        abbrev.reverse("OXyz")


@pytest.mark.parametrize("element, h, groups, key", [
    ("N", 1, ["Boc"], "NHBoc"), ("N", 0, ["Ts"], "NTs"), ("Sn", 0, ["Bu"] * 3, "SnBu3"),
    ("P", 0, ["t-Bu"] * 3, "P(t-Bu)3"), ("O", 0, ["TBS"], "OTBS"), ("N", 0, ["Me", "Me"], "NMe2"),
])
def test_compose_key(element, h, groups, key):
    assert abbrev.compose(element, h, groups) == key


def _mol(smiles, abbreviate, extra=None):
    reg = abbrev.default().with_entries(extra)
    return chem.from_smiles(smiles, "m", 14.4, abbreviate=abbreviate, registry=reg)


@pytest.mark.parametrize("smiles, abbreviate, key, degree", [
    ("CC(C)(C)OC(=O)Nc1ccccc1", ["Boc"], "NHBoc", 1),
    ("Cc1ccc(cc1)S(=O)(=O)N1CCCC1", ["Ts"], "NTs", 2),
    ("CCCC[Sn](CCCC)(CCCC)c1ccccc1", ["Bu"], "SnBu3", 1),
    ("CC(C)(C)P(C(C)(C)C)(C(C)(C)C)->[Pd]Cl", ["t-Bu"], "P(t-Bu)3", 1),
    ("CC(C)(C)[Si](C)(C)OCc1ccccc1", ["TBS"], "OTBS", 1),
])
def test_composition(smiles, abbreviate, key, degree):
    m = _mol(smiles, abbreviate)
    (a,) = [a for a in m.atoms if a.abbrev]
    assert a.abbrev == key and m.degree(a.idx) == degree and a.abbrev_info.kind == "composed"
    assert a.element == key[0] if key[:2] not in abbrev.ELEMENTS else key[:2]


def test_formula_entry_wins_over_composition():
    m = _mol("COC(=O)c1ccccc1", ["CO2Me"])
    assert [a.abbrev_info.kind for a in m.atoms if a.abbrev] == ["formula"]


def test_metal_keeps_its_groups():
    """Pd is not composed: Ph stays its own label on Pd (metals_x*)."""
    m = _mol("Br[Pd](c1ccccc1)(<-P(C)(C)C)<-P(C)(C)C", ["Ph", "PMe3"])
    assert sorted(a.abbrev for a in m.atoms if a.abbrev) == ["PMe3", "PMe3", "Ph"]


def test_spec_group_composes():
    m = _mol("COc1ccccc1", ["Xyz"], {"Xyz": {"smiles": "*C", "kind": "group"}})
    assert [a.abbrev for a in m.atoms if a.abbrev] == ["OXyz"]


def test_text_label_through_labels():
    from pathlib import Path

    from acs_chem_svg.spec import molecules
    spec = Spec.model_validate({"species": {"pd": {"smiles": "[*:1][Pd]([*:2])([*:3])[*:4]"}},
                                "labels": {1: "Ln", 2: "dppf", 3: "L", 4: "X"}})
    (m,) = molecules(spec, Path("."))
    got = {a.map: (a.abbrev, a.element) for a in m.atoms if a.map}
    assert got[1] == ("Ln", "*") and got[2] == ("dppf", "*") and got[3] == ("L", "*") and got[4] == (None, "X")
