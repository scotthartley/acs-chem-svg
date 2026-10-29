"""Carbenes and nitrenes (chem.nonbonding): two radicals in RDKit are a singlet's lone pair,
drawn always and usable as lp(n), unless the species lists the atom in triplet: (user decision
2026-10-29; placement unmeasured)."""

import math
from pathlib import Path

import pytest
from lxml import etree
from pydantic import ValidationError

from acs_chem_svg import chem, electrons
from acs_chem_svg.cdxml_writer import spec_cdxml
from acs_chem_svg.draw_mol import Drawer
from acs_chem_svg.geometry import dot, sub, unit
from acs_chem_svg.preset import load_preset
from acs_chem_svg.spec import CurvedArrow, Spec, load_spec, mark_triplets, molecules, parsed_species

ROOT = Path(__file__).resolve().parents[1]
P = load_preset("acs1996")


def mols(smiles: dict[str, str], triplet: dict[str, list] | None = None):
    out = {}
    for n, s in smiles.items():
        out[n] = chem.parse_smiles(s)
        mark_triplets(out[n], (triplet or {}).get(n, []))
    return out


def step(reactants, products, arrows, triplet=None):
    return electrons.check_step(1, mols(reactants, triplet), mols(products, triplet),
                                [CurvedArrow.model_validate(a) for a in arrows])


CYCLOPROPANATION = ({"carbene": "[Cl:2][C:1][Cl:3]", "alkene": "[CH2:7]=[CH2:8]"},
                    {"cp": "[Cl:2][C:1]1([Cl:3])[CH2:7][CH2:8]1"},
                    [{"from": "bond(7=8)", "to": "between(1,7)"}, {"from": "lp(1)", "to": "between(1,8)"}])
ALPHA = ({"anion": "[Cl:2][C-:1]([Cl:3])[Cl:4]"}, {"carbene": "[Cl:2][C:1][Cl:3]", "cl": "[Cl-:4]"},
         [{"from": "bond(1-4)", "to": "atom(4)"}])


# ---------------------------------------------------------------- bookkeeping


def test_singlet_starts_with_a_lone_pair():
    model = electrons.build_model(mols({"c": "Cl[C:1]Cl"}))
    a = model.atoms[model.by_map[1]]
    assert (a.lone_pairs, a.radicals, a.n) == (1, 0, 2)
    t = electrons.build_model(mols({"c": "Cl[C:1]Cl"}, {"c": [1]}))
    a = t.atoms[t.by_map[1]]
    assert (a.lone_pairs, a.radicals, a.n) == (0, 2, 2)


def test_nonbonding_helper():
    m = chem.parse_smiles("c1ccccc1[N:1]")
    n = next(a for a in m.GetAtoms() if a.GetAtomMapNum() == 1)
    assert chem.nonbonding(n) == (2, 0)
    mark_triplets(m, [1])
    assert chem.nonbonding(n) == (1, 2)


def test_cyclopropanation_passes():
    assert step(*CYCLOPROPANATION).problems == []


def test_alpha_elimination_passes():
    assert step(*ALPHA).problems == []


def test_lp_on_a_triplet_names_triplet():
    c = step(*CYCLOPROPANATION, triplet={"carbene": [1]})
    assert any("atom 1 (C) is a triplet (triplet: in carbene): start at electron(1), or drop it from triplet:" in p
               for p in c.problems), c.problems


def test_declared_triplet_product_is_a_spin_mismatch():
    c = step(*ALPHA, triplet={"carbene": [1]})
    assert c.result == c.expected
    assert any("the arrows leave a lone pair (singlet), but carbene lists it in triplet:; a step can't make a triplet"
               in p for p in c.problems), c.problems


def test_triplet_adds_to_an_alkene_with_fishhooks():
    c = step({"carbene": "[CH2:1]", "alkene": "[CH2:2]=[CH2:3]"}, {"diradical": "[CH2:1][CH2:2][CH2:3]"},
             [{"type": "fishhook", "from": "electron(1)", "to": "between(1,2)"},
              {"type": "fishhook", "from": "bond(2=3)", "to": "between(1,2)"},
              {"type": "fishhook", "from": "bond(2=3)", "to": "atom(3)"}], triplet={"carbene": [1]})
    assert c.problems == []


def test_cyclopropanation_example_validates():
    spec = load_spec(ROOT / "examples/m4_cyclopropanation_x1.0.yaml")
    assert [p for c in electrons.check_spec(spec, ROOT / "examples") for p in c.problems] == []


# ---------------------------------------------------------------- spec


def test_triplet_on_a_non_carbene_is_rejected():
    spec = Spec.model_validate({"species": {"m": {"smiles": "C[CH3]", "triplet": [1]}}})
    with pytest.raises(ValueError, match=r"species m: triplet: atom 1 \(C\) has 0 unpaired electrons; "
                                         r"triplet is for carbenes/nitrenes written \[C\] / \[CH2\] / \[N\]"):
        parsed_species(spec, ROOT)


def test_triplet_unknown_id_is_rejected():
    spec = Spec.model_validate({"species": {"m": {"smiles": "[CH2]", "triplet": [5]}}})
    with pytest.raises(ValueError, match="triplet: no atom with id 5"):
        molecules(spec, ROOT)


def test_triplet_is_a_list():
    with pytest.raises(ValidationError):
        Spec.model_validate({"species": {"m": {"smiles": "[CH2]", "triplet": 1}}})


# ---------------------------------------------------------------- drawing


def drawer(smiles, triplet=None, lewis=None):
    sp = {"smiles": smiles, **({"triplet": triplet} if triplet else {}), **({"lewis": lewis} if lewis else {})}
    (mol,) = molecules(Spec.model_validate({"species": {"m": sp}}), ROOT)
    return Drawer(mol, P)


def kinds(d):
    return sorted(ds.kind for sets in d.dot_sets().values() for ds in sets)


def test_methylene_draws_one_pair():
    d = drawer("[CH2]")
    assert kinds(d) == ["pair"] and d.mol.atoms[0].singlet


def test_triplet_draws_two_radicals():
    d = drawer("[CH2]", triplet=[1])
    assert kinds(d) == ["radical", "radical"] and not d.mol.atoms[0].singlet


def test_skeletal_triplet_dots_sit_on_either_side_of_the_bonds():
    """Ph₂C triplet: one dot in the wide gap, one between the bonds, so they don't read as a pair."""
    d = drawer("c1ccc(cc1)[C]c1ccccc1", triplet=[7])
    i = next(a.idx for a in d.mol.atoms if a.radicals == 2)
    first, second = d.dot_sets()[i]
    assert first.kind == second.kind == "radical"
    assert dot(first.direction, second.direction) == pytest.approx(-1)
    c = d.mol.atoms[i]
    to_bonds = [unit(sub(d.mol.atoms[j].p, c.p)) for j, _ in d.mol.neighbors(i)]
    assert all(dot(second.direction, u) > 0.4 for u in to_bonds)  # between the two bonds


@pytest.mark.parametrize("bonds, t0, want", [([], 0.0, 180), ([0.0], 90.0, 270), ([30.0, 150.0], 270.0, 90),
                                             ([0.0, 120.0, 240.0], 50.0, 180)])
def test_opposite_gap(bonds, t0, want):
    t = Drawer._opposite_gap(sorted(math.radians(b) for b in bonds), math.radians(t0))
    assert math.degrees(t) % 360 == pytest.approx(want)


def test_dichlorocarbene_pair_on_the_bisector_away_from_both_cl():
    d = drawer("Cl[C]Cl")
    i = next(a.idx for a in d.mol.atoms if a.element == "C")
    (ds,) = d.dot_sets()[i]
    assert ds.kind == "pair"
    c = d.mol.atoms[i]
    to_cl = [unit(sub(d.mol.atoms[j].p, c.p)) for j, _ in d.mol.neighbors(i)]
    assert all(dot(ds.direction, u) < -0.4 for u in to_cl)
    assert math.dist(ds.center, c.p) == pytest.approx(P.radical.distance)


def test_nhc_pair_points_out_of_the_ring():
    d = drawer("CN1C=CN(C)[C]1")
    i = next(a.idx for a in d.mol.atoms if a.singlet)
    (ds,) = d.dot_sets()[i]
    ring = next(r for r in d.mol.rings if i in r)
    cx = sum(d.mol.atoms[k].x for k in ring) / len(ring), sum(d.mol.atoms[k].y for k in ring) / len(ring)
    assert dot(ds.direction, unit(sub(d.mol.atoms[i].p, cx))) > 0.99


def test_nitrene_other_pair_follows_show_lone_pairs():
    d = drawer("c1ccccc1[N]")
    i = next(a.idx for a in d.mol.atoms if a.element == "N")
    assert d.mol.atoms[i].lone_pairs == 2 and kinds(d) == ["pair"]


def test_lewis_draws_the_carbene_pair_once():
    d = drawer("Cl[C]Cl", lewis=True)
    i = next(a.idx for a in d.mol.atoms if a.element == "C")
    assert [ds.kind for ds in d.dot_sets()[i]] == ["pair"]


def test_terminal_carbene_with_show_carbons_draws_ch():
    spec = Spec.model_validate({"species": {"m": {"smiles": "CC[CH]", "show_carbons": "terminal"}}})
    (mol,) = molecules(spec, ROOT)
    d = Drawer(mol, P)
    i = next(a.idx for a in mol.atoms if a.singlet)
    assert i in d.labels and [ds.kind for ds in d.dot_sets()[i]] == ["pair"]


# ---------------------------------------------------------------- CDXML


def test_cdxml_singlet_and_triplet():
    text = spec_cdxml(ROOT / "examples/m1_carbenes_x1.0.yaml", warn=lambda m: None)
    frags = dict(zip(load_spec(ROOT / "examples/m1_carbenes_x1.0.yaml").species,
                     etree.fromstring(text.encode()).find("page").findall("fragment")))

    def symbols(f, kind):
        return [g for g in f.findall("graphic") if g.get("SymbolType") == kind]
    for name in ("methylene", "dichlorocarbene", "diphenylcarbene", "nhc"):
        (n,) = [n for n in frags[name].findall("n") if n.get("Radical")]
        assert n.get("Radical") == "Singlet"
        (lp,) = symbols(frags[name], "LonePair")
        assert lp.find("represent").get("object") == n.get("id")
        assert not symbols(frags[name], "Electron")
    (n,) = [n for n in frags["diphenylcarbene_t"].findall("n") if n.get("Radical")]
    assert n.get("Radical") == "Singlet" and not symbols(frags["diphenylcarbene_t"], "LonePair")
    assert [g.find("represent").get("object") for g in symbols(frags["diphenylcarbene_t"], "Electron")] == [n.get("id")] * 2

