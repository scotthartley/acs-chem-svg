"""Dative bonds (P->[Pd]): drawn and laid out as single bonds, the donated pair
counted on the donor (chem.bond_order_sum), so a ligand shows no radical or extra lone pair,
and compared by the electron check as the P+–Pd− bond the arrows make."""

import math

import pytest
from lxml import etree
from rdkit import Chem

from acs_chem_svg import chem, electrons
from acs_chem_svg.figure import spec_cdxml
from acs_chem_svg.preset import load_preset
from acs_chem_svg.spec import CurvedArrow

P = load_preset("acs1996")
PD0 = "CP(C)(C)->[Pd]<-P(C)(C)C"
OXADD = "Br[Pd](c1ccccc1)(<-P(C)(C)C)<-P(C)(C)C"


@pytest.mark.parametrize("smiles, element, lp", [(PD0, "P", 0), ("[NH3]->[Pt](Cl)(Cl)<-[NH3]", "N", 0),
                                                 ("[NH3]->[Pt](Cl)(Cl)<-[NH3]", "Cl", 3)])
def test_donor_keeps_no_pair_and_no_radical(smiles, element, lp):
    mol = chem.parse_smiles(smiles)
    assert {chem.nonbonding(a) for a in mol.GetAtoms() if a.GetSymbol() == element} == {(lp, 0)}


def test_acceptor_keeps_its_own_electrons():
    """Pd(II) with Br and Ph: 10 − 2 = 8 nonbonding electrons; the two P->Pd pairs are P's."""
    mol = chem.parse_smiles(OXADD)
    pd = next(a for a in mol.GetAtoms() if a.GetSymbol() == "Pd")
    assert chem.nonbonding(pd) == (4, 0)


def test_layout_is_quiet(capfd):
    chem.from_smiles(OXADD, "m", P.lengths.bond_length)
    assert "unrecognized bond type" not in capfd.readouterr().err


def _single(smiles: str) -> str:
    """The same species with every dative bond written single: the layout to match."""
    return smiles.replace("<-", "").replace("->", "")


def _xy(mol: Chem.Mol) -> list[tuple[float, float]]:
    c = mol.GetConformer()
    return [(c.GetAtomPosition(i).x, c.GetAtomPosition(i).y) for i in range(mol.GetNumAtoms())]


@pytest.mark.parametrize("smiles", [PD0, OXADD])
def test_dative_bonds_drawn_single_and_laid_out_as_single(smiles):
    """The single-bond layout, which layout_2d lays out to keep CoordGen quiet (to 1e-4 pt:
    each Pd(II) ligand is laid out on its own, chem._square_planar)."""
    m, ref = (chem.from_smiles(s, "m", P.lengths.bond_length) for s in (smiles, _single(smiles)))
    pd = next(a.idx for a in m.atoms if a.element == "Pd")
    assert {b.order for b in m.bonds if pd in (b.a, b.b)} == {1}
    assert [(a.element, a.radicals, a.lone_pairs) for a in m.atoms if a.element == "P"] == [("P", 0, 0)] * 2
    assert [a.p for a in m.atoms] == [pytest.approx(a.p, abs=1e-4) for a in ref.atoms]


def test_aligned_dative_species_laid_out_as_single():
    """align_coords (CoordGen with a coordMap) too."""
    def aligned(pd0, oxadd):
        ref = chem.prepare(chem.parse_smiles(pd0), "pd0", P.lengths.bond_length)
        out = chem.prepare(chem.parse_smiles(oxadd), "oxadd", P.lengths.bond_length, refs=[ref])
        assert out.aligned_to == "pd0"
        return _xy(out.mol)
    flat = [v for p in aligned(PD0, OXADD) for v in p]
    assert flat == pytest.approx([v for p in aligned(_single(PD0), _single(OXADD)) for v in p], abs=1e-6)


def test_cdxml_writes_a_dative_bond_as_single(tmp_path):
    spec = tmp_path / "pd.yaml"
    spec.write_text(f'preset: acs1996\nspecies:\n  pd0: {{smiles: "{PD0}"}}\n')
    root = etree.fromstring(spec_cdxml(spec, warn=lambda m: None).encode())
    bonds = root.iter("b")
    assert [b.get("Order") for b in bonds] == [None] * 8
    (back,) = Chem.MolsFromCDXML(etree.tostring(root).decode())
    assert back.GetNumBonds() == 8


@pytest.mark.parametrize("reactants, products, arrow", [
    ({"p": "[P:1](C)(C)C", "m": "Cl[Pd:2]Cl"}, {"x": "[P:1](C)(C)(C)->[Pd:2](Cl)Cl"},
     {"from": "lp(1)", "to": "between(1,2)"}),
    ({"x": "[P:1](C)(C)(C)->[Pd:2](Cl)Cl"}, {"p": "[P:1](C)(C)C", "m": "Cl[Pd:2]Cl"},
     {"from": "bond(1-2)", "to": "atom(1)"}),
], ids=["association", "dissociation"])
def test_arrows_make_and_break_a_dative_bond(reactants, products, arrow):
    """The arrows give P+–Pd−; a declared P->[Pd] compares equal to it (electrons.canonical)."""
    check = electrons.check_step(1, {k: chem.parse_smiles(v) for k, v in reactants.items()},
                                 {k: chem.parse_smiles(v) for k, v in products.items()},
                                 [CurvedArrow.model_validate(arrow)])
    assert check.problems == []
