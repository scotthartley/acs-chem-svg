"""Charges (acs-style.md §2.2): clearance from lone pairs, circled signs, 2+ / override texts.
All placements here are unmeasured rules (user decision 2026-10-26)."""

from pathlib import Path

import pytest

from acs_chem_svg import chem
from acs_chem_svg.draw_mol import CHARGE_DOT_CLEARANCE, Drawer, parse_charge_text
from acs_chem_svg.geometry import segment_box_dist
from acs_chem_svg.preset import load_preset
from acs_chem_svg.scene import Circle, Polyline, TextRun
from acs_chem_svg.spec import Spec, molecules
from acs_chem_svg.svg_writer import to_svg

PRESETS = ["acs1996", "acs1996_x0.8", "acs1996_x0.6"]
HERE = Path(__file__).parent


def drawer(smiles, preset, lone_pairs=False, circled=False):
    mol = chem.from_smiles(smiles, "m", preset.lengths.bond_length)
    mol.circled_charges = circled
    show = {a.idx for a in mol.atoms if a.charge} if lone_pairs else None
    return Drawer(mol, preset, show_lone_pairs=show)


def charged(d):
    return next(a.idx for a in d.mol.atoms if a.charge)


@pytest.mark.parametrize("name", PRESETS)
@pytest.mark.parametrize("smiles", ["[OH-]", "[Br-]", "OS(=O)(=O)[O-]"])
def test_labeled_charge_clears_its_lone_pairs(smiles, name):
    """HO⁻, Br⁻, bisulfate O⁻ (m3_sn2, m4_nitration): the sign moves out along its corner's
    diagonal until it is CHARGE_DOT_CLEARANCE from every dot, and stays off the label."""
    P = load_preset(name)
    d = drawer(smiles, P, lone_pairs=True)
    i = charged(d)
    assert d.dot_sets()[i]
    segs = d.charge_sign(i).clash_segments(d.charge_center(i))
    assert d.charge_clearance(i, segs) >= CHARGE_DOT_CLEARANCE - 1e-6
    for run in d.labels[i].runs:
        assert min(segment_box_dist(a, b, run.bbox()) - hw for a, b, hw in segs) > 0


@pytest.mark.parametrize("name", PRESETS)
def test_charged_label_without_dots_keeps_its_corner(name):
    P = load_preset(name)
    c = P.charge
    d = drawer("C[N+](=O)[O-]", P)
    i = next(a.idx for a in d.mol.atoms if a.charge < 0)
    x0, y0, x1, y1 = d.labels[i].elem_ink_box()
    assert d.charge_center(i) == pytest.approx((x1 + c.labeled_offset_x, y0 + c.labeled_offset_y))


def spec(**kw):
    sp = {"smiles": "C[N+](=O)[O-]", **kw.pop("species", {})}
    return Spec.model_validate({"species": {"nitro": sp}, **kw})


def rings(mol, preset=load_preset("acs1996")):
    return [p for p in Drawer(mol, preset).draw("x").leaves() if isinstance(p, Circle) and p.stroke]


def test_circled_off_by_default_on_per_spec_and_overridden_per_species():
    assert not rings(molecules(spec(), HERE)[0])
    assert len(rings(molecules(spec(charges={"circled": True}), HERE)[0])) == 2
    off = spec(charges={"circled": True}, species={"circled_charges": False})
    assert not rings(molecules(off, HERE)[0])
    assert len(rings(molecules(spec(species={"circled_charges": True}), HERE)[0])) == 2


def test_circled_ring_is_written_unfilled_in_the_charges_group():
    mol = molecules(spec(charges={"circled": True}), HERE)[0]
    svg = to_svg(Drawer(mol, load_preset("acs1996")).draw("nitro"))[0]
    group = svg.split('id="species-nitro-charges"')[1].split("</g>")[0]
    assert group.count('fill="none" stroke="#000"') >= 2 and "<circle" in group


@pytest.mark.parametrize("name", PRESETS)
def test_circled_charge_on_unlabeled_carbon_clears_its_bonds(name):
    P = load_preset(name)
    d = drawer("C[CH+]C", P, circled=True)
    i = charged(d)
    ring = d.charge_sign(i)
    assert ring.circled
    assert d.charge_clearance(i, ring.clash_segments(d.charge_center(i))) >= CHARGE_DOT_CLEARANCE - 1e-6


def parts(d, i):
    return d.charge_parts(i)


@pytest.mark.parametrize("smiles, bars", [("[Fe+2]", 2), ("[O-2]", 1)])
def test_multiple_charge_draws_digits_then_sign(smiles, bars):
    d = drawer(smiles, load_preset("acs1996"), circled=True)  # 2+ is never circled
    i = charged(d)
    ps = parts(d, i)
    runs = [p for p in ps if isinstance(p, TextRun)]
    lines = [p for p in ps if isinstance(p, Polyline)]
    assert [r.text for r in runs] == ["2"] and len(lines) == bars
    assert not any(isinstance(p, Circle) for p in ps)
    assert runs[0].bbox()[2] < min(min(q[0] for q in ln.points) for ln in lines)  # digit before the sign
    assert runs[0].bbox()[0] > d.labels[i].elem_ink_box()[2]


def test_charge_label_overrides_the_text():
    s = Spec.model_validate({"species": {"fe": {"smiles": "[Fe+2:1]", "charge_labels": {1: "+3"}}},
                             "charges": {"circled": True}})
    mol = molecules(s, HERE)[0]
    assert mol.atoms[0].charge == 2 and mol.atoms[0].charge_label == "+3"
    d = Drawer(mol, load_preset("acs1996"))
    ps = parts(d, 0)
    (run,) = [p for p in ps if isinstance(p, TextRun)]
    lines = [p for p in ps if isinstance(p, Polyline)]
    assert run.text == "3" and len(lines) == 2 and not any(isinstance(p, Circle) for p in ps)
    assert run.bbox()[0] > max(max(q[0] for q in ln.points) for ln in lines)  # sign before the digit


def test_charge_label_on_a_neutral_atom_is_drawn():
    s = Spec.model_validate({"species": {"fe": {"smiles": "[Fe:1]", "charge_labels": {1: "3+"}}}})
    d = Drawer(molecules(s, HERE)[0], load_preset("acs1996"))
    assert [p.text for p in parts(d, 0) if isinstance(p, TextRun)] == ["3"]


@pytest.mark.parametrize("text, parsed", [("+", (1, "", False)), ("-", (-1, "", False)), ("2+", (1, "2", True)),
                                          ("+3", (1, "3", False)), ("2−", (-1, "2", True))])
def test_parse_charge_text(text, parsed):
    assert parse_charge_text(text) == parsed


@pytest.mark.parametrize("text", ["", "2", "+2+", "3x", "++"])
def test_bad_charge_text_is_rejected(text):
    with pytest.raises(ValueError, match="charge_labels"):
        Spec.model_validate({"species": {"fe": {"smiles": "[Fe:1]", "charge_labels": {1: text}}}})
