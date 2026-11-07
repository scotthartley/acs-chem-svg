"""Abbreviation labels vs ChemDraw (style-spec.md §3.5). abbrev_x*: the first 21 registry labels
on a benzene ring, bond leaving the label to the left (OMe) and to the right (MeO), plus ±30°
bonds for OMe, CO2Me and NO2. abbrev_2_x*: about 30 more groups and formulas the same way, ring
NTs turned in 60° steps (Ts beside, above and below N), OTBS / NHBoc up and down, OTf / NHBoc at
±30° and ±150°. ligands_x*: PMe3, PCy3, P(t-Bu)3, P(o-tol)3, OAc, Cl and Br round a Pd. Our
drawing of the templates (templates/abbrev/) is measured by the same code as ChemDraw's
(scripts/measure_refs.py): displayed text, every glyph origin relative to the atom, and each
bond end, within 0.1 pt, except the cases KNOWN lists."""

from functools import lru_cache

import pytest
from compare import render_reference

import measure_refs as mr
from chem_svg_cdxml.chem import ABBREVIATIONS, ABBREVIATIONS_REVERSED

pytestmark = pytest.mark.skipif(not (mr.REF_DIR / "abbrev_x1.0.cdxml").exists(),
                                reason="ChemDraw references (reference_renders/) are local only")

TOL = 0.1
REFS = ["abbrev_x1.0", "abbrev_x0.8", "abbrev_x0.6", "abbrev_2_x1.0", "abbrev_2_x0.8", "abbrev_2_x0.6",
        "ligands_x1.0", "ligands_x0.8", "ligands_x0.6"]
# abbrev_x*'s labels (the registry before abbrev_2_x*)
ABBREV_X = ["OMe", "OEt", "OBn", "OPh", "OAc", "SMe", "NHAc", "CO2Me", "CO2Et", "CO2H", "CF3", "CN", "NO2",
            "Ph", "Bn", "Et", "Ac", "TMS", "TES", "TBS", "TIPS"]
# deliberate differences (user decision 2026-11-03): one reversal rule, not ChemDraw's reading;
# right-justified only (the bond end differs with the anchor: ChemDraw's ₃, our P)
DECIDED_TEXT = {"P(o-tol)3": "(o-tol)3P", "IPr": "IPr"}
# bond ends ChemDraw puts elsewhere (style-spec.md §3.5):
# - a bond 3° or 9° off horizontal into a ring NTs: Open 28's cut, not emulated (as nh_x*)
NEAR_AXIS_KNOWN = 0.65
# - (reference, label as drawn[, bond direction]) -> allowed difference:
KNOWN = {
    # italic p: ChemDraw 0.12 nearer than the oblique outline's reach (x1.0 only has italics)
    ("abbrev_2_x1.0", "p-Tol"): 0.13,
    # abbrev_2_x0.8 / x0.6 set t-Bu, p-Tol, ... upright (CDXML face 96, no italic): ours are italic
    ("abbrev_2_x0.8", "t-Bu"): 0.75, ("abbrev_2_x0.6", "t-Bu"): 0.55,
    ("abbrev_2_x0.8", "p-Tol"): 0.35, ("abbrev_2_x0.6", "p-Tol"): 0.25,
    # ligands_x0.8 / x0.6 set Lₙ's n upright (CDXML face 32): ours is italic, its right side reaches farther
    ("ligands_x0.8", "Ln"): 0.5, ("ligands_x0.6", "Ln"): 0.4,
    # dppf above Pd: ChemDraw stops the bond below the p descenders, 1.9 pt farther (x0.8 / x0.6: no
    # bond; Open 31); ours clips against d alone
    ("ligands_x1.0", "dppf", 90.0): 1.95,
}


@lru_cache(maxsize=None)
def measured(ref):
    r = render_reference(ref)
    return r, mr.measure(r.ref), mr.measure(r.ours)


def _glyphs(doc, atom):
    _, whole = mr.label_glyphs(doc, atom)
    return [(g.ch, g.x - atom.p[0], g.base - atom.p[1], g.size) for g in whole]


@pytest.mark.parametrize("ref", REFS[:3])
def test_every_abbreviation_in_both_forms(ref):
    r, _, _ = measured(ref)
    shown = {(a.label, "".join(c for c, *_ in _glyphs(r.ref, a))) for a in r.ref.atoms.values() if a.nickname}
    assert shown == {(k, k) for k in ABBREV_X} | {(k, ABBREVIATIONS_REVERSED[k]) for k in ABBREV_X}
    assert all(k in ABBREVIATIONS for k in ABBREV_X)


@pytest.mark.parametrize("ref", REFS)
def test_label_glyphs_within_tolerance(ref):
    r, _, _ = measured(ref)
    fails = []
    for a in r.ref.atoms.values():
        if not a.nickname:
            continue
        if a.id not in r.ours.atoms:  # not drawn from a template (L, dppf: text labels, test_text_labels)
            continue
        cd, ours = _glyphs(r.ref, a), _glyphs(r.ours, r.ours.atoms[a.id])
        if a.label in DECIDED_TEXT and a.justification == "Right":
            assert "".join(g[0] for g in ours) == DECIDED_TEXT[a.label]
            continue
        assert [g[0] for g in ours] == [g[0] for g in cd], (a.id, a.label)
        for (ch, x, y, s), (_, xo, yo, so) in zip(cd, ours):
            if max(abs(xo - x), abs(yo - y), abs(so - s)) > TOL:
                fails.append(f"{a.label} atom {a.id} glyph {ch}: ChemDraw ({x:.3f}, {y:.3f}) ours ({xo:.3f}, {yo:.3f})")
    assert not fails, "\n".join(fails)


@pytest.mark.parametrize("ref", REFS)
def test_bond_ends_within_tolerance(ref):
    _, mc, mo = measured(ref)
    ours = {e["atom"]: e for e in mo["labels"]}
    fails, n = [], 0
    for e in mc["labels"]:
        if not e.get("nickname") or e["atom"] not in ours:
            continue
        assert "error" not in e and "error" not in ours[e["atom"]], e["atom"]
        if e["label"] in DECIDED_TEXT and ours[e["atom"]]["glyph_run"] == DECIDED_TEXT[e["label"]]:
            continue
        assert e["glyph_run"] == ours[e["atom"]]["glyph_run"]
        ob = {b["bond"]: b for b in ours[e["atom"]]["bonds"]}
        for b in e["bonds"]:
            if b["bond"] not in ob:
                continue
            n += 1
            d = ob[b["bond"]]["end_dist_from_atom"] - b["end_dist_from_atom"]
            off_axis = min(abs(b["dir_deg"]) % 90, 90 - abs(b["dir_deg"]) % 90)
            tol = NEAR_AXIS_KNOWN if 2 < off_axis < 12 else KNOWN.get(
                (ref, e["glyph_run"], b["dir_deg"]), KNOWN.get((ref, e["glyph_run"]), TOL))
            if abs(d) > tol:
                fails.append(f"{e['glyph_run']} atom {e['atom']} dir {b['dir_deg']}: diff {d:+.3f}")
    assert n >= N_ENDS[ref.rsplit("_x", 1)[0]] and not fails, "\n".join(fails)


# bond ends compared per reference (both drawings measured: ligands_x0.8 / x0.6 have fewer)
N_ENDS = {"abbrev": 48, "abbrev_2": 89, "ligands": 17}


@pytest.mark.parametrize("ref", REFS[3:6])
def test_ring_nts_tail_line(ref):
    """Ring NTs (abbrev_2_x*): ChemDraw's LabelAlignment and, with Ts above or below N, the Ts
    line's origin x and baseline pitch from N's, as ours draws them."""
    from chem_svg_cdxml.draw_mol import Drawer
    r, _, _ = measured(ref)
    summary = _forms_by_atom(ref)
    side = {"Left": "right", "Right": "left", "Above": "top", "Below": "bottom"}
    seen = set()
    for mol in r.molecules:
        d = Drawer(mol, r.preset)
        for a in mol.atoms:
            if a.abbrev != "NTs":
                continue
            cd = r.ref.atoms[a.id]
            lab = d.labels[a.idx]
            assert lab.h_side == side[cd.alignment or cd.justification], (a.id, lab.h_side, cd.alignment)
            seen.add(lab.h_side)
            if lab.h_side in ("top", "bottom"):
                tail = summary[a.id]["tail_line"]
                n_run, t_run = lab.runs[0], lab.runs[1]
                assert t_run.text == tail["text"] == "Ts"
                assert abs((t_run.x - n_run.x) - tail["origin_dx"]) < TOL
                assert abs(abs(t_run.y - n_run.y) - tail["pitch"]) < 0.01
    assert seen == {"left", "right", "top", "bottom"}


def _forms_by_atom(ref):
    return {f["atom"]: f for f in mr.abbrev_summary(ref)["forms"]}
