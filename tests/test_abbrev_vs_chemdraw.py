"""Abbreviation labels vs ChemDraw (abbrev_x*: every chem.ABBREVIATIONS key on a benzene ring,
bond leaving the label to the left (OMe) and to the right (MeO), plus ±30° bonds for OMe, CO2Me
and NO2; acs-style.md §3.3). Our drawing of the templates (templates/abbrev/) is measured by the
same code as ChemDraw's (scripts/measure_refs.py): displayed text, every glyph origin relative to
the atom, and each bond end, within 0.1 pt."""

from functools import lru_cache

import pytest
from compare import render_reference

import measure_refs as mr
from acs_chem_svg.chem import ABBREVIATIONS, ABBREVIATIONS_REVERSED

pytestmark = pytest.mark.skipif(not (mr.REF_DIR / "abbrev_x1.0.cdxml").exists(),
                                reason="ChemDraw references (reference_renders/) are local only")

TOL = 0.1
REFS = ["abbrev_x1.0", "abbrev_x0.8", "abbrev_x0.6"]


@lru_cache(maxsize=None)
def measured(ref):
    r = render_reference(ref)
    return r, mr.measure(r.ref), mr.measure(r.ours)


def _glyphs(doc, atom):
    _, whole = mr.label_glyphs(doc, atom)
    return [(g.ch, g.x - atom.p[0], g.base - atom.p[1], g.size) for g in whole]


@pytest.mark.parametrize("ref", REFS)
def test_every_abbreviation_in_both_forms(ref):
    r, _, _ = measured(ref)
    shown = {(a.label, "".join(c for c, *_ in _glyphs(r.ref, a))) for a in r.ref.atoms.values() if a.nickname}
    assert shown == {(k, k) for k in ABBREVIATIONS} | {(k, v) for k, v in ABBREVIATIONS_REVERSED.items()}


@pytest.mark.parametrize("ref", REFS)
def test_label_glyphs_within_tolerance(ref):
    r, _, _ = measured(ref)
    fails = []
    for a in r.ref.atoms.values():
        if not a.nickname:
            continue
        cd, ours = _glyphs(r.ref, a), _glyphs(r.ours, r.ours.atoms[a.id])
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
        if not e.get("nickname"):
            continue
        assert "error" not in e and "error" not in ours[e["atom"]], e["atom"]
        assert e["glyph_run"] == ours[e["atom"]]["glyph_run"]
        ob = {b["bond"]: b for b in ours[e["atom"]]["bonds"]}
        for b in e["bonds"]:
            n += 1
            d = ob[b["bond"]]["end_dist_from_atom"] - b["end_dist_from_atom"]
            if abs(d) > TOL:
                fails.append(f"{e['glyph_run']} atom {e['atom']} dir {b['dir_deg']}: diff {d:+.3f}")
    assert n == 48 and not fails, "\n".join(fails)
