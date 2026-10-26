"""M4 exact items vs ChemDraw (acs-style.md §4.1), within 0.1 pt, at all three presets.

Each resonance reference is redrawn on ChemDraw's page (compare.render_reference): structures
from templates/m4/, our double-headed arrows at ChemDraw's tips (arrows.resonance_arrow) and
our brackets on ChemDraw's bracket lines (brackets.square_bracket). So the differences come
from exact rules only: both heads, the shaft between the notches, the bracket stroke and lips.
"""

from functools import lru_cache

import pytest
from compare import render_reference

import measure_refs as mr

pytestmark = pytest.mark.skipif(not mr.REF_DIR.exists(),
                                reason="ChemDraw references (reference_renders/) are local only")

TOL = 0.1
REFS = ["resonance_x1.0", "resonance_x0.8", "resonance_x0.6"]


@lru_cache(maxsize=None)
def measured(ref):
    r = render_reference(ref)
    return r, mr.measure(r.ref), mr.measure(r.ours)


@pytest.mark.parametrize("ref", REFS)
def test_double_headed_arrows(ref):
    _, mc, mo = measured(ref)
    assert len(mc["arrows"]) == len(mo["arrows"]) == 2
    for c, o in zip(mc["arrows"], mo["arrows"]):
        for end in ("head", "tail"):
            hc, ho = c["heads"][end], o["heads"][end]
            assert ho["length"] == pytest.approx(hc["length"], abs=TOL)
            assert ho["inset_from_tip"] == pytest.approx(hc["inset_from_tip"], abs=TOL)
            for side in ("half_width_left", "half_width_right"):
                assert ho[side] == pytest.approx(hc[side], abs=TOL)
        for key in ("width", "u_start", "u_end_from_tip"):
            assert o["shaft"][key] == pytest.approx(c["shaft"][key], abs=TOL), key


@pytest.mark.parametrize("ref", REFS)
def test_brackets(ref):
    _, mc, mo = measured(ref)
    bc, bo = mc["brackets"], mo["brackets"]
    assert len(bc) == len(bo) == 2
    for c, o in zip(sorted(bc, key=lambda b: b["x"]), sorted(bo, key=lambda b: b["x"])):
        assert o["lips_point"] == c["lips_point"]
        for key in ("x", "y_top", "y_bottom", "stroke_width"):
            assert o[key] == pytest.approx(c[key], abs=TOL), key
        assert o["lip_length"]["median"] == pytest.approx(c["lip_length"]["median"], abs=TOL)
