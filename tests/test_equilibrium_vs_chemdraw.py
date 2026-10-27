"""Equilibrium arrows and the large arrowhead vs ChemDraw (acs-style.md §3.1, §3.6), within
0.1 pt, at all three presets.

Each reference is redrawn on ChemDraw's page (compare.render_reference): our arrows on
ChemDraw's (compare.ref_arrow, the head variant from each arrow's CDXML HeadSize). Only the
geometry ChemDraw generates is compared: shafts (offset from the axis, width, tip and butt),
half heads (length, barb, tip on the inner edge, inset, back-edge control), spacing, stagger
and the unbalanced short shaft. Conditions and gaps were placed by eye: not compared.
"""

from functools import lru_cache

import pytest
from compare import render_reference

import measure_refs as mr

pytestmark = pytest.mark.skipif(not mr.REF_DIR.exists(),
                                reason="ChemDraw references (reference_renders/) are local only")

TOL = 0.1
REFS = ["equilibrium_x1.0", "equilibrium_x0.8", "equilibrium_x0.6", "arrows_large_x1.0"]


@lru_cache(maxsize=None)
def measured(ref):
    r = render_reference(ref)
    return mr.measure(r.ref), mr.measure(r.ours)


def _pairs(ref):
    mc, mo = measured(ref)
    assert len(mc["arrows"]) == len(mo["arrows"]) > 0
    return list(zip(mc["arrows"], mo["arrows"]))


@pytest.mark.parametrize("ref", REFS)
def test_equilibrium_shafts_and_half_heads(ref):
    n = 0
    for c, o in _pairs(ref):
        if "equilibrium" not in c:
            continue
        ec, eo = c["equilibrium"], o["equilibrium"]
        assert eo["kind"] == ec["kind"]
        for key in ("spacing", "center_v", "stagger_tail", "stagger_head", "short_length", "long_length",
                    "short_center_minus_long_center"):
            if key in ec:
                assert eo[key] == pytest.approx(ec[key], abs=TOL), key
        assert len(ec["shafts"]) == len(eo["shafts"]) == 2
        for sc, so in zip(ec["shafts"], eo["shafts"]):
            assert (so["direction"], so["barb"]) == (sc["direction"], sc["barb"]) == (so["direction"], "outer")
            for key in ("v", "width", "tip_u", "butt_u", "u_end_from_tip"):
                assert so[key] == pytest.approx(sc[key], abs=TOL), (sc["direction"], key)
            for key in ("length", "barb_from_center", "tip_from_center", "inset_from_tip", "back_control_from_center"):
                assert so["head"][key] == pytest.approx(sc["head"][key], abs=TOL), (sc["direction"], key)
            n += 1
    assert n == 6  # three equilibrium arrows, two shafts each


def test_large_forward_and_resonance_heads():
    """arrows_large_x1.0's forward and ↔ arrows (HeadSize 1500) against arrow.heads.large."""
    pairs = [(c, o) for c, o in _pairs("arrows_large_x1.0") if "equilibrium" not in c]
    assert sorted(len(c["heads"]) for c, _ in pairs) == [1, 2]
    for c, o in pairs:
        assert c["cdxml"]["HeadSize"] == "1500"
        for end, hc in c["heads"].items():
            ho = o["heads"][end]
            for key in ("length", "inset_from_tip", "half_width_left", "half_width_right"):
                assert ho[key] == pytest.approx(hc[key], abs=TOL), (end, key)
        for key in ("width", "u_start", "u_end_from_tip"):
            assert o["shaft"][key] == pytest.approx(c["shaft"][key], abs=TOL), key
