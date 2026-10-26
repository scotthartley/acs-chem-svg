"""M3 exact items vs ChemDraw (acs-style.md §4.1), within 0.1 pt, at all three presets.

Each mech_basic reference is redrawn on ChemDraw's page (compare.render_reference): structures
from templates/m3/, our curved arrows on ChemDraw's CurvePoints (mechanism.curved_arrow), and
our dots at ChemDraw's lone-pair midpoints and radical centers. So the differences come from
exact rules only: head polygons on the end tangent (full and fishhook, with the fishhook's
offset spine), stroke width and end, dot radius and spacing.
"""

from functools import lru_cache

import pytest
from compare import render_reference

import measure_refs as mr

pytestmark = pytest.mark.skipif(not mr.REF_DIR.exists(),
                                reason="ChemDraw references (reference_renders/) are local only")
from acs_chem_svg.geometry import length, sub

TOL = 0.1
REFS = ["mech_basic_x1.0", "mech_basic_x0.8", "mech_basic_x0.6"]


@lru_cache(maxsize=None)
def measured(ref):
    r = render_reference(ref)
    return r, mr.measure(r.ref), mr.measure(r.ours)


@pytest.mark.parametrize("ref", REFS)
def test_curved_arrow_heads(ref):
    _, mc, mo = measured(ref)
    assert len(mc["curved_arrows"]) == len(mo["curved_arrows"]) == 4  # 2e⁻ O → C, C–Br → Br; two fishhooks
    for c, o in zip(mc["curved_arrows"], mo["curved_arrows"]):
        hc, ho = c["head"], o["head"]
        assert ho["n_vertices"] == hc["n_vertices"]
        for pc, po in zip(hc["pts_local"], ho["pts_local"]):  # same vertex order as ChemDraw
            assert po == pytest.approx(pc, abs=TOL)
        assert ho["tip_minus_stroke_end"] == pytest.approx(hc["tip_minus_stroke_end"], abs=TOL)


@pytest.mark.parametrize("ref", REFS)
def test_curved_arrow_strokes(ref):
    _, mc, mo = measured(ref)
    for c, o in zip(mc["curved_arrows"], mo["curved_arrows"]):
        assert o["stroke_width"] == pytest.approx(c["stroke_width"], abs=0.01)
        for pc, po in zip(c["svg_path_pts"][:4], o["svg_path_pts"][:4]):
            assert po == pytest.approx(pc, abs=TOL)


@pytest.mark.parametrize("ref", REFS)
def test_dots(ref):
    r, _, _ = measured(ref)
    pc, sc = mr._dot_sets(r.ref)
    po, so = mr._dot_sets(r.ours)
    assert (len(po), len(so)) == (len(pc), len(sc)) and pc
    for a in [c for p in po for c in p] + so:
        assert a.r == pytest.approx(pc[0][0].r, abs=TOL)
    for (a, b), (c, d) in zip(sorted(po, key=lambda p: p[0].center), sorted(pc, key=lambda p: p[0].center)):
        assert length(sub(a.center, b.center)) == pytest.approx(length(sub(c.center, d.center)), abs=TOL)
