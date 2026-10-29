"""Chairs and crossings vs ChemDraw (chair_x*: acs-style.md §14), within 0.1 pt.

Our drawing of each chair_x* template (templates/chair/, front bonds from ChemDraw's Z order) is
measured by the same code as ChemDraw's (scripts/measure_refs.py): the M1 items
(test_m1_vs_chemdraw.pairs: bond lengths, line widths and offsets, label positions, bond ends at
labels, the wedge outline) and, at every crossing, the rear bond's cut ends on each side of the
front bond (distance from the front axis, and the spread of the two corners: an oblique cut).
"""

import pytest
from test_m1_vs_chemdraw import _fails, measured, pairs

import measure_refs as mr

pytestmark = pytest.mark.skipif(not (mr.REF_DIR / "chair_x1.0.cdxml").exists(),
                                reason="ChemDraw references (reference_renders/) are local only")

TOL = 0.1
REFS = mr.CHAIR_REFS


@pytest.mark.parametrize("ref", REFS)
def test_exact_geometry_within_tolerance(ref):
    rows = pairs(ref)
    assert {"bond length", "line width", "label baseline", "bond end at label"} <= {p[0] for p in rows} \
        or ref != "chair_x1.0"
    fails = _fails(rows)
    assert not fails, "\n".join(fails)


@pytest.mark.parametrize("ref", REFS)
def test_crossing_cuts_within_tolerance(ref):
    r, _, _ = measured(ref)
    cd, ours = mr.crossing_summary(ref)["crossings"], mr.crossing_summary(ref, r.ours)["crossings"]
    assert len(cd) == len(ours) == 10
    fails = []
    for c, o in zip(cd, ours):
        assert (c["front"], c["rear"]) == (o["front"], o["rear"])
        assert c["rear_pieces"] <= o["rear_pieces"], c  # ours: stroke pieces plus oblique caps
        for side, e in c["cut"].items():
            f = o["cut"][side]
            for key in ("median", "min", "max"):
                d = f["beyond_axis"][key] - e["beyond_axis"][key]
                if abs(d) > TOL:
                    fails.append(f"front {c['front']} rear {c['rear']} ({c['angle_deg']} deg) side {side} {key}: "
                                 f"ChemDraw {e['beyond_axis'][key]} ours {f['beyond_axis'][key]}")
    assert not fails, "\n".join(fails)
