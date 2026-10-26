"""M2 exact items vs ChemDraw (acs-style.md §3), within 0.1 pt, at all three presets.

Each scheme reference is redrawn on ChemDraw's page (compare.render_reference): structures
from templates/m2/, our reaction arrows at ChemDraw's tail and head, and every caption line
set by richtext from its first glyph origin. So these differences come from exact rules only:
arrowhead and shaft geometry, text advances and kerning across ChemDraw's runs, sub- and
superscript size and offsets, the italic run, and the OMe/OBn bond ends.
"""

from functools import lru_cache

import pytest
from compare import label_bond_ends, render_reference

import measure_refs as mr

pytestmark = pytest.mark.skipif(not mr.REF_DIR.exists(),
                                reason="ChemDraw references (reference_renders/) are local only")
from acs_chem_svg import fontmetrics

TOL = 0.1
REFS = ["scheme_basic_x1.0", "scheme_basic_x0.8", "scheme_basic_x0.6",
        "scheme_multi_x1.0", "scheme_multi_x0.8", "scheme_multi_x0.6"]


@lru_cache(maxsize=None)
def measured(ref):
    r = render_reference(ref)
    return r, mr.measure(r.ref), mr.measure(r.ours)


@pytest.mark.parametrize("ref", REFS)
def test_arrowheads_and_shafts(ref):
    _, mc, mo = measured(ref)
    assert len(mc["arrows"]) == len(mo["arrows"]) > 0
    for c, o in zip(mc["arrows"], mo["arrows"]):
        hc, ho = c["heads"]["head"], o["heads"]["head"]
        assert ho["length"] == pytest.approx(hc["length"], abs=TOL)
        assert ho["inset_from_tip"] == pytest.approx(hc["inset_from_tip"], abs=TOL)
        for side in ("half_width_left", "half_width_right"):
            assert ho[side] == pytest.approx(hc[side], abs=TOL)
        assert o["shaft"]["width"] == pytest.approx(c["shaft"]["width"], abs=TOL)
        assert o["shaft"]["u_end_from_tip"] == pytest.approx(c["shaft"]["u_end_from_tip"], abs=TOL)
        assert o["shaft"]["u_start"] == pytest.approx(c["shaft"]["u_start"], abs=TOL)


@pytest.mark.parametrize("ref", REFS)
def test_caption_glyphs(ref):
    """Every glyph of every caption (conditions, plus, numbers) where ChemDraw put it: ChemDraw
    starts a new SVG run at each word, size and style change, so its run origins test our
    advances, kerning, sub/superscript offsets and italic advance one glyph at a time."""
    r, _, _ = measured(ref)
    ours = [g for g in r.ours.glyphs if g.box]
    n = 0
    for cap in r.ref.captions:
        for g in mr._caption_glyphs(r.ref, cap):
            o = min((x for x in ours if x.ch == g.ch and abs(x.size - g.size) < 0.01),
                    key=lambda x: abs(x.x - g.x) + abs(x.base - g.base))
            assert (o.x, o.base) == (pytest.approx(g.x, abs=TOL), pytest.approx(g.base, abs=TOL)), (cap["text"], g.ch)
            assert o.italic == g.italic and o.bold == g.bold, (cap["text"], g.ch)
            n += 1
    assert n > 20


@pytest.mark.parametrize("ref", REFS)
def test_abbreviation_and_label_bond_ends(ref):
    """OMe / OBn (scheme_multi) and the other labels: bond ends within 0.1 pt of ChemDraw."""
    _, mc, mo = measured(ref)
    ec, eo = label_bond_ends(mc), label_bond_ends(mo)
    assert ec.keys() == eo.keys()
    for k, e in ec.items():
        assert eo[k]["end_dist_from_atom"] == pytest.approx(e["end_dist_from_atom"], abs=TOL), (e["run"], k)
    if ref.startswith("scheme_multi"):
        assert {e["run"] for e in ec.values()} >= {"OMe", "OBn"}


@pytest.mark.skipif(not fontmetrics.SYSTEM_HELVETICA.exists(), reason="kerning is Helvetica's")
@pytest.mark.parametrize("ref, k", [("scheme_multi_x1.0", 1.0), ("scheme_multi_x0.8", 0.8), ("scheme_multi_x0.6", 0.6)])
def test_chemdraw_kerns_bold_11(ref, k):
    """Closes Open 6: ChemDraw's own PDF glyph origins step '1' -> '1' by the bold advance plus
    Helvetica Bold's -0.547 pt (x k) kern pair, as fontmetrics (and resvg) set it."""
    text, xs = mr.pdf_char_origins(ref)
    i = text.index("11")
    bold = fontmetrics.load_face("Bold")
    assert xs[i + 1] - xs[i] == pytest.approx(bold.advance("11", 10 * k) - bold.advance("1", 10 * k), abs=0.005)
    assert bold.kerning("1", "1", 10 * k) == pytest.approx(-0.547 * k, abs=0.002)
