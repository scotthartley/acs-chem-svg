"""Regression gallery: current renders vs the user-approved PNGs in tests/gallery/.

A failure writes the new render and a diff image to build/gallery/ (a CDXML case: the new
file). Show them to the user;
only after approval run `scripts/gallery.py --update` (CLAUDE.md).
"""

import pytest
from gallery import GALLERY, alpha, canonical, cases, cdxml_cases, png

from acs_chem_svg import fontmetrics

pytestmark = pytest.mark.skipif(not fontmetrics.SYSTEM_HELVETICA.exists(),
                                reason="gallery was rendered with system Helvetica")

CASES = cases()
CDXML_CASES = cdxml_cases()
MAX_DIFFERING = 1e-4  # fraction of pixels whose alpha changes by more than 0.25


@pytest.mark.parametrize("name", sorted(CASES))
def test_matches_approved_render(name):
    import numpy as np
    from PIL import Image

    approved = GALLERY / f"{name}.png"
    assert approved.exists(), f"no approved render {approved.name}"
    new = png(CASES[name])
    a, b = alpha(approved.read_bytes()), alpha(new)
    resized = a.shape != b.shape  # canvas changed: compare on the larger canvas
    shape = (max(a.shape[0], b.shape[0]), max(a.shape[1], b.shape[1]))
    a, b = (np.pad(x, [(0, shape[0] - x.shape[0]), (0, shape[1] - x.shape[1])]) for x in (a, b))
    d = np.abs(a - b) > 0.25
    if d.mean() > MAX_DIFFERING or resized:
        out = GALLERY.parents[1] / "build" / "gallery"
        out.mkdir(parents=True, exist_ok=True)
        (out / f"{name}_new.png").write_bytes(new)
        ov = np.dstack([255 * (1 - a), 255 * (1 - np.maximum(a, b)), 255 * (1 - b)]).astype(np.uint8)
        Image.fromarray(ov).save(out / f"{name}_diff.png")
        pytest.fail(f"{d.mean():.3%} of pixels differ from the approved render (canvas resized: {resized}); "
                    f"see build/gallery/{name}_diff.png (approved blue, new red)")


@pytest.mark.parametrize("name", sorted(CDXML_CASES))
def test_cdxml_matches_approved(name):
    approved = GALLERY / "cdxml" / f"{name}.cdxml"
    assert approved.exists(), f"no approved CDXML {approved.name}"
    new = CDXML_CASES[name]
    if canonical(approved.read_text()) != canonical(new):
        out = GALLERY.parents[1] / "build" / "gallery"
        out.mkdir(parents=True, exist_ok=True)
        (out / f"{name}_new.cdxml").write_text(new)
        pytest.fail(f"CDXML differs from the approved file; see build/gallery/{name}_new.cdxml")
