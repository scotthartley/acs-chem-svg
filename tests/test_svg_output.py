"""SVG output rules (design.md §4.4)."""

import re

import pytest
from lxml import etree

from acs_chem_svg import chem
from acs_chem_svg.draw_mol import draw_molecule
from acs_chem_svg.preset import load_preset
from acs_chem_svg.preview import rasterize
from acs_chem_svg.svg_writer import FONT_FAMILY, to_svg

NS = "{http://www.w3.org/2000/svg}"


@pytest.fixture(scope="module", params=["acs1996", "acs1996_x0.8", "acs1996_x0.6"])
def svg(request):
    p = load_preset(request.param)
    # labels, H sides, subscripts, charges, wedge, ring, triple
    mol = chem.from_smiles("C[C@@H](N)C(=O)[O-].c1ccccc1C#C[NH3+]", "test mol", p.lengths.bond_length)
    text, _ = to_svg(draw_molecule(mol, p))
    return text


def test_units_and_viewbox(svg):
    root = etree.fromstring(svg.encode())
    w, h = root.get("width"), root.get("height")
    assert w.endswith("pt") and h.endswith("pt")
    assert root.get("viewBox") == f"0 0 {w[:-2]} {h[:-2]}"


def test_no_transforms_clips_or_effects(svg):
    for attr in ("transform", "clip-path", "vector-effect"):
        assert attr not in svg
    assert "<clipPath" not in svg


def test_named_groups(svg):
    ids = [g.get("id") for g in etree.fromstring(svg.encode()).iter(NS + "g")]
    assert ids == ["species-test_mol", "species-test_mol-bonds", "species-test_mol-labels",
                   "species-test_mol-charges"]


def test_real_text(svg):
    texts = list(etree.fromstring(svg.encode()).iter(NS + "text"))
    assert {t.text for t in texts} >= {"O", "H", "2", "N"}
    assert all(t.get("font-family") == FONT_FAMILY for t in texts)
    assert "dominant-baseline" not in svg


def test_every_stroke_explicit(svg):
    for el in etree.fromstring(svg.encode()).iter(NS + "path"):
        if el.get("stroke", "none") != "none":
            assert float(el.get("stroke-width")) > 0
            assert el.get("stroke-linecap") and el.get("stroke-linejoin")
        else:
            assert el.get("fill") == "#000"


def test_coordinates_rounded(svg):
    nums = re.findall(r'(?:d|x|y)="([^"]*)"', svg)
    for n in re.findall(r"-?\d+\.\d+", " ".join(nums)):
        assert len(n.split(".")[1]) <= 2, n


def test_transparent_background(svg):
    assert "<rect" not in svg
    import io

    from PIL import Image
    img = Image.open(io.BytesIO(rasterize(svg, 2))).convert("RGBA")
    assert img.getpixel((0, 0))[3] == 0


def test_font_provenance_comment(svg):
    assert re.search(r"<!-- acs-chem-svg: text measured with .+\.(ttc|otf)", svg)


def test_bold_and_italic_by_editor():
    """Default (Inkscape, Chrome): the family list with font-weight / font-style; a face name
    is no family to them (user check, 2026-10-26). affinity=True: the bare face name, the only
    form Affinity imports as bold/italic (user checks, 2026-10-24 and 2026-10-26). Regular text
    keeps the family list either way; resvg renders both (rasterize maps face names back;
    scripts/check_rasterizer.py)."""
    from acs_chem_svg import richtext
    from acs_chem_svg.preset import load_preset
    from acs_chem_svg.scene import Group
    from acs_chem_svg.svg_writer import FONT_FAMILY, to_svg

    p = load_preset("acs1996")
    runs = (richtext.layout_line("*n*-BuLi", 0, 10, 10, p).runs
            + richtext.layout_line("11", 0, 30, 10, p, style="bold").runs)
    svg, _ = to_svg(Group("t", runs))
    assert 'font-family="Helvetica Bold"' not in svg and 'font-family="Helvetica Oblique"' not in svg
    assert f'font-family="{FONT_FAMILY}" font-size="10" font-weight="normal" font-style="italic"' in svg
    assert f'font-family="{FONT_FAMILY}" font-size="10" font-weight="bold"' in svg
    aff, _ = to_svg(Group("t", runs), affinity=True)
    assert 'font-family="Helvetica Oblique" font-size="10" font-weight="normal" font-style="italic"' in aff
    assert 'font-family="Helvetica Bold" font-size="10" font-weight="bold"' in aff
    for out in (svg, aff):
        assert f'font-family="{FONT_FAMILY}" font-size="10" font-weight="normal" fill="#000">-BuLi' in out
    assert rasterize(svg, 4) == rasterize(aff, 4)


# ---------------------------------------------------------------- arrowheads as markers

ARROW_SPECS = ["m4_nitration_x1.0", "m4_resonance_x0.8", "m3_radical_propagation_x1.0"]


def _old_reaction_head(tip, u, p):
    """The reaction head as arrows._head built it before markers (separate Path)."""
    from acs_chem_svg.geometry import add, mul, perp
    a, n = p.arrow, perp(u)

    def at(back, side):
        return add(tip, add(mul(u, -back), mul(n, side)))
    hw = a.head_half_width
    ctl = p.ratios.arrow_back_control_per_half_width * hw
    c1, c2 = at(a.head_length, hw), at(a.head_length, -hw)
    return [tip, c1, c1, at(a.head_inset, ctl), at(a.head_inset, 0.0), at(a.head_inset, -ctl), c2, c2]


def _close(ps, qs):
    return len(ps) == len(qs) and all(abs(p[0] - q[0]) < 1e-6 and abs(p[1] - q[1]) < 1e-6 for p, q in zip(ps, qs))


@pytest.mark.parametrize("preset", ["acs1996", "acs1996_x0.6"])
def test_reaction_and_resonance_markers_place_the_old_heads(preset):
    from acs_chem_svg.arrows import reaction_arrow, resonance_arrow
    from acs_chem_svg.geometry import mul, sub, unit
    p = load_preset(preset)
    tail, head = (3.0, 7.0), (40.0, -5.0)  # slanted, so a frame mix-up shows
    u = unit(sub(head, tail))
    (h,) = reaction_arrow(tail, head, p).children[0].head_leaves()
    assert h.stroke == p.arrow.head_outline and h.closed and h.fill
    assert _close(h.points, _old_reaction_head(head, u, p))
    h1, h2 = resonance_arrow(tail, head, p).children[0].head_leaves()
    assert _close(h1.points, _old_reaction_head(tail, mul(u, -1), p))
    assert _close(h2.points, _old_reaction_head(head, u, p))


@pytest.fixture(scope="module", params=ARROW_SPECS)
def arrow_case(request):
    from pathlib import Path

    from acs_chem_svg.layout import render_spec
    root = render_spec(Path(__file__).resolve().parents[1] / "examples" / f"{request.param}.yaml")
    return root, to_svg(root)[0], to_svg(root, affinity=True)[0]


def _arrow_groups(svg):
    return [g for g in etree.fromstring(svg.encode()).iter(NS + "g")
            if re.match(r"(arrow-step|resonance-arrow|curved-arrow-step)-", g.get("id")) and not g.findall(NS + "g")]


def test_markers_by_default(arrow_case):
    root, svg, _ = arrow_case
    tree = etree.fromstring(svg.encode())
    assert len(tree.findall(NS + "defs")) == 1 and svg.count("<defs") == 1
    ids = {m.get("id") for m in tree.iter(NS + "marker")}
    scale = root.id.rsplit("_", 1)[-1]
    preset = "acs1996" if scale == "x1.0" else f"acs1996_{scale}"
    assert ids and all(i.startswith(f"{preset}-head-") for i in ids)
    refs = re.findall(r'marker-(?:start|end)="url\(#([^)]*)\)"', svg)
    assert refs and set(refs) == ids
    for m in tree.iter(NS + "marker"):
        assert m.get("markerUnits") == "userSpaceOnUse" and m.get("orient") == "auto"
    groups = _arrow_groups(svg)
    assert groups
    for g in groups:
        paths = g.findall(NS + "path")
        assert all(p.get("fill") == "none" for p in paths)  # no head shapes outside <defs>
        assert all(p.get("marker-end") for p in paths)


def test_affinity_writes_separate_heads(arrow_case):
    root, _, svg = arrow_case
    assert "<defs" not in svg and "<marker" not in svg and "marker-" not in svg
    heads = sum(len(leaf.head_leaves()) for leaf in root.leaves() if hasattr(leaf, "head_leaves"))
    filled = sum(p.get("fill") == "#000" for g in _arrow_groups(svg) for p in g.findall(NS + "path"))
    assert heads and filled == heads


def test_markers_render_like_separate_heads(arrow_case):
    import io

    import numpy as np
    from PIL import Image

    _, svg, aff = arrow_case
    a, b = (np.asarray(Image.open(io.BytesIO(rasterize(s, 4))).convert("RGBA"))[:, :, 3] / 255 for s in (svg, aff))
    assert a.shape == b.shape
    assert (abs(a - b) > 0.3).mean() < 1e-4
