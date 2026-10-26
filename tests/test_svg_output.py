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


def test_bold_and_italic_name_their_face_for_affinity():
    """Affinity imports bold/italic only from a face name in font-family (user check,
    2026-10-24); regular text keeps the M1 family list. resvg still renders them (rasterize
    maps the face names back; scripts/check_rasterizer.py)."""
    from acs_chem_svg import richtext
    from acs_chem_svg.preset import load_preset
    from acs_chem_svg.scene import Group
    from acs_chem_svg.svg_writer import FONT_FAMILY, to_svg

    p = load_preset("acs1996")
    runs = (richtext.layout_line("*n*-BuLi", 0, 10, 10, p).runs
            + richtext.layout_line("11", 0, 30, 10, p, style="bold").runs)
    svg, _ = to_svg(Group("t", runs))
    assert 'font-family="Helvetica Oblique" font-size="10" font-weight="normal" font-style="italic"' in svg
    assert 'font-family="Helvetica Bold" font-size="10" font-weight="bold"' in svg
    assert f'font-family="{FONT_FAMILY}" font-size="10" font-weight="normal" fill="#000">-BuLi' in svg
