"""SVG -> PNG rasterization with resvg (resvg-py), restricted to explicit font files.

resvg-py 0.5.0 rejects root width/height given in pt ("SVG has an invalid
size"), so the root size is rewritten to unitless pixels before rendering.
Our SVGs use a viewBox in pt, so this sets the scale exactly to px_per_pt.

Bold and italic runs name their face in font-family ("Helvetica Bold", "Helvetica Oblique"),
as Affinity needs and ChemDraw writes. resvg matches family names only, so rasterize() maps
those to the Helvetica family list; font-weight / font-style then pick the face.
"""

from __future__ import annotations

import re
from pathlib import Path

import resvg_py
from lxml import etree

from .fontmetrics import FALLBACK_FILES, SYSTEM_HELVETICA

_LEN = re.compile(r"^\s*([0-9.]+)\s*(pt|px)?\s*$")
_FACE_FAMILIES = ("Helvetica Bold", "Helvetica Oblique")


def default_font_files() -> list[str]:
    files = [str(p) for p in FALLBACK_FILES.values()]
    if SYSTEM_HELVETICA.exists():
        files.insert(0, str(SYSTEM_HELVETICA))
    return files


def _pt(value: str) -> float:
    m = _LEN.match(value)
    if not m:
        raise ValueError(f"Unsupported root length {value!r}")
    v = float(m.group(1))
    return v if m.group(2) in (None, "pt") or m.group(2) == "pt" else v * 0.75


def rasterize(svg: str | bytes, px_per_pt: float = 4.0, font_files: list[str] | None = None,
              background: str | None = None) -> bytes:
    """Render an SVG whose user unit is 1 pt to PNG bytes at px_per_pt.

    Only `font_files` are visible to resvg (system fonts are skipped), so the
    raster uses the same files that fontmetrics measured.
    """
    root = etree.fromstring(svg.encode() if isinstance(svg, str) else svg)
    vb = root.get("viewBox")
    if vb:
        _, _, w, h = (float(v) for v in vb.replace(",", " ").split())
    else:
        w, h = _pt(root.get("width")), _pt(root.get("height"))
    from .svg_writer import FONT_FAMILY

    for el in root.iter("{http://www.w3.org/2000/svg}text"):
        if el.get("font-family") in _FACE_FAMILIES:
            el.set("font-family", FONT_FAMILY)
    root.set("width", f"{w * px_per_pt:.4f}")
    root.set("height", f"{h * px_per_pt:.4f}")
    png = resvg_py.svg_to_bytes(
        svg_string=etree.tostring(root, encoding="unicode"),
        skip_system_fonts=True,
        font_files=default_font_files() if font_files is None else font_files,
        background=background,
        text_rendering="geometric_precision",
    )
    return bytes(png)


def rasterize_file(svg_path: Path, png_path: Path, px_per_pt: float = 4.0, **kw) -> Path:
    png_path.write_bytes(rasterize(Path(svg_path).read_bytes(), px_per_pt, **kw))
    return png_path


# --------------------------------------------------------------------------- side-by-side comparisons


def alpha(svg_text: str, w_pt: float, h_pt: float, px_per_pt: float = 4.0, transform: str | None = None):
    """Render an SVG's content onto a w_pt x h_pt canvas (pt frame) and return its alpha
    channel as a float array in 0..1. Background rects are dropped. `transform` (an SVG
    transform string) is applied to the content, e.g. to map RDKit output onto ChemDraw's frame."""
    import io

    import numpy as np
    from PIL import Image

    root = etree.fromstring(svg_text.encode())
    for el in list(root.iter("{http://www.w3.org/2000/svg}rect")):
        el.getparent().remove(el)
    inner = "".join(etree.tostring(c, encoding="unicode") for c in root)
    g = f'<g transform="{transform}">{inner}</g>' if transform else inner
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w_pt}" height="{h_pt}" '
           f'viewBox="0 0 {w_pt} {h_pt}">{g}</svg>')
    png = rasterize(svg, px_per_pt)
    return np.asarray(Image.open(io.BytesIO(png)).convert("RGBA"))[:, :, 3].astype(float) / 255


def _gray(a):
    import numpy as np
    return np.dstack([255 * (1 - a)] * 3).astype(np.uint8)


def _overlay(a_ref, a_ours):
    """Reference in blue, ours in red, coincident ink dark."""
    import numpy as np
    return np.dstack([255 * (1 - a_ref), 255 * (1 - np.maximum(a_ref, a_ours)), 255 * (1 - a_ours)]).astype(np.uint8)


def compare_png(ref_svg: str, our_svg: str, w_pt: float, h_pt: float, out: Path,
                labels: tuple[str, str] = ("ChemDraw", "acs-chem-svg"), px_per_pt: float = 4.0,
                transform: str | None = None) -> Path:
    """Stack reference | ours | overlay (reference blue, ours red, coincident black) at one pt scale."""
    from PIL import Image, ImageDraw, ImageFont

    a_ref = alpha(ref_svg, w_pt, h_pt, px_per_pt)
    a_our = alpha(our_svg, w_pt, h_pt, px_per_pt, transform)
    h, w = a_ref.shape
    pad = 28
    canvas = Image.new("RGB", (w, 3 * (h + pad)), "white")
    d = ImageDraw.Draw(canvas)
    try:
        font = ImageFont.truetype(str(SYSTEM_HELVETICA), 20)
    except OSError:
        font = ImageFont.load_default()
    panels = ((_gray(a_ref), labels[0]), (_gray(a_our), labels[1]),
              (_overlay(a_ref, a_our), f"overlay: {labels[0]} blue, {labels[1]} red, coincident black"))
    for i, (img, lab) in enumerate(panels):
        y = i * (h + pad)
        d.text((6, y + 4), lab, fill="black", font=font)
        canvas.paste(Image.fromarray(img), (0, y + pad))
    canvas.save(out)
    return out


def zoom_crops(ref_svg: str, our_svg: str, w_pt: float, h_pt: float, boxes, out: Path,
               px_per_pt: float = 10.0, transform: str | None = None, margin_pt: float = 3.0) -> Path:
    """One row per box (pt, y down): reference | ours | overlay, side by side, at px_per_pt."""
    import numpy as np
    from PIL import Image

    a_ref = alpha(ref_svg, w_pt, h_pt, px_per_pt)
    a_our = alpha(our_svg, w_pt, h_pt, px_per_pt, transform)
    px, m = px_per_pt, margin_pt
    rows = []
    for (x0, y0, x1, y1) in boxes:
        sl = np.s_[max(0, int((y0 - m) * px)):int((y1 + m) * px), max(0, int((x0 - m) * px)):int((x1 + m) * px)]
        c, r = a_ref[sl], a_our[sl]
        gap = np.full((c.shape[0], 12, 3), 255, np.uint8)
        rows.append(np.hstack([_gray(c), gap, _gray(r), gap, _overlay(c, r)]))
    w = max(r.shape[1] for r in rows)
    rows = [np.pad(r, ((0, 12), (0, w - r.shape[1]), (0, 0)), constant_values=255) for r in rows]
    Image.fromarray(np.vstack(rows)).save(out)
    return out


# --------------------------------------------------------------------------- zoomed crops (design.md §6.5)


def crops(svg: str, boxes: dict[str, tuple], out_dir: Path, px_per_pt: float = 16.0,
          margin_pt: float = 4.0, prefix: str = "") -> list[Path]:
    """One PNG per named box (pt, the SVG's own frame) of our render alone, at px_per_pt on
    white: <out_dir>/<prefix><name>.png."""
    import io

    from PIL import Image

    img = Image.open(io.BytesIO(rasterize(svg, px_per_pt, background="white"))).convert("RGB")
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, (x0, y0, x1, y1) in boxes.items():
        m = margin_pt
        box = (max(0, int((x0 - m) * px_per_pt)), max(0, int((y0 - m) * px_per_pt)),
               min(img.width, int((x1 + m) * px_per_pt) + 1), min(img.height, int((y1 + m) * px_per_pt) + 1))
        p = out_dir / f"{prefix}{name}.png"
        img.crop(box).save(p)
        paths.append(p)
    return paths
