"""Scene graph -> SVG (design.md §4.4).

The root width, height and viewBox are in pt, so 1 user unit = 1 pt. There are no
transforms, clip paths or background. Every stroke has an explicit width, cap and join.
Coordinates are rounded to 0.01 pt. Text is real <text> so it stays editable. Groups
keep their scene ids, which Affinity shows as layer names.
"""

from __future__ import annotations

import re
from xml.sax.saxutils import escape

from .geometry import Point
from .scene import Circle, Group, Path, Polygon, Polyline, TextRun

FONT_FAMILY = "Helvetica, 'TeX Gyre Heros', 'Nimbus Sans', Arial, sans-serif"
# Affinity takes bold and italic only from a face name in font-family (user check, 2026-10-24:
# family list + font-weight="bold" imported regular; "Helvetica Bold" / "Helvetica Oblique",
# ChemDraw's italic name, imported correctly). preview.rasterize maps these back for resvg.
FACE_FAMILIES = {"bold": "Helvetica Bold", "italic": "Helvetica Oblique"}
PADDING = 2.0  # pt of empty border around the ink; an output choice, not a style value


def num(v: float) -> str:
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


def _pts(points: list[Point], dx: float, dy: float) -> str:
    return " L".join(f"{num(x + dx)} {num(y + dy)}" for x, y in points)


def safe_id(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", s)


def _element(node, dx: float, dy: float, indent: str) -> list[str]:
    if isinstance(node, Group):
        out = [f'{indent}<g id="{safe_id(node.id)}">']
        for c in node.children:
            out += _element(c, dx, dy, indent + "  ")
        return out + [f"{indent}</g>"]
    if isinstance(node, Polyline):
        d = "M" + _pts(node.points, dx, dy) + (" Z" if node.closed else "")
        return [f'{indent}<path d="{d}" fill="none" stroke="#000" stroke-width="{num(node.width)}" '
                f'stroke-linecap="{node.cap}" stroke-linejoin="{node.join}"/>']
    if isinstance(node, Path):
        d = " ".join(seg[0] + " ".join(f"{num(x + dx)} {num(y + dy)}" for x, y in seg[1:]) for seg in node.segments)
        d += " Z" if node.closed else ""
        fill = "#000" if node.fill else "none"
        stroke = (f'stroke="#000" stroke-width="{num(node.stroke)}" stroke-linejoin="miter"'
                  + ('' if node.closed else ' stroke-linecap="butt"') if node.stroke else 'stroke="none"')
        return [f'{indent}<path d="{d}" fill="{fill}" {stroke}/>']
    if isinstance(node, Circle):
        x, y = node.center
        return [f'{indent}<circle cx="{num(x + dx)}" cy="{num(y + dy)}" r="{num(node.r)}" fill="#000" stroke="none"/>']
    if isinstance(node, Polygon):
        return [f'{indent}<path d="M{_pts(node.points, dx, dy)} Z" fill="#000" stroke="none"/>']
    if isinstance(node, TextRun):
        style = ' font-style="italic"' if node.style == "italic" else ""
        family = FACE_FAMILIES["italic" if node.style == "italic" else node.weight] if (
            node.style == "italic" or node.weight == "bold") else FONT_FAMILY
        return [f'{indent}<text x="{num(node.x + dx)}" y="{num(node.y + dy)}" font-family="{family}" '
                f'font-size="{num(node.size)}" font-weight="{node.weight}"{style} fill="#000">{escape(node.text)}</text>']
    raise TypeError(f"cannot write {type(node).__name__}")


def faces_used(root: Group) -> list[str]:
    seen = {}
    for leaf in root.leaves():
        if isinstance(leaf, TextRun):
            seen.setdefault(leaf.face.describe(), leaf.face.is_fallback)
    return [f"{d}{' [FALLBACK]' if fb else ''}" for d, fb in seen.items()]


def to_svg(root: Group, padding: float = PADDING, size: tuple[float, float] | None = None,
           offset: Point | None = None) -> tuple[str, Point]:
    """Serialize `root`. By default the drawing is moved so its ink box sits `padding`
    from the top-left, and the canvas fits it. Pass `size` and `offset` to place it on a
    fixed canvas instead (the compare script uses the ChemDraw page). Returns the SVG text
    and the (dx, dy) that was added to every scene coordinate."""
    box = root.bbox() or (0.0, 0.0, 0.0, 0.0)
    if offset is None:
        offset = (padding - box[0], padding - box[1])
    dx, dy = offset
    if size is None:
        size = (box[2] - box[0] + 2 * padding, box[3] - box[1] + 2 * padding)
    w, h = (round(v, 2) for v in size)
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="{num(w)}pt" height="{num(h)}pt" '
        f'viewBox="0 0 {num(w)} {num(h)}">',
    ]
    faces = faces_used(root)
    if faces:
        lines.append(f"  <!-- acs-chem-svg: text measured with {'; '.join(faces)} -->")
    lines += _element(root, dx, dy, "  ")
    lines.append("</svg>")
    return "\n".join(lines) + "\n", (dx, dy)
