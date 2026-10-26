"""Scene graph -> SVG (design.md §4.4).

The root width, height and viewBox are in pt, so 1 user unit = 1 pt. There are no explicit
transforms, clip paths or background. Every stroke has an explicit width, cap and join.
Coordinates are rounded to 0.01 pt. Text is real <text> so it stays editable. Groups
keep their scene ids, which editors show as layer or object names.

Arrowheads are <marker orient="auto"> in one <defs>, the one implicit transform: the renderer
places each head on its path's end along the end tangent, so the head follows when the path is
edited (Inkscape, browsers, resvg).

Bold and italic text: the family list with font-weight / font-style, which Inkscape and Chrome
need (a face name such as "Helvetica Bold" is no family to them, so they substitute another font).

affinity=True (render.py --affinity) writes for Affinity instead, which needs the opposite on
both counts: it drops markers on import, so each head is a separate filled shape as the scene
places it, and it takes bold and italic only from a bare face name in font-family.
"""

from __future__ import annotations

import re
from xml.sax.saxutils import escape

from .geometry import Point
from .scene import Circle, Group, Marker, Path, Polygon, Polyline, TextRun

FONT_FAMILY = "Helvetica, 'TeX Gyre Heros', 'Nimbus Sans', Arial, sans-serif"
# Affinity takes bold and italic only from a bare face name in font-family (user checks:
# 2026-10-24, family list + font-weight="bold" imported regular; 2026-10-26, the face name
# first in the family list also imported regular; "Helvetica Bold" / "Helvetica Oblique",
# ChemDraw's italic name, alone imported correctly). Inkscape and Chrome can't resolve these
# (2026-10-26), so only affinity=True writes them. preview.rasterize maps them back for resvg.
FACE_FAMILIES = {"bold": "Helvetica Bold", "italic": "Helvetica Oblique"}
PADDING = 2.0  # pt of empty border around the ink; an output choice, not a style value
MARKER_MARGIN = 0.2  # pt of marker viewport around the head's control points, so its outline isn't clipped


def num(v: float) -> str:
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


def _pts(points: list[Point], dx: float, dy: float) -> str:
    return " L".join(f"{num(x + dx)} {num(y + dy)}" for x, y in points)


def safe_id(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", s)


def _d(segments, dx: float, dy: float, closed: bool) -> str:
    d = " ".join(seg[0] + " ".join(f"{num(x + dx)} {num(y + dy)}" for x, y in seg[1:]) for seg in segments)
    return d + (" Z" if closed else "")


def _markers(node, affinity: bool) -> str:
    if affinity:
        return ""
    return "".join(f' marker-{end}="url(#{safe_id(m.id)})"' for end, m in
                   (("start", node.marker_start), ("end", node.marker_end)) if m is not None)


def _stroke(width: float, closed: bool) -> str:
    return (f'stroke="#000" stroke-width="{num(width)}" stroke-linejoin="miter"'
            + ('' if closed else ' stroke-linecap="butt"'))


def _element(node, dx: float, dy: float, indent: str, affinity: bool = False) -> list[str]:
    if isinstance(node, Group):
        out = [f'{indent}<g id="{safe_id(node.id)}">']
        for c in node.children:
            out += _element(c, dx, dy, indent + "  ", affinity)
        return out + [f"{indent}</g>"]
    heads = [] if not affinity or not isinstance(node, (Polyline, Path)) else [
        line for h in node.head_leaves() for line in _element(h, dx, dy, indent, affinity)]
    if isinstance(node, Polyline):
        d = "M" + _pts(node.points, dx, dy) + (" Z" if node.closed else "")
        return [f'{indent}<path d="{d}" fill="none" stroke="#000" stroke-width="{num(node.width)}" '
                f'stroke-linecap="{node.cap}" stroke-linejoin="{node.join}"{_markers(node, affinity)}/>'] + heads
    if isinstance(node, Path):
        fill = "#000" if node.fill else "none"
        stroke = _stroke(node.stroke, node.closed) if node.stroke else 'stroke="none"'
        return [f'{indent}<path d="{_d(node.segments, dx, dy, node.closed)}" fill="{fill}" {stroke}'
                f'{_markers(node, affinity)}/>'] + heads
    if isinstance(node, Circle):
        x, y = node.center
        paint = (f'fill="none" stroke="#000" stroke-width="{num(node.stroke)}"' if node.stroke
                 else 'fill="#000" stroke="none"')
        return [f'{indent}<circle cx="{num(x + dx)}" cy="{num(y + dy)}" r="{num(node.r)}" {paint}/>']
    if isinstance(node, Polygon):
        return [f'{indent}<path d="M{_pts(node.points, dx, dy)} Z" fill="#000" stroke="none"/>']
    if isinstance(node, TextRun):
        style = ' font-style="italic"' if node.style == "italic" else ""
        family = FACE_FAMILIES["italic" if node.style == "italic" else node.weight] if affinity and (
            node.style == "italic" or node.weight == "bold") else FONT_FAMILY
        return [f'{indent}<text x="{num(node.x + dx)}" y="{num(node.y + dy)}" font-family="{family}" '
                f'font-size="{num(node.size)}" font-weight="{node.weight}"{style} fill="#000">{escape(node.text)}</text>']
    raise TypeError(f"cannot write {type(node).__name__}")


def markers_used(root: Group) -> list[Marker]:
    """The markers in the tree, one per id, in first-use order."""
    seen: dict[str, Marker] = {}
    for leaf in root.leaves():
        for m in (getattr(leaf, "marker_start", None), getattr(leaf, "marker_end", None)):
            if m is not None:
                seen.setdefault(m.id, m)
    return list(seen.values())


def _marker_def(m: Marker, indent: str) -> list[str]:
    """userSpaceOnUse, so the head keeps its pt size whatever the stroke width; the viewBox puts
    the shape at positive coordinates and refX/refY at the frame origin (the path's end)."""
    pts = [q for seg in m.segments for q in seg[1:]]
    x0, y0 = min(q[0] for q in pts) - MARKER_MARGIN, min(q[1] for q in pts) - MARKER_MARGIN
    w = max(q[0] for q in pts) + MARKER_MARGIN - x0
    h = max(q[1] for q in pts) + MARKER_MARGIN - y0
    fill = "#000" if m.fill else "none"
    stroke = _stroke(m.stroke, True) if m.stroke else 'stroke="none"'
    return [f'{indent}<marker id="{safe_id(m.id)}" markerUnits="userSpaceOnUse" orient="auto" '
            f'viewBox="0 0 {num(w)} {num(h)}" markerWidth="{num(w)}" markerHeight="{num(h)}" '
            f'refX="{num(-x0)}" refY="{num(-y0)}">',
            f'{indent}  <path d="{_d(m.segments, -x0, -y0, True)}" fill="{fill}" {stroke}/>',
            f'{indent}</marker>']


def faces_used(root: Group) -> list[str]:
    seen = {}
    for leaf in root.leaves():
        if isinstance(leaf, TextRun):
            seen.setdefault(leaf.face.describe(), leaf.face.is_fallback)
    return [f"{d}{' [FALLBACK]' if fb else ''}" for d, fb in seen.items()]


def to_svg(root: Group, padding: float = PADDING, size: tuple[float, float] | None = None,
           offset: Point | None = None, affinity: bool = False) -> tuple[str, Point]:
    """Serialize `root`. By default the drawing is moved so its ink box sits `padding`
    from the top-left, and the canvas fits it. Pass `size` and `offset` to place it on a
    fixed canvas instead (the compare script uses the ChemDraw page). `affinity`: write for
    Affinity (heads as separate shapes after their paths, bold/italic as bare face names).
    Returns the SVG text and the (dx, dy) that was added to every scene coordinate."""
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
    used = [] if affinity else markers_used(root)
    if used:
        lines.append("  <defs>")
        for m in used:
            lines += _marker_def(m, "    ")
        lines.append("  </defs>")
    lines += _element(root, dx, dy, "  ", affinity)
    lines.append("</svg>")
    return "\n".join(lines) + "\n", (dx, dy)
