"""Backend-agnostic scene graph in pt, y down (design.md §4.2).

draw_mol (and later layout/mechanism) build these; svg_writer serializes them and
geometry checks read them. Every leaf records its `source` ("atom:3", "bond:5",
"bond:1,2,3" for a chain), so checks and the compare loader can trace it back.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator, Union

import math

from .fontmetrics import FontFace
from .geometry import Point, add, bbox, dot, length, mul, perp, sub, union_box, unit


@dataclass(frozen=True)
class Marker:
    """An arrowhead the renderer places on a path's end (svg_writer: <marker orient="auto">),
    so the head follows the path when it is edited. `segments` ("M", p), ("L", p), ("C", c1,
    c2, p), closed, are in the marker frame: x along the path's tangent at that end, y along
    perp(x), origin at the path's end point. A start marker's shape is already turned to point
    backwards, since SVG 1.1 orients it along the path's own direction. `kind` is the leaf
    type place() returns ("polygon" | "path"), as the heads were drawn before markers."""
    id: str
    segments: tuple
    kind: str = "polygon"
    stroke: float | None = None  # hairline outline width in pt
    fill: bool = True

    def place(self, p: Point, t: Point, source: str = "") -> Union["Polygon", "Path"]:
        """The head with the frame origin at p and x along the unit vector t, in scene coordinates."""
        n = perp(t)

        def at(q: Point) -> Point:
            return add(p, add(mul(t, q[0]), mul(n, q[1])))
        segs = [(seg[0], *(at(q) for q in seg[1:])) for seg in self.segments]
        if self.kind == "polygon":
            return Polygon([seg[1] for seg in segs], source=source)
        return Path(segs, fill=self.fill, stroke=self.stroke, source=source)


def _heads(start: Point, start_dir: Point, end: Point, end_dir: Point, marker_start, marker_end, source: str) -> list:
    out = []
    if marker_start is not None:
        out.append(marker_start.place(start, unit(start_dir), source))
    if marker_end is not None:
        out.append(marker_end.place(end, unit(end_dir), source))
    return out


def _with_heads(box, heads: list, ink: bool = False):
    return union_box([box] + [h.ink_bbox() if ink and isinstance(h, Path) else h.bbox() for h in heads])


@dataclass
class Polyline:
    points: list[Point]
    width: float
    join: str = "miter"
    cap: str = "butt"
    closed: bool = False
    source: str = ""
    marker_start: Marker | None = None
    marker_end: Marker | None = None

    def head_leaves(self) -> list:
        """The marker heads placed in scene coordinates (start first), as svg_writer's
        markers=False writes them."""
        pts = self.points
        return _heads(pts[0], sub(pts[1], pts[0]), pts[-1], sub(pts[-1], pts[-2]),
                      self.marker_start, self.marker_end, self.source)

    def bbox(self):
        # miter tips can poke out further than width/2; 2x covers joins at >= 60 deg
        pad = self.width if self.join == "miter" and len(self.points) > 2 else self.width / 2
        xs, ys = [p[0] for p in self.points], [p[1] for p in self.points]
        return _with_heads((min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad), self.head_leaves())

    def ink_bbox(self):
        """Exact ink box of the stroke: butt-end corners and miter tips (SVG miter limit 4).
        bbox() stays conservative because M1 layouts were approved with it."""
        return _with_heads(self._stroke_ink_bbox(), self.head_leaves(), ink=True)

    def _stroke_ink_bbox(self):
        if self.cap != "butt" or self.join != "miter":
            return self.bbox()
        h = self.width / 2
        pts = self.points + ([self.points[0]] if self.closed else [])
        out = []
        for a, b in zip(pts, pts[1:]):
            u = unit(sub(b, a))
            n = mul(perp(u), h)
            out += [add(a, n), sub(a, n), add(b, n), sub(b, n)]
        joins = list(zip(pts, pts[1:], pts[2:])) + ([(pts[-2], pts[0], pts[1])] if self.closed else [])
        for a, v, b in joins:
            u1, u2 = unit(sub(v, a)), unit(sub(b, v))
            c = max(-1.0, min(1.0, dot(u1, u2)))
            half = math.acos(c) / 2  # half the turn angle
            if math.cos(half) < 1e-9 or 1 / math.cos(half) > 4:
                continue  # straight-through, or beyond the miter limit (bevel: corners cover it)
            bis = unit(sub(u1, u2)) if abs(c) < 1 - 1e-12 else perp(u1)
            out.append(add(v, mul(bis, h / math.cos(half))))
        return bbox(out)


@dataclass
class Polygon:
    points: list[Point]
    source: str = ""

    def bbox(self):
        xs, ys = [p[0] for p in self.points], [p[1] for p in self.points]
        return (min(xs), min(ys), max(xs), max(ys))


@dataclass
class Path:
    """Filled and/or stroked path with straight and cubic segments (reaction arrowheads; M3
    curved arrows). `segments`: ("M", p), ("L", p), ("C", c1, c2, p); `closed` adds Z."""
    segments: list[tuple]
    fill: bool = True
    stroke: float | None = None  # stroke width in pt, or None
    closed: bool = True
    source: str = ""
    marker_start: Marker | None = None
    marker_end: Marker | None = None

    def head_leaves(self) -> list:
        """The marker heads placed in scene coordinates (start first). Each end's tangent falls
        back to the next control point when a handle sits on that end (ChemDraw's C–Br arrow
        ends on its own handle), as SVG orients markers."""
        if self.marker_start is None and self.marker_end is None:
            return []
        start, end = self.segments[0][1], self.segments[-1][-1]
        first = [q for q in self.segments[1][1:]] if len(self.segments) > 1 else []
        last = [self.segments[-2][-1]] + list(self.segments[-1][1:-1]) if len(self.segments) > 1 else []
        ahead = next((q for q in first if length(sub(q, start)) > 1e-6), None)
        behind = next((q for q in reversed(last) if length(sub(end, q)) > 1e-6), None)
        return _heads(start, sub(ahead, start) if ahead else (1.0, 0.0), end, sub(end, behind) if behind else (1.0, 0.0),
                      self.marker_start, self.marker_end, self.source)

    @property
    def points(self) -> list[Point]:
        return [p for seg in self.segments for p in seg[1:]]

    @points.setter
    def points(self, pts: list[Point]) -> None:
        it = iter(pts)
        self.segments = [(seg[0], *(next(it) for _ in seg[1:])) for seg in self.segments]

    def bbox(self):
        # cubics stay inside their control hull, so the control points bound the curve
        pad = (self.stroke or 0.0) / 2
        xs, ys = [p[0] for p in self.points], [p[1] for p in self.points]
        return _with_heads((min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad), self.head_leaves())

    def ink_bbox(self):
        """Exact box of the segments (a cubic's extremes, not its control hull), padded by half
        the stroke: a curved arrow's handles reach well beyond its curve."""
        pts, cur = [], None
        for seg in self.segments:
            if seg[0] == "C":
                (c1, c2, p3) = seg[1:]
                pts += [cur, p3]
                for k in range(2):
                    a = -cur[k] + 3 * c1[k] - 3 * c2[k] + p3[k]
                    b = 2 * (cur[k] - 2 * c1[k] + c2[k])
                    c = c1[k] - cur[k]
                    ts = [-c / b] if abs(a) < 1e-12 and abs(b) > 1e-12 else []
                    if abs(a) >= 1e-12 and b * b - 4 * a * c >= 0:
                        r = math.sqrt(b * b - 4 * a * c)
                        ts = [(-b + r) / (2 * a), (-b - r) / (2 * a)]
                    pts += [_cubic(cur, c1, c2, p3, t) for t in ts if 0 < t < 1]
            else:
                pts.append(seg[1])
            cur = seg[-1]
        pad = (self.stroke or 0.0) / 2
        b = bbox(pts)
        return _with_heads((b[0] - pad, b[1] - pad, b[2] + pad, b[3] + pad), self.head_leaves(), ink=True)


def _cubic(p0: Point, c1: Point, c2: Point, p3: Point, t: float) -> Point:
    s = 1 - t
    return tuple(s ** 3 * p0[k] + 3 * s * s * t * c1[k] + 3 * s * t * t * c2[k] + t ** 3 * p3[k] for k in range(2))


@dataclass
class Circle:
    """Filled dot: lone pairs and radicals (M3), as ChemDraw draws them. With `stroke` (pt), a
    ring instead: the circle around a circled charge."""
    center: Point
    r: float
    source: str = ""
    stroke: float | None = None

    @property
    def points(self) -> list[Point]:
        return [self.center]

    @points.setter
    def points(self, pts: list[Point]) -> None:
        self.center = pts[0]

    def bbox(self):
        x, y = self.center
        r = self.r + (self.stroke or 0.0) / 2
        return (x - r, y - r, x + r, y + r)


@dataclass
class TextRun:
    x: float  # glyph origin of the first character
    y: float  # baseline
    size: float
    text: str
    face: FontFace
    weight: str = "normal"
    source: str = ""
    style: str = "normal"  # "italic" for oblique runs (n-BuLi)

    def bbox(self):
        b = self.face.text_ink_bounds(self.text, self.size)
        if b is None:
            return None
        return (self.x + b[0], self.y - b[3], self.x + b[2], self.y - b[1])

    def advance(self) -> float:
        return self.face.advance(self.text, self.size)


Leaf = Union[Polyline, Polygon, Path, Circle, TextRun]


@dataclass
class Group:
    id: str
    children: list[Union["Group", Leaf]] = field(default_factory=list)

    def leaves(self) -> Iterator[Leaf]:
        for c in self.children:
            if isinstance(c, Group):
                yield from c.leaves()
            else:
                yield c

    def bbox(self):
        return union_box(leaf.bbox() for leaf in self.leaves())

    def ink_bbox(self):
        """Exact ink box (Polyline.ink_bbox for strokes, Path.ink_bbox for cubics); the scheme
        layout measures gaps with it."""
        return union_box(leaf.ink_bbox() if isinstance(leaf, (Polyline, Path)) else leaf.bbox() for leaf in self.leaves())

    def translate(self, dx: float, dy: float) -> None:
        for leaf in self.leaves():
            if isinstance(leaf, TextRun):
                leaf.x += dx
                leaf.y += dy
            else:
                leaf.points = [(x + dx, y + dy) for x, y in leaf.points]
