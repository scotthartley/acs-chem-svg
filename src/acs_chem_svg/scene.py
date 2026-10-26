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
from .geometry import Point, add, bbox, dot, mul, perp, sub, union_box, unit


@dataclass
class Polyline:
    points: list[Point]
    width: float
    join: str = "miter"
    cap: str = "butt"
    closed: bool = False
    source: str = ""

    def bbox(self):
        # miter tips can poke out further than width/2; 2x covers joins at >= 60 deg
        pad = self.width if self.join == "miter" and len(self.points) > 2 else self.width / 2
        xs, ys = [p[0] for p in self.points], [p[1] for p in self.points]
        return (min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad)

    def ink_bbox(self):
        """Exact ink box of the stroke: butt-end corners and miter tips (SVG miter limit 4).
        bbox() stays conservative because M1 layouts were approved with it."""
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
        return (min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad)

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
        return (b[0] - pad, b[1] - pad, b[2] + pad, b[3] + pad)


def _cubic(p0: Point, c1: Point, c2: Point, p3: Point, t: float) -> Point:
    s = 1 - t
    return tuple(s ** 3 * p0[k] + 3 * s * s * t * c1[k] + 3 * s * t * t * c2[k] + t ** 3 * p3[k] for k in range(2))


@dataclass
class Circle:
    """Filled dot: lone pairs and radicals (M3), as ChemDraw draws them."""
    center: Point
    r: float
    source: str = ""

    @property
    def points(self) -> list[Point]:
        return [self.center]

    @points.setter
    def points(self, pts: list[Point]) -> None:
        self.center = pts[0]

    def bbox(self):
        x, y = self.center
        return (x - self.r, y - self.r, x + self.r, y + self.r)


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
