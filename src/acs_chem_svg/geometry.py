"""2D vector helpers on (x, y) tuples, in pt. Shared by the renderer and scripts/measure_refs.py."""

from __future__ import annotations

import math

Point = tuple[float, float]


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def mul(a, k):
    return (a[0] * k, a[1] * k)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def length(a):
    return math.hypot(*a)


def unit(a):
    n = length(a)
    return (a[0] / n, a[1] / n)


def perp(a):
    """a rotated +90° in the y-down frame: (x, y) -> (-y, x)."""
    return (-a[1], a[0])


def rotate(p, deg, about=(0.0, 0.0)):
    """Rotate p by `deg` (counter-clockwise on screen, i.e. y-down frame) about `about`."""
    t = math.radians(deg)
    c, s = math.cos(t), math.sin(t)
    x, y = sub(p, about)
    return add(about, (c * x + s * y, -s * x + c * y))


def angle_deg(a, b):
    """Unsigned angle between vectors a and b in degrees."""
    c = dot(a, b) / (length(a) * length(b))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def seg_point_dist(p, a, b):
    ab = sub(b, a)
    t = max(0.0, min(1.0, dot(sub(p, a), ab) / max(dot(ab, ab), 1e-12)))
    return length(sub(p, add(a, mul(ab, t))))


def seg_seg_dist(a, b, c, d):
    return min(seg_point_dist(a, c, d), seg_point_dist(b, c, d),
               seg_point_dist(c, a, b), seg_point_dist(d, a, b))


def segments_cross(a, b, c, d) -> bool:
    """Proper or touching intersection of segments ab and cd."""
    def orient(p, q, r):
        v = cross(sub(q, p), sub(r, p))
        return 0 if abs(v) < 1e-12 else (1 if v > 0 else -1)

    def on(p, q, r):  # r on segment pq, given collinear
        return min(p[0], q[0]) - 1e-12 <= r[0] <= max(p[0], q[0]) + 1e-12 and \
            min(p[1], q[1]) - 1e-12 <= r[1] <= max(p[1], q[1]) + 1e-12
    o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
    if o1 != o2 and o3 != o4:
        return True
    return any(o == 0 and on(p, q, r) for o, p, q, r in ((o1, a, b, c), (o2, a, b, d), (o3, c, d, a), (o4, c, d, b)))


def segment_distance(a, b, c, d) -> float:
    """True distance between segments ab and cd (0 when they cross)."""
    return 0.0 if segments_cross(a, b, c, d) else seg_seg_dist(a, b, c, d)


def point_box_dist(p, box) -> float:
    """Distance from p to an axis-aligned box (0 inside)."""
    dx = max(box[0] - p[0], 0.0, p[0] - box[2])
    dy = max(box[1] - p[1], 0.0, p[1] - box[3])
    return math.hypot(dx, dy)


def segment_box_dist(a, b, box) -> float:
    """Distance from segment ab to an axis-aligned box (0 when it enters the box)."""
    if point_box_dist(a, box) == 0 or point_box_dist(b, box) == 0:
        return 0.0
    x0, y0, x1, y1 = box
    edges = [((x0, y0), (x1, y0)), ((x1, y0), (x1, y1)), ((x1, y1), (x0, y1)), ((x0, y1), (x0, y0))]
    return min(segment_distance(a, b, c, d) for c, d in edges)


def ray_box_gap(p, direction, box):
    """Distance from p along `direction` to the first hit on an axis-aligned box
    (xmin, ymin, xmax, ymax); None if the ray misses."""
    tmin, tmax = -math.inf, math.inf
    for i in (0, 1):
        lo, hi = box[i], box[i + 2]
        if abs(direction[i]) < 1e-12:
            if not lo <= p[i] <= hi:
                return None
            continue
        t1, t2 = (lo - p[i]) / direction[i], (hi - p[i]) / direction[i]
        tmin, tmax = max(tmin, min(t1, t2)), min(tmax, max(t1, t2))
    if tmax < max(tmin, 0):
        return None
    return max(tmin, 0.0)


def ray_box_exit(p, direction, box):
    """Distance from p (inside the box) along `direction` to where the ray leaves the box."""
    t = math.inf
    for i in (0, 1):
        if direction[i] > 1e-12:
            t = min(t, (box[i + 2] - p[i]) / direction[i])
        elif direction[i] < -1e-12:
            t = min(t, (box[i] - p[i]) / direction[i])
    return t


def convex_hull(points):
    pts = sorted(set(points))

    def turn(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lo, up = [], []
    for p in pts:
        while len(lo) >= 2 and turn(lo[-2], lo[-1], p) <= 0:
            lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(up) >= 2 and turn(up[-2], up[-1], p) <= 0:
            up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]


def ray_polygon(p, d, poly):
    """Distance from p along d to the first crossing of closed polygon `poly`."""
    best = None
    for a, b in zip(poly, poly[1:] + poly[:1]):
        e = sub(b, a)
        den = -d[0] * e[1] + d[1] * e[0]
        if abs(den) < 1e-12:
            continue
        r = sub(a, p)
        t = (-r[0] * e[1] + r[1] * e[0]) / den
        s = (d[0] * r[1] - d[1] * r[0]) / den
        if 0 <= s <= 1 and t >= 0 and (best is None or t < best):
            best = t
    return best


def bbox(points):
    xs, ys = [p[0] for p in points], [p[1] for p in points]
    return (min(xs), min(ys), max(xs), max(ys))


def union_box(boxes):
    boxes = [b for b in boxes if b]
    if not boxes:
        return None
    return (min(b[0] for b in boxes), min(b[1] for b in boxes),
            max(b[2] for b in boxes), max(b[3] for b in boxes))
