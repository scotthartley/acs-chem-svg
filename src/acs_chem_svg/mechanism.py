"""Curved arrows (acs-style.md §4): two-electron arrows and fishhooks on a step's reactants.

Exact (§4.1, mech_basic_x*, the same at every preset):
- The stroke is a cubic Bézier, stroke_width wide, with no fill. It ends at the head's notch.
- The head sits on the stroke's end tangent. A full head is a straight 4-vertex polygon: the
  tip, the corners head_length back at ±head_half_width, and the notch head_inset back. A
  fishhook head has 3 vertices (notch, barb, tip) with the barb fishhook_barb out, on the
  outside of the curve's end; its axis runs fishhook_offset to the side of the stroke end,
  away from the barb.

Hand-placed (§4.2, round placement defaults):
- Tail:
  - lp(n): tail_gap beyond the outer edge of the chosen pair, outward. The pair is the one
    facing the bow side and the target, or `tail_side`.
  - electron(n): the same, from the radical dot.
  - bond(a-b): bond_tail_gap off the bond axis at its midpoint, on the bow side. Two arrows
    from one bond take opposite sides.
- Head: the tip points at the target, head_gap from it.
  - atom(n): the target is n's element ink box, or n's center when unlabeled. The tip keeps
    head_gap from the box and from n's dots, at the nearest point (clear_along), so it stops
    above a lone pair rather than between its dots; the planner also tries tips FARTHER out.
  - between(a,b), already bonded: head_gap off the bond axis at its midpoint, on the tail's side.
  - between(a,b), not bonded: for a two-electron arrow, the atom that is not part of the
    source (lp(1) -> between(1,2) points at 2). For a fishhook, and with neither or both, the
    midpoint of a and b, where the bond forms; when a plus sign sits there, head_gap past its
    ink box on the bow side (fishhooks meeting there share that side).
- Bow: away from the centroid of the species the arrow touches. When the chord passes near
  the centroid, it bows up. The planner also tries the other bow and the other lone pairs, and
  keeps whichever clashes least with bonds, labels, dots, charges and the step's earlier arrows
  (ties keep the default). `bow: left|right` (relative to the tail -> head direction) and
  `tail_side` override.
- Shape, in the chord frame (tail -> stroke end, along u and offset n on the bow side; ratios
  of the chord, from mech_basic_x*). A chord shorter than one bond is a loop
  (loop_curvature); longer ones are arcs (curvature). `curvature` overrides h.
  - arc: c1 = (0.26, h), c2 = (0.92, 0.48 h): leaves gently and lands steeply, like
    ChemDraw's 2e⁻ arrow;
  - loop: c1 = (1.19, h), c2 = (1.45, 0.59 h): ChemDraw's fishhooks, which curl back onto
    the atom.
  When the default clashes, the planner also tries h × 0.75 and × 0.5 (SOFTER), and tips
  FARTHER out. A head that doesn't aim (aims: into a labeled target's ink box past its dots,
  else within MISS_DEG of the center) counts as a clash. The tip, stroke end and controls
  are solved together (solve_end), since the end tangent decides where the tip stops.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

from .geometry import (Point, add, dot, length, mul, perp, point_box_dist, ray_box_exit, seg_point_dist,
                       segment_box_dist, segment_distance, sub, unit)
from .preset import Preset
from .draw_mol import DotSet
from .scene import Circle, Group, Marker, Path, Polygon, Polyline, TextRun
from .spec import ArrowOverride, CurvedArrow

ARC_HANDLES = ((0.26, 1.0), (0.92, 0.48))  # mech_basic_x*: (0.261, 0.531) and (0.921, 0.256) / 0.531
LOOP_HANDLES = ((1.19, 1.0), (1.45, 0.59))  # mech_basic_x*: (1.194, 1.311) and (1.446, 0.769) / 1.311
ITERATIONS = 12
SCAN = 360  # approach directions tried for an atom target (solve_end)
ROOTS = 3  # solve_end roots the planner tries per bow and tail set
FARTHER = (0.0, 1.5, 3.0)  # tip distances beyond head_gap the planner tries, so a head can clear the target's dots
MISS_DEG = 10.0  # a head pointing farther than this off an unlabeled target's center counts as a clash
SOFTER = (1.0, 0.75, 0.5)  # curvature multipliers the planner tries (flatter loops and arcs), in rank order
SAMPLES = 48  # Bézier segments for geometric checks


def bezier(p0: Point, c1: Point, c2: Point, p3: Point, t: float) -> Point:
    s = 1 - t
    return (s ** 3 * p0[0] + 3 * s * s * t * c1[0] + 3 * s * t * t * c2[0] + t ** 3 * p3[0],
            s ** 3 * p0[1] + 3 * s * s * t * c1[1] + 3 * s * t * t * c2[1] + t ** 3 * p3[1])


def head_marker(c1: Point, c2: Point, p3: Point, kind: str, preset: Preset) -> Marker:
    """The curved arrow's head in the marker frame (origin at stroke end p3, x along the end
    tangent); vertices as ChemDraw orders them. A fishhook's barb goes on the side the curve
    bends toward."""
    ca = preset.curved_arrow
    i, L = ca.head_inset, ca.head_length
    if kind == "two_electron":
        hw = ca.head_half_width
        return Marker(f"{preset.name}-head-2e", (("M", (i - L, hw)), ("L", (0.0, 0.0)), ("L", (i - L, -hw)), ("L", (i, 0.0))))
    t = unit(sub(p3, c2 if length(sub(p3, c2)) > 1e-6 else c1))
    bend = add(sub(p3, mul(c2, 2)), c1)  # B''(1) / 6: the curve bends toward this side
    s = -1.0 if dot(perp(t), bend) > 0 else 1.0  # +1: barb on perp(t), right of travel (y down)
    o, b = ca.fishhook_offset, ca.fishhook_barb  # the stroke ends beside the axis, toward the barb
    return Marker(f"{preset.name}-head-fishhook-{'right' if s > 0 else 'left'}",
                  (("M", (0.0, -s * o)), ("L", (i - L, -s * o + s * b)), ("L", (i, -s * o))))


def head_polygon(c1: Point, c2: Point, p3: Point, kind: str, preset: Preset) -> list[Point]:
    """Head at stroke end p3 on the end tangent (c2 -> p3); vertices as ChemDraw orders them."""
    path = Path([("M", c1), ("C", c1, c2, p3)], fill=False, closed=False,
                marker_end=head_marker(c1, c2, p3, kind, preset))
    return path.head_leaves()[0].points


def curved_arrow(p0: Point, c1: Point, c2: Point, p3: Point, kind: str, preset: Preset,
                 name: str = "curved-arrow") -> Group:
    """Stroke p0 -> p3 (a cubic with controls c1, c2) with its head at p3 as an end marker."""
    stroke = Path([("M", p0), ("C", c1, c2, p3)], fill=False, stroke=preset.curved_arrow.stroke_width,
                  closed=False, source=name, marker_end=head_marker(c1, c2, p3, kind, preset))
    return Group(name, [stroke])


# --------------------------------------------------------------------------- planning


@dataclass
class AtomRef:
    """An atom of a placed species, in scheme coordinates."""
    species: str
    idx: int
    p: Point
    box: tuple | None  # element ink box when labeled
    dots: list  # draw_mol.DotSet, moved into scheme coordinates


@dataclass
class MechArrow:
    step: int
    index: int
    spec: CurvedArrow
    p0: Point
    c1: Point
    c2: Point
    p3: Point
    group: Group
    head: list[Point] = field(default_factory=list)
    tail_source: str = ""  # "lp" | "electron" | "bond"
    target: tuple = ()  # ("box", box) | ("point", p) | ("bond", a, b)
    target_atom: tuple | None = None  # (species, atom index) the head lands on, for atom targets
    target_dots: list = field(default_factory=list)  # that atom's dot centers
    bow: str = ""  # the side chosen, "left" | "right" of tail -> head (as overrides.arrows write it)

    @property
    def tip(self) -> Point:
        return self.head[-1]

    def polyline(self, n: int = SAMPLES) -> list[Point]:
        return [bezier(self.p0, self.c1, self.c2, self.p3, k / n) for k in range(n + 1)]

    def shift(self, dx: float, dy: float) -> None:
        """Move the arrow's records (not its group: the layout moves groups with their row)."""
        def mv(p):
            return (p[0] + dx, p[1] + dy)
        self.p0, self.c1, self.c2, self.p3 = mv(self.p0), mv(self.c1), mv(self.c2), mv(self.p3)
        self.head = [mv(p) for p in self.head]
        self.target_dots = [mv(p) for p in self.target_dots]
        if self.target and self.target[0] == "box":
            b = self.target[1]
            self.target = ("box", (b[0] + dx, b[1] + dy, b[2] + dx, b[3] + dy))
        elif self.target:
            self.target = (self.target[0], *(mv(p) for p in self.target[1:]))


@dataclass
class Parts:
    """What a curved arrow must avoid in one species, in scheme coordinates."""
    bonds: list = field(default_factory=list)  # (a, b, half width, source)
    labels: list = field(default_factory=list)  # (ink box, text, source), one per text run
    charges: list = field(default_factory=list)  # (a, b, half width, source)
    dots: list = field(default_factory=list)  # (Circle, "lone pair" | "radical")


def _segments(points: list, closed: bool = False) -> list[tuple]:
    pts = points + ([points[0]] if closed and len(points) > 2 else [])
    return list(zip(pts, pts[1:]))


def species_parts(group: Group) -> Parts:
    """Bonds, labels, charges and dots of a drawn species (draw_mol.Drawer.draw groups)."""
    out = Parts()
    for g in group.children:
        if not isinstance(g, Group):
            continue
        kind = g.id.rsplit("-", 1)[-1]  # bonds | labels | charges | pairs (-lone-pairs) | radicals
        for leaf in g.leaves():
            if isinstance(leaf, Polyline):
                segs = [(a, b, leaf.width / 2, leaf.source) for a, b in _segments(leaf.points, leaf.closed)]
                (out.charges if kind == "charges" else out.bonds).extend(segs)
            elif isinstance(leaf, Polygon):
                out.bonds += [(a, b, 0.0, leaf.source) for a, b in _segments(leaf.points, True)]
            elif isinstance(leaf, TextRun) and leaf.bbox():
                out.labels.append((leaf.bbox(), leaf.text, leaf.source))
            elif isinstance(leaf, Circle):
                out.dots.append((leaf, "radical" if kind == "radicals" else "lone pair"))
    return out


def _seg_box(a, b) -> tuple:
    return (min(a[0], b[0]), min(a[1], b[1]), max(a[0], b[0]), max(a[1], b[1]))


def _box_gap(p: tuple, q: tuple) -> float:
    """A lower bound of every distance between the contents of boxes p and q."""
    return max(p[0] - q[2], q[0] - p[2], p[1] - q[3], q[1] - p[3], 0.0)


def clashes(arrow: "MechArrow", parts: dict[str, Parts], others: list["MechArrow"], preset: Preset) -> list[str]:
    """What the arrow's stroke and head run into (geometry_checks.check_mechanism). Pairs whose
    boxes are farther apart than the threshold are skipped: their distance can't be below it."""
    sw, margin = preset.curved_arrow.stroke_width, preset.lengths.margin_width
    stroke, head = _segments(arrow.polyline()), _segments(arrow.head, True)
    segs = [(a, b, _seg_box(a, b)) for a, b in stroke + head]
    whole = (min(s[2][0] for s in segs), min(s[2][1] for s in segs), max(s[2][2] for s in segs), max(s[2][3] for s in segs))

    def near(c, d, limit) -> bool:  # does any arrow segment come within `limit` of segment cd?
        box = _seg_box(c, d)
        if _box_gap(whole, box) >= limit:
            return False
        return any(_box_gap(sb, box) < limit and segment_distance(a, b, c, d) < limit for a, b, sb in segs)

    out = []
    for sp, pt in parts.items():
        hit = next((src for c, d, hw, src in pt.bonds if near(c, d, hw + sw / 2)), None)
        if hit:
            out.append(f"crosses a bond of {sp} ({hit})")
        for box, text, _ in pt.labels:
            if _box_gap(whole, box) >= max(sw / 2, margin):
                continue
            if any(segment_box_dist(a, b, box) < sw / 2 for a, b, _ in segs):
                out.append(f"crosses the {sp} label {text!r}")
                continue
            for end, d in (("head", min(segment_box_dist(a, b, box) for a, b in head)),
                           ("tail", point_box_dist(arrow.p0, box))):
                if d < margin:
                    out.append(f"{end} is {d:.2f} pt from the {sp} label {text!r} (margin_width {margin})")
        if any(near(c, d, hw + sw / 2) for c, d, hw, _ in pt.charges):
            out.append(f"crosses a charge in {sp}")
        for c, kind in pt.dots:
            lim = c.r + sw / 2
            x, y = c.center
            if _box_gap(whole, (x, y, x, y)) >= lim:
                continue
            if any(seg_point_dist(c.center, a, b) < lim for a, b, _ in segs):
                out.append(f"touches the {kind} on atom {c.source.split(':', 1)[-1]} of {sp}")
    for o in others:
        if any(near(c, d, sw) for c, d in _segments(o.polyline()) + _segments(o.head, True)):
            out.append(f"crosses {o.group.id}")
    return list(dict.fromkeys(out))


class Planner:
    """Lays out one step's curved arrows. `atoms` maps atom-map numbers of the step's
    reactants to AtomRef; `bonds` holds the bonded map-number pairs; `centroids` gives each
    species' atom centroid; `parts` what the arrows should avoid. Each arrow tries both bow
    sides and each candidate lone pair and keeps the one with the fewest clashes, preferring
    the default (away from the centroid, the pair facing the bow and the target). `pluses`
    holds the plus signs' ink boxes (fishhooks that meet between two species)."""

    def __init__(self, preset: Preset, atoms: dict[int, AtomRef], bonds: set[frozenset],
                 centroids: dict[str, Point], parts: dict[str, Parts] | None = None,
                 pluses: list[tuple] | None = None, strict: bool = True):
        self.preset, self.atoms, self.bonds, self.centroids = preset, atoms, bonds, centroids
        self.parts, self.pluses, self.strict = parts or {}, pluses or [], strict
        self.bond_sides: dict[frozenset, float] = {}
        self.meet_sides: dict[frozenset, float] = {}
        self.used_sets: set[int] = set()
        self.done: list[MechArrow] = []

    def ref(self, m: int, text: str) -> AtomRef:
        if m not in self.atoms:
            raise ValueError(f"arrow endpoint {text}: no drawn atom with map number {m} in the step's reactants "
                             "(atoms inside an abbreviation cannot carry arrows)")
        return self.atoms[m]

    def anchors(self, arr: CurvedArrow) -> tuple[Point, Point]:
        """Rough source and target points, to choose the bow and the lone pair."""
        src = arr.source
        if src.kind == "bond":
            a, b = (self.ref(m, src.text) for m in src.atoms)
            s0 = mul(add(a.p, b.p), 0.5)
        else:
            r = self.ref(src.atoms[0], src.text)
            rad = [d for d in r.dots if d.kind == "radical"]
            s0 = rad[0].center if src.kind == "electron" and rad else r.p
        kind, data = self.target(arr)
        if kind == "bond":
            return s0, mul(add(data[0].p, data[1].p), 0.5)
        return s0, data.p

    def target(self, arr: CurvedArrow):
        dst, src = arr.target, arr.source
        if dst.kind == "atom":
            return "atom", self.ref(dst.atoms[0], dst.text)
        a, b = (self.ref(m, dst.text) for m in dst.atoms)
        if frozenset(dst.atoms) in self.bonds:
            return "bond", (a, b)
        own = set(src.atoms)
        outside = [r for m, r in zip(dst.atoms, (a, b)) if m not in own]
        if len(outside) == 1 and arr.type == "two_electron":
            return "atom", outside[0]
        mid = mul(add(a.p, b.p), 0.5)
        gap = self.preset.curved_arrow.head_gap
        box = next((bx for bx in self.pluses if point_box_dist(mid, bx) < gap), None)
        return "atom", AtomRef(a.species, -1, mid, box, [])

    def default_bow(self, arr: CurvedArrow, s0: Point, p: Point) -> Point:
        u = unit(sub(p, s0))
        right = perp(u)  # y-down: perp() is the right-hand side of the travel direction
        species = {self.atoms[m].species for e in (arr.source, arr.target) for m in e.atoms if m in self.atoms}
        cs = [self.centroids[s] for s in species]
        c = (sum(x for x, _ in cs) / len(cs), sum(y for _, y in cs) / len(cs))
        d = dot(right, sub(mul(add(s0, p), 0.5), c))
        if abs(d) < 0.05 * length(sub(p, s0)):
            return right if right[1] < 0 or (abs(right[1]) < 1e-9 and right[0] > 0) else mul(right, -1)
        return right if d > 0 else mul(right, -1)

    def plan(self, step: int, index: int, arr: CurvedArrow, o: ArrowOverride, name: str,
             flip: bool = False) -> tuple[MechArrow, tuple] | None:
        """The arrow and its score (clashes, flip, rank). `flip` takes the other default bow; None
        when the bow left for it is ruled out (a bond's second arrow, a meeting side)."""
        src = arr.source
        s0, p_anchor = self.anchors(arr)
        u0 = unit(sub(p_anchor, s0))
        if o.bow:
            right = perp(u0)
            bows = [right if o.bow == "right" else mul(right, -1)]
        else:
            b0 = self.default_bow(arr, s0, p_anchor)
            bows = [mul(b0, -1)] if flip else [b0]
        key = frozenset(src.atoms) if src.kind == "bond" else None
        if key is not None and key in self.bond_sides and not o.bow and self.strict:
            nb = self.bond_normal(src)  # a second arrow from this bond: the other side
            bows = [b for b in bows if (1.0 if dot(nb, b) >= 0 else -1.0) != self.bond_sides[key]]
        sets: list = [None]
        if src.kind in ("lp", "electron"):
            ref = self.ref(src.atoms[0], src.text)
            kind = "pair" if src.kind == "lp" else "radical"
            sets = [d for d in ref.dots if d.kind == kind and id(d) not in self.used_sets] or \
                   [d for d in ref.dots if d.kind == kind]
            if not sets:  # the atom has none (electrons.py reports it): start beside the atom
                sets = [DotSet(ref.idx, kind, add(ref.p, mul(u0, self.preset.radical.distance)), u0,
                               "", [])]
            if o.tail_side:
                sets = [d for d in sets if d.side == o.tail_side]
                if not sets:
                    hint = ("" if ref.box else f"; on an unlabeled atom each {kind} sits on the largest bond gap's "
                            "bisector, so rotate the species (orientation.rotate) to move it")
                    raise ValueError(f"step {step} arrows[{index}]: tail_side {o.tail_side!r}: atom {src.atoms[0]} "
                                     f"has no {kind} there (has: {', '.join(d.side for d in ref.dots if d.kind == kind)})"
                                     f"{hint}")
        meet = self.meeting(arr)
        if meet is not None and meet[0] in self.meet_sides and not o.bow and self.strict:  # fishhooks meet on one side
            bows = [b for b in bows if (1.0 if dot(meet[1], b) >= 0 else -1.0) == self.meet_sides[meet[0]]]
        if not bows:
            return None
        best, r = None, self.preset.lone_pair.dot_radius
        sides = (0, 1) if self.target(arr)[0] == "bond" else (0,)  # a bond target: its tail's side, else the far one
        for rank_b, bow in enumerate(bows):
            ranked = sorted(sets, key=lambda d: 0.0 if d is None else
                            -(dot(d.direction, bow) + 0.5 * dot(d.direction, u0)))
            for rank_s, ds in enumerate(ranked):
                for rank_f, soft in enumerate(SOFTER if o.curvature is None else (1.0,)):
                    for rank_x, extra in enumerate(FARTHER):
                        for pick, far in itertools.product(range(ROOTS), sides):
                            cand = self.build(step, index, arr, o, name, bow, ds, s0, p_anchor, pick, soft, extra, far)
                            if cand is None:
                                if far == sides[-1]:
                                    break
                                continue
                            score = (len(clashes(cand, self.parts, self.done, self.preset)) + (not aims(cand, r)),
                                     flip, rank_b + rank_s + rank_f + rank_x + pick + far)
                            if best is None or score < best[0]:
                                best = (score, cand, bow, ds)
                        if best[0][0] == 0 and rank_x == 0:
                            break  # the head clears everything at head_gap
        score, arrow, bow, ds = best
        arrow.bow = "right" if dot(bow, perp(u0)) >= 0 else "left"
        if key is not None:
            self.bond_sides.setdefault(key, 1.0 if dot(self.bond_normal(src), bow) >= 0 else -1.0)
        if meet is not None:
            self.meet_sides.setdefault(meet[0], 1.0 if dot(meet[1], bow) >= 0 else -1.0)
        if ds is not None:
            self.used_sets.add(id(ds))
        self.done.append(arrow)
        return arrow, score

    def meeting(self, arr: CurvedArrow) -> tuple[frozenset, Point] | None:
        """For an arrow aimed at the midpoint of two unbonded atoms: its key and the unit normal
        of the a -> b line (a, b in map order), so arrows meeting there can share a side."""
        kind, ref = self.target(arr)
        if kind != "atom" or ref.idx != -1:
            return None
        a, b = (self.atoms[m] for m in sorted(arr.target.atoms))
        return frozenset(arr.target.atoms), perp(unit(sub(b.p, a.p)))

    def bond_normal(self, src) -> Point:
        """Unit normal of a source bond, oriented the same way for every arrow from it."""
        lo, hi = sorted(src.atoms)
        return perp(unit(sub(self.atoms[hi].p, self.atoms[lo].p)))

    def build(self, step, index, arr, o, name, bow, ds, s0, p_anchor, pick=0, soft=1.0,
              extra=0.0, far=0) -> MechArrow | None:
        """The arrow for one bow, tail set and solve_end root (`pick`), its curvature times
        `soft` and its tip `extra` beyond head_gap; None past the last root. A bond target's
        tip goes on the tail's side of the bond, or with `far` the other side (a ring bond next
        to the tail's bond, where the tail's side is inside the ring)."""
        P = self.preset
        ca, r = P.curved_arrow, P.lone_pair.dot_radius
        src = arr.source
        if src.kind == "bond":
            nb = self.bond_normal(src)
            side = 1.0 if dot(nb, bow) >= 0 else -1.0
            tail = add(s0, mul(nb, side * ca.bond_tail_gap))
        else:
            tail = add(ds.center, mul(ds.direction, r + ca.tail_gap))
        tkind, tdata = self.target(arr)
        if tkind == "bond":
            a, b = tdata
            m = mul(add(a.p, b.p), 0.5)
            nb = perp(unit(sub(b.p, a.p)))
            nb = nb if (dot(nb, sub(tail, m)) >= 0) != bool(far) else mul(nb, -1)
            fixed_tip, target = add(m, mul(nb, ca.head_gap)), ("bond", a.p, b.p)
        elif tdata.idx == -1 and tdata.box:  # a plus sign where fishhooks meet: just past it, on the bow side
            fixed_tip, (_, nm) = None, self.meeting(arr)
            d = nm if dot(nm, bow) >= 0 else mul(nm, -1)
            c = ((tdata.box[0] + tdata.box[2]) / 2, (tdata.box[1] + tdata.box[3]) / 2)
            target = ("point", add(c, mul(d, ray_box_exit(c, d, tdata.box) + ca.head_gap)))
        else:
            fixed_tip = None
            target = ("box", tdata.box) if tdata.box else ("point", tdata.p)
        dots = [q for d in (tdata.dots if tkind == "atom" else []) for q in d.dots]
        loop = length(sub(p_anchor, tail)) < P.lengths.bond_length
        handles = LOOP_HANDLES if loop else ARC_HANDLES
        h = o.curvature or soft * (ca.loop_curvature if loop else ca.curvature)

        def tip_for(t: Point) -> Point:
            if fixed_tip is not None:
                return fixed_tip
            box = target[1] if target[0] == "box" else None
            c = ((box[0] + box[2]) / 2, (box[1] + box[3]) / 2) if box else target[1]
            return add(c, mul(t, -clear_along(c, mul(t, -1), box, dots, r, ca.head_gap + extra)))

        def controls(end: Point) -> tuple[Point, Point]:
            ch = sub(end, tail)
            L, u = length(ch), unit(ch)
            n = perp(u) if dot(perp(u), bow) >= 0 else mul(perp(u), -1)
            return (add(tail, add(mul(u, L * handles[0][0]), mul(n, L * h * handles[0][1]))),
                    add(tail, add(mul(u, L * handles[1][0]), mul(n, L * h * handles[1][1]))))

        if fixed_tip is not None:
            end = fixed_tip
            for _ in range(ITERATIONS):
                c1, c2 = controls(end)
                t = unit(sub(end, c2))
                new_end = sub(fixed_tip, mul(t, ca.head_inset))
                done = length(sub(new_end, end)) < 1e-4
                end = new_end
                if done:
                    break
        else:
            ends = solve_end(tail, controls, tip_for, target, ca, dots, r)
            if pick >= len(ends):
                return None
            end = ends[pick]
        if pick and fixed_tip is not None:
            return None
        c1, c2 = controls(end)
        g = curved_arrow(tail, c1, c2, end, arr.type, P, name)
        hit = (tdata.species, tdata.idx) if tkind == "atom" and tdata.idx >= 0 else None
        return MechArrow(step, index, arr, tail, c1, c2, end, g, g.children[0].head_leaves()[0].points, src.kind, target, hit, dots)


def clear_along(c: Point, d: Point, box: tuple | None, dots: list[Point], r: float, gap: float) -> float:
    """How far from c along the unit direction d a point first keeps `gap` from the box (c
    inside it; c itself when None) and from the edges of the dots of radius r, and stays clear
    beyond: the nearest distance, not along the line through c, so a diagonal approach keeps
    its full gap and a head never pokes between a pair's dots."""
    s = gap if box is None else 0.0
    if box is not None:
        lo, hi = 0.0, ray_box_exit(c, d, box) + gap
        while point_box_dist(add(c, mul(d, hi)), box) < gap:  # a slanted ray needs more than gap
            hi += gap
        for _ in range(40):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if point_box_dist(add(c, mul(d, mid)), box) < gap else (lo, mid)
        s = hi
    for q in dots:
        w = sub(q, c)
        along, R = dot(w, d), r + gap
        off2 = dot(w, w) - along * along
        if off2 < R * R:
            s = max(s, along + math.sqrt(R * R - off2))
    return s


def miss(arrow: MechArrow) -> float:
    """Degrees between the head's direction and the line from its tip to the target's center
    (0 for a bond target, whose tip is fixed)."""
    kind, data = arrow.target[0], arrow.target[1]
    if kind == "bond":
        return 0.0
    c = data if kind == "point" else ((data[0] + data[2]) / 2, (data[1] + data[3]) / 2)
    t, w = sub(arrow.tip, arrow.p3), sub(c, arrow.tip)
    return math.degrees(math.acos(max(-1.0, min(1.0, dot(t, w) / (length(t) * length(w) or 1.0)))))


def aims(arrow: MechArrow, r: float) -> bool:
    """Whether the head points at its target: for a labeled atom, its line from the tip meets
    the element ink box without passing through one of the atom's dots (r their radius);
    otherwise within MISS_DEG of the target's center."""
    if arrow.target[0] != "box":
        return miss(arrow) <= MISS_DEG
    return _aims(arrow.tip, unit(sub(arrow.tip, arrow.p3)), arrow.target[1], arrow.target_dots, r)


def _aims(tip: Point, t: Point, box: tuple, dots: list[Point], r: float) -> bool:
    lo, hi = 0.0, math.inf  # slab test: where the ray tip + s t is inside the box
    for i in (0, 1):
        if abs(t[i]) < 1e-12:
            if not box[i] <= tip[i] <= box[i + 2]:
                return False
            continue
        a, b = (box[i] - tip[i]) / t[i], (box[i + 2] - tip[i]) / t[i]
        lo, hi = max(lo, min(a, b)), min(hi, max(a, b))
    if lo > hi:
        return False
    entry = add(tip, mul(t, lo))
    return all(seg_point_dist(q, tip, entry) >= r for q in dots)


def solve_end(tail: Point, controls, tip_for, target: tuple, ca, dots: list[Point] = (),
              r: float = 0.0) -> list[Point]:
    """Stroke ends for an atom target: the approach directions v (target -> tip) for which the
    end tangent points back along -v, scanned around the target and bisected. Several can
    exist (a loop over either side). Ordered for the planner to try: chords of at least two
    head lengths, shortest first, then the shorter ones down to one head length (a loop
    closing on itself beside the tail), longest first. Then, for a labeled target, the scanned ends whose heads still aim
    into its ink box past its dots (aims), the closest to the center first; with none of
    these, the closest miss."""
    c = target[1] if target[0] == "point" else ((target[1][0] + target[1][2]) / 2, (target[1][1] + target[1][3]) / 2)

    def at(th: float) -> tuple[float, float, Point, Point, Point]:
        v = (math.cos(th), math.sin(th))
        tip = tip_for(mul(v, -1))
        end = add(tip, mul(v, ca.head_inset))
        t = unit(sub(end, controls(end)[1]))
        return v[1] * t[0] - v[0] * t[1], -dot(t, v), end, tip, t  # cross(t, -v), dot(t, -v)

    n = SCAN
    ths = [2 * math.pi * k / n for k in range(n + 1)]
    vals = [at(th) for th in ths]
    roots = []
    for k in range(n):
        (fa, ga, *_), (fb, gb, *_) = vals[k], vals[k + 1]
        if fa * fb > 0 or min(ga, gb) <= 0:
            continue
        lo, hi = ths[k], ths[k + 1]
        for _ in range(40):
            mid = (lo + hi) / 2
            if at(lo)[0] * at(mid)[0] <= 0:
                hi = mid
            else:
                lo = mid
        f, g, end, _, _ = at((lo + hi) / 2)
        if abs(f) < 1e-3 and g > 0:  # not a jump (a dot pushing the tip out) posing as a root
            roots.append(end)
    chord = {id(e): length(sub(e, tail)) for e in roots}
    ok = sorted((e for e in roots if chord[id(e)] >= 2 * ca.head_length), key=lambda e: chord[id(e)])
    short = sorted((e for e in roots if ca.head_length <= chord[id(e)] < 2 * ca.head_length), key=lambda e: -chord[id(e)])

    def settle(tip: Point, end: Point) -> tuple[Point, Point] | None:
        """Off a root the head's own direction moves its tip: keep the tip and find the end
        (tip - w head_inset) whose tangent is w, the one nearest the scanned end's direction."""
        def h(ph: float) -> tuple[float, Point, Point]:
            w = (math.cos(ph), math.sin(ph))
            e = sub(tip, mul(w, ca.head_inset))
            t = unit(sub(e, controls(e)[1]))
            return (w[0] * t[1] - w[1] * t[0]) if dot(w, t) > 0 else math.nan, e, w
        w0 = unit(sub(tip, end))
        ph0 = math.atan2(w0[1], w0[0])
        phs = [ph0 + 2 * math.pi * k / 72 for k in range(73)]
        hs = [h(ph)[0] for ph in phs]
        found = []
        for k in range(72):
            if not (hs[k] * hs[k + 1] <= 0):  # also skips nan (tangent pointing back)
                continue
            lo, hi = phs[k], phs[k + 1]
            for _ in range(40):
                mid = (lo + hi) / 2
                lo, hi = (lo, mid) if h(lo)[0] * h(mid)[0] <= 0 else (mid, hi)
            f, e, w = h((lo + hi) / 2)
            if abs(f) < 1e-6:
                found.append((min(lo - ph0, 2 * math.pi - (lo - ph0)), e, w))
        return min(found)[1:] if found else None

    near = []
    if target[0] == "box":
        good = [k for k, (f, g, end, tip, t) in enumerate(vals[:n])
                if g > 0 and length(sub(end, tail)) >= ca.head_length and _aims(tip, t, target[1], dots, r)]
        runs = [k for k in good if (k - 1) % n not in good]  # one end per run of aiming directions
        for k0 in runs:
            run, k = [], k0
            while k in good and len(run) < n:
                run.append(k)
                k = (k + 1) % n
            near.append(min((vals[j] for j in run), key=lambda f: math.atan2(abs(f[0]), f[1])))
        near.sort(key=lambda f: math.atan2(abs(f[0]), f[1]))
        near = [x[0] for f, x in ((f, settle(f[3], f[2])) for f in near) if x and _aims(f[3], x[1], target[1], dots, r)]

    if ok or short or near:
        return ok + short + near
    f = min(vals, key=lambda f: (f[1] <= 0, abs(f[0])))
    x = settle(f[3], f[2])
    return [x[0] if x else f[2]]


def plan_step(step: int, arrows: list[CurvedArrow], overrides: list[ArrowOverride], preset: Preset,
              atoms: dict[int, AtomRef], bonds: set[frozenset], centroids: dict[str, Point],
              parts: dict[str, Parts] | None = None, pluses: list[tuple] | None = None,
              offset: int = 0) -> list[MechArrow]:
    """Each arrow is planned in turn against the ones before it; every combination of default
    or flipped bows is tried (arrows with a `bow` override take theirs), and the plan with the
    fewest clashes, then the fewest flips, then the lowest ranks, wins. So an early arrow gives way when its choice
    would force a later one into a bond (fishhooks meeting on one side, a bond's second arrow).
    If overrides leave no plan that keeps those side rules, they are dropped. `offset` numbers
    the arrows from there (a resonance step's later lists)."""
    best = None
    for strict, flips in itertools.product((True, False), itertools.product((False, True), repeat=len(arrows))):
        if any(f and overrides[k].bow for k, f in enumerate(flips)) or (best and not strict):
            continue
        pl = Planner(preset, atoms, bonds, centroids, parts, pluses, strict)
        out, total = [], (0, 0, 0)
        for k, a in enumerate(arrows):
            r = pl.plan(step, offset + k, a, overrides[k], f"curved-arrow-step-{step}-{offset + k}", flips[k])
            if r is None:
                break
            out.append(r[0])
            total = tuple(x + y for x, y in zip(total, r[1]))
        else:
            if best is None or total < best[0]:
                best = (total, out)
    return best[1]
