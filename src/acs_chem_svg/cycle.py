"""Catalytic cycles (acs-style.md §15; design.md §10): species on a circle, arc arrows between
them, and tails for species entering and leaving.

ChemDraw has no cycle layout, so every rule here is ours (the presets' `cycle:` keys,
unmeasured; user decision 2026-10-29 on the geometry spike, since retired). Only the arc
arrows are measured (arrows.arc_arrow, §15.1).

- Circle: centered at (0, 0); angles in radians counter-clockwise on screen from +x, so the
  start species sits at 90° (12 o'clock) and a clockwise cycle runs the angles down.
- Species: each species' box (its ink, plus its compound number) is centered on the circle. The
  circle leaves the box `clear` (Euclidean) from its ink at θk ∓ αk, exactly (bisection). The θk
  make every free arc between neighbouring boxes the same length (a fixed-point loop from
  θ0 = 90°), so a wide species takes more of the circle.
- Radius: the smallest one where the free arc reaches the longest step's need
  (scheme.arrow_min_length, plus tail_room per tail on that step's busier side), the boxes keep
  `clear` from each other, and the step labels keep label_clear from the species and the center
  text and label_pair_clear from each other. A label sits on the inward line through its arc's middle, as far out as its
  ink stays label_gap inside the circle. Tails never change the radius.
- Arc arrows: each step's arrow fills its free arc (arrows.arc_step_arrow).
- Tails: an arc of tail_radius, tangent to the circle from outside at its join, turning
  tail_turn away from it: an entering tail upstream (drawn from its free end to the join), a
  leaving one downstream (an arc arrow, its head on the chord as §15.1). Consecutive text entries
  on one side share a tail ("PhCHO, H+"); each species gets its own. The last entering tail joins
  at tail_join_in of the arc and earlier ones tail_room (along the arc) further back; the first
  leaving one at tail_join_out and later ones further on. On an equilibrium step a tail joins
  the outer shaft. Its text or structure sits tail_text_gap beyond the free end, along the tail. A tail that
  would touch a species, a label or another tail's item tries a larger radius, other turns and
  nearby joins (_TRIES), in that order; if none is clear it keeps the default and is reported
  (TailLayout.clash; validate.py).

Groups: species-<name>, numbers (number-<name>), cycle-arc-<n> (n from the start species round),
cycle-in-<n>-<i> / cycle-out-<n>-<i> (the tail, then its text runs or its species group),
cycle-label-<n>, cycle-center.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path as FsPath

from . import richtext
from .arrows import arc_arrow, arc_step_arrow
from .geometry import Point, add, arc_point, arc_segments, mul, point_box_dist, union_box
from .geometry_checks import check_labels
from .placed import Placed, place
from .preset import Preset, load_preset
from .scene import Group, Path, TextRun
from .spec import CycleStep, Spec, scheme_molecules

CENTER = (0.0, 0.0)
# a crowded tail's alternatives, tried in order after the default: radius factor, turn (deg) and
# join shift (fraction of the arc, away from the arc's middle first)
_TRIES = [(f, t, dj) for f in (1.0, 1.3, 1.6) for t in (None, 65.0, 100.0, 50.0, 120.0)
          for dj in (0.0, -0.06, 0.06, -0.12)]


@dataclass
class ArcLayout:
    n: int
    a0: float  # tail angle
    a1: float  # tip angle
    kind: str
    group: Group
    source: str
    target: str


@dataclass
class TailLayout:
    n: int  # step
    side: str  # "in" | "out"
    index: int  # 1-based, in spec order
    join: float  # angle where it meets the circle of radius r_join
    r_join: float
    radius: float  # its own arc's radius
    turn: float  # degrees
    free_end: Point
    group: Group
    text: str | None = None
    species: str | None = None
    clash: str | None = None  # what it still touches (None: clear)
    arc: tuple[Point, float, float] | None = None  # its own circle's center and angles, as drawn (from, to)


@dataclass
class CycleLayout:
    root: Group
    preset: Preset
    radius: float
    direction: float  # +1 counter-clockwise on screen, -1 clockwise
    free: float  # free arc between neighbouring boxes, pt along the circle
    species: dict[str, Placed] = field(default_factory=dict)
    theta: dict[str, float] = field(default_factory=dict)
    arcs: list[ArcLayout] = field(default_factory=list)
    tails: list[TailLayout] = field(default_factory=list)
    labels: dict[int, list[TextRun]] = field(default_factory=dict)
    numbers: dict[str, list[TextRun]] = field(default_factory=dict)
    center_text: list[TextRun] = field(default_factory=list)

    @property
    def center(self) -> Point:
        return CENTER


# --------------------------------------------------------------------------- solver


def _box_at(half: Point, p: Point):
    return (p[0] - half[0], p[1] - half[1], p[0] + half[0], p[1] + half[1])


def _gap(a, b) -> float:
    """Distance between two boxes (0 if they overlap)."""
    return math.hypot(max(a[0] - b[2], b[0] - a[2], 0.0), max(a[1] - b[3], b[1] - a[3], 0.0))


def half_angles(R: float, theta: float, half: Point, clear: float) -> tuple[float, float] | None:
    """(α counter-clockwise, α clockwise): how far round the circle of radius R, each way from θ,
    a point first lies `clear` from the box of half-extents `half` centered at θ. None if it
    never does within a half turn (R too small for the box)."""
    box = _box_at(half, arc_point(CENTER, R, theta))
    out = []
    step = min(0.05, max(clear, 0.5) / R)
    for s in (1.0, -1.0):
        lo, d, hi = 0.0, step, None
        while d < math.pi:
            if point_box_dist(arc_point(CENTER, R, theta + s * d), box) >= clear:
                hi = d
                break
            lo, d = d, d + step
        if hi is None:
            return None
        for _ in range(40):
            mid = (lo + hi) / 2
            lo, hi = (lo, mid) if point_box_dist(arc_point(CENTER, R, theta + s * mid), box) >= clear else (mid, hi)
        out.append(hi)
    return out[0], out[1]


def _lead_trail(al: tuple[float, float], sign: float) -> tuple[float, float]:
    """(α toward the travel, α against it)."""
    return (al[0], al[1]) if sign > 0 else (al[1], al[0])


def solve_thetas(R: float, halves: list[Point], clear: float, sign: float):
    """(θk with θ0 = 90° and equal free arcs, (α ccw, α cw) per box, free angle per step), or
    None if a box swallows the circle or the boxes fill it."""
    n = len(halves)
    thetas = [math.pi / 2 + sign * 2 * math.pi * k / n for k in range(n)]
    new = thetas
    for _ in range(100):
        al = [half_angles(R, t, h, clear) for t, h in zip(thetas, halves)]
        if any(a is None for a in al):
            return None
        lt = [_lead_trail(a, sign) for a in al]
        free = 2 * math.pi - sum(a + b for a, b in lt)
        if free <= 0:
            return None
        new = [math.pi / 2]
        for k in range(1, n):
            new.append(new[-1] + sign * (lt[k - 1][0] + free / n + lt[k][1]))
        if max(abs(a - b) for a, b in zip(new, thetas)) < 1e-9:
            break
        thetas = new
    return new, al, free / n


def arc_ends(thetas: list[float], al: list, sign: float) -> list[tuple[float, float]]:
    """(tail angle, tip angle) of each step's free arc, step k from box k to box k + 1 (the last
    back to the first, unwrapped so the angles run the travel's way)."""
    n = len(thetas)
    out = []
    for k in range(n):
        nk = (k + 1) % n
        a0 = thetas[k] + sign * _lead_trail(al[k], sign)[0]
        a1 = thetas[nk] - sign * _lead_trail(al[nk], sign)[1] + (sign * 2 * math.pi if nk == 0 else 0.0)
        out.append((a0, a1))
    return out


def _text_beyond(line: richtext.Line, anchor: Point, direction: Point, gap: float) -> None:
    """Move `line` so its ink box sits `gap` beyond `anchor` along the unit `direction`."""
    b = line.ink()
    _move_box_beyond(lambda dx, dy: line.shift(dx, dy), b, anchor, direction, gap)


def _move_box_beyond(shift, box, anchor: Point, direction: Point, gap: float) -> None:
    w, h = box[2] - box[0], box[3] - box[1]
    support = abs(direction[0]) * w / 2 + abs(direction[1]) * h / 2
    cx, cy = add(anchor, mul(direction, gap + support))
    shift(cx - (box[0] + box[2]) / 2, cy - (box[1] + box[3]) / 2)


def _label_line(text: str, R: float, a0: float, a1: float, preset: Preset, n: int) -> richtext.Line:
    """A step label inside the circle, on the inward line through the arc's middle: as far out
    as it goes with its ink at least label_gap inside the circle (every corner of its box)."""
    mid = (a0 + a1) / 2
    gap = preset.cycle.label_gap
    ln = richtext.layout_line(text, 0.0, 0.0, preset.text.caption_size, preset, f"cycle-label:{n}")
    b = ln.ink()
    hw, hh = (b[2] - b[0]) / 2, (b[3] - b[1]) / 2
    u = (math.cos(mid), -math.sin(mid))  # outward

    def reach(d: float) -> float:  # farthest corner from the center, box center d out along u
        cx, cy = mul(u, d)
        return max(math.hypot(cx + sx * hw, cy + sy * hh) for sx in (-1, 1) for sy in (-1, 1))
    lo, hi = 0.0, R - gap
    if reach(lo) > R - gap:  # too big for the circle even at its center: centered (the radius grows)
        d = 0.0
    else:
        for _ in range(50):
            m = (lo + hi) / 2
            lo, hi = (m, hi) if reach(m) <= R - gap else (lo, m)
        d = lo
    cx, cy = mul(u, d)
    ln.shift(cx - (b[0] + b[2]) / 2, cy - (b[1] + b[3]) / 2)
    return ln


def solve(halves: list[Point], need: float, preset: Preset, sign: float, labels: list[str | None],
          center_box=None):
    """The radius (module docstring) and its (θk, α per box, free angle)."""
    c = preset.cycle

    def ok(R: float) -> bool:
        s = solve_thetas(R, halves, c.clear, sign)
        if s is None or R * s[2] < need:
            return False
        boxes = [_box_at(h, arc_point(CENTER, R, t)) for t, h in zip(s[0], halves)]
        if any(_gap(boxes[i], boxes[j]) < c.clear for i in range(len(boxes)) for j in range(i + 1, len(boxes))):
            return False
        lb = [_label_line(t, R, a0, a1, preset, k).ink() for k, (t, (a0, a1)) in
              enumerate(zip(labels, arc_ends(s[0], s[1], sign))) if t]
        others = boxes + ([center_box] if center_box else [])
        return all(_gap(a, b) >= c.label_clear for a in lb for b in others) and \
            all(_gap(lb[i], lb[j]) >= c.label_pair_clear for i in range(len(lb)) for j in range(i + 1, len(lb)))

    lo, hi = 1.0, 40.0
    while not ok(hi):
        lo, hi = hi, hi * 1.5
        if hi > 1e5:
            raise ValueError("cycle: no radius fits these species")
    while hi - lo > 0.005:
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if ok(mid) else (mid, hi)
    return hi, solve_thetas(hi, halves, c.clear, sign)


# --------------------------------------------------------------------------- tails


def tail_path(R: float, join: float, radius: float, turn_deg: float, sign: float, incoming: bool,
              preset: Preset, name: str, variant: str | None = None) -> tuple[Path, Point, Point, tuple]:
    """A tail on the circle of radius R at angle `join` (module docstring): its Path, its free
    end, the unit direction leaving the tail there (away from the join), and its arc as drawn
    (center, from angle, to angle: free end -> join entering, join -> free end leaving)."""
    J = arc_point(CENTER, R, join)
    o = add(J, mul((math.cos(join), -math.sin(join)), radius))  # the tail circle's center, outside
    beta = join + math.pi  # J seen from o
    # going the travel's way at J, the tail circle runs the other way round
    turn = math.radians(turn_deg) * (1.0 if incoming else -1.0) * sign
    end = beta + turn
    F = arc_point(o, radius, end)
    d = 1.0 if end > beta else -1.0
    away = (-math.sin(end) * d, -math.cos(end) * d)
    w = preset.arrow.shaft_width
    if incoming:
        return (Path(arc_segments(o, radius, end, beta), fill=False, stroke=w, closed=False, source=name), F, away,
                (o, end, beta))
    # the head on its chord (§15.1)
    return arc_arrow(o, radius, beta, end, preset, name, variant).children[0], F, away, (o, beta, end)


def _samples(path: Path, n: int = 8) -> list[Point]:
    out, cur = [], None
    for seg in path.segments:
        if seg[0] == "C":
            c1, c2, p3 = seg[1:]
            for k in range(1, n + 1):
                t, s = k / n, 1 - k / n
                out.append(tuple(s ** 3 * cur[j] + 3 * s * s * t * c1[j] + 3 * s * t * t * c2[j] + t ** 3 * p3[j]
                                 for j in range(2)))
        else:
            out.append(seg[1])
        cur = seg[-1]
    return out + [h for leaf in path.head_leaves() for h in leaf.points]


# --------------------------------------------------------------------------- build


def _drawers(spec: Spec, base_dir: FsPath, preset: Preset) -> dict[str, Placed]:
    return {n: place(n, mol, preset, spec.species[n].show_lone_pairs)
            for n, mol in scheme_molecules(spec, base_dir).items()}


def _number(spec: Spec, p: Placed, preset: Preset) -> list[TextRun]:
    """The species' bold compound number, ink-centered under its ink, number_gap below."""
    if p.name not in (spec.numbers or {}):
        return []
    ln = richtext.layout_line(str(spec.numbers[p.name]), 0.0, 0.0, preset.text.caption_size, preset,
                              f"number:{p.name}", style="bold")
    b, ink = p.group.ink_bbox(), ln.ink()
    ln.shift((b[0] + b[2]) / 2 - (ink[0] + ink[2]) / 2, b[3] + preset.scheme.number_gap - ink[1])
    return ln.runs


def _shift_runs(runs: list[TextRun], dx: float, dy: float) -> None:
    for r in runs:
        r.x += dx
        r.y += dy


def build_cycle(spec: Spec, base_dir: FsPath, name: str = "cycle") -> CycleLayout:
    cyc = spec.cycle
    preset = load_preset(spec.preset)
    c = preset.cycle
    sign = -1.0 if cyc.direction == "clockwise" else 1.0
    placed = _drawers(spec, base_dir, preset)
    members = cyc.members
    steps: list[CycleStep] = cyc.ordered_steps()

    numbers = {n: _number(spec, placed[n], preset) for n in members}
    halves = []
    for n in members:  # each member's box (ink and number) centered on the origin
        p = placed[n]
        b = union_box([p.group.ink_bbox()] + [r.bbox() for r in numbers[n]])
        dx, dy = -(b[0] + b[2]) / 2, -(b[1] + b[3]) / 2
        p.translate(dx, dy)
        _shift_runs(numbers[n], dx, dy)
        halves.append(((b[2] - b[0]) / 2, (b[3] - b[1]) / 2))

    center_line = None
    if cyc.center:
        center_line = richtext.layout_line(cyc.center, 0.0, 0.0, preset.text.caption_size, preset, "cycle-center")
        b = center_line.ink()
        center_line.shift(-(b[0] + b[2]) / 2, -(b[1] + b[3]) / 2)

    need = max(preset.scheme.arrow_min_length + c.tail_room * max(len(tail_items(st.in_)), len(tail_items(st.out)))
               for st in steps)
    R, (thetas, al, free) = solve(halves, need, preset, sign, [st.label for st in steps],
                                  center_line.ink() if center_line else None)
    out = CycleLayout(Group(name), preset, R, sign, R * free)
    boxes = {}
    for n, t in zip(members, thetas):
        p = placed[n]
        x, y = arc_point(CENTER, R, t)
        p.translate(x, y)
        _shift_runs(numbers[n], x, y)
        out.species[n], out.theta[n] = p, t
        out.root.children.append(p.group)
        boxes[n] = union_box([p.group.ink_bbox()] + [r.bbox() for r in numbers[n]])

    ends = arc_ends(thetas, al, sign)
    label_boxes = []
    for k, (st, (a0, a1)) in enumerate(zip(steps, ends), 1):
        g = arc_step_arrow(st.arrow_type, CENTER, R, a0, a1, preset, f"cycle-arc-{k}", spec.arrow_head)
        out.arcs.append(ArcLayout(k, a0, a1, st.arrow_type, g, st.from_, st.to))
        out.root.children.append(g)
        if st.label:
            ln = _label_line(st.label, R, a0, a1, preset, k)
            out.labels[k] = ln.runs
            label_boxes.append(ln.ink())
            out.root.children.append(Group(f"cycle-label-{k}", ln.runs))
    if center_line:
        out.center_text = center_line.runs
        out.root.children.append(Group("cycle-center", center_line.runs))
        label_boxes.append(center_line.ink())

    obstacles = list(boxes.values()) + label_boxes
    for k, (st, (a0, a1)) in enumerate(zip(steps, ends), 1):
        r_join = R + (preset.equilibrium.shaft_spacing / 2 if st.arrow_type != "forward" else 0.0)
        spacing = c.tail_room / (R * abs(a1 - a0))  # fraction of the arc between tails on one side
        for side, tails in (("in", st.in_), ("out", st.out)):
            items = tail_items(tails)
            for i, (text, species) in enumerate(items):
                frac = (c.tail_join_in - spacing * (len(items) - 1 - i) if side == "in"
                        else c.tail_join_out + spacing * i)
                lay = _place_tail(spec, out, placed, k, side, i + 1, text, species, a0, a1, frac, r_join, obstacles)
                obstacles.append(lay.group.ink_bbox())
                out.tails.append(lay)
                out.root.children.append(lay.group)

    if any(numbers.values()):
        out.numbers = {n: runs for n, runs in numbers.items() if runs}
        out.root.children.append(Group("numbers", [Group(f"number-{n}", runs) for n, runs in out.numbers.items()]))
    return out


def tail_items(tails) -> list[tuple[str | None, str | None]]:
    """(text, species) per tail drawn on one side of a step: consecutive text entries share one
    tail, joined by ", " (reagents written together); each species gets its own."""
    out: list[list] = []
    for t in tails:
        if t.text is not None and out and out[-1][0] is not None:
            out[-1][0] += ", " + t.text
        else:
            out.append([t.text, t.species])
    return [tuple(x) for x in out]


def _place_tail(spec: Spec, out: CycleLayout, placed: dict[str, Placed], k: int, side: str, i: int,
                text: str | None, species: str | None, a0: float, a1: float, frac: float, r_join: float,
                obstacles: list) -> TailLayout:
    """One tail, with the first of _TRIES whose path and item stay clear (module docstring)."""
    preset, c, sign = out.preset, out.preset.cycle, out.direction
    incoming = side == "in"
    name = f"cycle-{side}-{k}-{i}"
    away_mid = -1.0 if incoming else 1.0  # join shifts away from the arc's middle first
    p = placed.get(species) if species else None
    first = None
    for f, turn, dj in _TRIES:
        radius, turn = c.tail_radius * f, turn if turn is not None else c.tail_turn
        fr = min(max(frac + away_mid * dj, 0.05), 0.95)
        join = a0 + fr * (a1 - a0)
        path, F, away, arc = tail_path(r_join, join, radius, turn, sign, incoming, preset, name, spec.arrow_head)
        if p is not None:
            b = p.group.ink_bbox()
            _move_box_beyond(p.translate, b, F, away, c.tail_text_gap)
            item, item_box = p.group, p.group.ink_bbox()
        else:
            ln = richtext.layout_line(text, 0.0, 0.0, preset.text.caption_size, preset, f"cycle-{side}:{k}")
            _text_beyond(ln, F, away, c.tail_text_gap)
            item, item_box = ln.runs, ln.ink()
        clash = _tail_clash(path, item_box, obstacles, out, c.label_clear)
        lay = TailLayout(k, side, i, join, r_join, radius, turn, F,
                         Group(name, [path] + ([item] if p is not None else item)), text, species, clash, arc)
        if first is None:
            first = lay
        if clash is None:
            break
    else:
        lay = first
        if p is not None:  # the species was moved by the later tries: put it back at the default
            path0 = lay.group.children[0]
            _, F, away, _ = tail_path(r_join, lay.join, lay.radius, lay.turn, sign, incoming, preset, name, spec.arrow_head)
            _move_box_beyond(p.translate, p.group.ink_bbox(), F, away, c.tail_text_gap)
            lay.group.children = [path0, p.group]
    if p is not None:
        out.species[species] = p
    return lay


def _tail_clash(path: Path, item_box, obstacles: list, out: CycleLayout, clear: float) -> str | None:
    """What the tail touches: its path within `clear` of an obstacle box (species, labels, the
    center text, earlier tails), or its item within `clear` of one or of the circle's arrows."""
    pts = _samples(path)
    for b in obstacles:
        if _gap(item_box, b) < clear:
            return "its text or structure is too close to another item"
        if min(point_box_dist(q, b) for q in pts) < clear:
            return "its path runs too close to another item"
    for a in out.arcs:  # the circle's arrows, sampled about every 2 pt
        n = max(2, int(out.radius * abs(a.a1 - a.a0) / 2))
        if any(point_box_dist(arc_point(CENTER, out.radius, a.a0 + (a.a1 - a.a0) * j / n), item_box) < clear
               for j in range(n + 1)):
            return f"its text or structure is too close to the arrow of step {a.n}"
    return None


def check_cycle(lay: CycleLayout) -> list[str]:
    """Problems: each species' label checks (geometry_checks.check_labels) and every tail still
    crowded after its tries."""
    out = [p for n, pl in lay.species.items() for p in check_labels(n, pl.group, pl.mol)]
    return out + [f"cycle step {t.n}: {'entering' if t.side == 'in' else 'leaving'} tail {t.index} "
                  f"({t.text or t.species}): {t.clash}; shorten its text (tails never move the species)"
                  for t in lay.tails if t.clash]


def render_cycle(spec: Spec, base_dir: FsPath, name: str) -> Group:
    return build_cycle(spec, base_dir, name).root
