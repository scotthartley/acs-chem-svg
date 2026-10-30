"""Catalytic-cycle geometry spike (plan phase 1; no schema yet): species on a circle, the first
at 12 o'clock, equal free arcs between their ink boxes, arc arrows, and branch tails for species
entering and leaving.

Solver: species k's box sits centered on the circle at angle θk (y down, θ counterclockwise
from +x, so 12 o'clock is 90° and clockwise runs θ down). The circle leaves the box's
`clear`-offset outline (Euclidean distance `clear` from its ink) at θk ∓ αk; the αk follow from
R and θk exactly (bisection on the distance). The θk are set so every free arc
(θk − αk⁻) − (θk+1 + αk+1⁺) is the same, with θ0 = 90°, by a fixed-point loop; R is the smallest
radius whose free arc reaches the longest step minimum (arrow_min_length; more for a step with
tails) and whose boxes keep `clear` from each other. Tails don't change R.

Every number here is a spike value (unmeasured), to settle the rules by eye.

Usage: uv run python scripts/cycle_spike.py [--out spike/cycle] [--case NAME]
Writes <case>.svg and <case>.png (4 px/pt) per case and prints R, the θk and the free arc.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, field
from pathlib import Path

from acs_chem_svg import richtext
from acs_chem_svg.arrows import head_marker
from acs_chem_svg.geometry import Point, add, mul, point_box_dist, sub
from acs_chem_svg.layout import m1_placed
from acs_chem_svg.preset import Preset, load_preset
from acs_chem_svg.preview import rasterize
from acs_chem_svg.scene import Group, Path as SPath
from acs_chem_svg.spec import Spec
from acs_chem_svg.svg_writer import to_svg

ROOT = Path(__file__).resolve().parent.parent

# spike values (pt, unmeasured)
CLEAR = 9.0           # species ink -> arc arrow tail / tip
TAIL_ROOM = 24.0      # extra free arc a step with tails needs
TAIL_JOIN = (0.38, 0.62)  # where an incoming / outgoing tail meets the arc (fraction of the arc)
TAIL_RADIUS = 20.0    # radius of a tail's arc (tangent to the circle from outside)
TAIL_TURN = 80.0      # degrees a tail turns from its join to its free end
TAIL_TEXT_GAP = 2.5   # tail end -> its text
LABEL_GAP = 6.0       # arc -> step label (inside the circle)
LABEL_CLEAR = 4.0     # step label ink <-> species ink, other labels, center text


@dataclass
class CycleStep:
    label: str = ""
    kind: str = "forward"  # forward | equilibrium
    tail_in: list[str] = field(default_factory=list)  # text
    tail_out: list[str] = field(default_factory=list)


@dataclass
class Case:
    name: str
    species: dict[str, dict]  # name -> Species fields
    steps: list[CycleStep]  # steps[k]: species k -> k+1 (last closes)
    center: str = ""
    clockwise: bool = True
    preset: str = "acs1996"


L2 = "<-P(C)(C)C"  # PMe3 in full, bound by a dative bond
CASES = [
    Case("suzuki", {
        "pd0": {"smiles": f"[Pd]({L2}){L2}"},
        "oxadd": {"smiles": f"Br[Pd](c1ccccc1)({L2}){L2}"},
        "trans": {"smiles": f"c1ccc(cc1)[Pd](c1ccc(C)cc1)({L2}){L2}"},
    }, [
        CycleStep("oxidative addition", tail_in=["PhBr"]),
        CycleStep("transmetalation", tail_in=["ArB(OH)2, base"], tail_out=["B(OH)3, Br–"]),
        CycleStep("reductive elimination", tail_out=["Ph–Ar"]),
    ], center="Pd(0)"),
    Case("enamine4", {
        "cat": {"smiles": "OC(=O)[C@@H]1CCCN1"},
        "iminium": {"smiles": "OC(=O)[C@@H]1CCC[N+]1=C1CCCCC1"},
        "enamine": {"smiles": "OC(=O)[C@@H]1CCCN1C1=CCCCC1"},
        "adduct": {"smiles": "OC(=O)[C@@H]1CCC[N+]1=C1CCCCC1C(O)c1ccccc1"},
    }, [
        CycleStep("condensation", "equilibrium", tail_in=["cyclohexanone"], tail_out=["H2O"]),
        CycleStep("tautomerization", "equilibrium"),
        CycleStep("C–C bond formation", tail_in=["PhCHO"]),
        CycleStep("hydrolysis", tail_in=["H2O"], tail_out=["aldol product"]),
    ]),
    Case("six_mixed", {
        "a": {"smiles": "c1ccccc1"},
        "b": {"smiles": "CC(=O)Oc1ccccc1C(=O)O"},
        "c": {"smiles": "O"},
        "d": {"smiles": "c1ccc2ccccc2c1"},
        "e": {"smiles": "CCCCCCCC"},
        "f": {"smiles": "CN1CCC[C@H]1c1cccnc1"},
    }, [CycleStep(f"step {k}") for k in range(1, 7)], center="six species"),
    Case("five_ccw", {
        "a": {"smiles": "CC(=O)O"},
        "b": {"smiles": "c1ccc2ccccc2c1"},
        "c": {"smiles": "OCCO"},
        "d": {"smiles": "CC(C)(C)OC(=O)N"},
        "e": {"smiles": "C1CCCCC1"},
    }, [CycleStep(kind="equilibrium"), CycleStep(tail_in=["X"]), CycleStep(), CycleStep(tail_out=["Y"]), CycleStep()],
        clockwise=False),
]


def at(c: Point, r: float, phi: float) -> Point:
    return (c[0] + r * math.cos(phi), c[1] - r * math.sin(phi))


def outward(phi: float) -> Point:
    return (math.cos(phi), -math.sin(phi))


def travel(phi: float, sign: float) -> Point:
    """Unit tangent of travel at phi (sign -1: clockwise, θ falling)."""
    return (-math.sin(phi) * sign, -math.cos(phi) * sign)


def box_at(half: tuple[float, float], p: Point):
    return (p[0] - half[0], p[1] - half[1], p[0] + half[0], p[1] + half[1])


def half_angles(R: float, theta: float, half: tuple[float, float], clear: float) -> tuple[float, float] | None:
    """(α ccw, α cw): where the circle reaches `clear` from the box centered at θ, each side.
    None if it never does within a half turn (R too small)."""
    box = box_at(half, at((0.0, 0.0), R, theta))
    out = []
    for s in (1.0, -1.0):
        lo, hi, step = 0.0, None, 0.01
        d = step
        while d < math.pi:
            if point_box_dist(at((0.0, 0.0), R, theta + s * d), box) >= clear:
                hi = d
                break
            lo, d = d, d + step
        if hi is None:
            return None
        for _ in range(50):
            mid = (lo + hi) / 2
            if point_box_dist(at((0.0, 0.0), R, theta + s * mid), box) >= clear:
                hi = mid
            else:
                lo = mid
        out.append(hi)
    return out[0], out[1]


def solve_thetas(R: float, halves: list, clear: float, sign: float):
    """θk (θ0 = 90°) with equal free arcs, and the (α ccw, α cw) per species; None if a box
    swallows the circle or the boxes fill it."""
    n = len(halves)
    thetas = [math.pi / 2 + sign * 2 * math.pi * k / n for k in range(n)]
    for _ in range(100):
        al = [half_angles(R, t, h, clear) for t, h in zip(thetas, halves)]
        if any(a is None for a in al):
            return None
        # leaving side / arriving side in the direction of travel
        lead = [a[1] if sign < 0 else a[0] for a in al]  # toward travel
        trail = [a[0] if sign < 0 else a[1] for a in al]
        free = 2 * math.pi - sum(lead) - sum(trail)
        if free <= 0:
            return None
        new = [math.pi / 2]
        for k in range(1, n):
            new.append(new[-1] + sign * (lead[k - 1] + free / n + trail[k]))
        if max(abs(a - b) for a, b in zip(new, thetas)) < 1e-10:
            break
        thetas = new
    return new, al, free / n


def boxes_clear(R: float, thetas, halves, clear: float) -> bool:
    boxes = [box_at(h, at((0.0, 0.0), R, t)) for t, h in zip(thetas, halves)]
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i], boxes[j]
            dx = max(a[0] - b[2], b[0] - a[2], 0.0)
            dy = max(a[1] - b[3], b[1] - a[3], 0.0)
            if math.hypot(dx, dy) < clear:
                return False
    return True


def _gap(a, b) -> float:
    dx = max(a[0] - b[2], b[0] - a[2], 0.0)
    dy = max(a[1] - b[3], b[1] - a[3], 0.0)
    return math.hypot(dx, dy)


def solve(halves, need: float, clear: float, sign: float, labels: list[str], center: tuple | None, preset: Preset):
    """Smallest R whose free arc is >= need (bisection), whose boxes keep `clear` apart and
    whose inner step labels keep LABEL_CLEAR from the species, each other and the center
    text; with its θk, α and free arc length."""
    def ok(R):
        s = solve_thetas(R, halves, clear, sign)
        if s is None or R * s[2] < need or not boxes_clear(R, s[0], halves, clear):
            return False
        lb = [b for b in label_boxes(R, s[0], s[1], halves, labels, sign, preset) if b]
        sp = [box_at(h, at((0.0, 0.0), R, t)) for t, h in zip(s[0], halves)] + ([center] if center else [])
        return all(_gap(a, b) >= LABEL_CLEAR for a in lb for b in sp) and \
            all(_gap(lb[i], lb[j]) >= LABEL_CLEAR for i in range(len(lb)) for j in range(i + 1, len(lb)))
    lo, hi = 1.0, 50.0
    while not ok(hi):
        lo, hi = hi, hi * 1.5
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if ok(mid) else (mid, hi)
    th, al, free = solve_thetas(hi, halves, clear, sign)
    return hi, th, al, free


def arc_segments(c: Point, r: float, a0: float, a1: float) -> list[tuple]:
    """Cubics for the arc a0 -> a1 (radians, either direction), each at most 45°."""
    n = max(1, math.ceil(abs(a1 - a0) / (math.pi / 4)))
    d = (a1 - a0) / n
    k = 4 / 3 * math.tan(d / 4)
    segs = [("M", at(c, r, a0))]
    for i in range(n):
        p, q = a0 + i * d, a0 + (i + 1) * d
        P, Q = at(c, r, p), at(c, r, q)
        tp = (-math.sin(p), -math.cos(p))  # d/dphi of at(), unit
        tq = (-math.sin(q), -math.cos(q))
        segs.append(("C", add(P, mul(tp, k * r)), sub(Q, mul(tq, k * r)), Q))
    return segs


def arc_arrow(c, R, a0, a1, preset: Preset, name: str) -> Group:
    """Forward arrow along the circle, tail at a0, tip at a1; the shaft stops
    shaft_end_from_tip short of the tip (measured along the arc), the head on its end tangent."""
    h = preset.head()
    s = 1.0 if a1 > a0 else -1.0
    end = a1 - s * h.shaft_end_from_tip / R
    return Group(name, [SPath(arc_segments(c, R, a0, end), fill=False, stroke=preset.arrow.shaft_width, closed=False,
                              source=name, marker_end=head_marker(preset))])


def arc_equilibrium(c, R, a0, a1, preset: Preset, name: str, clockwise: bool) -> Group:
    """Two concentric half-headed shafts (equilibrium_arrow's balanced rule, lengths along the
    arc): forward on the left of travel (outside for clockwise travel on screen)."""
    h, e = preset.head(), preset.equilibrium
    s = 1.0 if a1 > a0 else -1.0
    r_fwd = R + (e.shaft_spacing / 2 if clockwise else -e.shaft_spacing / 2)
    r_back = 2 * R - r_fwd
    st, se = e.stagger / R, h.shaft_end_from_tip / R

    def shaft(r, b, t, half):
        d = 1.0 if t > b else -1.0
        return SPath(arc_segments(c, r, b, t - d * se), fill=False, stroke=preset.arrow.shaft_width, closed=False,
                     source=name, marker_end=head_marker(preset, half=half))
    return Group(name, [shaft(r_fwd, a0 + s * st, a1, "forward"), shaft(r_back, a1 - s * st, a0, "back")])


def text_at(text: str, anchor: Point, direction: Point, gap: float, preset: Preset, src: str):
    """A caption line whose ink box sits `gap` beyond `anchor` along unit `direction`."""
    ln = richtext.layout_line(text, 0.0, 0.0, preset.text.caption_size, preset, src)
    b = ln.ink()
    w, hgt = b[2] - b[0], b[3] - b[1]
    support = abs(direction[0]) * w / 2 + abs(direction[1]) * hgt / 2
    cx, cy = add(anchor, mul(direction, gap + support))
    ln.shift(cx - (b[0] + b[2]) / 2, cy - (b[1] + b[3]) / 2)
    return ln


def tail(c, R, join: float, preset: Preset, name: str, sign: float, incoming: bool) -> tuple[SPath, Point, Point]:
    """A branch: an arc of radius TAIL_RADIUS tangent to the circle at `join` from outside,
    turning TAIL_TURN away from it: upstream for an incoming tail (drawn free end -> join),
    downstream for an outgoing one (join -> free end, with a head). Returns the path, its free
    end and the unit direction pointing away from the tail there (for the text)."""
    J = at(c, R, join)
    o = add(J, mul(outward(join), TAIL_RADIUS))  # the tail circle's center
    beta = join + math.pi  # J as seen from o
    # walking forward along the main circle (sign) at J, the tail circle runs the other way round
    turn = math.radians(TAIL_TURN) * (1.0 if incoming else -1.0) * sign
    end = beta + turn
    F = at(o, TAIL_RADIUS, end)
    if incoming:
        segs = arc_segments(o, TAIL_RADIUS, end, beta)
        return SPath(segs, fill=False, stroke=preset.arrow.shaft_width, closed=False, source=name), F, \
            _away(o, end, beta)
    hd = preset.head()
    stop = end - (turn / abs(turn)) * hd.shaft_end_from_tip / TAIL_RADIUS
    segs = arc_segments(o, TAIL_RADIUS, beta, stop)
    return SPath(segs, fill=False, stroke=preset.arrow.shaft_width, closed=False, source=name,
                 marker_end=head_marker(preset)), F, _away(o, end, beta)


def _away(o: Point, end: float, beta: float) -> Point:
    """Unit tangent at the free end, pointing off the tail (away from the join)."""
    d = 1.0 if end > beta else -1.0
    return (-math.sin(end) * d, -math.cos(end) * d)


def label_boxes(R: float, thetas, al, halves, labels: list[str], sign: float, preset: Preset):
    """Step label ink boxes inside the circle (arc midpoints, LABEL_GAP in) with the circle at
    (0, 0); None where a step has no label."""
    out = []
    n = len(thetas)
    for k, txt in enumerate(labels):
        if not txt:
            out.append(None)
            continue
        nk = (k + 1) % n
        a0 = thetas[k] + sign * (al[k][1] if sign < 0 else al[k][0])
        a1 = thetas[nk] - sign * (al[nk][0] if sign < 0 else al[nk][1]) + (sign * 2 * math.pi if nk == 0 else 0.0)
        mid = (a0 + a1) / 2
        ln = text_at(txt, at((0.0, 0.0), R - LABEL_GAP, mid), mul(outward(mid), -1.0), 0.0, preset, "")
        out.append(ln.ink())
    return out


def build(case: Case) -> tuple[Group, dict]:
    preset = load_preset(case.preset)
    spec = Spec.model_validate({"preset": case.preset, "species": case.species})
    placed = m1_placed(spec, ROOT)
    halves = []
    for p in placed:
        b = p.group.ink_bbox()
        p.translate(-(b[0] + b[2]) / 2, -(b[1] + b[3]) / 2)
        halves.append(((b[2] - b[0]) / 2, (b[3] - b[1]) / 2))
    sign = -1.0 if case.clockwise else 1.0
    s = preset.scheme
    need = max(s.arrow_min_length + (TAIL_ROOM if st.tail_in or st.tail_out else 0.0) for st in case.steps)
    cbox = None
    if case.center:
        cl = text_at(case.center, (0.0, 0.0), (0.0, 0.0), 0.0, preset, "center")
        cbox = cl.ink()
    R, th, al, free = solve(halves, need, CLEAR, sign, [st.label for st in case.steps], cbox, preset)
    c = (0.0, 0.0)
    root = Group(case.name)
    for p, t in zip(placed, th):
        p.translate(*at(c, R, t))
        root.children.append(p.group)
    n = len(placed)
    info = {"R": R, "theta": [math.degrees(t) for t in th], "free_pt": R * free, "arcs": []}
    for k, st in enumerate(case.steps):
        nk = (k + 1) % n
        lead = al[k][1] if sign < 0 else al[k][0]
        trail = al[nk][0] if sign < 0 else al[nk][1]
        a0 = th[k] + sign * lead
        a1 = th[nk] - sign * trail
        if nk == 0:
            a1 += sign * 2 * math.pi
        nm = f"cycle-arc-{k + 1}"
        g = arc_arrow(c, R, a0, a1, preset, nm) if st.kind == "forward" else \
            arc_equilibrium(c, R, a0, a1, preset, nm, case.clockwise)
        root.children.append(g)
        info["arcs"].append(round(R * abs(a1 - a0), 3))
        span = a1 - a0
        rj = R + (preset.equilibrium.shaft_spacing / 2 if st.kind != "forward" else 0.0)  # an outer shaft
        for i, txt in enumerate(st.tail_in):
            tp, F, away = tail(c, rj, a0 + TAIL_JOIN[0] * span, preset, f"cycle-in-{k + 1}-{i + 1}", sign, True)
            ln = text_at(txt, F, away, TAIL_TEXT_GAP, preset, f"cycle-in-{k + 1}")
            root.children.append(Group(tp.source, [tp, *ln.runs]))
        for i, txt in enumerate(st.tail_out):
            tp, F, away = tail(c, rj, a0 + TAIL_JOIN[1] * span, preset, f"cycle-out-{k + 1}-{i + 1}", sign, False)
            ln = text_at(txt, F, away, TAIL_TEXT_GAP , preset, f"cycle-out-{k + 1}")
            root.children.append(Group(tp.source, [tp, *ln.runs]))
        if st.label:
            mid = (a0 + a1) / 2
            ln = text_at(st.label, at(c, R - LABEL_GAP, mid), mul(outward(mid), -1.0), 0.0, preset, f"label-{k + 1}")
            root.children.append(Group(f"cycle-label-{k + 1}", ln.runs))
    if case.center:
        ln = text_at(case.center, c, (0.0, 0.0), 0.0, preset, "center")
        root.children.append(Group("cycle-center", ln.runs))
    return root, info


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=ROOT / "spike" / "cycle")
    ap.add_argument("--case", action="append")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for case in CASES:
        if args.case and case.name not in args.case:
            continue
        root, info = build(case)
        svg, _ = to_svg(root)
        (args.out / f"{case.name}.svg").write_text(svg)
        (args.out / f"{case.name}.png").write_bytes(rasterize(svg, 4.0, background="white"))
        print(f"{case.name}: R {info['R']:.2f} pt, free arc {info['free_pt']:.2f} pt, arcs {info['arcs']}, "
              f"θ {[round(t, 1) for t in info['theta']]}")


if __name__ == "__main__":
    main()
