"""Reaction arrows (style-spec.md §3.1; exact, measured on every arrow in scheme_basic_x*,
scheme_multi_x* and arrows_large_x1.0).

- Shaft: a butt-capped stroke of shaft_width from the tail to shaft_end_from_tip before the
  tip, so it runs into the notch.
- Head: tip, two back corners head_length behind the tip at ±head_half_width, and a concave
  back edge: a cubic from each corner to the notch (head_inset behind the tip) with its
  controls at the corner and at (notch, ±arrow_back_control_per_half_width × half-width).
  Filled, with a head_outline hairline stroke as ChemDraw draws it. ChemDraw's head is
  slightly asymmetric (1.55 / 1.60); ours uses the mean half-width on both sides.
- Two head sizes (preset arrow.heads): standard (ChemDraw's default) and large; each arrow
  takes the preset's arrow.default_head unless `variant` names one (Preset.head).

Resonance arrows (style-spec.md §4.1; resonance_x*): the same head at both ends, and the shaft
from shaft_end_from_tip after one tip to shaft_end_from_tip before the other.

Equilibrium arrows (style-spec.md §3.6; equilibrium_x*, arrows_large_x1.0): two shafts
equilibrium.shaft_spacing apart, centered on the tail -> head axis. The forward shaft (left of
the direction of travel: above for a left-to-right arrow) points at `head`, the back shaft at
`tail`; each ends in a half head, the outer half of the full head with its tip and notch on
the shaft's inner edge (shaft_width / 2 from its center line) and its barb outward. The back
shaft's half head has its back-edge control back_control_shift nearer the center line.
`tail` and `head` are the arrow's ends: balanced, each shaft's butt sits equilibrium.stagger
behind the other's tip; unbalanced (favoring products or reactants), the long shaft runs end to
end and the short one is centered on it, its butt -> notch run 1 / ratios.equilibrium_long_per_short
of the long one's.

Arc arrows (style-spec.md §15.1; arcs_x*): a forward arrow along a circle, tail and tip on it.
The notch is on the circle too, head_inset (the chord) back from the tip, and the head is the
straight head on that chord, so it does not follow the circle. The shaft is the arc from the tail
to the notch. Its <marker> is turned in its frame by the angle between the arc's end tangent and
the chord (SVG orients a marker by the tangent), one marker id per turn.
Arc equilibrium arrows (style-spec.md §15.1; arcs_eq_x*): two concentric shafts
equilibrium.shaft_spacing apart, the forward one on the left of travel. Unlike a straight one,
no stagger: each shaft runs from end to end of the arc (its butt on one end, its tip on the
other). Unbalanced: the long shaft end to end, the short one centered, L / ratio long (butt to
tip), lengths along the center circle. Each half head is on its own shaft's chord, its tip and
notch on the shaft's center line (not its inner edge, as straight ones).

Heads are scene.Markers on the shaft (written as SVG <marker>s, or as separate shapes with
svg_writer.to_svg(affinity=True)); Polyline.head_leaves() gives them in scheme coordinates.
"""

from __future__ import annotations

import math

from .geometry import Point, add, arc_point, arc_segments, arc_tangent, dot, length, mul, perp, sub, unit
from .preset import Preset
from .scene import Group, Marker, Path, Polyline

ARROW_TYPES = ("forward", "equilibrium", "equilibrium_products", "equilibrium_reactants")
FAVORS = {"equilibrium": None, "equilibrium_products": "products", "equilibrium_reactants": "reactants"}


def head_marker(preset: Preset, start: bool = False, variant: str | None = None, half: str | None = None,
                arc_turn: float | None = None) -> Marker:
    """The arrowhead in the marker frame: origin at the shaft's end, the tip shaft_end_from_tip
    ahead of it. `start`: turned to point back along the shaft, for a resonance arrow's first
    end (SVG 1.1 orients a start marker along the path). `half` ("forward" | "back"): an
    equilibrium shaft's half head, barb on the left of the direction of travel. `arc_turn`
    (degrees, from the frame's x toward its y): an arc arrow's head, origin at the notch
    (head_inset behind the tip) and turned onto the chord."""
    h, a = preset.head(variant), preset.arrow
    k = -1.0 if start else 1.0
    origin = h.shaft_end_from_tip if arc_turn is None else h.head_inset
    cs, sn = (1.0, 0.0) if arc_turn is None else (math.cos(math.radians(arc_turn)), math.sin(math.radians(arc_turn)))

    def at(back: float, side: float) -> Point:
        x, y = k * (origin - back), k * side
        return (x * cs - y * sn, x * sn + y * cs)

    hw = h.head_half_width
    ctl = preset.ratios.arrow_back_control_per_half_width * hw
    tip, notch = at(0.0, 0.0), at(h.head_inset, 0.0)
    variant = variant or a.default_head
    mid = f"{preset.name}-head-{variant}" + ("" if arc_turn is None else f"-arc{arc_turn:+.3f}")
    if half:
        if half == "back":
            ctl -= preset.equilibrium.back_control_shift
        # +y: the right of travel, toward the other shaft; an arc's tip and notch are on the shaft's
        # center line (arcs_eq_x*: tip offset = shaft offset, ±0.90)
        edge = a.shaft_width / 2 if arc_turn is None else 0.0
        tip, notch, corner = at(0.0, edge), at(h.head_inset, edge), at(h.head_length, -hw)
        return Marker(f"{mid}-half-{half}", (("M", tip), ("L", notch), ("C", at(h.head_inset, -ctl), corner, corner)),
                      kind="path", stroke=a.head_outline)
    c1, c2 = at(h.head_length, hw), at(h.head_length, -hw)
    return Marker(f"{mid}{'-start' if start else ''}",
                  (("M", tip), ("L", c1), ("C", c1, at(h.head_inset, ctl), notch), ("C", at(h.head_inset, -ctl), c2, c2)),
                  kind="path", stroke=a.head_outline)


def reaction_arrow(tail: Point, head: Point, preset: Preset, name: str = "arrow", variant: str | None = None) -> Group:
    h = preset.head(variant)
    u = unit(sub(head, tail))
    if length(sub(head, tail)) <= h.head_length:
        raise ValueError(f"arrow {name} is shorter than its head ({length(sub(head, tail)):.2f} pt)")
    shaft = Polyline([tail, add(head, mul(u, -h.shaft_end_from_tip))], preset.arrow.shaft_width, source=name,
                     marker_end=head_marker(preset, variant=variant))
    return Group(name, [shaft])


def resonance_arrow(tail: Point, head: Point, preset: Preset, name: str = "resonance-arrow",
                    variant: str | None = None) -> Group:
    """Double-headed arrow, tip to tip from `tail` to `head`."""
    h = preset.head(variant)
    u = unit(sub(head, tail))
    if length(sub(head, tail)) <= 2 * h.head_length:
        raise ValueError(f"arrow {name} is shorter than its two heads ({length(sub(head, tail)):.2f} pt)")
    shaft = Polyline([add(tail, mul(u, h.shaft_end_from_tip)), add(head, mul(u, -h.shaft_end_from_tip))],
                     preset.arrow.shaft_width, source=name, marker_start=head_marker(preset, start=True, variant=variant),
                     marker_end=head_marker(preset, variant=variant))
    return Group(name, [shaft])


def equilibrium_arrow(tail: Point, head: Point, preset: Preset, name: str = "equilibrium-arrow",
                      favors: str | None = None, variant: str | None = None) -> Group:
    """Two half-headed shafts from `tail` to `head` (module docstring): balanced, or `favors`
    "products" (the forward shaft long) or "reactants" (the back shaft long). One group, the
    forward shaft first."""
    h, e = preset.head(variant), preset.equilibrium
    L = length(sub(head, tail))
    u = unit(sub(head, tail))
    side = mul(perp(u), e.shaft_spacing / 2)  # the back shaft's side: right of travel
    if favors not in (None, "products", "reactants"):
        raise ValueError(f"arrow {name}: favors must be 'products' or 'reactants', not {favors!r}")
    if favors:
        short = h.head_inset + (L - h.head_inset) / preset.ratios.equilibrium_long_per_short
        if short <= h.head_length:
            raise ValueError(f"arrow {name}: its short shaft is shorter than its head ({short:.2f} pt)")
        lo, hi = (L - short) / 2, (L + short) / 2
        fwd, back = ((0.0, L), (hi, lo)) if favors == "products" else ((lo, hi), (L, 0.0))
    else:
        if L <= e.stagger + h.head_length:
            raise ValueError(f"arrow {name} is shorter than its heads ({L:.2f} pt)")
        fwd, back = (e.stagger, L), (L - e.stagger, 0.0)

    def shaft(butt: float, tip: float, offset: Point, half: str) -> Polyline:
        d = 1.0 if tip > butt else -1.0
        return Polyline([add(add(tail, mul(u, butt)), offset), add(add(tail, mul(u, tip - d * h.shaft_end_from_tip)), offset)],
                        preset.arrow.shaft_width, source=name, marker_end=head_marker(preset, variant=variant, half=half))
    return Group(name, [shaft(*fwd, mul(side, -1.0), "forward"), shaft(*back, side, "back")])


def _arc_head(center: Point, r: float, tip_angle: float, direction: float, preset: Preset, variant: str | None,
              half: str | None = None) -> tuple[float, Marker]:
    """(notch angle, marker) for a head whose tip is at tip_angle on the circle, travelling
    `direction` (+1 counter-clockwise on screen): the notch head_inset back along the chord."""
    inset = preset.head(variant).head_inset
    if inset >= 2 * r:
        raise ValueError(f"arc radius {r:.2f} pt is too small for the head ({inset:.2f} pt inset)")
    notch = tip_angle - direction * 2 * math.asin(inset / (2 * r))
    t = arc_tangent(notch, direction)
    chord = unit(sub(arc_point(center, r, tip_angle), arc_point(center, r, notch)))
    turn = round(math.degrees(math.atan2(dot(chord, perp(t)), dot(chord, t))), 3)
    return notch, head_marker(preset, variant=variant, half=half, arc_turn=turn)


def arc_shaft(center: Point, r: float, a0: float, a1: float, preset: Preset, name: str, variant: str | None,
               half: str | None = None) -> Path:
    direction = 1.0 if a1 > a0 else -1.0
    notch, marker = _arc_head(center, r, a1, direction, preset, variant, half)
    if direction * (notch - a0) <= 0:
        raise ValueError(f"arrow {name} is shorter than its head ({r * abs(a1 - a0):.2f} pt along the arc)")
    return Path(arc_segments(center, r, a0, notch), fill=False, stroke=preset.arrow.shaft_width, closed=False,
                source=name, marker_end=marker)


def arc_arrow(center: Point, r: float, a0: float, a1: float, preset: Preset, name: str = "arc-arrow",
              variant: str | None = None) -> Group:
    """Forward arrow along the circle (center, r) from angle a0 (tail) to a1 (tip), radians,
    counter-clockwise on screen from +x; a1 < a0 runs clockwise (module docstring)."""
    if preset.head(variant).head_length >= r * abs(a1 - a0):
        raise ValueError(f"arrow {name} is shorter than its head ({r * abs(a1 - a0):.2f} pt along the arc)")
    return Group(name, [arc_shaft(center, r, a0, a1, preset, name, variant)])


def arc_equilibrium_shafts(r: float, a0: float, a1: float, preset: Preset, name: str = "arc-equilibrium-arrow",
                           favors: str | None = None, variant: str | None = None) -> list[tuple[float, float, float]]:
    """The two shafts of arc_equilibrium_arrow as (radius, tail angle, tip angle), the forward
    one first (module docstring; arcs_eq_x*)."""
    h, e = preset.head(variant), preset.equilibrium
    direction = 1.0 if a1 > a0 else -1.0
    L = r * abs(a1 - a0)
    if favors not in (None, "products", "reactants"):
        raise ValueError(f"arrow {name}: favors must be 'products' or 'reactants', not {favors!r}")
    if L <= h.head_length:
        raise ValueError(f"arrow {name} is shorter than its heads ({L:.2f} pt)")
    if favors:
        short = L / preset.ratios.equilibrium_long_per_short  # butt to tip (arcs_eq_x*: L / 3)
        if short <= h.head_length:
            raise ValueError(f"arrow {name}: its short shaft is shorter than its head ({short:.2f} pt)")
        lo, hi = (L - short) / 2, (L + short) / 2
        fwd, back = ((0.0, L), (hi, lo)) if favors == "products" else ((lo, hi), (L, 0.0))
    else:  # no stagger on an arc: each shaft runs end to end
        fwd, back = (0.0, L), (L, 0.0)
    left = -direction  # left of travel: inward (-1) for counter-clockwise travel on screen
    r_fwd, r_back = r + left * e.shaft_spacing / 2, r - left * e.shaft_spacing / 2

    def angle(u: float) -> float:
        return a0 + direction * u / r

    return [(r_fwd, angle(fwd[0]), angle(fwd[1])), (r_back, angle(back[0]), angle(back[1]))]


def arc_equilibrium_arrow(center: Point, r: float, a0: float, a1: float, preset: Preset,
                          name: str = "arc-equilibrium-arrow", favors: str | None = None,
                          variant: str | None = None) -> Group:
    """equilibrium_arrow along the circle from a0 to a1 (module docstring; unmeasured): the
    forward shaft on the left of travel (outside for a clockwise arrow), lengths along the
    circle of radius r. One group, the forward shaft first."""
    shafts = arc_equilibrium_shafts(r, a0, a1, preset, name, favors, variant)
    return Group(name, [arc_shaft(center, rr, b0, b1, preset, name, variant, half)
                        for (rr, b0, b1), half in zip(shafts, ("forward", "back"))])


def arc_step_arrow(kind: str, center: Point, r: float, a0: float, a1: float, preset: Preset, name: str,
                   variant: str | None = None) -> Group:
    """A cycle step's arc arrow by its spec arrow_type (ARROW_TYPES)."""
    if kind == "forward":
        return arc_arrow(center, r, a0, a1, preset, name, variant)
    return arc_equilibrium_arrow(center, r, a0, a1, preset, name, FAVORS[kind], variant)


def step_arrow(kind: str, tail: Point, head: Point, preset: Preset, name: str, variant: str | None = None) -> Group:
    """A scheme step's arrow by its spec arrow_type (ARROW_TYPES)."""
    if kind == "forward":
        return reaction_arrow(tail, head, preset, name, variant)
    favors = {"equilibrium": None, "equilibrium_products": "products", "equilibrium_reactants": "reactants"}[kind]
    return equilibrium_arrow(tail, head, preset, name, favors, variant)


def heads_along(kind: str) -> int:
    """Notches a step arrow's conditions must clear: the forward arrow's one, an equilibrium
    arrow's two (one at each end)."""
    return 1 if kind == "forward" else 2


def notch_x(tail: Point, head: Point, preset: Preset, variant: str | None = None) -> float:
    """x of the head's notch for a horizontal arrow."""
    return head[0] - preset.head(variant).head_inset * (1 if head[0] >= tail[0] else -1)


def conditions_mid(tail: Point, head: Point, preset: Preset, kind: str = "forward", variant: str | None = None) -> float:
    """x the conditions are centered on: tail -> notch for a forward arrow, the middle of an
    equilibrium arrow (symmetric; unmeasured rule, style-spec.md §3.6)."""
    if kind == "forward":
        return (tail[0] + notch_x(tail, head, preset, variant)) / 2
    return (tail[0] + head[0]) / 2
