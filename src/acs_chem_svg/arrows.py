"""Reaction arrows (acs-style.md §3.1; exact, measured on every arrow in scheme_basic_x*,
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

Resonance arrows (acs-style.md §4.1; resonance_x*): the same head at both ends, and the shaft
from shaft_end_from_tip after one tip to shaft_end_from_tip before the other.

Equilibrium arrows (acs-style.md §3.6; equilibrium_x*, arrows_large_x1.0): two shafts
equilibrium.shaft_spacing apart, centered on the tail -> head axis. The forward shaft (left of
the direction of travel: above for a left-to-right arrow) points at `head`, the back shaft at
`tail`; each ends in a half head, the outer half of the full head with its tip and notch on
the shaft's inner edge (shaft_width / 2 from its center line) and its barb outward. The back
shaft's half head has its back-edge control back_control_shift nearer the center line.
`tail` and `head` are the arrow's ends: balanced, each shaft's butt sits equilibrium.stagger
behind the other's tip; unbalanced (favoring products or reactants), the long shaft runs end to
end and the short one is centered on it, its butt -> notch run 1 / ratios.equilibrium_long_per_short
of the long one's.

Heads are scene.Markers on the shaft (written as SVG <marker>s, or as separate shapes with
svg_writer.to_svg(affinity=True)); Polyline.head_leaves() gives them in scheme coordinates.
"""

from __future__ import annotations

from .geometry import Point, add, length, mul, perp, sub, unit
from .preset import Preset
from .scene import Group, Marker, Polyline

ARROW_TYPES = ("forward", "equilibrium", "equilibrium_products", "equilibrium_reactants")


def head_marker(preset: Preset, start: bool = False, variant: str | None = None, half: str | None = None) -> Marker:
    """The arrowhead in the marker frame: origin at the shaft's end, the tip shaft_end_from_tip
    ahead of it. `start`: turned to point back along the shaft, for a resonance arrow's first
    end (SVG 1.1 orients a start marker along the path). `half` ("forward" | "back"): an
    equilibrium shaft's half head, barb on the left of the direction of travel."""
    h, a = preset.head(variant), preset.arrow
    k = -1.0 if start else 1.0

    def at(back: float, side: float) -> Point:
        return (k * (h.shaft_end_from_tip - back), k * side)

    hw = h.head_half_width
    ctl = preset.ratios.arrow_back_control_per_half_width * hw
    tip, notch = at(0.0, 0.0), at(h.head_inset, 0.0)
    variant = variant or a.default_head
    mid = f"{preset.name}-head-{variant}"
    if half:
        if half == "back":
            ctl -= preset.equilibrium.back_control_shift
        edge = a.shaft_width / 2  # +y: the right of travel, toward the other shaft
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
    equilibrium arrow (symmetric; unmeasured rule, acs-style.md §3.6)."""
    if kind == "forward":
        return (tail[0] + notch_x(tail, head, preset, variant)) / 2
    return (tail[0] + head[0]) / 2
