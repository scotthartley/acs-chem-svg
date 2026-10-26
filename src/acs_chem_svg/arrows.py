"""Reaction arrows (acs-style.md §3; exact, measured on every arrow in scheme_basic_x* and
scheme_multi_x*).

- Shaft: a butt-capped stroke of shaft_width from the tail to shaft_end_from_tip before the
  tip, so it runs into the notch.
- Head: tip, two back corners head_length behind the tip at ±head_half_width, and a concave
  back edge: a cubic from each corner to the notch (head_inset behind the tip) with its
  controls at the corner and at (notch, ±arrow_back_control_per_half_width × half-width).
  Filled, with a head_outline hairline stroke as ChemDraw draws it. ChemDraw's head is
  slightly asymmetric (1.55 / 1.60); ours uses the mean half-width on both sides.

Resonance arrows (acs-style.md §4.1; resonance_x*): the same head at both ends, and the shaft
from shaft_end_from_tip after one tip to shaft_end_from_tip before the other.
"""

from __future__ import annotations

from .geometry import Point, add, length, mul, perp, sub, unit
from .preset import Preset
from .scene import Group, Path, Polyline


def _head(tip: Point, u: Point, preset: Preset, name: str) -> Path:
    """The reaction arrowhead with its tip at `tip`, pointing along unit vector u."""
    a = preset.arrow
    n = perp(u)

    def at(back: float, side: float) -> Point:
        return add(tip, add(mul(u, -back), mul(n, side)))

    hw = a.head_half_width
    ctl = preset.ratios.arrow_back_control_per_half_width * hw
    notch = at(a.head_inset, 0.0)
    c1, c2 = at(a.head_length, hw), at(a.head_length, -hw)
    return Path([("M", tip), ("L", c1), ("C", c1, at(a.head_inset, ctl), notch),
                 ("C", at(a.head_inset, -ctl), c2, c2)], fill=True, stroke=a.head_outline, source=name)


def reaction_arrow(tail: Point, head: Point, preset: Preset, name: str = "arrow") -> Group:
    a = preset.arrow
    u = unit(sub(head, tail))
    if length(sub(head, tail)) <= a.head_length:
        raise ValueError(f"arrow {name} is shorter than its head ({length(sub(head, tail)):.2f} pt)")
    shaft = Polyline([tail, add(head, mul(u, -a.shaft_end_from_tip))], a.shaft_width, source=name)
    return Group(name, [shaft, _head(head, u, preset, name)])


def resonance_arrow(tail: Point, head: Point, preset: Preset, name: str = "resonance-arrow") -> Group:
    """Double-headed arrow, tip to tip from `tail` to `head`."""
    a = preset.arrow
    u = unit(sub(head, tail))
    if length(sub(head, tail)) <= 2 * a.head_length:
        raise ValueError(f"arrow {name} is shorter than its two heads ({length(sub(head, tail)):.2f} pt)")
    shaft = Polyline([add(tail, mul(u, a.shaft_end_from_tip)), add(head, mul(u, -a.shaft_end_from_tip))],
                     a.shaft_width, source=name)
    return Group(name, [shaft, _head(tail, mul(u, -1), preset, name), _head(head, u, preset, name)])


def notch_x(tail: Point, head: Point, preset: Preset) -> float:
    """x of the head's notch for a horizontal arrow: conditions are centered on tail -> notch."""
    return head[0] - preset.arrow.head_inset * (1 if head[0] >= tail[0] else -1)
