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

Heads are scene.Markers on the shaft (written as SVG <marker>s, or as separate shapes with
svg_writer.to_svg(affinity=True)); Polyline.head_leaves() gives them in scheme coordinates.
"""

from __future__ import annotations

from .geometry import Point, add, length, mul, sub, unit
from .preset import Preset
from .scene import Group, Marker, Polyline


def head_marker(preset: Preset, start: bool = False) -> Marker:
    """The reaction arrowhead in the marker frame: origin at the shaft's end, the tip
    shaft_end_from_tip ahead of it. `start`: turned to point back along the shaft, for a
    resonance arrow's first end (SVG 1.1 orients a start marker along the path)."""
    a = preset.arrow
    k = -1.0 if start else 1.0

    def at(back: float, side: float) -> Point:
        return (k * (a.shaft_end_from_tip - back), k * side)

    hw = a.head_half_width
    ctl = preset.ratios.arrow_back_control_per_half_width * hw
    tip, notch = at(0.0, 0.0), at(a.head_inset, 0.0)
    c1, c2 = at(a.head_length, hw), at(a.head_length, -hw)
    return Marker(f"{preset.name}-head-reaction{'-start' if start else ''}",
                  (("M", tip), ("L", c1), ("C", c1, at(a.head_inset, ctl), notch), ("C", at(a.head_inset, -ctl), c2, c2)),
                  kind="path", stroke=a.head_outline)


def reaction_arrow(tail: Point, head: Point, preset: Preset, name: str = "arrow") -> Group:
    a = preset.arrow
    u = unit(sub(head, tail))
    if length(sub(head, tail)) <= a.head_length:
        raise ValueError(f"arrow {name} is shorter than its head ({length(sub(head, tail)):.2f} pt)")
    shaft = Polyline([tail, add(head, mul(u, -a.shaft_end_from_tip))], a.shaft_width, source=name,
                     marker_end=head_marker(preset))
    return Group(name, [shaft])


def resonance_arrow(tail: Point, head: Point, preset: Preset, name: str = "resonance-arrow") -> Group:
    """Double-headed arrow, tip to tip from `tail` to `head`."""
    a = preset.arrow
    u = unit(sub(head, tail))
    if length(sub(head, tail)) <= 2 * a.head_length:
        raise ValueError(f"arrow {name} is shorter than its two heads ({length(sub(head, tail)):.2f} pt)")
    shaft = Polyline([add(tail, mul(u, a.shaft_end_from_tip)), add(head, mul(u, -a.shaft_end_from_tip))],
                     a.shaft_width, source=name, marker_start=head_marker(preset, start=True),
                     marker_end=head_marker(preset))
    return Group(name, [shaft])


def notch_x(tail: Point, head: Point, preset: Preset) -> float:
    """x of the head's notch for a horizontal arrow: conditions are centered on tail -> notch."""
    return head[0] - preset.arrow.head_inset * (1 if head[0] >= tail[0] else -1)
