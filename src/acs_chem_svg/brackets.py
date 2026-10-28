"""Square brackets around resonance contributors (acs-style.md §4.1, §4.3; resonance_x*).

ChemDraw strokes each bracket as three segments: a vertical bar and two lips of lip_length
pointing toward the content, stroke_width wide (exempt). Ours is one stroked polyline with
miter joins, so the outer corners are square (ChemDraw's butt-capped segments leave a
stroke_width / 2 notch there). The layout places the bar gap_x from the content ink and the
lips gap_above / gap_below beyond it (stroke centers).
"""

from __future__ import annotations

from .preset import Preset
from .scene import Group, Polyline


def square_bracket(x: float, top: float, bottom: float, side: str, preset: Preset, name: str = "bracket") -> Group:
    """'[' (side "left": lips point right) or ']' (side "right") with its bar's stroke center at
    x and its lips' stroke centers at top and bottom."""
    b = preset.bracket
    lip = x + (b.lip_length if side == "left" else -b.lip_length)
    pts = [(lip, top), (x, top), (x, bottom), (lip, bottom)]
    return Group(name, [Polyline(pts, b.stroke_width, source=name)])


def bracket_ys(box: tuple, preset: Preset) -> tuple[float, float]:
    """The lips' stroke centers (top, bottom) of brackets around a content ink box."""
    return box[1] - preset.bracket.gap_above, box[3] + preset.bracket.gap_below


def bracket_pair(box: tuple, preset: Preset, name: str) -> tuple[Group, Group]:
    """Brackets around a content ink box (x0, y0, x1, y1): `<name>-left` and `<name>-right`."""
    b = preset.bracket
    top, bottom = bracket_ys(box, preset)
    return (square_bracket(box[0] - b.gap_x, top, bottom, "left", preset, f"{name}-left"),
            square_bracket(box[2] + b.gap_x, top, bottom, "right", preset, f"{name}-right"))
