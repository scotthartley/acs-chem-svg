"""A spec file as a figure: the one place that decides what kind of spec it is.

- cycle: a catalytic cycle (cycle.build_cycle), when the spec has `cycle:`;
- scheme: an M2 scheme or M3/M4 mechanism (layout.build_scheme), when it has `steps:`;
- species: M1 species side by side (layout.m1_placed) otherwise.

The SVG is the figure's root (svg_writer.to_svg); its CDXML (Figure.cdxml) comes from the same
layout, in the SVG's frame.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from .cdxml_writer import cycle_records, scheme_records, species_records, to_cdxml, warn_stderr
from .cycle import CycleLayout, build_cycle
from .geometry import Point
from .layout import SchemeLayout, build_scheme, m1_placed
from .placed import Placed
from .scene import Group
from .spec import Spec, load_spec
from .svg_writer import to_svg


@dataclass
class Figure:
    spec: Spec
    kind: str  # "cycle" | "scheme" | "species"
    name: str  # the spec file's stem
    root: Group
    scheme: SchemeLayout | None = None
    cycle: CycleLayout | None = None
    placed: list[Placed] | None = None  # M1 species

    def cdxml(self, svg_offset: Point, warn: Callable[[str], None] = warn_stderr) -> str:
        """CDXML for the figure, in the frame of the SVG that to_svg wrote with `svg_offset`."""
        if self.cycle is not None:
            records = cycle_records(self.cycle, self.spec, svg_offset)
        elif self.scheme is not None:
            records = scheme_records(self.scheme, self.spec, svg_offset)
        else:
            records = species_records(self.placed, svg_offset)
        return to_cdxml(records, self.spec.preset, f"{self.name}.cdxml", warn)


def build_figure(path: Path) -> Figure:
    """Lay out the spec at `path` (module docstring)."""
    path = Path(path)
    spec = load_spec(path)
    if spec.cycle is not None:
        cyc = build_cycle(spec, path.parent, path.stem)
        return Figure(spec, "cycle", path.stem, cyc.root, cycle=cyc)
    if spec.steps:
        lay = build_scheme(spec, path.parent, path.stem)
        return Figure(spec, "scheme", path.stem, lay.root, scheme=lay)
    placed = m1_placed(spec, path.parent)
    return Figure(spec, "species", path.stem, Group(path.stem, [p.group for p in placed]), placed=placed)


def render_spec(path: Path) -> Group:
    """Scene for a spec file (build_figure)."""
    return build_figure(path).root


def spec_cdxml(path: Path, warn: Callable[[str], None] = warn_stderr) -> str:
    """What render.py writes as <stem>.cdxml for a spec (build_figure, Figure.cdxml)."""
    f = build_figure(path)
    _, offset = to_svg(f.root)
    return f.cdxml(offset, warn)
