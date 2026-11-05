"""A species drawn and placed by a layout (layout.py, cycle.py): its group, molecule, Drawer
and the offset the layout added to its atom coordinates; and its compound number's line."""

from __future__ import annotations

from collections.abc import Collection
from dataclasses import dataclass

from . import richtext
from .chem import Molecule
from .draw_mol import Drawer
from .geometry import Point
from .preset import Preset
from .scene import Group
from .spec import Spec


@dataclass
class Placed:
    name: str
    group: Group
    mol: Molecule
    ring_centers: list[Point]
    offset: Point = (0.0, 0.0)  # added to mol atom coordinates by the layout
    drawer: Drawer | None = None  # labels and dots, in mol coordinates

    def translate(self, dx: float, dy: float) -> None:
        self.group.translate(dx, dy)
        self.ring_centers = [(x + dx, y + dy) for x, y in self.ring_centers]
        self.offset = (self.offset[0] + dx, self.offset[1] + dy)


def ring_centers(mol: Molecule) -> list[Point]:
    return [(sum(mol.atoms[i].x for i in r) / len(r), sum(mol.atoms[i].y for i in r) / len(r)) for r in mol.rings]


def drawer(mol: Molecule, preset: Preset, lone_pair_maps: Collection[int],
           tails: dict[int, dict[str, Point]] | None = None) -> Drawer:
    """A Drawer showing the lone pairs of the atoms whose map numbers are in lone_pair_maps."""
    return Drawer(mol, preset, show_lone_pairs={a.idx for a in mol.atoms if a.map in lone_pair_maps}, tails=tails)


def place(name: str, mol: Molecule, preset: Preset, lone_pair_maps: Collection[int],
          tails: dict[int, dict[str, Point]] | None = None) -> Placed:
    """The species drawn (group `name`) at its own coordinates."""
    d = drawer(mol, preset, lone_pair_maps, tails)
    return Placed(name, d.draw(name), mol, ring_centers(mol), drawer=d)


def number_line(spec: Spec, name: str, preset: Preset) -> richtext.Line:
    """The species' compound number (bold, caption size) at the origin; each layout places it."""
    return richtext.layout_line(str(spec.numbers[name]), 0.0, 0.0, preset.text.caption_size, preset,
                                f"number:{name}", style="bold")
