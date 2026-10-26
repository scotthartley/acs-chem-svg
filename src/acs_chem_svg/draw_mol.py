"""ACS molecule renderer: Molecule -> scene Group (design.md §3.3, references/acs-style.md).

Every rule cites its section of acs-style.md and the references it was measured on.
Unmeasured choices say so in their docstrings (acs-style.md §6 Open).

Frames: a bond runs from atom a to atom b. u is the distance along a->b, and v is
the offset along n = perp(u) (y-down frame). Lines are stored as (v, u0, u1).
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

from . import richtext
from .chem import Bond, Molecule
from .fontmetrics import FontFace, load_face
from .geometry import Point, add, angle_deg, dot, length, mul, perp, ray_box_exit, seg_point_dist, sub, unit
from .preset import Preset
from .scene import Circle, Group, Polygon, Polyline, TextRun

COLLINEAR_DEG = 175.0  # a substituent this close to straight-through counts as collinear
# hash count step: 2 hashes at D/hash_spacing = 1.997 (m1_set_x0.8), 3 at 2.05 (m1_set_x1.0)
HASH_COUNT_SHIFT = 0.025


# --------------------------------------------------------------------------- labels (§2)


@dataclass
class Label:
    atom: int
    element: str
    num_h: int
    h_side: str  # "right" | "left"; also the justification (left = bonds point right)
    runs: list[TextRun]
    elem_origin: Point  # glyph origin of the element symbol (x, baseline)
    size: float
    face: FontFace
    outline: list[list[Point]] = field(default_factory=list)  # element glyph contours, relative to the atom

    def elem_ink_box(self):
        b = self.face.text_ink_bounds(self.element, self.size)
        x, base = self.elem_origin
        return (x + b[0], base - b[3], x + b[2], base - b[1])

    def shift(self, d: Point) -> None:
        """Move the whole label (runs, origin and outline) by d."""
        for run in self.runs:
            run.x, run.y = run.x + d[0], run.y + d[1]
        self.elem_origin = add(self.elem_origin, d)
        self.outline = [[add(p, d) for p in c] for c in self.outline]


def is_labeled(mol: Molecule, i: int) -> bool:
    """Heteroatoms and explicit H atoms get labels. Carbon is labeled when isolated or asked
    for (Atom.show_label: CH3 in stereo_x*). ChemDraw's ShowTerminalCarbonLabels is off, and
    charged carbons stay unlabeled with the charge on the bisector (resonance_x1.0)."""
    a = mol.atoms[i]
    return a.element != "C" or a.show_label or a.abbrev is not None or mol.degree(i) == 0


def h_side(mol: Molecule, i: int) -> str:
    """H goes opposite the mean bond direction in x; right on a tie (§2, H₂N vs NH₂)."""
    a = mol.atoms[i]
    sx = sum(unit(sub(mol.atoms[j].p, a.p))[0] for j, _ in mol.neighbors(i))
    return "left" if sx > 1e-6 else "right"


def element_x(face: FontFace, element: str, size: float, side: str) -> float:
    """Glyph origin of the element symbol relative to the atom x (§2): the first letter is
    centered on the atom, or the last letter when the label is right-justified (bonds point
    right, H on the left): Cl -> C centered, but -Cl -> l centered (clearance_x*)."""
    if side == "left":
        return -face.advance(element[:-1], size) - face.advance(element[-1], size) / 2
    return -face.advance(element[0], size) / 2


def element_outline(face: FontFace, element: str, size: float, x0: float, baseline: float) -> list[list[Point]]:
    """Glyph contours of the element symbol, y down, with its origin at (x0, baseline)."""
    out, x = [], x0
    for ch in element:
        out += [[(x + px, baseline - py) for px, py in c] for c in face.outline(ch, size)]
        x += face.advance(ch, size)
    return out


def layout_label(mol: Molecule, i: int, preset: Preset, face: FontFace, side: str | None = None) -> Label:
    """Text runs for one atom label (§2). The anchor letter of the element (element_x) is
    centered on the atom; baseline = atom y + 0.3585 x size; H count as a subscript (0.75 x
    size, dropped 0.23 x size). Left-hand forms are set right to left. Runs split at size changes."""
    atom = mol.atoms[i]
    r = preset.ratios
    size = preset.text.label_size
    sub_size = r.subscript_scale * size
    base = atom.y + r.label_baseline_per_size * size
    sub_base = base + r.subscript_drop_per_size * size
    side = side or h_side(mol, i)
    src = f"atom:{atom.id}"
    if atom.abbrev is not None:
        return layout_abbreviation(mol, i, preset, face, side, base)
    ex = atom.x + element_x(face, atom.element, size, side)
    pieces = [(atom.element, size, base)]  # (text, size, baseline), in reading order
    if atom.num_h:
        hs = [("H", size, base)] + ([(str(atom.num_h), sub_size, sub_base)] if atom.num_h > 1 else [])
        pieces = pieces + hs if side == "right" else hs + pieces
    # glyph origins: right-hand pieces follow the element; left-hand ones precede it
    xs, x = [], ex
    k = next(n for n, p in enumerate(pieces) if p[0] == atom.element)
    for text, sz, _ in pieces[k:]:
        xs.append(x)
        x += face.advance(text, sz)
    x = ex
    left = []
    for text, sz, _ in reversed(pieces[:k]):
        x -= face.advance(text, sz)
        left.append(x)
    xs = list(reversed(left)) + xs
    runs: list[TextRun] = []
    for (text, sz, b), x0 in zip(pieces, xs):
        prev = runs[-1] if runs else None
        if (prev and prev.size == sz and prev.y == b and abs(prev.x + prev.advance() - x0) < 1e-9
                and face.kerning(prev.text[-1], text[0], sz) == 0):
            prev.text += text  # same size, contiguous, unkerned: one run, as ChemDraw writes it
        else:
            runs.append(TextRun(x0, b, sz, text, face, source=src))
    outline = element_outline(face, atom.element, size, ex - atom.x, base - atom.y)
    return Label(i, atom.element, atom.num_h, side, runs, (ex, base), size, face, outline)


def anchor_symbol(text: str, element: str) -> str:
    """The glyphs of an abbreviation that stand for its attachment atom: the element symbol when
    the label starts with it (OMe: O, CO2Me: C), else the first character (Ph, Bn)."""
    return element if text.startswith(element) else text[:1]


def layout_abbreviation(mol: Molecule, i: int, preset: Preset, face: FontFace, side: str, base: float) -> Label:
    """Abbreviation label (OMe, OBn: scheme_multi_x*): set left-justified as one rich-text run
    with the first letter centered on the atom, on the label baseline. The bond clips against
    the anchor glyphs (anchor_symbol) by the §2.1 band model (OMe/OBn bond ends within 0.08 pt
    of ChemDraw at all presets). The right-justified form (MeO–, bonds leaving to the right)
    is unmeasured (acs-style.md §6 Open)."""
    atom = mol.atoms[i]
    if side == "left":
        raise NotImplementedError(f"atom {atom.id}: {atom.abbrev} with its bond on the right needs the "
                                  "right-justified form (MeO–), which is unmeasured (acs-style.md §6 Open); "
                                  "turn the structure (orientation) so the bond comes from the left or below")
    size = preset.text.label_size
    parts = richtext.parse(atom.abbrev)
    anchor = anchor_symbol("".join(p.text for p in parts), atom.element)
    first = richtext.face_for(parts[0].style)
    ex = atom.x - first.advance(parts[0].text[0], size) / 2
    line = richtext.layout_line(parts, ex, base, size, preset, source=f"atom:{atom.id}")
    outline = element_outline(first, anchor, size, ex - atom.x, base - atom.y)
    return Label(i, anchor, 0, "right", line.runs, (ex, base), size, first, outline)


# --------------------------------------------------------------------------- clearance (§2.1)


def band_extent(contours: list[list[Point]], direction: Point, band: float) -> float:
    """Furthest reach along `direction` (unit) of the outline parts within ±band of the axis
    through the origin. Edges are clipped to the band, so long straight sides count."""
    dx, dy = direction
    best = 0.0
    for c in contours:
        for p, q in zip(c, c[1:] + c[:1]):
            u0, v0 = p[0] * dx + p[1] * dy, p[0] * dy - p[1] * dx
            u1, v1 = q[0] * dx + q[1] * dy, q[0] * dy - q[1] * dx
            lo, hi, dv = 0.0, 1.0, v1 - v0
            if abs(dv) < 1e-12:
                if abs(v0) > band:
                    continue
            else:
                ta, tb = (band - v0) / dv, (-band - v0) / dv
                lo, hi = max(lo, min(ta, tb)), min(hi, max(ta, tb))
            if lo <= hi:
                best = max(best, u0 + lo * (u1 - u0), u0 + hi * (u1 - u0))
    return best


def clip_distance_band(contours: list[list[Point]], direction: Point, band: float, margin: float,
                       per_margin: float, offset: Point = (0.0, 0.0)) -> float:
    """§2.1: a bond stops per_margin x margin beyond the furthest point, along the bond, of the
    element outline lying within ±band of the bond axis. The outline is shifted by `offset`
    (fitted). `direction` is the unit vector out of the label toward the bonded atom."""
    shifted = [[(x + offset[0], y + offset[1]) for x, y in c] for c in contours]
    return band_extent(shifted, direction, band) + per_margin * margin


def clip_distance(label: Label, direction: Point, preset: Preset) -> float:
    """Where a bond leaving `label` along `direction` must stop, from the atom center."""
    L, r = preset.lengths, preset.ratios
    return clip_distance_band(label.outline, direction, r.label_band_per_size * label.size, L.margin_width,
                              r.label_clearance_per_margin, (L.label_outline_offset_x, L.label_outline_offset_y))


# --------------------------------------------------------------------------- bond geometry (§1)


@dataclass
class Frame:
    o: Point
    u: Point
    n: Point
    L: float

    def at(self, u: float, v: float) -> Point:
        return add(self.o, add(mul(self.u, u), mul(self.n, v)))


@dataclass
class BondGeometry:
    bond: int
    frame: Frame
    lines: list[tuple[float, float, float]] = field(default_factory=list)  # (v, u0, u1)
    main: tuple[float, float] | None = None  # (u0, u1) of the axis line, joinable at unlabeled atoms
    placement: str = ""  # double bonds: ring | centered | offset
    polygon: list[Point] | None = None  # wedge
    hashes: list[tuple[float, float]] = field(default_factory=list)  # (u center, half length)


# --------------------------------------------------------------------------- dots (§4, M3)

SIDES = {"top": (0.0, -1.0), "right": (1.0, 0.0), "bottom": (0.0, 1.0), "left": (-1.0, 0.0)}
SIDE_BLOCK_DEG = 45.0  # a bond (or the H run) within this angle of a side's direction takes the side
CHARGE_CORNER_DEG = 60.0  # a bond within this angle of a corner's diagonal moves a labeled charge on: upper right, upper left, lower right, lower left
CHARGE_CORNER_MIN_DEG = 45.0  # ... unless every corner has one: then the first corner with none within this (m1_set_x*: N⁺ with bonds at 0° and 180° keeps the upper right), else the largest gap


@dataclass
class DotSet:
    """A lone pair (two dots) or a radical (one dot) on atom `atom`: `center` is the pair
    midpoint or the dot, `direction` points out from the atom, `side` names it (top / right /
    bottom / left; the nearest one on unlabeled atoms)."""
    atom: int
    kind: str  # "pair" | "radical"
    center: Point
    direction: Point
    side: str
    dots: list[Point]


def _side_name(d: Point) -> str:
    return max(SIDES, key=lambda k: dot(SIDES[k], d))


class Drawer:
    def __init__(self, mol: Molecule, preset: Preset, face: FontFace | None = None,
                 show_lone_pairs: set[int] | None = None,
                 tails: dict[int, dict[str, Point]] | None = None):
        self.mol, self.preset = mol, preset
        self.face = face or load_face("Regular")
        self.labels = {i: layout_label(mol, i, preset, self.face)
                       for i in range(len(mol.atoms)) if is_labeled(mol, i)}
        self.center_terminal_labels()
        self.show_lone_pairs = set(show_lone_pairs or ())  # atom indices whose lone pairs are drawn
        # atom index -> {"pair" | "radical": direction toward its arrow's target}: that set faces it
        self.tails = tails or {}
        self._dots: dict[int, list[DotSet]] | None = None

    def center_terminal_labels(self) -> None:
        """A label on a degree-1 atom at the end of an offset double bond sits centered
        between the two lines: moved bond_spacing/2 toward the second line (double_bonds_2,
        acetaldehyde O, 1.30 pt at both orientations). The bond lines don't move."""
        for b in self.mol.bonds:
            if b.order != 2 or b.ring or b.stereo:
                continue
            ends = [i for i in (b.a, b.b) if i in self.labels and self.mol.degree(i) == 1]
            if not ends:
                continue
            fr = self.frame(b)
            kind, side = self.double_placement(b, fr)
            if kind == "offset":
                for i in ends:
                    self.labels[i].shift(mul(fr.n, side * self.preset.lengths.bond_spacing / 2))

    # ---------------------------------------------------------------- helpers

    def frame(self, b: Bond) -> Frame:
        pa, pb = self.mol.atoms[b.a].p, self.mol.atoms[b.b].p
        u = unit(sub(pb, pa))
        return Frame(pa, u, perp(u), length(sub(pb, pa)))

    def end_clip(self, b: Bond, i: int) -> float:
        """Distance from atom i to where bond b's lines stop (0 at unlabeled atoms)."""
        if i not in self.labels:
            return 0.0
        d = unit(sub(self.mol.atoms[b.other(i)].p, self.mol.atoms[i].p))
        return clip_distance(self.labels[i], d, self.preset)

    def substituents(self, b: Bond, i: int, fr: Frame, plain_only: bool = False) -> list[tuple[int, float]]:
        """(side, angle φ) of atom i's other bonds: side is the sign of v in fr; φ is the angle
        between the bond (from i toward the other end) and the substituent. `plain_only`
        skips wedge and hash bonds (they draw nothing at the atom center)."""
        pi = self.mol.atoms[i].p
        into = unit(sub(self.mol.atoms[b.other(i)].p, pi))
        out = []
        for j, nb in self.mol.neighbors(i):
            if j == b.other(i) or (plain_only and nb.stereo):
                continue
            w = unit(sub(self.mol.atoms[j].p, pi))
            phi = angle_deg(into, w)
            side = 0 if phi >= COLLINEAR_DEG or phi < 1e-6 else (1 if dot(w, fr.n) > 0 else -1)
            out.append((side, phi))
        return out

    def nearest_on_side(self, b: Bond, i: int, fr: Frame, side: int, plain_only: bool = False) -> float | None:
        phis = [phi for s, phi in self.substituents(b, i, fr, plain_only) if s == side]
        return min(phis) if phis else None

    def is_open(self, b: Bond, i: int, fr: Frame) -> bool:
        """Terminal, or with only collinear substituents."""
        return all(s == 0 for s, _ in self.substituents(b, i, fr))

    def ring_doubles(self, ring: list[int]) -> int:
        members = set(ring)
        return sum(1 for b in self.mol.bonds if b.order == 2 and b.a in members and b.b in members)

    def ring_side(self, b: Bond, fr: Frame) -> int:
        """Side of the ring the second line goes in. A bond shared by two rings uses the ring
        with more double bonds, then the larger ring (indene: the benzene ring), then the
        v < 0 side (naphthalene) (double_bonds_2)."""
        mid = fr.at(fr.L / 2, 0)

        def side(ring):
            c = (sum(self.mol.atoms[i].x for i in ring) / len(ring), sum(self.mol.atoms[i].y for i in ring) / len(ring))
            return 1 if dot(sub(c, mid), fr.n) > 0 else -1
        rings = [r for r in self.mol.rings if b.a in r and b.b in r]
        best = max(rings, key=lambda r: (self.ring_doubles(r), len(r), -side(r)))
        return side(best)

    def double_placement(self, b: Bond, fr: Frame) -> tuple[str, int]:
        """§1.1. In a ring: second line inside. Acyclic: centered when the terminal (or
        collinear) end faces a degree-3 atom or a label, or both ends are open; otherwise
        offset toward the substituents."""
        if b.ring:
            return "ring", self.ring_side(b, fr)
        open_a, open_b = self.is_open(b, b.a, fr), self.is_open(b, b.b, fr)
        if open_a and open_b:
            return "centered", 0
        if open_a or open_b:
            facing = b.b if open_a else b.a
            n_sub = sum(1 for s, _ in self.substituents(b, facing, fr) if s != 0)
            if n_sub >= 2 or facing in self.labels:
                return "centered", 0
        count = sum(s for i in (b.a, b.b) for s, _ in self.substituents(b, i, fr))
        if count:
            return "offset", 1 if count > 0 else -1
        # tie (E-alkenes): the v < 0 side of a->b, at every rotation (double_bonds_2, four
        # orientations; double_bonds_x1.0). It follows the bond's drawn direction, not the screen.
        return "offset", -1

    def inner_inset(self, phi: float) -> float:
        """Inset of an offset line at an end with a substituent at angle φ on its side: the
        line stops on the angle bisector: ring_inner_inset at 120°, scaled by cot(φ/2)
        (measured at 108° too: 1.88 in indene's five-membered ring, double_bonds_2)."""
        return self.preset.lengths.ring_inner_inset * (1 / math.tan(math.radians(phi / 2))) / (1 / math.tan(math.radians(60)))

    # ---------------------------------------------------------------- per bond

    def bond_geometry(self, b: Bond) -> BondGeometry:
        L, r = self.preset.lengths, self.preset.ratios
        fr = self.frame(b)
        da, db = self.end_clip(b, b.a), self.end_clip(b, b.b)
        u0, u1 = da, fr.L - db
        g = BondGeometry(b.idx, fr)
        s = L.bond_spacing
        if b.stereo == "wedge":
            g.polygon = self.wedge_polygon(b, fr, u0, u1)
            return g
        if b.stereo == "hash":
            g.hashes = self.hash_lines(u0, u1)
            return g
        if b.order == 1:
            g.main = (u0 - self.behind_wedge(b, b.a), u1 + self.behind_wedge(b, b.b))
        elif b.order == 3:  # three lines at -s, 0, +s, all full length
            g.main = (u0, u1)
            g.lines = [(-s, u0, u1), (s, u0, u1)]
        elif b.order == 2:
            kind, side = self.double_placement(b, fr)
            g.placement = kind
            if kind == "centered":
                h = s / 2
                for sd in (1, -1):
                    ends = []
                    for i, d in ((b.a, da), (b.b, db)):
                        if i in self.labels:
                            ends.append(d)
                            continue
                        # runs on until it meets the neighbor bond's centerline on its side
                        phi = self.nearest_on_side(b, i, fr, sd)
                        ends.append(h / math.tan(math.radians(phi)) if phi is not None else 0.0)
                    g.lines.append((sd * h, ends[0], fr.L - ends[1]))
            else:
                g.main = (u0, u1)
                ends = []
                for i, d in ((b.a, da), (b.b, db)):
                    if i in self.labels:
                        ends.append(d)
                        continue
                    phi = self.nearest_on_side(b, i, fr, side)
                    ends.append(self.inner_inset(phi) if phi is not None else 0.0)
                g.lines.append((side * s, ends[0], fr.L - ends[1]))
        else:
            raise ValueError(f"bond order {b.order} not supported")
        return g

    def wedge_polygon(self, b: Bond, fr: Frame, u0: float, u1: float) -> list[Point]:
        """§1.2 wedge: sides run linearly from half-width line_width/2 at the stereocenter
        to wedge_wide_end/2 at the clipped end or the wide atom (measured, all presets).
        - Narrow end at an unlabeled stereocenter: cut along the bisector with the neighbor
          bond on each side (wedges and hashes ignored), meeting at the atom center
          (alanine: cut toward CH3, not toward the hashed NH2).
        - Wide end at a labeled or degree-1 atom: flat (stereo_x*).
        - Wide end at an unlabeled atom with other bonds: wedge_wide_end()."""
        L = self.preset.lengths
        h0, h1 = L.line_width / 2, L.wedge_wide_end / 2
        k = (h1 - h0) / (u1 - u0)
        sides = []
        for sd in (1, -1):
            phi = self.nearest_on_side(b, b.a, fr, sd, plain_only=True) if b.a not in self.labels else None
            t = None
            if phi is not None:
                c, s = math.cos(math.radians(phi / 2)), math.sin(math.radians(phi / 2))
                if s - k * c > 1e-9:
                    t = h0 / (s - k * c)  # bisector ray meets the side line
            u = u0 + (t * c if t is not None else 0.0)
            sides.append(fr.at(u, sd * (h0 + k * (u - u0))))
        wide = self.wedge_wide_end(b, fr, u0, h0, k) if b.b not in self.labels else None
        if wide is None:
            wide = [fr.at(u1, h1), fr.at(u1, -h1)]
        pts = [sides[0], *wide, sides[1]]
        return ([fr.at(u0, 0)] if b.a not in self.labels else []) + pts

    def wedge_wide_end(self, b: Bond, fr: Frame, u0: float, h0: float, k: float) -> list[Point] | None:
        """Wide end at an unlabeled atom with other bonds (stereo_x*, within 0.01 pt). The
        sides run past the atom until they meet the nearer edge of the neighbor bond's
        stroke on their side; the end is then cut back to the atom center (a V notch). A side
        with no neighbor stops on the line through the other corner and the center."""
        lw = self.preset.lengths.line_width
        nbrs = [j for j, nb in self.mol.neighbors(b.b) if j != b.a and not nb.stereo]
        if not nbrs:
            return None
        pb = self.mol.atoms[b.b].p
        corners = {}
        for sd in (1, -1):
            best = None
            for j in nbrs:
                w = unit(sub(self.mol.atoms[j].p, pb))
                wl = (dot(w, fr.u), dot(w, fr.n))
                if (wl[1] > 0) != (sd > 0) or abs(wl[1]) < 1e-9:
                    continue
                for off in (lw / 2, -lw / 2):
                    q = (fr.L - wl[1] * off, wl[0] * off)  # point on the edge line, local frame
                    # side line (u, sd*(h0 + k(u - u0))) meets q + t*wl:
                    #   u - t*wl0 = q0,  sd*k*u - t*wl1 = q1 - sd*(h0 - k*u0)
                    a11, a12, a21, a22 = 1.0, -wl[0], sd * k, -wl[1]
                    det = a11 * a22 - a12 * a21
                    if abs(det) < 1e-12:
                        continue
                    r1, r2 = q[0], q[1] - sd * (h0 - k * u0)
                    u = (r1 * a22 - a12 * r2) / det
                    pt = (u, sd * (h0 + k * (u - u0)))
                    dist = math.dist(pt, (fr.L, 0.0))
                    if best is None or dist < best[0]:
                        best = (dist, pt)
            if best:
                corners[sd] = best[1]
        if not corners:
            return None
        for sd in (1, -1):
            if sd not in corners:  # continue the cut from the other corner through the center
                cu, cv = corners[-sd]
                du, dv = fr.L - cu, -cv
                # side line v = sd*(h0 + k(u-u0)) meets (cu, cv) + t*(du, dv)
                t = (sd * (h0 + k * (cu - u0)) - cv) / (dv - sd * k * du)
                corners[sd] = (cu + t * du, cv + t * dv)
        return [fr.at(*corners[1]), fr.at(fr.L, 0.0), fr.at(*corners[-1])]

    def behind_wedge(self, b: Bond, i: int) -> float:
        """How far plain bond b's line starts behind unlabeled atom i when it is the only plain
        bond there and a wedge's wide end sits on i: back to the wedge corner on its far edge,
        so it fills the notch (stereo_x1.0 butane: ChemDraw's bond polygon starts 1.51 pt behind
        the atom, at the wedge end). Where a chain passes through (rings), the miter fills it."""
        if i in self.labels:
            return 0.0
        nbs = self.mol.neighbors(i)
        if any(nb.idx != b.idx and not nb.stereo for _, nb in nbs):
            return 0.0
        pi = self.mol.atoms[i].p
        out = unit(sub(self.mol.atoms[b.other(i)].p, pi))
        back = 0.0
        for _, wb in nbs:
            if wb.stereo == "wedge" and wb.b == i:
                fr = self.frame(wb)
                poly = self.wedge_polygon(wb, fr, self.end_clip(wb, wb.a), fr.L)
                wide = [p for p in poly if math.dist(p, pi) < fr.L / 2]  # corners at i, not the narrow end
                back = max(back, *(-dot(sub(p, pi), out) for p in wide))
        return back

    def hash_lines(self, u0: float, u1: float) -> list[tuple[float, float]]:
        """§1.2 hashed wedge. The outline's half-width runs linearly from line_width/2 at
        u0 + hash_spacing to wedge_wide_end/2 at u1 (measured, all presets). The first
        hash starts at u0 + hash_spacing and the last ends flush with u1. Count:
        floor(D / hash_spacing - HASH_COUNT_SHIFT) + 1 evenly spaced hashes, with
        D = u1 - lw - (u0 + hash_spacing) (27 of 27 hashes in m1_set_x*, stereo_x*, lengths
        0.5-2.2 x bond length). The data bracket the step from 2 to 3 hashes between
        D / hash_spacing = 1.997 and 2.05, so the cut sits mid-bracket rather than at 2.0.
        At an unlabeled atom with other bonds the last hash still ends at the atom (unmeasured)."""
        L = self.preset.lengths
        lw, hs = L.hash_line_width, L.hash_spacing
        first, last = u0 + hs, u1 - lw
        span = last - first
        if span <= 0:
            nears = [max(u0, last)]
        else:
            n = math.floor(span / hs - HASH_COUNT_SHIFT) + 1
            nears = [first + span * k / (n - 1) for k in range(n)] if n > 1 else [last]
        hw0, hw1 = lw / 2, L.wedge_wide_end / 2

        def half(u):
            return hw1 if span <= 0 else hw0 + (hw1 - hw0) * (u - first) / (u1 - first)
        return [(near + lw / 2, half(near + lw)) for near in nears]

    # ---------------------------------------------------------------- chains

    def chains(self, geoms: dict[int, BondGeometry]) -> list[tuple[list[Point], bool, list[int]]]:
        """Join main lines into polylines through unlabeled atoms (miter joins). At each
        atom, ends pair up most-collinear first; an odd end stops at the atom center."""
        ends: dict[int, list[tuple[int, int]]] = {}  # atom -> [(bond idx, which end 0/1)]
        for bi, g in geoms.items():
            if g.main is None:
                continue
            b = self.mol.bonds[bi]
            for which, i in ((0, b.a), (1, b.b)):
                if i not in self.labels:
                    ends.setdefault(i, []).append((bi, which))
        partner: dict[tuple[int, int], tuple[int, int]] = {}
        for i, lst in ends.items():
            pi = self.mol.atoms[i].p
            dirs = {e: unit(sub(self.mol.atoms[self.mol.bonds[e[0]].other(i)].p, pi)) for e in lst}
            free = list(lst)
            while len(free) >= 2:
                e1, e2 = max(((x, y) for k, x in enumerate(free) for y in free[k + 1:]),
                             key=lambda xy: angle_deg(dirs[xy[0]], dirs[xy[1]]))
                partner[e1], partner[e2] = e2, e1
                free.remove(e1)
                free.remove(e2)

        def pts(bi, start_end):
            g = geoms[bi]
            p0, p1 = g.frame.at(g.main[0], 0), g.frame.at(g.main[1], 0)
            return (p0, p1) if start_end == 0 else (p1, p0)

        seen, out = set(), []
        order = sorted(geoms)
        starts = [(bi, w) for bi in order if geoms[bi].main for w in (0, 1) if (bi, w) not in partner]
        for bi, w in starts + [(bi, 0) for bi in order if geoms[bi].main]:
            if bi in seen:
                continue
            cur, which = bi, w
            p0, p1 = pts(cur, which)
            line, bonds, closed = [p0, p1], [cur], False
            seen.add(cur)
            while True:
                nxt = partner.get((cur, 1 - which))
                if nxt is None:
                    break
                if nxt[0] in seen:
                    closed = nxt[0] == bi
                    break
                cur, which = nxt
                seen.add(cur)
                bonds.append(cur)
                line.append(pts(cur, which)[1])
            if closed:
                line.pop()
            out.append((line, closed, bonds))
        return out

    # ---------------------------------------------------------------- charges (§2.2)

    def charge_center(self, i: int) -> Point:
        a = self.mol.atoms[i]
        c = self.preset.charge
        if i in self.labels:
            lab = self.labels[i]
            if lab.num_h and lab.h_side == "right":  # NH₃⁺: centered above the element
                return (a.x, a.y - c.above_label)
            x0, y0, x1, y1 = lab.elem_ink_box()
            bonds = [unit(sub(self.mol.atoms[j].p, a.p)) for j, _ in self.mol.neighbors(i)]

            def spare(corner: Point) -> float:  # angle to the nearest bond
                return min((angle_deg(corner, b) for b in bonds), default=180.0)
            corners = [((1.0, -1.0), (x1 + c.labeled_offset_x, y0 + c.labeled_offset_y)),  # upper right
                       ((-1.0, -1.0), (x0 - c.labeled_offset_x, y0 + c.labeled_offset_y)),  # HO⁻ drawn O-H: upper left
                       ((1.0, 1.0), (x1 + c.labeled_offset_x, y1 - c.labeled_offset_y)),  # nitro N⁺: a lower corner
                       ((-1.0, 1.0), (x0 - c.labeled_offset_x, y1 - c.labeled_offset_y))]
            for limit in (CHARGE_CORNER_DEG, CHARGE_CORNER_MIN_DEG):
                free = [p for d, p in corners if spare(unit(d)) >= limit - 1e-9]
                if free:
                    return free[0]
            # a bond close to every corner (nitro N⁺ in a ring compound): on the largest gap's
            # bisector, as far beyond the element ink box as the corner offset (unmeasured rule,
            # acs-style.md §2.2)
            ctr = ((x0 + x1) / 2, (y0 + y1) / 2)
            d = self._largest_gap_direction(i)
            dist = ray_box_exit(ctr, d, (x0, y0, x1, y1)) + math.hypot(c.labeled_offset_x, c.labeled_offset_y)
            return add(ctr, mul(d, dist))
        d = self._largest_gap_direction(i)
        return (a.x + c.unlabeled_distance * d[0], a.y + c.unlabeled_distance * d[1])

    def _largest_gap_direction(self, i: int) -> Point:
        """Unit vector along the bisector of the largest angular gap between atom i's bonds
        (+x for an atom without bonds)."""
        a = self.mol.atoms[i]
        angs = sorted(math.atan2(self.mol.atoms[j].y - a.y, self.mol.atoms[j].x - a.x) for j, _ in self.mol.neighbors(i))
        if not angs:
            return (1.0, 0.0)
        gaps = [((angs[(k + 1) % len(angs)] - angs[k]) % (2 * math.pi) or 2 * math.pi, angs[k]) for k in range(len(angs))]
        gap, start = max(gaps)
        t = start + gap / 2
        return (math.cos(t), math.sin(t))

    def charge_bars(self, i: int) -> list[Polyline]:
        a = self.mol.atoms[i]
        if abs(a.charge) != 1:
            raise NotImplementedError(f"atom {a.id}: only charges of ±1 are drawn in M1 (got {a.charge})")
        c = self.preset.charge
        cx, cy = self.charge_center(i)
        h = c.bar_length / 2
        src = f"atom:{a.id}"
        bars = [Polyline([(cx - h, cy), (cx + h, cy)], c.bar_thickness, source=src)]
        if a.charge > 0:
            bars.append(Polyline([(cx, cy - h), (cx, cy + h)], c.bar_thickness, source=src))
        return bars

    # ---------------------------------------------------------------- lone pairs and radicals (§4)

    def dot_sets(self) -> dict[int, list[DotSet]]:
        """atom index -> its drawn lone pairs and radicals. Radicals are always drawn; lone
        pairs on the atoms in show_lone_pairs (all of the atom's pairs)."""
        if self._dots is None:
            self._dots = {}
            for i, a in enumerate(self.mol.atoms):
                pairs = a.lone_pairs if i in self.show_lone_pairs else 0
                if pairs or a.radicals:
                    sets = (self.label_dots(i, pairs, a.radicals) if i in self.labels
                            else self.bisector_dots(i, pairs, a.radicals))
                    self._dots[i] = sets
        return self._dots

    def label_dots(self, i: int, pairs: int, radicals: int) -> list[DotSet]:
        """Labeled atom (§4.2, mech_basic_x*: O⁻, Br⁻, Cl·): each lone pair takes a side of
        the element ink box (top / right / bottom / left), its dots lone_pair.gap from the box
        and dot_spacing apart along it, centered on the box. A side is blocked by a bond or the
        H run within 45° of it, and by the charge when it sits above the label (NH₃⁺). A set
        an arrow starts from (`tails`) takes the free side facing the arrow's target. Arrowheads
        never move dots: the arrows work around them (mechanism.Planner). The other
        pairs take the free sides that keep the most clearance from bonds, then spread out,
        then top, right, bottom, left (O⁻ with its bond on the left: top, right, bottom, as in
        ChemDraw). Radicals take the next free sides (Cl·: left), a single dot at the side
        center."""
        a, lab = self.mol.atoms[i], self.labels[i]
        lp = self.preset.lone_pair
        blocked = [unit(sub(self.mol.atoms[j].p, a.p)) for j, _ in self.mol.neighbors(i)]
        if lab.num_h:
            blocked.append(SIDES[lab.h_side])
            if a.charge and lab.h_side == "right":
                blocked.append(SIDES["top"])  # NH₃⁺: the sign sits above the element
        clearance = {k: min((angle_deg(d, b) for b in blocked), default=180.0) for k, d in SIDES.items()}
        free = [k for k in SIDES if clearance[k] > SIDE_BLOCK_DEG]
        order = list(SIDES)

        def choose(n: int, avail: list[str]) -> list[str]:
            if n <= 0:
                return []
            if len(avail) < n:  # not enough free sides: the least blocked others too (checks report clashes)
                avail = avail + sorted((k for k in order if k not in avail), key=lambda k: -clearance[k])[:n - len(avail)]
            best = max(itertools.combinations(avail, n), key=lambda ks: (
                round(min(clearance[k] for k in ks)),
                min((angle_deg(SIDES[x], SIDES[y]) for x, y in itertools.combinations(ks, 2)), default=0),
                [-order.index(k) for k in sorted(ks, key=order.index)]))
            return sorted(best, key=order.index)
        # an arrow's tail set first: the free side facing its target (Br· abstracting an H to its right)
        first: dict[str, str] = {}
        for kind, n in (("radical", radicals), ("pair", pairs)):
            d = self.tails.get(i, {}).get(kind)
            avail = [k for k in free if k not in first.values()]
            if d is not None and n and avail:
                first[kind] = max(avail, key=lambda k: dot(SIDES[k], d))
        rest = [k for k in free if k not in first.values()]
        pair_sides = ([first["pair"]] if "pair" in first else []) + choose(min(pairs, 4) - ("pair" in first), rest)
        rest = [k for k in rest if k not in pair_sides]
        rad_sides = ([first["radical"]] if "radical" in first else []) + choose(
            min(radicals, 4 - len(pair_sides)) - ("radical" in first), rest)
        pair_sides.sort(key=order.index)
        x0, y0, x1, y1 = lab.elem_ink_box()
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        out = []
        off = lp.gap + lp.dot_radius
        for kind, sides in (("pair", pair_sides), ("radical", rad_sides)):
            for k in sides:
                d = SIDES[k]
                mid = {"top": (cx, y0 - off), "bottom": (cx, y1 + off),
                       "left": (x0 - off, cy), "right": (x1 + off, cy)}[k]
                along = perp(d)
                if kind == "pair":
                    mid = self._clear_of_bonds(i, mid, along, (x0, y0, x1, y1))
                dots = ([add(mid, mul(along, lp.dot_spacing / 2)), add(mid, mul(along, -lp.dot_spacing / 2))]
                        if kind == "pair" else [mid])
                out.append(DotSet(i, kind, mid, d, k, dots))
        return out

    def _clear_of_bonds(self, i: int, mid: Point, along: Point, box: tuple) -> Point:
        """A pair on a crowded side (×0.6: a bond 50° off the side) slides along it, at most
        half the side, until neither dot touches a bond (the bond line beyond the element ink
        box); otherwise it stays centered."""
        a, lp = self.mol.atoms[i], self.preset.lone_pair
        reach = lp.dot_radius + self.preset.lengths.line_width / 2
        segs = []
        for j, _ in self.mol.neighbors(i):
            u = unit(sub(self.mol.atoms[j].p, a.p))
            segs.append((add(a.p, mul(u, ray_box_exit(a.p, u, box))), self.mol.atoms[j].p))

        def clear(m: Point) -> bool:
            return all(seg_point_dist(add(m, mul(along, sgn * lp.dot_spacing / 2)), p, q) >= reach
                       for sgn in (1, -1) for p, q in segs)
        if clear(mid):
            return mid
        half = abs(dot(along, (box[2] - box[0], box[3] - box[1]))) / 2
        for k in range(1, int(half / 0.1) + 1):
            for t in (0.1 * k, -0.1 * k):
                if clear(add(mid, mul(along, t))):
                    return add(mid, mul(along, t))
        return mid

    def bisector_dots(self, i: int, pairs: int, radicals: int) -> list[DotSet]:
        """Unlabeled atom (§4.2, mech_basic_x*: tert-butyl radical): each dot set goes on the
        bisector of the largest angular gap between the bonds, the charge (charge_center) and
        the sets already placed, as charges do; ties go the way charge_center breaks them. A set
        an arrow starts from (`tails`) goes first, into the widest gap facing its target.
        A radical sits radical.distance from the atom; a lone pair's midpoint does too
        (unmeasured: no reference has one), its dots dot_spacing apart across the bisector."""
        a = self.mol.atoms[i]
        lp, dist = self.preset.lone_pair, self.preset.radical.distance
        angs = [math.atan2(self.mol.atoms[j].y - a.y, self.mol.atoms[j].x - a.x) for j, _ in self.mol.neighbors(i)]
        if a.charge:
            cx, cy = self.charge_center(i)
            angs.append(math.atan2(cy - a.y, cx - a.x))
        out = []
        kinds = ["pair"] * pairs + ["radical"] * radicals
        hint = self.tails.get(i, {})
        kinds.sort(key=lambda k: k not in hint)  # an arrow's tail set is placed first
        for n, kind in enumerate(kinds):
            srt = sorted(angs)
            if srt:
                gaps = [((srt[(k + 1) % len(srt)] - srt[k]) % (2 * math.pi) or 2 * math.pi, srt[k])
                        for k in range(len(srt))]
                gap, start = max(gaps)
                if n == 0 and kind in hint:  # among the widest gaps (within 1°), the one facing the target
                    wide = [g for g in gaps if g[0] > gap - math.radians(1)]
                    d = hint[kind]
                    gap, start = max(wide, key=lambda g: math.cos(g[1] + g[0] / 2 - math.atan2(d[1], d[0])))
                t = start + gap / 2
            else:
                t = math.atan2(hint[kind][1], hint[kind][0]) if kind in hint else 0.0
            angs.append(t)
            d = (math.cos(t), math.sin(t))
            mid = add(a.p, mul(d, dist))
            along = perp(d)
            dots = ([add(mid, mul(along, lp.dot_spacing / 2)), add(mid, mul(along, -lp.dot_spacing / 2))]
                    if kind == "pair" else [mid])
            out.append(DotSet(i, kind, mid, d, _side_name(d), dots))
        return out

    # ---------------------------------------------------------------- assemble

    def draw(self, name: str) -> Group:
        lw = self.preset.lengths.line_width
        geoms = {b.idx: self.bond_geometry(b) for b in self.mol.bonds}
        bonds = Group(f"species-{name}-bonds")
        for line, closed, bis in self.chains(geoms):
            src = "bond:" + ",".join(self.mol.bonds[bi].id for bi in bis)
            bonds.children.append(Polyline(line, lw, closed=closed, source=src))
        for bi, g in geoms.items():
            src = f"bond:{self.mol.bonds[bi].id}"
            for v, u0, u1 in g.lines:
                bonds.children.append(Polyline([g.frame.at(u0, v), g.frame.at(u1, v)], lw, source=src))
            if g.polygon:
                bonds.children.append(Polygon(g.polygon, source=src))
            for uc, hl in g.hashes:
                bonds.children.append(Polyline([g.frame.at(uc, hl), g.frame.at(uc, -hl)],
                                                self.preset.lengths.hash_line_width, source=src))
        labels = Group(f"species-{name}-labels", [r for i in sorted(self.labels) for r in self.labels[i].runs])
        charges = Group(f"species-{name}-charges", [p for i, a in enumerate(self.mol.atoms) if a.charge
                                                    for p in self.charge_bars(i)])
        r = self.preset.lone_pair.dot_radius
        sets = [ds for i in sorted(self.dot_sets()) for ds in self.dot_sets()[i]]
        dots = [Group(f"species-{name}-{kind}", [Circle(p, r, source=f"atom:{self.mol.atoms[ds.atom].id}")
                                                 for ds in sets if ds.kind == k for p in ds.dots])
                for k, kind in (("pair", "lone-pairs"), ("radical", "radicals"))]
        return Group(f"species-{name}", [g for g in (bonds, labels, charges, *dots) if g.children])


def draw_molecule(mol: Molecule, preset: Preset, name: str | None = None, face: FontFace | None = None) -> Group:
    return Drawer(mol, preset, face).draw(name or mol.name)

