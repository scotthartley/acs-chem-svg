"""Molecules and schemes -> CDXML, so a figure can be edited in ChemDraw (acs-style.md §13).
Phase 1: M1 species. Phase 2 (acs-style.md §13.1): M2 schemes, i.e. reaction and equilibrium
arrows, conditions, plus signs and compound numbers as captions, and <scheme><step>. Phase 3
(acs-style.md §13.2): M3 mechanisms, i.e. curved arrows and lone-pair and radical dots. Phase 4a
(acs-style.md §13.3): abbreviations as nickname nodes, resonance ↔ arrows and brackets, and so
the wrapped and multi-row specs that use them. Phase 4b (acs-style.md §13.4, charges_lewis_x*):
circled signs, 2+ / 2- and charge_labels texts, two unpaired electrons, Lewis structures.

The header (root attributes, color and font tables, page) is presets/cdxml/<preset>.xml, and
each abbreviation's group presets/cdxml/<preset>.nicknames.xml, both copied from the ChemDraw
references by scripts/cdxml_headers.py; this module types no style numbers. Per species, one
<fragment>:
- <n> per atom at its drawn position: Element (atomic number; carbon has none), and on a labeled
  atom NumHydrogens and NeedsClean="yes", as ChemDraw writes its own labeled atoms; Charge (the
  SMILES charge, whatever is drawn); Radical="Doublet", or "Singlet" for two unpaired electrons
  (ChemDraw's own for two Electron dots, charges_lewis_x*' CH2; RDKit reads two radicals).
- A labeled atom's <t>: the text in CDXML order (NH2 even when drawn H2N), one <s> in the
  header's label font, size and face. p: x is the element symbol's glyph origin, or for a
  right-justified label (H on the left) the end of its advance (labels_x1.0, clearance_x*: -Cl
  sits 1.11 right of the atom); y is the label's center line plus half the ink height of the
  element's first glyph, as ChemDraw stores it (all references: 0.3684 x size for O, C, S, Cl;
  0.3586 for N, H, F, I, Br), though it draws every label's baseline at 0.3585 x size: half
  the ink top of the label's full-size glyphs (NO2 and TMS: 3.68, for the O and S).
  BoundingBox: the union of the drawn runs' ink boxes, moved by the same stored - drawn offset.
  LabelJustification / LabelAlignment Left or Right from Label.h_side, plus
  Justification="Right"; a one-character label (-F, -H) is always Left, and only labels of more
  than one character get LabelAlignment, as in the references.
- An abbreviation (phase 4a, abbrev_x*): <n NodeType="Fragment"> (OMe, CO2Me, NO2, ...) or
  "Nickname" (Ph, Bn, Et, Ac, silyl groups; IsNickname="yes" on the <t>), no Element,
  NeedsClean="yes". Inside, ChemDraw's group from the template, turned and moved rigidly so its
  attachment atom is on the atom and its ExternalConnectionPoint points at the one neighbor. The
  <t> is the key text (OMe) in both forms, placed as a label: a right-justified one (drawn MeO)
  has p.x at the end of the drawn text's advance.
- <b>: Order 2 / 3; wedges Display WedgeBegin / WedgedHashBegin with B at the narrow end
  (Bond.a). No DoublePosition: ChemDraw places double bonds itself. A Lewis structure drawn with
  dots for bonds: Visible="no" (ChemDraw hides it; the chemistry stays).
- A charge (_Writer.charge): <graphic GraphicType="Symbol" SymbolType="Plus|Minus"> with
  <represent attribute="Charge">, its BoundingBox the bars' center (Drawer.charge_center) and
  (center x - cdxml.charge_box_dx, center y). Circled: CirclePlus / CircleMinus, dx
  cdxml.circled_box_dx (ChemDraw sizes the drawn ring from the box). ChemDraw binds every sign
  symbol to the nearest atom on open and resets that atom's charge from it (phase 4b gate), so
  the drawing keeps each sign nearest its own atom (Drawer._nearest_atom). Digits (2+, a
  charge_labels text such as +3, a Lewis oxidation state) are one free <t> instead with
  InterpretChemically="no" ("2+", "+3", "+"; digit size, its ink centered where the drawing
  centers digits and bars), written after the fragments (user decisions 2026-10-28/29; ChemDraw
  has no multiple-charge symbol, a free sign symbol gets bound, and a free "+" text is read as
  a charge). ChemDraw drops the Charge of a node with no sign symbol (Fe2+ opens as Fe): the
  CDXML is for appearance, not chemistry.
- Dots (phase 3, mech_basic_x*), each bound to its atom by <represent attribute="Radical">: a
  lone pair is <graphic SymbolType="LonePair"> with BoundingBox the two dot centers,
  counterclockwise on screen around the atom; a radical dot SymbolType="Electron", BoundingBox
  (dot center, center x - cdxml.electron_box_dx, center y), and the node keeps Radical.
- Lewis structures (phase 4b): every atom labeled. A bond's electrons: each column of two dots
  across the bond a free LonePair (shared pairs; assigned ones counterclockwise around their
  owner; ChemDraw binds the off-center ones to the nearer atom, which only warns), an odd
  assigned electron a filled oval (<graphic GraphicType="Oval" OvalType="Circle Filled">,
  radius lone_pair.dot_radius; unmeasured form, user decision 2026-10-28), because ChemDraw
  makes any Electron near a node a radical of it (charges_lewis_x*' HCN, phase 4b gate). Bond
  lines (lewis_lines): plain bonds and no shared pairs.

Scheme objects (phase 2), from the neutral records CdxScheme / CdxArrow / CdxText / CdxStep
(scheme_records builds them from a SchemeLayout; tests build them from a reference):
- A caption <t>: every conditions line set, plus sign and compound number. A stacked caption
  with a sub- or superscript is one <t> per line instead (ChemDraw's line pitch for explicit
  script runs is not ours; user decision 2026-10-28), each placed as below. p is the first line's
  first glyph origin on its baseline (Left), or that line's advance center (Center, with
  CaptionJustification and Justification); LineHeight="auto"; lines joined with "\\n" and, for
  more than one, LineStarts (each line's end offset in UTF-8 bytes, its newline included:
  scheme_multi_x*'s "90 °C, 12 h" ends at 34). One <s> per richtext run at the nominal caption
  size: Helvetica (no face), Helvetica Oblique face 2, Helvetica Bold face 1, plus 32 for a
  subscript and 64 for a superscript. BoundingBox: the ink box, except that a text of more than
  one line reaches the advance end of each line that ends in a full-size run (scheme_*:
  ChemDraw's box of a stacked caption). No Warning attribute (ChemDraw's
  interpretation note; it writes its own).
- An arrow: <graphic GraphicType="Line" SupersededBy> (ArrowType FullHead, Equilibrium or Resonance,
  BoundingBox Head3D then Tail3D) plus the <arrow> it points to, with the same Z. Head
  attributes per head variant from cdxml.arrow_heads; ↔ arrows Full at both ends (phase 4a);
  equilibrium arrows HalfLeft at both ends,
  ArrowShaftSpacing = equilibrium.shaft_spacing and ArrowEquilibriumRatio =
  ratios.equilibrium_long_per_short, x 100 / LineWidth and x 100. Head3D / Tail3D: a forward
  arrow's tip and tail; a balanced equilibrium arrow's ends moved in by equilibrium.stagger; one
  favoring products, its ends; one favoring reactants, reversed (Head3D on the left, as ChemDraw
  stores it, so the step's reactants are the species on its right). Center3D is the midpoint,
  MajorAxisEnd3D = Center3D + (Head3D - Tail3D), MinorAxisEnd3D = Center3D + that turned 90°.
  BoundingBox: along the axis the drawn ends; across it the half width (ArrowheadWidth x
  LineWidth / 100, plus half the shaft spacing), with cdxml.arrow_box_pad more on the top of a
  forward arrow and on both sides of an equilibrium arrow whose Head3D is on the right; a ↔
  arrow's, cdxml.resonance_box_pad on both sides.
- <scheme> per drawn sequence (a new_row step starts one; the rows of a wrapped sequence stay
  one), <step> with ReactionStepReactants / Products (fragment ids), ReactionStepArrows (the
  <graphic> id, as ChemDraw cites it) and ReactionStepObjectsAbove / BelowArrow (<t> ids, and
  the fragments of species drawn above the arrow: a reactant a mechanism step adds, or `over`;
  user decision 2026-10-28). A resonance step: one <step> per ↔, contributor k to k + 1, in
  its sequence's <scheme>.
- A bracket (phase 4a, CdxBracket): a free <graphic GraphicType="Bracket" BracketType="Square"
  LipSize=cdxml.bracket_lip_size>, no <bracketedgroup> (user decision 2026-10-28); BoundingBox
  the bar's x and its lips' stroke centers, (x, bottom, x, top) for "[" and (x, top, x, bottom)
  for "]", as resonance_x* stores them.
- A curved arrow (phase 3, CdxCurve): a free <curve> on the page (no <group>; no step cites it;
  user decision 2026-10-28), no BoundingBox. Two-electron: CurveType 8, ArrowheadHead Full,
  CurvePoints (2p0 - c1, p0, c1, c2, p3, 2p3 - c2). Fishhook: CurveType 32, the same six points
  (the references repeat the last anchor: not written, see _Writer.curve), ArrowheadHead HalfRight when mechanism.fishhook_side is +1 (the
  references' fishhooks), HalfLeft otherwise. p3 is the stroke end, not the head tip: every
  reference curve's last anchor is its stroke end.
Page order as in the references: fragments, captions, graphics (arrows, then brackets), curves,
schemes, arrows.

Not written: atom maps. They warn, except in a mechanism (CdxScheme.warn_maps): its arrows
need them, and the figure never shows them.

Fragment BoundingBox: the drawn ink box (Drawer.draw: stroke widths, second lines, charge signs)
with the stored label boxes and each dot center plus cdxml.lone_pair_box_reach /
electron_box_reach, as ChemDraw's; the root's also covers the captions, arrows, brackets and
curves' ink.
Coordinates are in the caller's frame: each species' offset is added to its atom positions, and
CdxScheme.offset to the arrows, captions, curves and brackets (render.py passes to_svg's offset, so the
file lines up with the SVG).
Numbers are written with two decimals, as ChemDraw writes them. CreationProgram names this
package, and the page grows by whole pages to fit the drawing.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Callable
from functools import lru_cache
from importlib.metadata import version
from dataclasses import dataclass, field

from lxml import etree
from rdkit import Chem

from . import richtext
from .brackets import square_bracket
from .chem import Molecule
from .draw_mol import Drawer
from .geometry import Point, add, mul, sub, union_box, unit
from .layout import Placed, SchemeLayout, build_scheme, m1_placed
from .mechanism import curved_arrow, fishhook_side
from .preset import PRESET_DIR, load_preset
from .scene import Group, TextRun
from .spec import ResonanceStep, Spec, load_spec
from .svg_writer import to_svg

HEADER_DIR = PRESET_DIR / "cdxml"
_PT = Chem.GetPeriodicTable()


class CdxmlUnsupported(ValueError):
    """The species can't be written as correct chemistry in this phase."""


@dataclass
class CdxArrow:
    """A step's arrow as drawn: tail and head as the layout passes them to arrows.step_arrow."""
    tail: Point
    head: Point
    kind: str = "forward"  # forward | equilibrium | equilibrium_products | equilibrium_reactants | resonance
    variant: str | None = None  # arrow.heads variant; None: the preset's default_head


@dataclass
class CdxText:
    """A caption: laid-out lines (richtext.Line with its parts and runs) and its justification."""
    lines: list[richtext.Line]
    justify: str = "Left"  # Left | Center


@dataclass
class CdxBracket:
    """A square bracket as drawn: its bar's stroke center x, its lips' stroke centers top and
    bottom; side "left" ([) or "right" (])."""
    x: float
    top: float
    bottom: float
    side: str


@dataclass
class CdxCurve:
    """A curved arrow as drawn: the cubic p0 -> p3 with controls c1, c2 (mechanism.MechArrow;
    p3 is the stroke end, not the head tip)."""
    p0: Point
    c1: Point
    c2: Point
    p3: Point
    kind: str = "two_electron"  # two_electron | fishhook


@dataclass
class CdxStep:
    """Indices into CdxScheme.species, .arrows and .texts; above_species / below_species index
    CdxScheme.species drawn above / below the arrow (a reactant a mechanism step adds)."""
    reactants: list[int]
    products: list[int]
    arrow: int
    above: list[int] = field(default_factory=list)
    below: list[int] = field(default_factory=list)
    above_species: list[int] = field(default_factory=list)
    below_species: list[int] = field(default_factory=list)


@dataclass
class CdxScheme:
    species: list[tuple[Molecule, Drawer, Point]]  # (molecule, its Drawer, offset added to its atoms)
    arrows: list[CdxArrow] = field(default_factory=list)
    texts: list[CdxText] = field(default_factory=list)
    schemes: list[list[CdxStep]] = field(default_factory=list)  # one <scheme> per drawn sequence
    curves: list[CdxCurve] = field(default_factory=list)
    brackets: list[CdxBracket] = field(default_factory=list)
    offset: Point = (0.0, 0.0)  # added to the arrows, captions, curves and brackets
    warn_maps: bool = True  # False for a mechanism: every species is mapped (user decision 2026-10-28)


def f2(v: float) -> str:
    s = f"{v:.2f}"
    return "0.00" if s == "-0.00" else s


def pts(*xy: float) -> str:
    return " ".join(f2(v) for v in xy)


@dataclass
class Ids:
    """Sequential object ids and Z order, as ChemDraw numbers them."""
    next_id: int = 1
    z: int = 1

    def new(self) -> str:
        self.next_id += 1
        return str(self.next_id - 1)

    def next_z(self) -> str:
        self.z += 1
        return str(self.z - 1)


@lru_cache(maxsize=None)
def load_nicknames(preset_name: str) -> dict[str, etree._Element]:
    """ABBREVIATIONS key -> its <nickname> template (presets/cdxml/<preset>.nicknames.xml)."""
    path = HEADER_DIR / f"{preset_name}.nicknames.xml"
    if not path.exists():
        raise FileNotFoundError(f"no CDXML nickname templates for preset {preset_name!r} ({path}); "
                                "run scripts/cdxml_headers.py")
    root = etree.parse(str(path), etree.XMLParser(load_dtd=False, no_network=True, remove_blank_text=True)).getroot()
    return {el.get("key"): el for el in root}


def load_header(preset_name: str) -> etree._ElementTree:
    path = HEADER_DIR / f"{preset_name}.xml"
    if not path.exists():
        raise FileNotFoundError(f"no CDXML header for preset {preset_name!r} ({path}); run scripts/cdxml_headers.py")
    return etree.parse(str(path), etree.XMLParser(load_dtd=False, no_network=True, remove_blank_text=True))


def label_text(element: str, num_h: int) -> str:
    return element + ("H" + (str(num_h) if num_h > 1 else "") if num_h else "")


class _Writer:
    def __init__(self, header: etree._ElementTree, preset_name: str, warn: Callable[[str], None]):
        self.root = header.getroot()
        self.preset = load_preset(preset_name)
        self.warn = warn
        self.warn_maps = True
        self.ids = Ids()
        font_name = self.preset.font_family
        font = next((f for f in self.root.iter("font") if f.get("name") == font_name), None)
        if font is None:
            raise ValueError(f"header for {preset_name} has no font {font_name!r}")
        self.s_attrs = {"font": font.get("id"), "size": self.root.get("LabelSize"), "color": self.root.get("color", "0"),
                        "face": self.root.get("LabelFace")}
        fonts = {f.get("name"): f.get("id") for f in self.root.iter("font")}
        self.caption_fonts = {style: fonts.get(name) for style, name in (
            ("regular", font_name), ("italic", f"{font_name} Oblique"), ("bold", f"{font_name} Bold"))}
        self.line_width = float(self.root.get("LineWidth"))
        self.preset_name = preset_name
        self.digit_texts: list[etree._Element] = []  # charge digits' <t>s, written after the fragments

    def check(self, mol: Molecule) -> None:
        if self.warn_maps and any(a.map for a in mol.atoms):
            self.warn(f"{mol.name}: atom maps not written")

    def fragment(self, page, mol: Molecule, d: Drawer, off: Point) -> tuple:
        self.check(mol)
        frag = etree.SubElement(page, "fragment", {"id": self.ids.new()})
        node_ids, label_boxes = {}, []
        for a in mol.atoms:
            x, y = a.x + off[0], a.y + off[1]
            node_ids[a.idx] = nid = self.ids.new()
            attrs = {"id": nid, "p": pts(x, y), "Z": self.ids.next_z()}
            if a.abbrev is not None:  # a nickname node: no Element; its atoms are in the inner fragment
                tpl = self.nickname_template(mol, a)
                n = etree.SubElement(frag, "n", {**attrs, "NodeType": tpl.get("NodeType"), "NeedsClean": "yes"})
                self.nickname(n, tpl, mol, a, off)
                label_boxes.append(self.label(n, d.labels[a.idx], a.abbrev, off, tpl.get("IsNickname")))
                continue
            if a.element != "C":
                attrs["Element"] = str(_PT.GetAtomicNumber(a.element))
            lab = d.labels.get(a.idx)
            if lab is not None:
                attrs["NumHydrogens"] = str(a.num_h)
            if a.charge:
                attrs["Charge"] = str(a.charge)
            if a.radicals:  # two Electron symbols: ChemDraw writes Singlet (charges_lewis_x*' CH2)
                attrs["Radical"] = "Doublet" if a.radicals == 1 else "Singlet"
            if lab is not None:
                attrs["NeedsClean"] = "yes"
            n = etree.SubElement(frag, "n", attrs)
            if lab is not None:
                label_boxes.append(self.label(n, lab, label_text(a.element, a.num_h), off))
        zs = [self.ids.next_z() for _ in mol.bonds]
        rank = z_order(mol, d)
        for b in mol.bonds:
            attrs = {"id": self.ids.new(), "Z": zs[rank[b.idx]], "B": node_ids[b.a], "E": node_ids[b.b]}
            if b.order != 1:
                attrs["Order"] = str(b.order)
            if b.stereo:
                attrs["Display"] = {"wedge": "WedgeBegin", "hash": "WedgedHashBegin"}[b.stereo]
            if mol.lewis and not mol.lewis_lines:  # the dots show the bond; the bond keeps the chemistry
                attrs["Visible"] = "no"
            etree.SubElement(frag, "b", attrs)
        for a in mol.atoms:
            if d.charge_sign(a.idx) is not None:
                self.charge(frag, d, a, off, node_ids[a.idx])
        dot_boxes = self.dots(frag, mol, d, off, node_ids)
        # ChemDraw's fragment box: the drawn ink (strokes, second lines, signs), the stored label
        # boxes, and each dot center plus cdxml.*_box_reach
        drawn = Group("ink", [g for g in d.draw(mol.name).children
                              if not g.id.endswith(("-labels", "-lone-pairs", "-radicals", "-bonding-pairs", "-bond-electrons"))])
        b = drawn.ink_bbox()
        box = union_box(([(b[0] + off[0], b[1] + off[1], b[2] + off[0], b[3] + off[1])] if b else [])
                        + label_boxes + dot_boxes)
        frag.set("BoundingBox", pts(*box))
        frag.set("Z", self.ids.next_z())
        return box

    def digit_text(self, d: Drawer, a, off: Point) -> tuple[str, Point, tuple]:
        """Atom a's drawn charge digits and sign ("2+", "+3"): the text, its origin with its ink
        centered where the drawing centers it, and its ink box."""
        sign, center = d.charge_sign(a.idx), d.charge_center(a.idx)
        mark = "+" if sign.sign > 0 else "\u2212"
        text = sign.digits + mark if sign.digits_first else mark + sign.digits
        b = sign.face.text_ink_bounds(text, sign.size)
        x = center[0] - (b[0] + b[2]) / 2 + off[0]
        y = center[1] + (b[1] + b[3]) / 2 + off[1]
        return text, (x, y), (x + b[0], y - b[3], x + b[2], y - b[1])

    def charge(self, frag, d: Drawer, a, off: Point, nid: str) -> None:
        """Atom a's drawn charge (§13.4). A plain ±1 or circled sign is ChemDraw's symbol bound to
        the node (<represent attribute="Charge">), box (center, center x - charge_box_dx /
        circled_box_dx, center y). Digits (2+, charge_labels texts such as +3, oxidation states)
        are one free <t> ("2+", "+3", "+"), its ink centered where the drawing centers digits and
        bars (user decision 2026-10-28): ChemDraw has no symbol for them and rebinds any free
        sign to the nearest atom. The node's Charge stays the SMILES charge."""
        sign, cx = d.charge_sign(a.idx), self.preset.cdxml
        center = d.charge_center(a.idx)
        if sign.digits or a.charge_label is not None:
            # free text: ChemDraw binds any sign symbol to the nearest atom and resets its charge
            # to the symbol's (phase 4b gate: Fe2+ came back Fe+, every oxidation-state H a H+),
            # and reads a "+" text beside an atom as a charge unless told not to (gate round 2)
            text, (x, y), box = self.digit_text(d, a, off)
            t = etree.Element("t", {"id": self.ids.new(), "p": pts(x, y), "BoundingBox": pts(*box),
                                    "Z": self.ids.next_z(), "InterpretChemically": "no", "LineHeight": "auto"})
            etree.SubElement(t, "s", {"font": self.s_attrs["font"], "size": f"{sign.size:g}",
                                      "color": self.s_attrs["color"]}).text = text
            self.digit_texts.append(t)
            return
        segs, _ = sign._layout(center)
        p, q = segs[0]
        c = add(mul(add(p, q), 0.5), off)
        if sign.circled:
            sym, dx = ("CirclePlus" if sign.sign > 0 else "CircleMinus"), cx.circled_box_dx
        else:
            sym, dx = ("Plus" if sign.sign > 0 else "Minus"), cx.charge_box_dx
        g = etree.SubElement(frag, "graphic", {"id": self.ids.new(), "BoundingBox": pts(c[0], c[1], c[0] - dx, c[1]),
                                               "Z": self.ids.next_z(), "GraphicType": "Symbol", "SymbolType": sym})
        etree.SubElement(g, "represent", {"attribute": "Charge", "object": nid})

    def dots(self, frag, mol: Molecule, d: Drawer, off: Point, node_ids: dict) -> list[tuple]:
        """Each atom's lone pairs (LonePair, the two dot centers counterclockwise on screen around
        the atom) and radical dots (Electron, (center, center x - cdxml.electron_box_dx, center
        y)), bound to the atom by <represent attribute="Radical"> (mech_basic_x*). The boxes
        they add to the fragment's."""
        cx = self.preset.cdxml
        boxes = []

        def symbol(sym: str, bb: tuple, centers: list[Point], owner: int | None) -> None:
            g = etree.SubElement(frag, "graphic", {"id": self.ids.new(), "BoundingBox": pts(*bb), "Z": self.ids.next_z(),
                                                   "GraphicType": "Symbol", "SymbolType": sym})
            if owner is not None:
                etree.SubElement(g, "represent", {"attribute": "Radical", "object": node_ids[owner]})
            r = cx.lone_pair_box_reach if sym == "LonePair" else cx.electron_box_reach
            boxes.extend((c[0] - r, c[1] - r, c[0] + r, c[1] + r) for c in centers)

        def pair(p: Point, q: Point, around: Point, owner: int | None) -> None:
            u, v = sub(p, around), sub(q, around)
            if u[0] * v[1] - u[1] * v[0] > 0:  # clockwise on screen (y down): swap
                p, q = q, p
            symbol("LonePair", (*p, *q), [p, q], owner)

        def oval(p: Point) -> None:  # a dot that is only a drawing (user decision 2026-10-28)
            r = self.preset.lone_pair.dot_radius
            etree.SubElement(frag, "graphic", {
                "id": self.ids.new(), "BoundingBox": pts(p[0], p[1], p[0] + r, p[1]), "Z": self.ids.next_z(),
                "GraphicType": "Oval", "OvalType": "Circle Filled", "Center3D": _p3(p),
                "MajorAxisEnd3D": _p3((p[0] + r, p[1])), "MinorAxisEnd3D": _p3((p[0], p[1] + r))})
            boxes.append((p[0] - r, p[1] - r, p[0] + r, p[1] + r))

        def electron(p: Point, owner: int | None) -> None:
            symbol("Electron", (p[0], p[1], p[0] - cx.electron_box_dx, p[1]), [p], owner)

        sets = [ds for i in sorted(d.dot_sets()) for ds in d.dot_sets()[i]]
        for ds in sets + (d.lewis_dot_sets() if mol.lewis else []):
            i, a = ds.atom, mol.atoms[ds.atom]
            dots = [add(t, off) for t in ds.dots]
            if ds.kind == "pair":
                pair(*dots, add(a.p, off), i)
            elif ds.kind == "radical":
                electron(dots[0], i)
            elif ds.kind in ("shared", "assigned"):
                # Lewis (§13.4): each column of two dots across the bond a free LonePair (ChemDraw
                # binds an off-center one to the nearer atom: only a warning); an odd last
                # electron a filled oval, since ChemDraw makes any Electron near a node a radical
                b = mol.bonds[ds.bond]
                mid = add(mul(add(mol.atoms[b.a].p, mol.atoms[b.b].p), 0.5), off)
                around = add(a.p, off) if ds.kind == "assigned" else mid
                for k in range(0, len(dots) - 1, 2):
                    pair(dots[k], dots[k + 1], around, None)
                if len(dots) % 2:
                    oval(dots[-1])
            else:
                raise CdxmlUnsupported(f"{mol.name}: {ds.kind} dots are not written to CDXML")
        return boxes

    def nickname_template(self, mol: Molecule, a) -> etree._Element:
        tpl = load_nicknames(self.preset_name).get(a.abbrev)
        if tpl is None:
            raise CdxmlUnsupported(f"{mol.name}: no CDXML nickname template for {a.abbrev!r}")
        if mol.degree(a.idx) != 1:
            raise CdxmlUnsupported(f"{mol.name}: abbreviation {a.abbrev} with {mol.degree(a.idx)} bonds")
        return tpl

    def nickname(self, n, tpl: etree._Element, mol: Molecule, a, off: Point) -> None:
        """The nickname's inner <fragment> (abbrev_x*): ChemDraw's group from the template, turned
        and moved rigidly so its attachment atom is on `a` and its ExternalConnectionPoint (on -x
        in the template) points at a's one neighbor."""
        (j, _), = mol.neighbors(a.idx)
        u = unit(sub(mol.atoms[j].p, a.p))
        c, s = -u[0], -u[1]  # the template's +x
        src = tpl.find("fragment")
        frag = etree.SubElement(n, "fragment", {"id": self.ids.new()})
        ids = {m.get("id"): self.ids.new() for m in src.findall("n")}
        for m in src.findall("n"):
            vx, vy = (float(v) for v in m.get("p").split())
            p = (a.x + off[0] + vx * c - vy * s, a.y + off[1] + vx * s + vy * c)
            etree.SubElement(frag, "n", {"id": ids[m.get("id")], "p": pts(*p),
                                         **{k: v for k, v in m.attrib.items() if k not in ("id", "p")}})
        for b in src.findall("b"):
            etree.SubElement(frag, "b", {"id": self.ids.new(), "B": ids[b.get("B")], "E": ids[b.get("E")],
                                         **{k: v for k, v in b.attrib.items() if k not in ("id", "B", "E")}})

    def label(self, n, lab, text: str, off: Point, is_nickname: str | None = None) -> tuple:
        """A label's <t>: `text` is the CDXML text (NH2, or an abbreviation's key, OMe, in both
        forms)."""
        # bonds to the right: right-justified (H2N, HO, -Cl, MeO); one glyph (-F, -H) stays Left
        right = lab.h_side == "left" and len(text) > 1
        x, base = lab.elem_origin  # the element's (an abbreviation's anchor glyph's) origin
        if right:  # the drawn text ends with the element or anchor glyph
            x += lab.face.advance(lab.element, lab.size)
        # ChemDraw stores y half the ink top of the label's full-size glyphs below the label's
        # center line (atom y, or where center_terminal_labels moved it), not the drawn baseline:
        # 3.68 for O, C, S (OH, Cl, NO2, TMS), 3.59 for N, H, F, T (NH2, NHAc), 3.60 for Ph
        center = base - self.preset.ratios.label_baseline_per_size * lab.size
        y = center + max(r.face.text_ink_bounds(r.text, r.size)[3] for r in lab.runs if r.size == lab.size) / 2
        box = union_box(r.bbox() for r in lab.runs)  # in the stored frame too: moved by y - base
        box = (box[0] + off[0], box[1] + y - base + off[1], box[2] + off[0], box[3] + y - base + off[1])
        just = "Right" if right else "Left"
        attrs = {"p": pts(x + off[0], y + off[1]), "BoundingBox": pts(*box), "LabelJustification": just}
        if right:
            attrs["Justification"] = "Right"
        if len(text) > 1:
            attrs["LabelAlignment"] = just
        if is_nickname:
            attrs["IsNickname"] = is_nickname
        t = etree.SubElement(n, "t", attrs)
        etree.SubElement(t, "s", self.s_attrs).text = text
        return box


    # ------------------------------------------------------------------ captions

    def _runs(self, text: CdxText) -> list[tuple[dict, str]]:
        """(<s> attributes, text) per run, line breaks on the last run of each line, adjacent
        runs with the same attributes merged."""
        size = f"{self.preset.text.caption_size:g}"
        out: list[tuple[dict, str]] = []
        for k, ln in enumerate(text.lines):
            for j, part in enumerate(ln.parts):
                font = self.caption_fonts[part.style]
                if font is None:
                    raise ValueError(f"CDXML header has no font for {part.style} captions; run scripts/cdxml_headers.py")
                face = {"regular": 0, "italic": 2, "bold": 1}[part.style] | {"base": 0, "sub": 32, "sup": 64}[part.script]
                attrs = {"font": font, "size": size, "color": self.root.get("color", "0")}
                if face:
                    attrs["face"] = str(face)
                t = part.text + ("\n" if j == len(ln.parts) - 1 and k < len(text.lines) - 1 else "")
                if out and out[-1][0] == attrs:
                    out[-1] = (attrs, out[-1][1] + t)
                else:
                    out.append((attrs, t))
        return out

    def captions(self, page, text: CdxText, off: Point) -> list[tuple[str, tuple]]:
        """The <t>s of a caption: one, or one per line for a stacked caption with a sub- or
        superscript (ChemDraw sets those lines 11.775 apart at x1.0, not our 10.825 / 10.975: gate
        2026-10-28, acs-style.md §13.1). [(id, BoundingBox)]."""
        if len(text.lines) > 1 and any(p.script != "base" for ln in text.lines for p in ln.parts):
            return [self.caption(page, CdxText([ln], text.justify), off) for ln in text.lines]
        return [self.caption(page, text, off)]

    def caption(self, page, text: CdxText, off: Point) -> tuple[str, tuple]:
        """A <t> on the page; (its id, its BoundingBox)."""
        first = text.lines[0]
        x0, base = _origin(first, self.preset)
        x = x0 + (first.advance / 2 if text.justify == "Center" else 0.0)
        box = list(union_box(ln.ink() for ln in text.lines))
        if len(text.lines) > 1:  # ChemDraw's stacked caption box reaches each full-size line end
            for ln in text.lines:
                if ln.parts and ln.parts[-1].script == "base":
                    box[2] = max(box[2], _origin(ln, self.preset)[0] + ln.advance)
        box = (box[0] + off[0], box[1] + off[1], box[2] + off[0], box[3] + off[1])
        tid = self.ids.new()
        attrs = {"id": tid, "p": pts(x + off[0], base + off[1]), "BoundingBox": pts(*box), "Z": self.ids.next_z()}
        if text.justify == "Center":
            attrs["CaptionJustification"] = attrs["Justification"] = "Center"
        elif text.justify != "Left":
            raise ValueError(f"caption justification {text.justify!r}: Left or Center")
        attrs["LineHeight"] = "auto"
        if len(text.lines) > 1:
            ends, n = [], 0
            for k, ln in enumerate(text.lines):
                n += len(ln.text.encode()) + (1 if k < len(text.lines) - 1 else 0)  # UTF-8 bytes: "°" counts 2
                ends.append(str(n))
            attrs["LineStarts"] = " ".join(ends)
        t = etree.SubElement(page, "t", attrs)
        for a, txt in self._runs(text):
            etree.SubElement(t, "s", a).text = txt
        return tid, box

    # ------------------------------------------------------------------ arrows

    def arrow_ends(self, a: CdxArrow) -> tuple[Point, Point]:
        """(Tail3D, Head3D) for a drawn arrow (module docstring)."""
        if a.kind == "equilibrium":
            d = mul(unit(sub(a.head, a.tail)), self.preset.equilibrium.stagger)
            return add(a.tail, d), sub(a.head, d)
        if a.kind == "equilibrium_reactants":
            return a.head, a.tail
        if a.kind in ("forward", "equilibrium_products", "resonance"):
            return a.tail, a.head
        raise ValueError(f"unknown arrow kind {a.kind!r}")

    def arrow_attrs(self, a: CdxArrow) -> dict:
        """The <arrow> head attributes, in ChemDraw's order."""
        variant = a.variant or self.preset.arrow.default_head
        if variant not in self.preset.cdxml.arrow_heads:
            raise ValueError(f"no CDXML head attributes for arrow head {variant!r} (preset cdxml.arrow_heads)")
        size, center, width = self.preset.cdxml.arrow_heads[variant]
        eq = a.kind.startswith("equilibrium")
        attrs = {"FillType": "None", "ArrowheadHead": "HalfLeft" if eq else "Full"}
        if eq or a.kind == "resonance":
            attrs["ArrowheadTail"] = "HalfLeft" if eq else "Full"
        attrs.update({"ArrowheadType": "Solid", "HeadSize": str(size), "ArrowheadCenterSize": str(center),
                      "ArrowheadWidth": str(width)})
        if eq:
            attrs["ArrowShaftSpacing"] = str(round(self.preset.equilibrium.shaft_spacing / self.line_width * 100))
            if a.kind != "equilibrium":
                attrs["ArrowEquilibriumRatio"] = str(round(self.preset.ratios.equilibrium_long_per_short * 100))
        return attrs

    def arrow_box(self, a: CdxArrow, attrs: dict, tail3d: Point, head3d: Point) -> tuple:
        """<arrow> BoundingBox of a horizontal arrow (module docstring)."""
        if abs(a.head[1] - a.tail[1]) > 1e-6:
            raise CdxmlUnsupported("CDXML arrow boxes are measured for horizontal arrows only")
        lw = self.line_width
        half = int(attrs["ArrowheadWidth"]) * lw / 100 + int(attrs.get("ArrowShaftSpacing", 0)) * lw / 200
        pad = self.preset.cdxml.arrow_box_pad
        if a.kind == "forward":
            top, bottom = pad, 0.0
        elif a.kind == "resonance":
            top = bottom = self.preset.cdxml.resonance_box_pad
        else:
            top = bottom = pad if head3d[0] > tail3d[0] else 0.0
        y = a.tail[1]
        return (min(a.tail[0], a.head[0]), y - half - top, max(a.tail[0], a.head[0]), y + half + bottom)

    def arrow(self, page, a: CdxArrow, off: Point) -> tuple[str, str, tuple, etree._Element]:
        """The legacy <graphic> on the page and its <arrow> (returned, appended later); (graphic
        id, arrow id, box, arrow element)."""
        attrs = self.arrow_attrs(a)
        tail3d, head3d = (add(p, off) for p in self.arrow_ends(a))
        a = CdxArrow(add(a.tail, off), add(a.head, off), a.kind, a.variant)
        box = self.arrow_box(a, attrs, tail3d, head3d)
        gid, aid, z = self.ids.new(), self.ids.new(), self.ids.next_z()
        etree.SubElement(page, "graphic", {
            "id": gid, "SupersededBy": aid, "BoundingBox": pts(*head3d, *tail3d), "Z": z, "GraphicType": "Line",
            "ArrowType": {"forward": "FullHead", "resonance": "Resonance"}.get(a.kind, "Equilibrium"),
            "HeadSize": attrs["HeadSize"]})
        c = mul(add(head3d, tail3d), 0.5)
        d = sub(head3d, tail3d)
        el = etree.Element("arrow", {"id": aid, "BoundingBox": pts(*box), "Z": z, **attrs,
                                     "Head3D": _p3(head3d), "Tail3D": _p3(tail3d), "Center3D": _p3(c),
                                     "MajorAxisEnd3D": _p3(add(c, d)), "MinorAxisEnd3D": _p3(add(c, (-d[1], d[0])))})
        return gid, aid, box, el


    # ------------------------------------------------------------------ brackets

    def bracket(self, page, br: CdxBracket, off: Point) -> tuple:
        """A free bracket <graphic> (resonance_x*; no <bracketedgroup>, user decision
        2026-10-28): BoundingBox the bar's ends, the bottom first for "[" and the top first for
        "]", as ChemDraw writes them; its drawn ink box."""
        x, top, bottom = br.x + off[0], br.top + off[1], br.bottom + off[1]
        ys = (bottom, top) if br.side == "left" else (top, bottom)
        etree.SubElement(page, "graphic", {
            "id": self.ids.new(), "BoundingBox": pts(x, ys[0], x, ys[1]), "Z": self.ids.next_z(),
            "GraphicType": "Bracket", "BracketType": "Square", "LipSize": str(self.preset.cdxml.bracket_lip_size)})
        return square_bracket(x, top, bottom, br.side, self.preset).ink_bbox()

    # ------------------------------------------------------------------ curved arrows

    def curve(self, page, c: CdxCurve, off: Point) -> tuple:
        """A free <curve> on the page (§13.2); its ink box. CurvePoints: (in-handle, anchor,
        out-handle) per anchor, the end handles mirrored (2p0 - c1, 2p3 - c2). The last anchor
        is the stroke end (mech_basic_x*). A fishhook has the same six points: the references'
        third anchor (a repeat of the last, with its handles) is a zero-length loop that ChemDraw
        draws and edits as an extra node, as long as our end handle (gate 2026-10-28)."""
        p0, c1, c2, p3 = (add(p, off) for p in (c.p0, c.c1, c.c2, c.p3))
        cps = [sub(mul(p0, 2), c1), p0, c1, c2, p3, sub(mul(p3, 2), c2)]
        attrs = {"id": self.ids.new(), "Z": self.ids.next_z()}
        if c.kind == "two_electron":
            attrs.update({"CurveType": "8", "ArrowheadHead": "Full"})
        elif c.kind == "fishhook":
            side = fishhook_side(c1, c2, p3)
            attrs.update({"CurveType": "32", "ArrowheadHead": "HalfRight" if side > 0 else "HalfLeft"})
        else:
            raise ValueError(f"unknown curved arrow kind {c.kind!r}")
        attrs.update({"ArrowheadType": "Solid", "CurvePoints": pts(*(v for q in cps for v in q))})
        etree.SubElement(page, "curve", attrs)
        return curved_arrow(p0, c1, c2, p3, c.kind, self.preset).ink_bbox()


def _p3(p: Point) -> str:
    return f"{pts(*p)} 0"


def _origin(ln: richtext.Line, preset) -> Point:
    """A line's first glyph origin and baseline, from its runs (they move with their groups;
    Line.x0 / baseline may not)."""
    r = ln.runs[0]
    return r.x, r.y - richtext.run_metrics(ln.parts[0], preset.text.caption_size, preset)[1]


def _stderr(msg: str) -> None:
    print(f"warning: {msg}", file=sys.stderr)


def to_cdxml(scheme: CdxScheme | list[tuple[Molecule, Drawer, Point]], preset_name: str, name: str = "",
             warn: Callable[[str], None] = _stderr) -> str:
    """CDXML text for a CdxScheme, or for a list of species: (molecule, its Drawer, the offset
    added to its coordinates). The Drawer's labels must be in molecule coordinates (Drawer.draw
    leaves them so)."""
    if not isinstance(scheme, CdxScheme):
        scheme = CdxScheme(list(scheme))
    header = load_header(preset_name)
    w = _Writer(header, preset_name, warn)
    w.warn_maps = scheme.warn_maps
    root = w.root
    page = root.find("page")
    page.set("id", w.ids.new())
    frag_ids, boxes = [], []
    for mol, d, off in scheme.species:
        boxes.append(w.fragment(page, mol, d, off))
        frag_ids.append(page[-1].get("id"))
    page.extend(w.digit_texts)
    text_ids = []  # per CdxText, its <t> ids
    for t in scheme.texts:
        written = w.captions(page, t, scheme.offset)
        text_ids.append([tid for tid, _ in written])
        boxes += [b for _, b in written]
    graphic_ids, arrows = [], []
    for a in scheme.arrows:
        gid, _, b, el = w.arrow(page, a, scheme.offset)
        graphic_ids.append(gid)
        arrows.append(el)
        boxes.append(b)
    for br in scheme.brackets:
        boxes.append(w.bracket(page, br, scheme.offset))
    for c in scheme.curves:
        boxes.append(w.curve(page, c, scheme.offset))
    for steps in scheme.schemes:
        sch = etree.SubElement(page, "scheme", {"id": w.ids.new()})
        for st in steps:
            attrs = {"id": w.ids.new()}
            for key, ids in (("ReactionStepReactants", [frag_ids[i] for i in st.reactants]),
                             ("ReactionStepProducts", [frag_ids[i] for i in st.products]),
                             ("ReactionStepArrows", [graphic_ids[st.arrow]]),
                             ("ReactionStepObjectsAboveArrow", [frag_ids[i] for i in st.above_species]
                              + [j for i in st.above for j in text_ids[i]]),
                             ("ReactionStepObjectsBelowArrow", [frag_ids[i] for i in st.below_species]
                              + [j for i in st.below for j in text_ids[i]])):
                if ids:
                    attrs[key] = " ".join(ids)
            etree.SubElement(sch, "step", attrs)
    page.extend(arrows)
    box = union_box(boxes) if boxes else (0.0, 0.0, 0.0, 0.0)
    attrs = dict(root.attrib)
    root.attrib.clear()
    # CreationProgram, Name and BoundingBox first, as ChemDraw orders them
    root.set("CreationProgram", f"acs-chem-svg {version('acs-chem-svg')}")
    root.set("Name", name)
    root.set("BoundingBox", pts(*box))
    for k, v in attrs.items():
        root.set(k, v)
    # Whole pages of the template's size, enough for the drawing: ChemDraw offers to expand a
    # page it doesn't fit in, and writes BoundingBox="0 0 1080 720" WidthPages="2" (re-save of
    # m1_smiles_x1.0, 2026-10-28).
    x0, y0, pw, ph = (float(v) for v in page.get("BoundingBox").split())
    pw, ph = (pw - x0) / int(page.get("WidthPages")), (ph - y0) / int(page.get("HeightPages"))
    wide, tall = max(1, math.ceil((box[2] - x0) / pw)), max(1, math.ceil((box[3] - y0) / ph))
    page.set("BoundingBox", f"{x0:g} {y0:g} {x0 + wide * pw:g} {y0 + tall * ph:g}")
    page.set("WidthPages", str(wide))
    page.set("HeightPages", str(tall))
    # page id first, as in the references
    pattrs = dict(page.attrib)
    page.attrib.clear()
    page.set("id", pattrs.pop("id"))
    for k, v in pattrs.items():
        page.set(k, v)
    return etree.tostring(header, xml_declaration=True, encoding="UTF-8", doctype=header.docinfo.doctype,
                          pretty_print=True).decode()


def placed_cdxml(placed: list[Placed], svg_offset: Point, preset_name: str, name: str = "",
                 warn: Callable[[str], None] = _stderr) -> str:
    """CDXML for laid-out species (layout.m1_placed), in the frame of the SVG that to_svg wrote
    with `svg_offset`."""
    return to_cdxml([(p.mol, p.drawer, add(p.offset, svg_offset)) for p in placed], preset_name, name, warn)


def _line(runs: list[TextRun], parts: list[richtext.Run]) -> richtext.Line:
    """A richtext.Line for runs placed by the layout (pluses, numbers): positions from the runs."""
    return richtext.Line(runs, runs[0].x, runs[0].y, sum(r.advance() for r in runs),
                         richtext.has_sub(parts), "".join(p.text for p in parts), parts)


def scheme_records(lay: SchemeLayout, spec: Spec, svg_offset: Point = (0.0, 0.0)) -> CdxScheme:
    """The CdxScheme of a laid-out scheme, in the frame of the SVG that to_svg wrote with
    `svg_offset`. No layout here: positions are the layout's records."""
    names = list(lay.species)
    out = CdxScheme([(p.mol, p.drawer, add(p.offset, svg_offset)) for p in lay.species.values()], offset=svg_offset)
    for run in lay.pluses:
        out.texts.append(CdxText([_line([run], [richtext.Run("+")])]))
    by_n = {st.n: st for st in lay.steps}
    by_res = {res.n: res for res in lay.resonance}
    for i, st in enumerate(spec.steps, 1):
        if i == 1 or getattr(st, "new_row", False):
            out.schemes.append([])
        if isinstance(st, ResonanceStep):  # a <step> per ↔, left to right; the brackets free
            res = by_res[i]
            for k, (tail, head, _) in enumerate(res.arrows):
                out.arrows.append(CdxArrow(tail, head, "resonance", spec.arrow_head))
                out.schemes[-1].append(CdxStep([names.index(res.contributors[k])],
                                               [names.index(res.contributors[k + 1])], len(out.arrows) - 1))
            out.brackets += [CdxBracket(x, *res.bar_y, side) for x, side in zip(res.bar_x, ("left", "right"))]
            continue
        sl = by_n[i]
        out.arrows.append(CdxArrow(sl.tail, sl.head, sl.kind, sl.variant))
        sides = []
        for lines in (sl.above, sl.below):
            k0 = len(out.texts)
            if sl.staged:  # one Left caption per stage, as ChemDraw's "1. ..." / "2. ..."
                out.texts += [CdxText([ln]) for ln in lines]
            elif lines:
                out.texts.append(CdxText(list(lines), "Center"))
            sides.append(list(range(k0, len(out.texts))))
        over = spec.over_names(i)  # drawn above the arrow: cited there, not as reactants (user decision 2026-10-28)
        src = [names.index(n) for n in st.sources if n in lay.species and n not in over]
        dst = [names.index(n) for n in st.drawn_products]
        if sl.kind == "equilibrium_reactants":  # Head3D on the left: ChemDraw reads the step right to left
            src, dst = dst, src
        out.schemes[-1].append(CdxStep(src, dst, len(out.arrows) - 1, *sides,
                                       above_species=[names.index(n) for n in over if n in lay.species]))
    for n, runs in lay.numbers.items():
        out.texts.append(CdxText([_line(runs, richtext.parse(str(spec.numbers[n]), "bold"))]))
    out.curves = [CdxCurve(m.p0, m.c1, m.c2, m.p3, m.spec.type) for m in lay.curved]
    out.warn_maps = not lay.curved
    return out


def spec_cdxml(spec_path, warn: Callable[[str], None] = _stderr) -> str:
    """What render.py --cdxml writes for a spec: M1 species, or an M2 scheme or M3 mechanism."""
    spec = load_spec(spec_path)
    if spec.steps:
        lay = build_scheme(spec, spec_path.parent, spec_path.stem)
        _, offset = to_svg(lay.root)
        return to_cdxml(scheme_records(lay, spec, offset), spec.preset, f"{spec_path.stem}.cdxml", warn)
    placed = m1_placed(spec, spec_path.parent)
    _, offset = to_svg(Group(spec_path.stem, [p.group for p in placed]))
    return placed_cdxml(placed, offset, spec.preset, f"{spec_path.stem}.cdxml", warn)


def z_order(mol, d) -> dict[int, int]:
    """Bond idx -> its rank in Z (acs-style.md §14.2): bond order, except that a bond drawn behind
    another where they cross (draw_mol.Drawer.crossings) comes before it, so ChemDraw, which cuts
    the lower-Z bond of a crossing (chair_x*), draws the same gap."""
    after: dict[int, set[int]] = {b.idx: set() for b in mol.bonds}  # rear -> fronts
    need = {b.idx: 0 for b in mol.bonds}
    for f, r in d.crossings():
        if f not in after[r]:
            after[r].add(f)
            need[f] += 1
    ready = sorted(i for i, n in need.items() if n == 0)
    out = []
    while ready:
        i = ready.pop(0)
        out.append(i)
        for f in sorted(after[i]):
            need[f] -= 1
            if need[f] == 0:
                ready.append(f)
                ready.sort()
    out += [i for i in sorted(need) if i not in out]  # cycles (can't happen with a strict order)
    return {i: k for k, i in enumerate(out)}
