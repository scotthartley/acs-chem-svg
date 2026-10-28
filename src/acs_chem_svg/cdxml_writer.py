"""Molecules and schemes -> CDXML, so a figure can be edited in ChemDraw (acs-style.md §13).
Phase 1: M1 species. Phase 2 (acs-style.md §13.1): M2 schemes, i.e. reaction and equilibrium
arrows, conditions, plus signs and compound numbers as captions, and <scheme><step>.

The header (root attributes, color and font tables, page) is presets/cdxml/<preset>.xml,
copied from the ChemDraw references by scripts/cdxml_headers.py; this module types no style
numbers. Per species, one <fragment>:
- <n> per atom at its drawn position: Element (atomic number; carbon has none), and on a labeled
  atom NumHydrogens and NeedsClean="yes", as ChemDraw writes its own labeled atoms; Charge;
  Radical="Doublet".
- A labeled atom's <t>: the text in CDXML order (NH2 even when drawn H2N), one <s> in the
  header's label font, size and face. p: x is the element symbol's glyph origin, or for a
  right-justified label (H on the left) the end of its advance (labels_x1.0, clearance_x*: -Cl
  sits 1.11 right of the atom); y is the label's center line plus half the ink height of the
  element's first glyph, as ChemDraw stores it (all references: 0.3684 x size for O, C, S, Cl;
  0.3586 for N, H, F, I, Br), though it draws every label's baseline at 0.3585 x size.
  BoundingBox: the union of the drawn runs' ink boxes, moved by the same stored - drawn offset.
  LabelJustification / LabelAlignment Left or Right from Label.h_side, plus
  Justification="Right"; a one-character label (-F, -H) is always Left, and only labels of more
  than one character get LabelAlignment, as in the references.
- <b>: Order 2 / 3; wedges Display WedgeBegin / WedgedHashBegin with B at the narrow end
  (Bond.a). No DoublePosition: ChemDraw places double bonds itself.
- A charge: <graphic GraphicType="Symbol" SymbolType="Plus|Minus"> with
  <represent attribute="Charge">, its BoundingBox the sign's center (Drawer.charge_center) and
  (center x - cdxml.charge_box_dx, center y).

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
- An arrow: <graphic GraphicType="Line" SupersededBy> (ArrowType FullHead or Equilibrium,
  BoundingBox Head3D then Tail3D) plus the <arrow> it points to, with the same Z. Head
  attributes per head variant from cdxml.arrow_heads; equilibrium arrows HalfLeft at both ends,
  ArrowShaftSpacing = equilibrium.shaft_spacing and ArrowEquilibriumRatio =
  ratios.equilibrium_long_per_short, x 100 / LineWidth and x 100. Head3D / Tail3D: a forward
  arrow's tip and tail; a balanced equilibrium arrow's ends moved in by equilibrium.stagger; one
  favoring products, its ends; one favoring reactants, reversed (Head3D on the left, as ChemDraw
  stores it, so the step's reactants are the species on its right). Center3D is the midpoint,
  MajorAxisEnd3D = Center3D + (Head3D - Tail3D), MinorAxisEnd3D = Center3D + that turned 90°.
  BoundingBox: along the axis the drawn ends; across it the half width (ArrowheadWidth x
  LineWidth / 100, plus half the shaft spacing), with cdxml.arrow_box_pad more on the top of a
  forward arrow and on both sides of an equilibrium arrow whose Head3D is on the right.
- <scheme> per drawn sequence (a new_row step starts one; the rows of a wrapped sequence stay
  one), <step> with ReactionStepReactants / Products (fragment ids), ReactionStepArrows (the
  <graphic> id, as ChemDraw cites it) and ReactionStepObjectsAbove / BelowArrow (<t> ids).
Page order as in the references: fragments, captions, graphics, schemes, arrows.

Not written: abbreviations and Lewis structures (CdxmlUnsupported), curved arrows, resonance
steps and brackets (CdxmlUnsupported, phases 3-4), lone-pair and radical dots (the node keeps
Radical), circled signs and charge_labels texts (a plain sign with the formal charge), atom
maps. Each of the last ones warns.

Fragment BoundingBox: the drawn ink box (Drawer.draw: stroke widths, second lines, charge signs)
with the stored label boxes, as ChemDraw's; the root's also covers the captions and arrows.
Coordinates are in the caller's frame: each species' offset is added to its atom positions, and
CdxScheme.offset to the arrows and captions (render.py passes to_svg's offset, so the file
lines up with the SVG).
Numbers are written with two decimals, as ChemDraw writes them. CreationProgram names this
package, and the page grows by whole pages to fit the drawing.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Callable
from importlib.metadata import version
from dataclasses import dataclass, field

from lxml import etree
from rdkit import Chem

from . import richtext
from .chem import Molecule
from .draw_mol import Drawer
from .geometry import Point, add, mul, sub, union_box, unit
from .layout import Placed, SchemeLayout, build_scheme, m1_placed
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
    kind: str = "forward"  # forward | equilibrium | equilibrium_products | equilibrium_reactants
    variant: str | None = None  # arrow.heads variant; None: the preset's default_head


@dataclass
class CdxText:
    """A caption: laid-out lines (richtext.Line with its parts and runs) and its justification."""
    lines: list[richtext.Line]
    justify: str = "Left"  # Left | Center


@dataclass
class CdxStep:
    """Indices into CdxScheme.species, .arrows and .texts."""
    reactants: list[int]
    products: list[int]
    arrow: int
    above: list[int] = field(default_factory=list)
    below: list[int] = field(default_factory=list)


@dataclass
class CdxScheme:
    species: list[tuple[Molecule, Drawer, Point]]  # (molecule, its Drawer, offset added to its atoms)
    arrows: list[CdxArrow] = field(default_factory=list)
    texts: list[CdxText] = field(default_factory=list)
    schemes: list[list[CdxStep]] = field(default_factory=list)  # one <scheme> per drawn sequence
    offset: Point = (0.0, 0.0)  # added to the arrows and captions


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

    def check(self, mol: Molecule, d: Drawer) -> None:
        name, msgs = mol.name, []
        if mol.lewis:
            raise CdxmlUnsupported(f"{name}: Lewis structures are not written to CDXML (phases 1-2)")
        if any(a.abbrev for a in mol.atoms):
            raise CdxmlUnsupported(f"{name}: abbreviations (nicknames) are not written to CDXML yet (phases 1-2)")
        dots = d.dot_sets()
        if any(ds.kind == "pair" for sets in dots.values() for ds in sets):
            msgs.append(f"{name}: lone-pair dots omitted from CDXML")
        if any(a.radicals for a in mol.atoms):
            msgs.append(f"{name}: radical dots omitted from CDXML (the atom keeps Radical)")
        if any(a.radicals > 1 for a in mol.atoms):
            msgs.append(f"{name}: more than one unpaired electron on an atom is not written")
        for a in mol.atoms:
            sign = d.charge_sign(a.idx)
            if sign is None:
                continue
            if a.charge_label is not None:
                msgs.append(f"{name}: charge label {a.charge_label!r} written as the formal charge {a.charge:+d}")
            elif sign.circled:
                msgs.append(f"{name}: circled charge written as a plain sign")
            elif sign.digits:
                msgs.append(f"{name}: charge {a.charge:+d} written as a plain sign with Charge={a.charge}")
        if any(a.map for a in mol.atoms):
            msgs.append(f"{name}: atom maps not written")
        for m in dict.fromkeys(msgs):  # once per species
            self.warn(m)

    def fragment(self, page, mol: Molecule, d: Drawer, off: Point) -> tuple:
        self.check(mol, d)
        frag = etree.SubElement(page, "fragment", {"id": self.ids.new()})
        node_ids, label_boxes = {}, []
        for a in mol.atoms:
            x, y = a.x + off[0], a.y + off[1]
            node_ids[a.idx] = nid = self.ids.new()
            attrs = {"id": nid, "p": pts(x, y), "Z": self.ids.next_z()}
            if a.element != "C":
                attrs["Element"] = str(_PT.GetAtomicNumber(a.element))
            lab = d.labels.get(a.idx)
            if lab is not None:
                attrs["NumHydrogens"] = str(a.num_h)
            if a.charge:
                attrs["Charge"] = str(a.charge)
            if a.radicals == 1:
                attrs["Radical"] = "Doublet"
            if lab is not None:
                attrs["NeedsClean"] = "yes"
            n = etree.SubElement(frag, "n", attrs)
            if lab is not None:
                label_boxes.append(self.label(n, lab, a.element, a.num_h, off))
        for b in mol.bonds:
            attrs = {"id": self.ids.new(), "Z": self.ids.next_z(), "B": node_ids[b.a], "E": node_ids[b.b]}
            if b.order != 1:
                attrs["Order"] = str(b.order)
            if b.stereo:
                attrs["Display"] = {"wedge": "WedgeBegin", "hash": "WedgedHashBegin"}[b.stereo]
            etree.SubElement(frag, "b", attrs)
        for a in mol.atoms:
            sign = d.charge_sign(a.idx)
            if sign is None:
                continue
            cx, cy = d.charge_center(a.idx)
            cx, cy = cx + off[0], cy + off[1]
            g = etree.SubElement(frag, "graphic", {
                "id": self.ids.new(), "BoundingBox": pts(cx, cy, cx - self.preset.cdxml.charge_box_dx, cy),
                "Z": self.ids.next_z(), "GraphicType": "Symbol", "SymbolType": "Plus" if sign.sign > 0 else "Minus"})
            etree.SubElement(g, "represent", {"attribute": "Charge", "object": node_ids[a.idx]})
        # ChemDraw's fragment box: the drawn ink (strokes, second lines, signs) and the stored label boxes
        drawn = Group("ink", [g for g in d.draw(mol.name).children if not g.id.endswith("-labels")])
        b = drawn.ink_bbox()
        box = union_box(([(b[0] + off[0], b[1] + off[1], b[2] + off[0], b[3] + off[1])] if b else []) + label_boxes)
        frag.set("BoundingBox", pts(*box))
        frag.set("Z", self.ids.next_z())
        return box

    def label(self, n, lab, element: str, num_h: int, off: Point) -> tuple:
        text = label_text(element, num_h)
        # bonds to the right: right-justified (H2N, HO, -Cl); one glyph (-F, -H) stays Left
        right = lab.h_side == "left" and len(text) > 1
        x, base = lab.elem_origin
        if right:
            x += lab.face.advance(element, lab.size)
        # ChemDraw stores y half the first glyph's ink height below the label's center line
        # (atom y, or where center_terminal_labels moved it), not the drawn baseline: 3.68 for O
        center = base - self.preset.ratios.label_baseline_per_size * lab.size
        y = center + lab.face.text_ink_bounds(element[0], lab.size)[3] / 2
        box = union_box(r.bbox() for r in lab.runs)  # in the stored frame too: moved by y - base
        box = (box[0] + off[0], box[1] + y - base + off[1], box[2] + off[0], box[3] + y - base + off[1])
        just = "Right" if right else "Left"
        attrs = {"p": pts(x + off[0], y + off[1]), "BoundingBox": pts(*box), "LabelJustification": just}
        if right:
            attrs["Justification"] = "Right"
        if len(text) > 1:
            attrs["LabelAlignment"] = just
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
        if a.kind in ("forward", "equilibrium_products"):
            return a.tail, a.head
        raise ValueError(f"unknown arrow kind {a.kind!r}")

    def arrow_attrs(self, a: CdxArrow) -> dict:
        """The <arrow> head attributes, in ChemDraw's order."""
        variant = a.variant or self.preset.arrow.default_head
        if variant not in self.preset.cdxml.arrow_heads:
            raise ValueError(f"no CDXML head attributes for arrow head {variant!r} (preset cdxml.arrow_heads)")
        size, center, width = self.preset.cdxml.arrow_heads[variant]
        eq = a.kind != "forward"
        attrs = {"FillType": "None", "ArrowheadHead": "HalfLeft" if eq else "Full"}
        if eq:
            attrs["ArrowheadTail"] = "HalfLeft"
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
            "ArrowType": "FullHead" if a.kind == "forward" else "Equilibrium", "HeadSize": attrs["HeadSize"]})
        c = mul(add(head3d, tail3d), 0.5)
        d = sub(head3d, tail3d)
        el = etree.Element("arrow", {"id": aid, "BoundingBox": pts(*box), "Z": z, **attrs,
                                     "Head3D": _p3(head3d), "Tail3D": _p3(tail3d), "Center3D": _p3(c),
                                     "MajorAxisEnd3D": _p3(add(c, d)), "MinorAxisEnd3D": _p3(add(c, (-d[1], d[0])))})
        return gid, aid, box, el


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
    root = w.root
    page = root.find("page")
    page.set("id", w.ids.new())
    frag_ids, boxes = [], []
    for mol, d, off in scheme.species:
        boxes.append(w.fragment(page, mol, d, off))
        frag_ids.append(page[-1].get("id"))
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
    for steps in scheme.schemes:
        sch = etree.SubElement(page, "scheme", {"id": w.ids.new()})
        for st in steps:
            attrs = {"id": w.ids.new()}
            for key, ids in (("ReactionStepReactants", [frag_ids[i] for i in st.reactants]),
                             ("ReactionStepProducts", [frag_ids[i] for i in st.products]),
                             ("ReactionStepArrows", [graphic_ids[st.arrow]]),
                             ("ReactionStepObjectsAboveArrow", [j for i in st.above for j in text_ids[i]]),
                             ("ReactionStepObjectsBelowArrow", [j for i in st.below for j in text_ids[i]])):
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


def check_scheme(spec: Spec) -> None:
    """Refuse what phase 2 doesn't write: resonance steps (with their brackets) and curved arrows."""
    for i, st in enumerate(spec.steps or [], 1):
        if isinstance(st, ResonanceStep):
            raise CdxmlUnsupported(f"step {i}: resonance steps and brackets are not written to CDXML yet (phase 2)")
        if st.arrows:
            raise CdxmlUnsupported(f"step {i}: curved arrows are not written to CDXML yet (phase 2)")


def _line(runs: list[TextRun], parts: list[richtext.Run]) -> richtext.Line:
    """A richtext.Line for runs placed by the layout (pluses, numbers): positions from the runs."""
    return richtext.Line(runs, runs[0].x, runs[0].y, sum(r.advance() for r in runs),
                         richtext.has_sub(parts), "".join(p.text for p in parts), parts)


def scheme_records(lay: SchemeLayout, spec: Spec, svg_offset: Point = (0.0, 0.0)) -> CdxScheme:
    """The CdxScheme of a laid-out scheme, in the frame of the SVG that to_svg wrote with
    `svg_offset`. No layout here: positions are the layout's records."""
    check_scheme(spec)
    if lay.resonance or lay.curved:
        raise CdxmlUnsupported("resonance brackets and curved arrows are not written to CDXML yet (phase 2)")
    names = list(lay.species)
    out = CdxScheme([(p.mol, p.drawer, add(p.offset, svg_offset)) for p in lay.species.values()], offset=svg_offset)
    for run in lay.pluses:
        out.texts.append(CdxText([_line([run], [richtext.Run("+")])]))
    by_n = {st.n: st for st in lay.steps}
    for i, st in enumerate(spec.steps, 1):
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
        src = [names.index(n) for n in st.sources if n in lay.species]
        dst = [names.index(n) for n in st.drawn_products]
        if sl.kind == "equilibrium_reactants":  # Head3D on the left: ChemDraw reads the step right to left
            src, dst = dst, src
        if i == 1 or st.new_row:
            out.schemes.append([])
        out.schemes[-1].append(CdxStep(src, dst, len(out.arrows) - 1, *sides))
    for n, runs in lay.numbers.items():
        out.texts.append(CdxText([_line(runs, richtext.parse(str(spec.numbers[n]), "bold"))]))
    return out


def spec_cdxml(spec_path, warn: Callable[[str], None] = _stderr) -> str:
    """What render.py --cdxml writes for a spec: M1 species, or an M2 scheme."""
    spec = load_spec(spec_path)
    if spec.steps:
        check_scheme(spec)
        lay = build_scheme(spec, spec_path.parent, spec_path.stem)
        _, offset = to_svg(lay.root)
        return to_cdxml(scheme_records(lay, spec, offset), spec.preset, f"{spec_path.stem}.cdxml", warn)
    placed = m1_placed(spec, spec_path.parent)
    _, offset = to_svg(Group(spec_path.stem, [p.group for p in placed]))
    return placed_cdxml(placed, offset, spec.preset, f"{spec_path.stem}.cdxml", warn)
