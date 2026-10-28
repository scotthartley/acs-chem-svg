"""Molecules -> CDXML, so a figure can be edited in ChemDraw (acs-style.md §13). Phase 1: M1
species only.

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

Not written in phase 1: abbreviations and Lewis structures (CdxmlUnsupported), lone-pair and
radical dots (the node keeps Radical), circled signs and charge_labels texts (a plain sign with
the formal charge), atom maps. Each of those warns.

Fragment and root BoundingBox: the drawn ink box (Drawer.draw: stroke widths, second lines,
charge signs) with the stored label boxes, as ChemDraw's. Coordinates are in the caller's
frame: each species' offset is added to its atom positions (render.py passes the layout offset
plus to_svg's, so the file lines up with the SVG).
Numbers are written with two decimals, as ChemDraw writes them. CreationProgram names this
package, and the page grows by whole pages to fit the drawing.
"""

from __future__ import annotations

import math
import sys
from collections.abc import Callable
from importlib.metadata import version
from dataclasses import dataclass

from lxml import etree
from rdkit import Chem

from .chem import Molecule
from .draw_mol import Drawer
from .geometry import Point, add, union_box
from .layout import Placed, m1_placed
from .preset import PRESET_DIR, load_preset
from .scene import Group
from .spec import load_spec
from .svg_writer import to_svg

HEADER_DIR = PRESET_DIR / "cdxml"
_PT = Chem.GetPeriodicTable()


class CdxmlUnsupported(ValueError):
    """The species can't be written as correct chemistry in this phase."""


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

    def check(self, mol: Molecule, d: Drawer) -> None:
        name, msgs = mol.name, []
        if mol.lewis:
            raise CdxmlUnsupported(f"{name}: Lewis structures are not written to CDXML (phase 1)")
        if any(a.abbrev for a in mol.atoms):
            raise CdxmlUnsupported(f"{name}: abbreviations (nicknames) are not written to CDXML yet (phase 1)")
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


def _stderr(msg: str) -> None:
    print(f"warning: {msg}", file=sys.stderr)


def to_cdxml(species: list[tuple[Molecule, Drawer, Point]], preset_name: str, name: str = "",
             warn: Callable[[str], None] = _stderr) -> str:
    """CDXML text for `species`: (molecule, its Drawer, the offset added to its coordinates).
    The Drawer's labels must be in molecule coordinates (Drawer.draw leaves them so)."""
    header = load_header(preset_name)
    w = _Writer(header, preset_name, warn)
    root = w.root
    page = root.find("page")
    page.set("id", w.ids.new())
    boxes = [w.fragment(page, mol, d, off) for mol, d, off in species]
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


def spec_cdxml(spec_path, warn: Callable[[str], None] = _stderr) -> str:
    """What render.py --cdxml writes for an M1 spec."""
    spec = load_spec(spec_path)
    if spec.steps:
        raise CdxmlUnsupported("CDXML export covers M1 species only (phase 1)")
    placed = m1_placed(spec, spec_path.parent)
    _, offset = to_svg(Group(spec_path.stem, [p.group for p in placed]))
    return placed_cdxml(placed, offset, spec.preset, f"{spec_path.stem}.cdxml", warn)
