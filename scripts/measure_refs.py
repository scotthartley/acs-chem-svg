"""Measure ChemDraw reference renders into numbers for references/acs-style.md.

Atom positions, bond orders, stereo flags, and object types come from the
.cdxml saved next to each export (exact, no clustering heuristics). Drawn
geometry comes from the .svg, mapped into the CDXML frame by the path
transform (svg = cdxml + (e, f)). SVG user units are checked against the PDF
page size, so every number is in pt.

ChemDraw draws every bond line, wedge, hash, charge sign, and arrow shaft as a
filled polygon; lone pairs as filled circles; curved arrows as stroked cubics.

Usage:
    uv run python scripts/measure_refs.py                 # all references
    uv run python scripts/measure_refs.py m1_set_x1.0     # one file
    uv run python scripts/measure_refs.py --json out.json
    uv run python scripts/measure_refs.py --scheme        # scheme parameters per preset (§3)
    uv run python scripts/measure_refs.py --mech          # mechanism parameters per preset (§4)
    uv run python scripts/measure_refs.py --resonance     # resonance arrows and brackets per preset (§4.3)
    uv run python scripts/measure_refs.py --wrap          # multi-row scheme parameters per preset (§3.4)
    uv run python scripts/measure_refs.py --abbrev        # abbreviation labels per preset (§3.5)
    uv run python scripts/measure_refs.py --equilibrium   # equilibrium arrows, head variants (§3.1, §3.6)
    uv run python scripts/measure_refs.py --cdxml-symbols # charge-symbol <graphic> boxes (§13, CDXML export)
    uv run python scripts/measure_refs.py --cdxml-scheme  # CDXML captions, arrows and steps (§13.1)
    uv run python scripts/measure_refs.py --charges-lewis # circled signs, symbols, two-electron radicals (§2.2, §13.4)
    uv run python scripts/measure_refs.py --chair         # chair templates, axial / equatorial bonds (§14.1)
    uv run python scripts/measure_refs.py --crossing      # bond crossings: front by Z, rear cuts (§14.2)
"""

from __future__ import annotations

import argparse
import copy
import itertools
import json
import math
import re
import statistics
import sys
from dataclasses import dataclass, field
from pathlib import Path

import pypdfium2
from lxml import etree

from acs_chem_svg.fontmetrics import load_face
from acs_chem_svg.geometry import (  # shared with the renderer
    add, convex_hull, dot, length, mul, ray_box_gap, ray_polygon, seg_point_dist, seg_seg_dist, sub, unit,
)

REPO = Path(__file__).resolve().parents[1]
REF_DIR = REPO / "reference_renders" / "chemdraw"
ELEMENTS = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 14: "Si", 15: "P", 16: "S", 17: "Cl", 35: "Br", 53: "I"}


def stats(values):
    vals = [v for v in values if v is not None]
    if not vals:
        return None
    return {"median": round(statistics.median(vals), 3), "min": round(min(vals), 3),
            "max": round(max(vals), 3), "n": len(vals)}


# --------------------------------------------------------------------------- model


@dataclass
class Atom:
    id: str
    p: tuple[float, float]
    element: str = "C"
    label: str | None = None
    charge: int = 0
    num_h: int | None = None
    frag: str = ""
    nickname: bool = False
    radicals: int = 0  # CDXML Radical="Doublet"
    justification: str = ""  # CDXML LabelJustification of the node's <t> (Left / Right; nicknames: MeO is Right)


@dataclass
class Bond:
    id: str
    a: str
    b: str
    order: int = 1
    display: str = ""
    ring: bool = False
    frag: str = ""


@dataclass
class Shape:
    kind: str  # poly | circle | stroke | outline
    pts: list
    curved: bool = False
    width: float = 0.0  # stroke width (stroke/outline)
    center: tuple | None = None
    r: float = 0.0
    owner: str | None = None  # what claimed it
    group: str = ""  # our SVG (load_rendered): the id of the <g> it is drawn in

    @property
    def bbox(self):
        xs, ys = [p[0] for p in self.pts], [p[1] for p in self.pts]
        if self.kind == "circle":
            return (self.center[0] - self.r, self.center[1] - self.r,
                    self.center[0] + self.r, self.center[1] + self.r)
        return (min(xs), min(ys), max(xs), max(ys))

    @property
    def centroid(self):
        if self.kind == "circle":
            return self.center
        return (statistics.fmean(p[0] for p in self.pts), statistics.fmean(p[1] for p in self.pts))


@dataclass
class Glyph:
    ch: str
    x: float  # origin
    base: float  # baseline y (y down)
    size: float
    bold: bool
    run: int
    box: tuple | None = None  # ink box, y down
    adv: float = 0.0
    italic: bool = False

    def outline(self):
        face = load_face("Bold" if self.bold else "Oblique" if self.italic else "Regular")
        return [[(self.x + x, self.base - y) for x, y in c] for c in face.outline(self.ch, self.size)]


@dataclass
class Doc:
    name: str
    atoms: dict = field(default_factory=dict)
    bonds: list = field(default_factory=list)
    shapes: list = field(default_factory=list)
    glyphs: list = field(default_factory=list)
    runs: list = field(default_factory=list)  # (text, x, base, size, bold)
    arrows: list = field(default_factory=list)  # dicts in svg frame
    curves: list = field(default_factory=list)
    captions: list = field(default_factory=list)
    graphics: list = field(default_factory=list)
    settings: dict = field(default_factory=dict)
    page_pt: tuple | None = None
    svg_size: tuple | None = None
    offset: tuple = (0.0, 0.0)  # svg = cdxml + offset
    fragments: dict = field(default_factory=dict)  # frag id -> bbox (cdxml)
    text_frags: set = field(default_factory=set)  # formula fragments (K2CO3, LiAlH4): text, no atoms


# --------------------------------------------------------------------------- loading

_MATRIX = re.compile(r"matrix\(([^)]*)\)")
_TOK = re.compile(r"[MLCAZ]|-?\d*\.?\d+(?:e-?\d+)?")


def _matrix(el):
    m = _MATRIX.search(el.get("transform", ""))
    if not m:
        return (1, 0, 0, 1, 0, 0)
    return tuple(float(v) for v in m.group(1).replace(",", " ").split())


def _apply(m, x, y):
    a, b, c, d, e, f = m
    return (a * x + c * y + e, b * x + d * y + f)


def parse_path(d: str, m) -> Shape:
    toks = _TOK.findall(d)
    pts, curved, arcs = [], False, []
    i, cmd = 0, None
    while i < len(toks):
        t = toks[i]
        if t in "MLCAZ":
            cmd = t
            i += 1
            continue
        if cmd in ("M", "L"):
            pts.append(_apply(m, float(toks[i]), float(toks[i + 1])))
            i += 2
        elif cmd == "C":
            curved = True
            for k in range(3):
                pts.append(_apply(m, float(toks[i + 2 * k]), float(toks[i + 2 * k + 1])))
            i += 6
        elif cmd == "A":
            rx = float(toks[i])
            end = _apply(m, float(toks[i + 5]), float(toks[i + 6]))
            arcs.append((rx * m[0], end))
            pts.append(end)
            i += 7
        else:
            i += 1
    if arcs:  # ChemDraw dots: two semicircles
        p0, p1 = pts[0], arcs[0][1]
        return Shape("circle", pts, center=mul(add(p0, p1), 0.5), r=arcs[0][0])
    # drop the closing duplicate vertex
    while len(pts) > 1 and length(sub(pts[0], pts[-1])) < 1e-6:
        pts.pop()
    return Shape("poly", pts, curved=curved)


def load_svg(doc: Doc, path: Path, faces):
    root = etree.parse(str(path)).getroot()
    ns = "{http://www.w3.org/2000/svg}"
    vb = [float(v) for v in root.get("viewBox").split()]
    doc.svg_size = (vb[2], vb[3])
    offset = None
    for el in root.iter(ns + "path", ns + "text"):
        m = _matrix(el)
        if el.tag == ns + "path":
            offset = offset or (m[4], m[5])
            sh = parse_path(el.get("d"), m)
            stroke = el.get("stroke", "none")
            if stroke != "none" and el.get("fill") == "none":
                w = float(el.get("stroke-width", "1")) * m[0]
                sh.kind = "stroke" if w >= 0.2 else "outline"
                sh.width = w
            doc.shapes.append(sh)
        else:
            _add_text(doc, el, m, faces)
    return offset or (0.0, 0.0)


def text_faces() -> dict:
    return {k: load_face(k) for k in ("Regular", "Bold", "Oblique")}


def _add_text(doc: Doc, el, m, faces):
    """Glyphs of one <text> run, positioned by its transform (ChemDraw) or x/y (ours).
    Italic runs: ChemDraw writes font-family="Helvetica Oblique" with font-style="italic"."""
    text = (el.text or "").strip("\n")
    size = float(el.get("font-size").rstrip("px")) * m[0]
    is_bold = el.get("font-weight") == "bold"
    italic = el.get("font-style") == "italic" or "Oblique" in el.get("font-family", "")
    face = faces["Bold" if is_bold else "Oblique" if italic else "Regular"]
    x, base = m[4] + float(el.get("x", "0")), m[5] + float(el.get("y", "0"))
    ri = len(doc.runs)
    doc.runs.append((text, x, base, size, is_bold))
    for k, ch in enumerate(text):
        g = Glyph(ch, x, base, size, is_bold, ri, italic=italic)
        if ch.strip():
            b = face.ink_bounds(ch, size)
            if b:
                g.box = (x + b[0], base - b[3], x + b[2], base - b[1])
        g.adv = face.advance(ch, size)
        doc.glyphs.append(g)
        x += g.adv
        if k + 1 < len(text):
            x += face.kerning(ch, text[k + 1], size)


def stroke_rects(pts, width, closed=False) -> list[Shape]:
    """A stroked polyline as one butt-capped rectangle per segment (ChemDraw-comparable polys)."""
    out = []
    segs = list(zip(pts, pts[1:])) + ([(pts[-1], pts[0])] if closed and len(pts) > 2 else [])
    for p, q in segs:
        u = unit(sub(q, p))
        n = mul((-u[1], u[0]), width / 2)
        out.append(Shape("poly", [add(p, n), add(q, n), sub(q, n), sub(p, n)]))
    return out


_OWN_PATH = re.compile(r"([MLZ])\s*(-?[\d.]+)?\s+(-?[\d.]+)?")


def load_rendered(svg_text: str, name: str, atoms: dict, bonds: list, settings: dict,
                  arrows: list | None = None, curves: list | None = None) -> Doc:
    """Our own SVG (svg_writer: M/L/C/Z paths, x/y text, no transforms) as a Doc, with atoms,
    bonds and arrows (dicts as load_cdxml makes them) supplied by the caller. Stroked paths
    become butt-capped rectangles per segment, so assign() and the label measurements treat
    both renderers alike. Text outside species groups becomes captions, one per innermost
    group (conditions-step-<n>-above, number-<name>, plus-<i>), like ChemDraw's <t> objects.
    Circles (dots) and stroked cubics (curved arrows, matched to `curves` as load_cdxml makes
    them) load as ChemDraw's do, and bracket polylines as one stroked segment each."""
    doc = Doc(name, atoms=atoms, bonds=bonds, settings=settings, arrows=list(arrows or []),
              curves=list(curves or []))
    faces = text_faces()
    root = etree.fromstring(svg_text.encode())
    ns = "{http://www.w3.org/2000/svg}"
    vb = [float(v) for v in root.get("viewBox").split()]
    doc.svg_size = (vb[2], vb[3])
    captions: dict = {}
    for el in root.iter(ns + "path", ns + "text", ns + "circle"):
        if el.tag == ns + "circle":
            c = (float(el.get("cx")), float(el.get("cy")))
            doc.shapes.append(Shape("circle", [c], center=c, r=float(el.get("r"))))
            continue
        if el.tag == ns + "text":
            first = len(doc.glyphs)
            _add_text(doc, el, (1, 0, 0, 1, 0, 0), faces)
            gid = el.getparent().get("id", "")
            if not gid.startswith(("species-", "frag", "compare-")) and el.getparent() is not root:
                cap = captions.setdefault(gid, {"text": "", "glyphs": [], "bold": False, "formula": False})
                cap["glyphs"] += doc.glyphs[first:]
                cap["text"] += el.text or ""
                cap["bold"] |= el.get("font-weight") == "bold"
            continue
        d = el.get("d")
        if "C" in d and el.get("fill") == "none":  # curved-arrow stroke
            sh = parse_path(d, (1, 0, 0, 1, 0, 0))
            sh.kind, sh.width = "stroke", float(el.get("stroke-width"))
            doc.shapes.append(sh)
            continue
        if "C" in d:  # arrowheads: filled cubic outline, as ChemDraw's
            doc.shapes.append(parse_path(d.replace("Z", ""), (1, 0, 0, 1, 0, 0)))
            continue
        pts = [(float(x), float(y)) for c, x, y in _OWN_PATH.findall(d) if c in "ML"]
        if "bracket" in el.getparent().get("id", ""):  # one stroked segment each, as ChemDraw's
            doc.shapes += [Shape("stroke", [p, q], width=float(el.get("stroke-width"))) for p, q in zip(pts, pts[1:])]
            continue
        if el.get("fill", "none") == "none":
            doc.shapes += stroke_rects(pts, float(el.get("stroke-width")), closed="Z" in d)
        else:
            doc.shapes.append(Shape("poly", pts, group=el.getparent().get("id", "")))
    doc.captions = list(captions.values())
    _mark_rings(doc)
    return doc


def _xy(s, off):
    v = [float(t) for t in s.split()]
    return (v[0] + off[0], v[1] + off[1])


def load_cdxml(doc: Doc, path: Path, off):
    root = etree.parse(str(path)).getroot()
    for k in ("BondLength", "LineWidth", "BoldWidth", "MarginWidth", "HashSpacing", "BondSpacing",
              "LabelSize", "CaptionSize", "ChainAngle", "CreationProgram"):
        doc.settings[k] = root.get(k)
    for frag in root.iter("fragment"):
        if frag.getparent().tag == "n":
            continue  # nickname internals (e.g. CH3Cl), not drawn as atoms
        fid = frag.get("id")
        bb = [float(v) for v in frag.get("BoundingBox").split()]
        doc.fragments[fid] = (bb[0] + off[0], bb[1] + off[1], bb[2] + off[0], bb[3] + off[1])
        nodes = frag.findall("n")
        if (len(nodes) == 1 and not frag.findall("b") and nodes[0].find("t") is not None
                and nodes[0].get("NodeType")):
            # a formula typed as text (K2CO3, LiAlH4, PCC): ChemDraw parses it into a fragment,
            # but it is drawn and placed as a caption. Kept in doc.fragments so template numbering holds.
            # A lone atom (Br⁻, Cl· in mech_basic) has an Element and no NodeType: it is an atom.
            doc.text_frags.add(fid)
            _add_caption(doc, nodes[0].find("t"), off, formula=True)
            continue
        for n in nodes:
            t = n.find("t")
            label = "".join(s.text or "" for s in t.iter("s")) if t is not None else None
            el = ELEMENTS.get(int(n.get("Element", "6")), "?")
            nick = n.get("NodeType") in ("Fragment", "Nickname")
            if nick:  # the attachment atom of the nickname sits at the node (OMe: O)
                el = next((ELEMENTS.get(int(m.get("Element", "6")), "?") for m in n.iter("n")
                           if m is not n and m.get("p") == n.get("p")), "X")
            doc.atoms[n.get("id")] = Atom(
                n.get("id"), _xy(n.get("p"), off), el, label,
                int(n.get("Charge", "0")), int(n.get("NumHydrogens")) if n.get("NumHydrogens") else None,
                fid, nick, 1 if n.get("Radical") == "Doublet" else 0,
                t.get("LabelJustification", "") if t is not None else "")
        for b in frag.findall("b"):
            order = b.get("Order", "1")
            doc.bonds.append(Bond(b.get("id"), b.get("B"), b.get("E"),
                                  int(float(order)) if order.replace(".", "").isdigit() else 1,
                                  b.get("Display", ""), frag=fid))
    _mark_rings(doc)
    for a in root.iter("arrow"):
        doc.arrows.append({
            "head": _xy(a.get("Head3D").rsplit(" ", 1)[0], off),
            "tail": _xy(a.get("Tail3D").rsplit(" ", 1)[0], off),
            "both": a.get("ArrowheadTail") == "Full",
            # two shafts with a half head each (ChemDraw's equilibrium tools)
            "equilibrium": (a.get("ArrowheadHead") or "").startswith("Half") and bool(a.get("ArrowShaftSpacing")),
            "attrs": {k: a.get(k) for k in ("HeadSize", "ArrowheadCenterSize", "ArrowheadWidth",
                                             "ArrowheadHead", "ArrowheadTail", "ArrowheadType",
                                             "ArrowShaftSpacing", "ArrowEquilibriumRatio")},
        })
    for c in root.iter("curve"):
        v = [float(t) for t in c.get("CurvePoints").split()]
        pts = [(v[i] + off[0], v[i + 1] + off[1]) for i in range(0, len(v), 2)]
        doc.curves.append({"pts": pts, "head": c.get("ArrowheadHead"), "type": c.get("CurveType")})
    for g in root.iter("graphic"):
        if g.get("SupersededBy"):
            continue
        doc.graphics.append({k: g.get(k) for k in ("GraphicType", "SymbolType", "BracketType",
                                                    "LipSize", "BoundingBox")})
    for t in root.iter("t"):
        if t.getparent().tag == "n":
            continue
        _add_caption(doc, t, off)


def _add_caption(doc: Doc, t, off, formula: bool = False):
    """CDXML face bits: 1 bold, 2 italic; 96 = formula (ChemDraw subscripts the digits)."""
    runs = [(s.text or "", float(s.get("size")), s.get("face")) for s in t.iter("s")]
    bb = [float(v) for v in t.get("BoundingBox").split()]
    doc.captions.append({
        "text": "".join(r[0] for r in runs), "runs": runs, "p": _xy(t.get("p"), off),
        "bbox": (bb[0] + off[0], bb[1] + off[1], bb[2] + off[0], bb[3] + off[1]),
        "justification": t.get("Justification"),
        "bold": any(int(r[2] or 0) & 1 for r in runs),
        "formula": formula,
    })


def _mark_rings(doc: Doc):
    adj = {}
    for b in doc.bonds:
        adj.setdefault(b.a, set()).add(b.b)
        adj.setdefault(b.b, set()).add(b.a)
    for b in doc.bonds:
        seen, stack = {b.a}, [b.a]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if (u, w) in ((b.a, b.b), (b.b, b.a)) or w in seen:
                    continue
                seen.add(w)
                stack.append(w)
        b.ring = b.b in seen


def load(name: str, ref_dir: Path = REF_DIR) -> Doc:
    doc = Doc(name)
    off = load_svg(doc, ref_dir / f"{name}.svg", text_faces())
    doc.offset = off
    pdf = ref_dir / f"{name}.pdf"
    if pdf.exists():
        doc.page_pt = tuple(pypdfium2.PdfDocument(str(pdf))[0].get_size())
    load_cdxml(doc, ref_dir / f"{name}.cdxml", off)
    return doc


# --------------------------------------------------------------------------- bond frames


@dataclass
class Frame:
    o: tuple
    u: tuple
    n: tuple
    L: float

    def local(self, p):
        d = sub(p, self.o)
        return (dot(d, self.u), dot(d, self.n))

    def world(self, u, v):
        return add(self.o, add(mul(self.u, u), mul(self.n, v)))


def frame(doc, a_id, b_id):
    a, b = doc.atoms[a_id].p, doc.atoms[b_id].p
    u = unit(sub(b, a))
    return Frame(a, u, (-u[1], u[0]), length(sub(b, a)))


def neighbors(doc, atom_id):
    out = []
    for b in doc.bonds:
        if b.a == atom_id:
            out.append(b.b)
        elif b.b == atom_id:
            out.append(b.a)
    return out


# --------------------------------------------------------------------------- shape -> bond assignment


@dataclass
class Line:
    shape: Shape
    v: float  # centerline offset from bond axis (+ = left of a->b in y-down frame)
    width: float
    u0: float
    u1: float


def _line_from_shape(sh, fr) -> Line:
    loc = [fr.local(p) for p in sh.pts]
    vs = [p[1] for p in loc]
    us = [p[0] for p in loc]
    width = max(vs) - min(vs)
    tol = max(width * 1.5, 0.5)
    u0 = statistics.fmean(u for u in us if u <= min(us) + tol)
    u1 = statistics.fmean(u for u in us if u >= max(us) - tol)
    return Line(sh, (max(vs) + min(vs)) / 2, width, u0, u1)


def assign(doc: Doc):
    s_bl = float(doc.settings["BondLength"])
    spacing = s_bl * float(doc.settings["BondSpacing"]) / 100
    bold = float(doc.settings["BoldWidth"])
    polys = [s for s in doc.shapes if s.kind == "poly" and not s.curved]
    result = {}
    # stereo bonds first
    for b in doc.bonds:
        if "Wedge" not in b.display:
            continue
        a_id, e_id = (b.a, b.b) if b.display.endswith("Begin") else (b.b, b.a)
        fr = frame(doc, a_id, e_id)
        mine = []
        for sh in polys:
            if sh.owner:
                continue
            loc = [fr.local(p) for p in sh.pts]
            if all(-1.0 <= u <= fr.L + 1.0 and abs(v) <= bold for u, v in loc):
                mine.append((sh, loc))
        for sh, _ in mine:
            sh.owner = b.id
        result[b.id] = {"frame": fr, "stereo": b.display, "shapes": mine}
    # plain lines
    for sh in polys:
        if sh.owner:
            continue
        best = None
        for b in doc.bonds:
            if "Wedge" in b.display:
                continue
            fr = frame(doc, b.a, b.b)
            loc = [fr.local(p) for p in sh.pts]
            us, vs = [p[0] for p in loc], [p[1] for p in loc]
            # lines may run up to 2 pt past an atom: the bond beside a notched wedge end is cut
            # along the wedge's end line (stereo_x*: 1.5 pt)
            if min(us) < -2.0 or max(us) > fr.L + 2.0 or max(abs(v) for v in vs) > spacing * 1.3 + 0.5:
                continue
            if (max(us) - min(us)) < 3 * (max(vs) - min(vs)):
                continue
            score = abs(statistics.fmean(vs))
            if best is None or score < best[0]:
                best = (score, b, fr)
        if best:
            _, b, fr = best
            sh.owner = b.id
            result.setdefault(b.id, {"frame": fr, "lines": []}).setdefault("lines", []).append(
                _line_from_shape(sh, fr))
    for r in result.values():
        if "lines" in r:
            r["lines"].sort(key=lambda ln: ln.v)
    return result


# --------------------------------------------------------------------------- labels


def label_glyphs(doc: Doc, atom: Atom):
    """(element glyphs, whole-label glyphs) for a labeled atom. The element
    symbol is the matching glyph run closest to the atom; the whole label is
    the chain of glyphs set contiguously with it (origin of each = previous
    origin + advance), which picks up H, subscripts, and left-side H2N forms.
    A nickname's anchor is found by position, since ChemDraw displays a right-justified one
    reversed (CDXML OMe, drawn MeO: abbrev_x*): the glyph whose advance center is on the atom
    starts or ends the chain, and draw_mol.anchor_symbol picks the anchor glyphs from that end."""
    if not atom.label:
        return [], []
    if atom.nickname:
        return _nickname_glyphs(doc, atom)
    sym = atom.element
    best = None
    for i, g in enumerate(doc.glyphs):
        if not g.box or g.ch != sym[0]:
            continue
        if len(sym) == 2 and (i + 1 >= len(doc.glyphs) or doc.glyphs[i + 1].ch != sym[1]):
            continue
        d = length(sub(_box_center(g.box), atom.p))
        if best is None or d < best[0]:
            best = (d, i)
    if best is None or best[0] > 6:
        return [], []
    i = best[1]
    elem = doc.glyphs[i:i + len(sym)]
    return elem, [g for g in _glyph_chain(doc, elem) if g.box]


def _glyph_chain(doc: Doc, seed: list) -> list:
    """`seed` extended both ways by glyphs set contiguously with it on (about) its baseline."""
    whole = list(seed)
    tol = 0.05
    changed = True
    while changed:
        changed = False
        first, last = whole[0], whole[-1]
        for g in doc.glyphs:
            if g in whole or abs(g.base - first.base) > 0.35 * first.size:
                continue
            if abs(g.x - (last.x + last.adv)) < tol:
                whole.append(g)
                changed = True
                break
            if abs(g.x + g.adv - first.x) < tol:
                whole.insert(0, g)
                changed = True
                break
    return whole


def _nickname_glyphs(doc: Doc, atom: Atom):
    from acs_chem_svg.draw_mol import anchor_symbol

    best = None
    for g in doc.glyphs:
        if not g.box:
            continue
        d = length(sub((g.x + g.adv / 2, (g.box[1] + g.box[3]) / 2), atom.p))
        if best is None or d < best[0]:
            best = (d, g)
    if best is None or best[0] > 6:
        return [], []
    whole = _glyph_chain(doc, [best[1]])
    text = "".join(g.ch for g in whole)
    side = "left" if best[1] is whole[-1] and len(whole) > 1 else "right"
    k = len(anchor_symbol(text, atom.element, side))
    elem = whole[-k:] if side == "left" else whole[:k]
    return elem, [g for g in whole if g.box]


def _box_center(b):
    return ((b[0] + b[2]) / 2, (b[1] + b[3]) / 2)


def union(boxes):
    boxes = [b for b in boxes if b]
    return (min(b[0] for b in boxes), min(b[1] for b in boxes),
            max(b[2] for b in boxes), max(b[3] for b in boxes)) if boxes else None


# --------------------------------------------------------------------------- measurements


def measure(doc: Doc) -> dict:
    out: dict = {"file": doc.name, "settings": doc.settings}
    out["units"] = {"svg_viewbox": doc.svg_size, "pdf_page_pt": doc.page_pt,
                    "svg_unit_is_pt": bool(doc.page_pt) and all(
                        abs(a - b) < 0.01 for a, b in zip(doc.svg_size, doc.page_pt))}
    asg = assign(doc)
    bl = float(doc.settings["BondLength"])
    sp_nominal = bl * float(doc.settings["BondSpacing"]) / 100

    # bond lengths (CDXML atom coordinates)
    out["bond_length"] = stats(frame(doc, b.a, b.b).L for b in doc.bonds)

    # line widths of all plain lines
    widths = [ln.width for r in asg.values() for ln in r.get("lines", [])]
    out["line_width"] = stats(widths)

    # double / triple bonds
    ring_inner, ring_spacing, acyclic, triples = [], [], [], []
    for b in doc.bonds:
        r = asg.get(b.id)
        if not r or "lines" not in r:
            continue
        lines, fr = r["lines"], r["frame"]
        A, B = doc.atoms[b.a], doc.atoms[b.b]
        if b.order == 2 and len(lines) == 2:
            v = [ln.v for ln in lines]
            spacing = v[1] - v[0]
            if b.ring:
                main = min(lines, key=lambda ln: abs(ln.v))
                inner = max(lines, key=lambda ln: abs(ln.v))
                ring_spacing.append(spacing)
                ring_inner.append({
                    "bond": b.id, "spacing": round(spacing, 3),
                    "inset_a": round(inner.u0, 3), "inset_b": round(fr.L - inner.u1, 3),
                    "main_u": (round(main.u0, 3), round(main.u1, 3)), "L": round(fr.L, 3),
                    "a_labeled": bool(A.label), "b_labeled": bool(B.label),
                })
            else:
                acyclic.append(_acyclic_entry(doc, b, fr, lines, spacing))
        elif b.order == 3 and len(lines) == 3:
            triples.append({"bond": b.id, "v": [round(ln.v, 3) for ln in lines],
                            "u": [(round(ln.u0, 3), round(fr.L - ln.u1, 3)) for ln in lines]})
    out["bond_spacing_ring"] = stats(ring_spacing)
    out["bond_spacing_nominal"] = round(sp_nominal, 4)
    out["ring_inner_line"] = {
        "inset_unlabeled_ends": stats([e["inset_a"] for e in ring_inner if not e["a_labeled"]]
                                      + [e["inset_b"] for e in ring_inner if not e["b_labeled"]]),
        "inset_formula_s_tan30": round(sp_nominal * math.tan(math.radians(30)), 3),
        "entries": ring_inner,
    }
    out["acyclic_double"] = acyclic
    out["triple"] = triples

    # stereo
    out["stereo"] = [_stereo_entry(doc, b, asg[b.id]) for b in doc.bonds if b.id in asg and "stereo" in asg[b.id]]

    # labels
    out["labels"] = [e for a in doc.atoms.values() if a.label for e in [_label_entry(doc, a, asg)] if e]
    out["charges"] = _charges(doc)
    out["dots"] = _dots(doc)
    out["arrows"] = [_arrow_entry(doc, a) for a in doc.arrows]
    out["curved_arrows"] = [_curve_entry(doc, c) for c in doc.curves]
    out["scheme"] = _r(_scheme(doc, out["arrows"]))
    out["brackets"] = _brackets(doc)
    unclaimed = [s for s in doc.shapes if s.kind == "poly" and not s.owner]
    out["unclaimed_polys"] = [{"bbox": tuple(round(v, 2) for v in s.bbox), "n": len(s.pts),
                               "curved": s.curved} for s in unclaimed]
    return out


def _side(doc, fr, atom_id, exclude):
    """Signed perpendicular sides of an atom's other neighbours in frame fr."""
    return [fr.local(doc.atoms[n].p)[1] for n in neighbors(doc, atom_id) if n != exclude]


def _acyclic_entry(doc, b, fr, lines, spacing):
    A, B = doc.atoms[b.a], doc.atoms[b.b]
    sides_a, sides_b = _side(doc, fr, b.a, b.b), _side(doc, fr, b.b, b.a)
    v = [ln.v for ln in lines]
    if abs(v[0] + v[1]) < 0.15 * spacing:
        placement = "centered"
    else:
        off = max(lines, key=lambda ln: abs(ln.v)).v
        subs = sides_a + sides_b
        toward = sum(1 for s in subs if s * off > 0)
        away = sum(1 for s in subs if s * off < 0)
        placement = f"offset (second line toward {toward} / away from {away} substituents)"
    return {
        "bond": b.id,
        "a": f"{A.label or 'C'}(deg {len(neighbors(doc, b.a))})",
        "b": f"{B.label or 'C'}(deg {len(neighbors(doc, b.b))})",
        "v": [round(x, 3) for x in v], "spacing": round(spacing, 3), "placement": placement,
        "subst_sides_a": [round(s, 2) for s in sides_a], "subst_sides_b": [round(s, 2) for s in sides_b],
        "ends": [(round(ln.u0, 3), round(fr.L - ln.u1, 3)) for ln in lines], "L": round(fr.L, 3),
    }


def _stereo_entry(doc, b, r):
    fr = r["frame"]
    shapes = r["shapes"]
    entry = {"bond": b.id, "display": b.display, "L": round(fr.L, 3),
             "narrow_atom": doc.atoms[fr_id(doc, fr)].label or "C"}
    end_atom = doc.atoms[b.b if b.display.endswith("Begin") else b.a]
    entry["wide_atom"] = end_atom.label or "C"
    if "Hash" in b.display:
        hashes = []
        for sh, loc in shapes:
            us, vs = [p[0] for p in loc], [p[1] for p in loc]
            hashes.append({"u": round(statistics.fmean(us), 3),
                           "thickness_u": round(max(us) - min(us), 3),
                           "len_v": round(max(vs) - min(vs), 3)})
        hashes.sort(key=lambda h: h["u"])
        gaps = [h2["u"] - h1["u"] for h1, h2 in zip(hashes, hashes[1:])]
        # each hash line is a parallelogram: its true thickness is measured
        # perpendicular to its long edge, not along u
        thick = []
        for sh, _ in shapes:
            e = max(((sh.pts[i], sh.pts[(i + 1) % len(sh.pts)]) for i in range(len(sh.pts))),
                    key=lambda ab: length(sub(ab[1], ab[0])))
            d = unit(sub(e[1], e[0]))
            nn = (-d[1], d[0])
            proj = [dot(p, nn) for p in sh.pts]
            thick.append(max(proj) - min(proj))
        # extrapolate hash length to u = L (linear in u)
        if len(hashes) >= 2:
            k = (hashes[-1]["len_v"] - hashes[0]["len_v"]) / (hashes[-1]["u"] - hashes[0]["u"])
            at_L = hashes[-1]["len_v"] + k * (fr.L - hashes[-1]["u"])
            at_0 = hashes[0]["len_v"] - k * hashes[0]["u"]
        else:
            at_L = at_0 = None
        entry.update({"count": len(hashes), "hashes": hashes, "pitch": stats(gaps),
                      "line_thickness": stats(thick),
                      "len_extrapolated_at_narrow_atom": round(at_0, 3) if at_0 is not None else None,
                      "len_extrapolated_at_wide_atom": round(at_L, 3) if at_L is not None else None})
    else:
        pts = [p for sh, loc in shapes for p in loc]
        umax, umin = max(p[0] for p in pts), min(p[0] for p in pts)
        wide = [p[1] for p in pts if p[0] > umax - 0.3]
        narrow = [p[1] for p in pts if p[0] < umin + 0.3]
        w_end = max(wide) - min(wide)
        entry.update({"u_range": (round(umin, 3), round(umax, 3)),
                      "wide_end_width": round(w_end, 3),
                      "narrow_end_width": round(max(narrow) - min(narrow), 3),
                      "width_extrapolated_at_wide_atom": round(w_end * fr.L / umax, 3) if umax else None,
                      "n_vertices": len(pts)})
    return entry


def fr_id(doc, fr):
    for aid, a in doc.atoms.items():
        if length(sub(a.p, fr.o)) < 1e-6:
            return aid
    return None


def _label_entry(doc, atom, asg):
    elem, whole = label_glyphs(doc, atom)
    if not elem:
        return {"atom": atom.id, "label": atom.label, "error": "element glyph not found"}
    face = load_face("Regular")
    ebox = union(g.box for g in elem)
    eadv_c = elem[0].x + sum(g.adv for g in elem) / 2
    entry = {
        "atom": atom.id, "label": atom.label, "element": atom.element, "nickname": atom.nickname,
        "glyph_run": "".join(g.ch for g in whole),
        "anchor_ink_center_minus_atom": tuple(round(v, 3) for v in sub(_box_center(ebox), atom.p)),
        "anchor_advance_center_x_minus_atom": round(eadv_c - atom.p[0], 3),
        "baseline_minus_atom_y": round(elem[0].base - atom.p[1], 3),
        "cap_height": round(face.cap_height(elem[0].size), 3),
        "size": elem[0].size,
    }
    # H, subscripts, charges relative to the element
    others = []
    for g in whole:
        if g in elem:
            continue
        others.append({"ch": g.ch, "size": round(g.size, 3),
                       "dx_origin_from_elem_origin": round(g.x - elem[0].x, 3),
                       "dbaseline": round(g.base - elem[0].base, 3)})
    entry["other_glyphs"] = others
    # bonds into this label
    bonds = []
    outlines = [c for g in elem for c in g.outline()]
    hull = convex_hull([p for c in outlines for p in c])
    margin = float(doc.settings["MarginWidth"])
    for b in doc.bonds:
        if atom.id not in (b.a, b.b) or b.id not in asg:
            continue
        r = asg[b.id]
        other = b.b if b.a == atom.id else b.a
        fr = frame(doc, other, atom.id)  # u runs toward the label
        lines = r.get("lines") or []
        if "stereo" in r:
            continue
        for ln in lines:
            ln2 = _line_from_shape(ln.shape, fr)
            end = fr.world(ln2.u1, ln2.v)
            c1, c2 = fr.world(ln2.u1, ln2.v - ln2.width / 2), fr.world(ln2.u1, ln2.v + ln2.width / 2)
            gap_box = ray_box_gap(end, fr.u, ebox)
            gap_outline = min(seg_seg_dist(c1, c2, p, q) for c in outlines
                              for p, q in zip(c, c[1:] + c[:1]))
            gap_hull = ray_polygon(fr.world(ln2.u1, 0.0), fr.u, hull)
            bonds.append({"bond": b.id, "order": b.order, "line_v": round(ln2.v, 3),
                          "dir_deg": round(math.degrees(math.atan2(-fr.u[1], fr.u[0])), 1),
                          "end_dist_from_atom": round(fr.L - ln2.u1, 3),
                          "gap_to_elem_inkbox_along_bond": round(gap_box, 3) if gap_box is not None else None,
                          "gap_to_elem_outline": round(gap_outline, 3),
                          "axial_gap_to_elem_hull": round(gap_hull, 3) if gap_hull is not None else None,
                          "axial_gap_to_hull_over_margin": round(gap_hull / margin, 3) if gap_hull is not None else None,
                          "other_end_labeled": bool(doc.atoms[other].label)})
    entry["bonds"] = bonds
    return entry


def _small_bars(doc):
    bars = []
    for s in doc.shapes:
        if s.kind != "poly" or s.owner or s.curved or len(s.pts) != 4 or s.group.endswith("-bonds"):
            continue  # ours: a crossing's oblique cut end (draw_mol.cut_caps) is a short polygon among the bonds
        x0, y0, x1, y1 = s.bbox
        w, h = x1 - x0, y1 - y0
        if 1.0 < max(w, h) < 6.0 and min(w, h) < 1.0:
            bars.append((s, "h" if w > h else "v"))
    return bars


def _charges(doc):
    bars = _small_bars(doc)
    used, out = set(), []
    for i, (s, o) in enumerate(bars):
        if i in used:
            continue
        mate = None
        if o == "h":
            for j, (t, o2) in enumerate(bars):
                if o2 == "v" and j not in used and length(sub(s.centroid, t.centroid)) < 0.3:
                    mate = j
        if o == "v":
            continue
        used.add(i)
        sign = "+" if mate is not None else "-"
        s.owner = "charge"
        if mate is not None:
            used.add(mate)
            bars[mate][0].owner = "charge"
        c = s.centroid
        x0, y0, x1, y1 = s.bbox
        atom = min(doc.atoms.values(), key=lambda a: length(sub(a.p, c)))
        rec = {"sign": sign, "center": (round(c[0], 3), round(c[1], 3)),
               "bar_length": round(x1 - x0, 3), "bar_thickness": round(y1 - y0, 3),
               "atom_label": atom.label or "C", "center_minus_atom": tuple(round(v, 3) for v in sub(c, atom.p))}
        if atom.label:
            elem, whole = label_glyphs(doc, atom)
            if elem:
                eb = union(g.box for g in elem)
                lb = union(g.box for g in whole)
                rec["center_minus_elem_inkbox_topright"] = (round(c[0] - eb[2], 3), round(c[1] - eb[1], 3))
                rec["center_minus_elem_inkbox_center"] = tuple(round(v, 3) for v in sub(c, _box_center(eb)))
                rec["label_run"] = "".join(g.ch for g in whole)
                rec["label_inkbox"] = tuple(round(v, 2) for v in lb)
        else:
            nb = neighbors(doc, atom.id)
            vecs = [unit(sub(doc.atoms[n].p, atom.p)) for n in nb]
            rec["bond_dirs_deg"] = [round(math.degrees(math.atan2(-v[1], v[0])), 1) for v in vecs]
            d = sub(c, atom.p)
            rec["charge_dir_deg"] = round(math.degrees(math.atan2(-d[1], d[0])), 1)
            rec["charge_dist"] = round(length(d), 3)
        out.append(rec)
    return out


def _dots(doc):
    circles = [s for s in doc.shapes if s.kind == "circle"]
    if not circles:
        return {}
    pairs, singles, used = [], [], set()
    for i, a in enumerate(circles):
        if i in used:
            continue
        mate = None
        for j, b in enumerate(circles):
            if j != i and j not in used and length(sub(a.center, b.center)) < 4.5:
                mate = j
                break
        atom = min(doc.atoms.values(), key=lambda at: length(sub(at.p, a.center)))
        if mate is None:
            used.add(i)
            singles.append({"center_minus_atom": tuple(round(v, 3) for v in sub(a.center, atom.p)),
                            "atom": atom.label or "C", "dist": round(length(sub(a.center, atom.p)), 3)})
            continue
        used |= {i, mate}
        b = circles[mate]
        mid = mul(add(a.center, b.center), 0.5)
        atom = min(doc.atoms.values(), key=lambda at: length(sub(at.p, mid)))
        rec = {"atom": atom.label or "C", "dot_spacing": round(length(sub(a.center, b.center)), 3),
               "mid_minus_atom": tuple(round(v, 3) for v in sub(mid, atom.p)),
               "mid_dist_from_atom": round(length(sub(mid, atom.p)), 3)}
        elem, _ = label_glyphs(doc, atom) if atom.label else ([], [])
        if elem:
            eb = union(g.box for g in elem)
            # gap from the nearer dot edge to the element ink box, measured on the side the pair sits
            d = sub(mid, _box_center(eb))
            if abs(d[0]) * (eb[3] - eb[1]) > abs(d[1]) * (eb[2] - eb[0]):  # left/right side
                gap = (min(a.center[0], b.center[0]) - eb[2]) if d[0] > 0 else (eb[0] - max(a.center[0], b.center[0]))
                side = "right" if d[0] > 0 else "left"
            else:
                gap = (min(a.center[1], b.center[1]) - eb[3]) if d[1] > 0 else (eb[1] - max(a.center[1], b.center[1]))
                side = "below" if d[1] > 0 else "above"
            rec.update({"side": side, "gap_dot_edge_to_elem_inkbox": round(gap - a.r, 3),
                        "mid_minus_elem_inkbox_center": tuple(round(v, 3) for v in d)})
        pairs.append(rec)
    return {"radius": stats([c.r for c in circles]), "pairs": pairs, "singles": singles}


def _head_geometry(sh, tip, axis):
    """Arrowhead polygon geometry relative to its tip and axis direction (unit,
    pointing toward the tip)."""
    n = (-axis[1], axis[0])
    loc = [(dot(sub(p, tip), axis), dot(sub(p, tip), n)) for p in sh.pts]
    back = min(u for u, _ in loc)
    vs = [v for _, v in loc]
    # inset point: control/vertex on the axis (|v| small) farthest from the tip, excluding the tip
    on_axis = [u for u, v in loc if abs(v) < 0.05 and u < -0.5]
    return {"length": round(-back, 3), "half_width_left": round(max(vs), 3),
            "half_width_right": round(-min(vs), 3),
            "inset_from_tip": round(-min(on_axis), 3) if on_axis else None}


def _arrow_entry(doc, arr):
    if arr.get("equilibrium"):
        return _equilibrium_entry(doc, arr)
    head, tail = arr["head"], arr["tail"]
    axis = unit(sub(head, tail))
    n = (-axis[1], axis[0])
    rec = {"tail": tuple(round(v, 2) for v in tail), "head": tuple(round(v, 2) for v in head),
           "length": round(length(sub(head, tail)), 3), "cdxml": arr["attrs"]}
    heads = []
    for sh in doc.shapes:
        if sh.kind == "poly" and sh.curved and not sh.owner:
            for tip, ax, nm in ((head, axis, "head"), (tail, mul(axis, -1), "tail")):
                if any(length(sub(p, tip)) < 0.05 for p in sh.pts):
                    sh.owner = f"arrow-{nm}"
                    heads.append((nm, _head_geometry(sh, tip, ax)))
    rec["heads"] = dict(heads)
    outlines = [s.width for s in doc.shapes if s.kind == "outline"]
    rec["head_outline_stroke"] = stats(outlines)
    # shaft: an unclaimed straight polygon along the axis
    for sh in doc.shapes:
        if sh.kind != "poly" or sh.curved or sh.owner:
            continue
        loc = [(dot(sub(p, tail), axis), dot(sub(p, tail), n)) for p in sh.pts]
        us, vs = [p[0] for p in loc], [p[1] for p in loc]
        if max(abs(v) for v in vs) < 0.6 and min(us) > -1 and max(us) < rec["length"] + 1 and max(us) - min(us) > 3:
            sh.owner = "arrow-shaft"
            rec["shaft"] = {"width": round(max(vs) - min(vs), 3), "u_start": round(min(us), 3),
                            "u_end_from_tip": round(rec["length"] - max(us), 3)}
    return rec


def _half_head_geometry(lp, tip, d, vc, o):
    """Half arrowhead (one barb) in its shaft's frame: `lp` the polygon's (u, v) in the arrow frame,
    `tip` its (u, v), d = +1 / -1 the shaft's direction along u, vc the shaft center's v and
    o = +1 / -1 the side away from the other shaft. back = distance behind the tip, side =
    offset from the shaft's center line, positive away from the other shaft (where the barb is)."""
    loc = [(d * (tip[0] - u), o * (v - vc)) for u, v in lp]
    tip_side = o * (tip[1] - vc)
    flat = [b for b, s in loc if abs(s - tip_side) < 0.05]
    barb = max(s for _, s in loc)
    inset = max(flat)
    ctl = [s for b, s in loc if abs(b - inset) < 0.05 and tip_side + 0.05 < s < barb - 0.05]
    return {"length": round(max(b for b, _ in loc), 3), "barb_from_center": round(barb, 3),
            "tip_from_center": round(tip_side, 3), "inset_from_tip": round(inset, 3),
            "back_control_from_center": round(ctl[0], 3) if ctl else None,
            "back_control_per_barb": round(ctl[0] / barb, 4) if ctl else None}


def _equilibrium_entry(doc, arr):
    """Equilibrium arrow (acs-style.md §3.6): two straight shafts parallel to Tail3D -> Head3D,
    each with a filled half head (curved polygon) at one end. Per shaft: center offset v from the
    CDXML axis, width, tip and butt along the axis (u from Tail3D), direction (forward: tip
    toward Head3D) and its half head (_half_head_geometry). Per arrow: shaft spacing (center to
    center), balanced (equal shafts: the stagger of each butt behind the other shaft's tip) or
    unbalanced (the short shaft's center against the long one's, and long / short butt-to-notch
    ratio), and the ink extent along the axis."""
    head, tail = arr["head"], arr["tail"]
    axis = unit(sub(head, tail))
    n = (-axis[1], axis[0])
    L = length(sub(head, tail))
    shafts, heads = [], []
    for sh in doc.shapes:
        if sh.kind != "poly" or sh.owner:
            continue
        lp = [(dot(sub(p, tail), axis), dot(sub(p, tail), n)) for p in sh.pts]
        us, vs = [q[0] for q in lp], [q[1] for q in lp]
        if max(abs(v) for v in vs) > 4 or min(us) < -3 or max(us) > L + 3:
            continue
        if sh.curved:
            heads.append((sh, lp))
        elif max(vs) - min(vs) < 1.0 and max(us) - min(us) > 3:
            shafts.append((sh, lp))
    recs = []
    for sh, lp in shafts:
        us, vs = [q[0] for q in lp], [q[1] for q in lp]
        vc, u0, u1 = (max(vs) + min(vs)) / 2, min(us), max(us)
        o = 1 if vc > 0 else -1
        for hs, hp in heads:
            hu = [q[0] for q in hp]
            if hs.owner or o * statistics.fmean(q[1] for q in hp) < 0:
                continue
            if min(hu) < u1 < max(hu) + 0.5:
                d, tip = 1, max(hp, key=lambda q: q[0])
            elif min(hu) - 0.5 < u0 < max(hu):
                d, tip = -1, min(hp, key=lambda q: q[0])
            else:
                continue
            sh.owner, hs.owner = "arrow-shaft", "arrow-head"
            butt = u0 if d == 1 else u1
            hg = _half_head_geometry(hp, tip, d, vc, o)
            recs.append({"v": round(vc, 3), "width": round(max(vs) - min(vs), 3),
                         "direction": "forward" if d == 1 else "back",
                         "tip_u": round(tip[0], 3), "butt_u": round(butt, 3),
                         "length": round(abs(tip[0] - butt), 3),
                         "u_end_from_tip": round(d * (tip[0] - (u1 if d == 1 else u0)), 3),
                         "head": hg, "barb": "outer" if hg["barb_from_center"] > 0 else "inner",
                         "y": add(tail, mul(n, vc))[1]})
            break
    recs.sort(key=lambda r: r["v"])
    rec = {"tail": tuple(round(v, 2) for v in tail), "head": tuple(round(v, 2) for v in head),
           "length": round(L, 3), "cdxml": arr["attrs"], "equilibrium": {"shafts": recs}, "heads": {}}
    if len(recs) != 2:
        rec["equilibrium"]["error"] = f"{len(recs)} shafts with heads"
        return rec
    e = rec["equilibrium"]
    a, b = recs
    lo = min(min(r["tip_u"], r["butt_u"]) for r in recs)
    hi = max(max(r["tip_u"], r["butt_u"]) for r in recs)
    rec["ends"] = (add(tail, mul(axis, lo)), add(tail, mul(axis, hi)))
    rec["outer_y"] = tuple(sorted(r["y"] for r in recs))
    e.update({"spacing": round(b["v"] - a["v"], 3), "center_v": round((a["v"] + b["v"]) / 2, 3),
              "extent_u": (round(lo, 3), round(hi, 3))})
    fwd = next((r for r in recs if r["direction"] == "forward"), None)
    back = next((r for r in recs if r["direction"] == "back"), None)
    if abs(a["length"] - b["length"]) < 1.0:
        e["kind"] = "balanced"
        if fwd and back:
            e["stagger_tail"] = round(fwd["butt_u"] - back["tip_u"], 3)  # forward butt behind the back tip
            e["stagger_head"] = round(fwd["tip_u"] - back["butt_u"], 3)
    else:
        long, short = (a, b) if a["length"] > b["length"] else (b, a)
        e["kind"] = "unbalanced"
        e["long"] = long["direction"]
        e["short_length"], e["long_length"] = short["length"], long["length"]
        e["short_center_minus_long_center"] = round((short["tip_u"] + short["butt_u"]) / 2
                                                    - (long["tip_u"] + long["butt_u"]) / 2, 3)
        e["long_over_short_butt_to_notch"] = round((long["length"] - long["head"]["inset_from_tip"])
                                                   / (short["length"] - short["head"]["inset_from_tip"]), 4)
    return rec


def _curve_entry(doc, c):
    pts = c["pts"]
    stroke = None
    end = pts[-2] if len(pts) % 3 == 0 else pts[-1]
    # ChemDraw CurvePoints: (in-handle, anchor, out-handle) per anchor
    anchors = pts[1::3]
    best = None
    for sh in doc.shapes:
        if sh.kind == "stroke" and sh.curved:
            d = min(length(sub(sh.pts[0], anchors[0])), length(sub(sh.pts[-1], anchors[0])))
            if best is None or d < best[0]:
                best = (d, sh)
    rec = {"head_type": c["head"], "curve_type": c["type"],
           "anchors": [tuple(round(v, 2) for v in a) for a in anchors],
           "handles": [tuple(round(v, 2) for v in p) for i, p in enumerate(pts) if i % 3 != 1]}
    if best:
        stroke = best[1]
        rec["stroke_width"] = round(stroke.width, 3)
        rec["svg_path_pts"] = [tuple(round(v, 2) for v in p) for p in stroke.pts]
        p_end = stroke.pts[3] if len(stroke.pts) >= 4 else stroke.pts[-1]
        # end tangent from the last control point that isn't the end (ChemDraw's C–Br arrow has c2 = p3)
        c_last = next(p for p in reversed(stroke.pts[:3]) if length(sub(p_end, p)) > 1e-6)
        tangent = unit(sub(p_end, c_last))
        # the curve's last anchor is the stroke end (the head's notch), not the tip (§13.2)
        rec["last_anchor_minus_stroke_end"] = tuple(round(v, 3) for v in sub(anchors[-1], p_end))
        # head polygon: unclaimed straight poly with a vertex near the stroke end
        heads = [s for s in doc.shapes if s.kind == "poly" and not s.owner and not s.curved
                 and min(length(sub(p, p_end)) for p in s.pts) < 3.0 and len(s.pts) in (3, 4)]
        if heads:
            h = min(heads, key=lambda s: length(sub(s.centroid, p_end)))
            h.owner = "curve-head"
            tip = max(h.pts, key=lambda p: dot(sub(p, p_end), tangent))
            ax = tangent
            n = (-ax[1], ax[0])
            loc = [(dot(sub(p, tip), ax), dot(sub(p, tip), n)) for p in h.pts]
            pe = (dot(sub(p_end, tip), ax), dot(sub(p_end, tip), n))
            rec["head"] = {"n_vertices": len(h.pts), "length": round(-min(u for u, _ in loc), 3),
                           # stroke end in the head frame: fishhooks end beside the head's axis, toward the barb
                           "stroke_end_local": (round(pe[0], 3), round(pe[1], 3)),
                           "extent_left": round(max(v for _, v in loc), 3),
                           "extent_right": round(-min(v for _, v in loc), 3),
                           "tip_minus_stroke_end": tuple(round(v, 3) for v in sub(tip, p_end)),
                           "pts_local": [(round(u, 3), round(v, 3)) for u, v in loc]}
        # chord / bow (approximate: hand-shaped)
        p0, p3 = stroke.pts[0], p_end
        ch = sub(p3, p0)
        L = length(ch)
        nn = (-ch[1] / L, ch[0] / L)
        bow = max(abs(dot(sub(p, p0), nn)) for p in stroke.pts[:4])
        rec["chord"] = round(L, 3)
        rec["control_offset_max_over_chord"] = round(bow / L, 3)
    return rec


MULTISTEP = re.compile(r"^\d+\.")  # "1. n-BuLi, ..." : a numbered stage, set left-justified


def _caption_glyphs(doc, cap):
    """Glyphs of one caption: given directly (our renders) or found inside the CDXML text box."""
    if "glyphs" in cap:
        return [g for g in cap["glyphs"] if g.box]
    bb = cap["bbox"]
    return [g for g in doc.glyphs if g.box and bb[0] - 1 <= _box_center(g.box)[0] <= bb[2] + 1
            and bb[1] - 1 <= _box_center(g.box)[1] <= bb[3] + 1]


def text_lines(glyphs) -> list[dict]:
    """Group a caption's glyphs into lines, top to bottom. Full-size glyphs set the baselines;
    smaller ones (sub- and superscripts) join the nearest line. Advances: x0 is the first
    glyph origin, x1 the last glyph's origin + advance."""
    main = max(g.size for g in glyphs)
    bases = []
    for b in sorted(g.base for g in glyphs if abs(g.size - main) < 0.01):
        if not bases or b - bases[-1] > 0.05:
            bases.append(b)
    lines = [{"base": b, "glyphs": []} for b in bases]
    for g in glyphs:
        min(lines, key=lambda ln: abs(ln["base"] - g.base))["glyphs"].append(g)
    out = []
    for ln in lines:
        gl = sorted(ln["glyphs"], key=lambda g: g.x)
        base = ln["base"]
        small = [g for g in gl if g.size < main - 0.01]
        subs = [g for g in small if g.base > base + 0.01]
        sups = [g for g in small if g.base < base - 0.01]
        rec = {"text": "".join(g.ch for g in gl), "base": base, "size": main,
               "x0": gl[0].x, "x1": gl[-1].x + gl[-1].adv, "ink": union(g.box for g in gl),
               "has_sub": bool(subs), "has_sup": bool(sups), "bold": any(g.bold for g in gl),
               "italic": "".join(g.ch for g in gl if g.italic)}
        if subs:
            rec["sub"] = {"size_ratio": subs[0].size / main, "drop_over_size": (subs[0].base - base) / main}
        if sups:
            rec["sup"] = {"size_ratio": sups[0].size / main, "raise_over_size": (base - sups[0].base) / main}
        out.append(rec)
    return out


def _rings(doc, atom_ids) -> list[list[str]]:
    """Smallest cycle through each ring bond, as atom-id lists (no duplicates)."""
    adj = {}
    for b in doc.bonds:
        if b.a in atom_ids:
            adj.setdefault(b.a, set()).add(b.b)
            adj.setdefault(b.b, set()).add(b.a)
    rings = {}
    for b in doc.bonds:
        if not b.ring or b.a not in atom_ids:
            continue
        prev, queue = {b.a: None}, [b.a]
        while queue and b.b not in prev:
            u = queue.pop(0)
            for w in adj[u]:
                if w not in prev and {u, w} != {b.a, b.b}:
                    prev[w] = u
                    queue.append(w)
        if b.b in prev:
            path, u = [], b.b
            while u is not None:
                path.append(u)
                u = prev[u]
            rings.setdefault(frozenset(path), path)
    return list(rings.values())


def _structures(doc) -> dict:
    """fragment id -> {ink, rings} for fragments drawn with atoms (formula text excluded)."""
    out = {}
    for fid in dict.fromkeys(a.frag for a in doc.atoms.values()):
        if fid in doc.text_frags:
            continue
        ids = {a.id for a in doc.atoms.values() if a.frag == fid}
        bond_ids = {b.id for b in doc.bonds if b.frag == fid}
        boxes = [s.bbox for s in doc.shapes if s.owner in bond_ids]
        for aid in ids:
            if doc.atoms[aid].label:
                boxes += [g.box for g in label_glyphs(doc, doc.atoms[aid])[1]]
        if not boxes:
            continue
        centers = [(statistics.fmean(doc.atoms[i].p[0] for i in r), statistics.fmean(doc.atoms[i].p[1] for i in r))
                   for r in _rings(doc, ids)]
        out[fid] = {"ink": union(boxes), "ring_centers": centers}
    return out


def _rows(structs: dict) -> list[set]:
    """Structures grouped into rows by overlapping vertical ink extents."""
    rows: list[set] = []
    for fid, s in sorted(structs.items(), key=lambda kv: kv[1]["ink"][1]):
        y0, y1 = s["ink"][1], s["ink"][3]
        hit = [r for r in rows if any(structs[o]["ink"][1] < y1 and y0 < structs[o]["ink"][3] for o in r)]
        merged = set().union(*hit) | {fid} if hit else {fid}
        rows = [r for r in rows if r not in hit] + [merged]
    return rows


def _scheme(doc, arrows: list[dict]):
    """Scheme layout (acs-style.md §3): each caption is assigned to an arrow (horizontal overlap,
    then nearest in y), a plus sign or a compound number; structures are grouped in rows.
    `arrows` are the _arrow_entry records (they carry the measured head inset)."""
    structs = _structures(doc)
    rows = _rows(structs)

    def row_of(y):
        return next((r for r in rows if any(structs[f]["ink"][1] <= y <= structs[f]["ink"][3] for f in r)), set())

    arrow_recs = []
    for a in arrows:
        tail, head = a["tail"], a["head"]
        inset = (a["heads"].get("head") or {}).get("inset_from_tip") or 0.0
        if "ends" in a:  # equilibrium: its ink extent, left to right, centered (no notch)
            tail, head = sorted(a["ends"])
        sgn = 1 if head[0] >= tail[0] else -1
        row = row_of(tail[1])
        left = [structs[f]["ink"] for f in row if structs[f]["ink"][2] <= min(tail[0], head[0])]
        right = [structs[f]["ink"] for f in row if structs[f]["ink"][0] >= max(tail[0], head[0])]
        arrow_recs.append({
            "tail": tail, "head": head, "axis_y": tail[1], "length": a["length"], "notch": head[0] - sgn * inset,
            "mid": (tail[0] + head[0]) / 2 if "ends" in a else None,
            # conditions are measured from the outer shafts of an equilibrium arrow
            "axis_above": a["outer_y"][0] if "outer_y" in a else tail[1],
            "axis_below": a["outer_y"][1] if "outer_y" in a else tail[1],
            "row": row, "above": [], "below": [],
            "gap_left_struct_to_tail": min(tail[0], head[0]) - max(b[2] for b in left) if left else None,
            "gap_head_to_right_struct": min(b[0] for b in right) - max(tail[0], head[0]) if right else None,
            "structures": [{"ink_center_minus_axis": (structs[f]["ink"][1] + structs[f]["ink"][3]) / 2 - tail[1],
                            "ring_centers_minus_axis": sorted((c[1] - tail[1] for c in structs[f]["ring_centers"]),
                                                              key=abs)}
                           for f in sorted(row, key=lambda f: structs[f]["ink"][0])],
        })
    pluses, numbers, other = [], [], []
    for cap in doc.captions:
        glyphs = _caption_glyphs(doc, cap)
        if not glyphs:
            continue
        lines = text_lines(glyphs)
        ink = union(ln["ink"] for ln in lines)
        cx, cy = _box_center(ink)
        if cap["text"].strip() == "+":
            row = row_of(cy)
            left = [structs[f]["ink"] for f in row if structs[f]["ink"][2] < ink[0]]
            right = [structs[f]["ink"] for f in row if structs[f]["ink"][0] > ink[2]]
            axis = next((r["axis_y"] for r in arrow_recs if r["row"] == row), None)
            pluses.append({"gap_left_struct": ink[0] - max(b[2] for b in left) if left else None,
                           "gap_right_struct": min(b[0] for b in right) - ink[2] if right else None,
                           "ink_center_minus_axis": cy - axis if axis is not None else None})
        elif all(ln["bold"] for ln in lines):
            above = [f for f in structs if structs[f]["ink"][3] < ink[1] and structs[f]["ink"][0] < cx < structs[f]["ink"][2]]
            if not above:
                continue
            fid = max(above, key=lambda f: structs[f]["ink"][3])
            sb = structs[fid]["ink"]
            row = next(r for r in rows if fid in r)
            numbers.append({"text": cap["text"], "row": row, "base": lines[0]["base"], "ink": ink,
                            "advance": lines[0]["x1"] - lines[0]["x0"],
                            "gap_own_struct_to_ink_top": ink[1] - sb[3],
                            "gap_row_bottom_to_ink_top": ink[1] - max(structs[f]["ink"][3] for f in row),
                            "center_minus_struct_center_x": cx - (sb[0] + sb[2]) / 2})
        else:
            cands = [r for r in arrow_recs if ink[0] < max(r["tail"][0], r["head"][0])
                     and ink[2] > min(r["tail"][0], r["head"][0])]
            if not cands:
                other.append(cap["text"])
                continue
            r = min(cands, key=lambda r: abs(cy - r["axis_y"]))
            for ln in lines:
                ln["caption"] = cap["text"]
                ln["formula"] = cap.get("formula", False)
                ln["left_justified"] = bool(MULTISTEP.match(lines[0]["text"]))
            (r["above"] if cy < r["axis_y"] else r["below"]).append(lines)
    for r in arrow_recs:
        r["above"] = [ln for lines in sorted(r["above"], key=lambda ls: ls[0]["base"]) for ln in lines]
        r["below"] = [ln for lines in sorted(r["below"], key=lambda ls: ls[0]["base"]) for ln in lines]
        tail_x, notch = r["tail"][0], r["notch"]
        mid = r.pop("mid") or (tail_x + notch) / 2
        for side in ("above", "below"):
            for ln in r[side]:
                ln["center_minus_tail_notch_mid"] = (ln["x0"] + ln["x1"]) / 2 - mid
                ln["ink_center_minus_tail_notch_mid"] = (ln["ink"][0] + ln["ink"][2]) / 2 - mid
                ln["x0_minus_tail"] = ln["x0"] - tail_x
                ln["ink_left_minus_tail"] = ln["ink"][0] - tail_x
                ln["notch_minus_x1"] = notch - ln["x1"]
        if r["above"]:
            last = r["above"][-1]
            r["above_axis_minus_last_base"] = r["axis_above"] - last["base"]
            r["above_axis_minus_ink_bottom"] = r["axis_above"] - max(ln["ink"][3] for ln in r["above"])
        if r["below"]:
            first = r["below"][0]
            r["below_first_base_minus_axis"] = first["base"] - r["axis_below"]
            r["below_ink_top_minus_axis"] = min(ln["ink"][1] for ln in r["below"]) - r["axis_below"]
        lines = r["above"] + r["below"]
        if lines:
            w = max(lines, key=lambda ln: ln["x1"] - ln["x0"])
            r["widest"] = {"text": w["text"], "advance": w["x1"] - w["x0"], "ink_width": w["ink"][2] - w["ink"][0],
                           "pad_tail": w["x0"] - tail_x, "pad_notch": notch - w["x1"],
                           "left_justified": w["left_justified"]}
        r["pitches"] = [{"pitch": b["base"] - a["base"], "size": a["size"], "sub_above": a["has_sub"],
                         "sub_below": b["has_sub"], "texts": (a["text"], b["text"])}
                        for side in ("above", "below") for a, b in zip(r[side], r[side][1:])
                        if a["caption"] == b["caption"]]
    by_row = {}
    for n in numbers:
        by_row.setdefault(frozenset(n["row"]), []).append(n)
    number_rows = [{"texts": [n["text"] for n in ns], "bases": [n["base"] for n in ns],
                    "gap_row_bottom_to_ink_top": [n["gap_row_bottom_to_ink_top"] for n in ns],
                    "center_minus_struct_center_x": [n["center_minus_struct_center_x"] for n in ns],
                    "numbers": ns}
                   for ns in by_row.values()]
    for r in arrow_recs:
        r["row"] = sorted(r["row"])
    for n in numbers:
        n["row"] = sorted(n["row"])
    return {"arrows": arrow_recs, "pluses": pluses, "number_rows": number_rows, "unassigned": other}


def _r(v, nd=3):
    """Round floats (recursively) for JSON output."""
    if isinstance(v, float):
        return round(v, nd)
    if isinstance(v, (list, tuple)):
        return type(v)(_r(x, nd) for x in v)
    if isinstance(v, dict):
        return {k: _r(x, nd) for k, x in v.items()}
    return v


def _brackets(doc):
    """Square brackets: ChemDraw draws each as three stroked segments."""
    segs = [s for s in doc.shapes if s.kind == "stroke" and len(s.pts) == 2]
    out = []
    for v in segs:
        (x0, y0), (x1, y1) = v.pts
        if abs(x1 - x0) > 0.01 or abs(y1 - y0) < 20:
            continue
        lips = [h for h in segs if h is not v and abs(h.pts[0][1] - h.pts[1][1]) < 0.01
                and min(length(sub(p, q)) for p in h.pts for q in v.pts) < 0.01]
        for sh in [v, *lips]:
            sh.owner = "bracket"
        lip_len = [abs(h.pts[0][0] - h.pts[1][0]) for h in lips]
        lip_dir = "right" if lips and max(p[0] for h in lips for p in h.pts) > x0 + 0.01 else "left"
        out.append({"x": round(x0, 3), "y_top": round(min(y0, y1), 3), "y_bottom": round(max(y0, y1), 3),
                    "height": round(abs(y1 - y0), 3), "stroke_width": round(v.width, 3),
                    "lip_length": stats(lip_len), "lips_point": lip_dir})
    if out:
        content = [s.bbox for s in doc.shapes if s.owner and s.owner != "bracket"]
        content += [g.box for g in doc.glyphs if g.box]
        cb = union(content)
        for b in out:
            b["content_ink"] = tuple(round(v, 2) for v in cb)
            b["extends_above_content"] = round(cb[1] - b["y_top"], 3)
            b["extends_below_content"] = round(b["y_bottom"] - cb[3], 3)
            b["gap_to_content_x"] = round(cb[0] - b["x"], 3) if b["lips_point"] == "right" else round(b["x"] - cb[2], 3)
    return out


# --------------------------------------------------------------------------- clip-model fit


def clip_rows(names=None) -> list[dict]:
    """One row per drawn bond end at a label: the measured distance from the label atom to
    where the bond stops (plain lines: line end; wedges and hashes: the wide end), with the
    direction from the label atom toward the bond's other atom (unit, y down) and the element
    glyph outline as ChemDraw placed it, relative to the atom (y down)."""
    rows = []
    for n in names or all_names():
        doc = load(n)
        asg = assign(doc)
        margin = float(doc.settings["MarginWidth"])
        for b in doc.bonds:
            r = asg.get(b.id)
            if not r:
                continue
            for atom_id in (b.a, b.b):
                atom = doc.atoms[atom_id]
                if not atom.label:
                    continue
                elem, _ = label_glyphs(doc, atom)
                if not elem:
                    continue
                other = b.b if atom_id == b.a else b.a
                fr = frame(doc, other, atom_id)  # u runs toward the label
                if "stereo" in r:
                    narrow = b.a if b.display.endswith("Begin") else b.b
                    if narrow == atom_id or not r["shapes"]:
                        continue
                    end = fr.L - max(u for _, loc in r["shapes"] for u, _ in loc)
                    kind = "hash" if "Hash" in b.display else "wedge"
                elif r.get("lines"):
                    end = statistics.median(fr.L - _line_from_shape(ln.shape, fr).u1 for ln in r["lines"])
                    kind = "line"
                else:
                    continue
                outline = [[(x - atom.p[0], y - atom.p[1]) for x, y in c] for g in elem for c in g.outline()]
                rows.append({"file": n, "atom": atom_id, "bond": b.id, "kind": kind, "outline": outline,
                             "element": "".join(g.ch for g in elem), "size": elem[0].size, "margin": margin,
                             "dir": (-fr.u[0], -fr.u[1]), "end": end,
                             "dir_deg": round(math.degrees(math.atan2(fr.u[1], -fr.u[0])), 1)})
    return rows


def fit_label_clip(rows=None) -> dict:
    """Fit the bond-to-label clip model used by draw_mol.clip_distance (acs-style.md §2.1).

    A bond stops c x MarginWidth beyond the furthest point, along the bond, of the element
    glyph outline (as ChemDraw placed it) within ±band of the bond axis, with
    band = beta x label size and the outline shifted by (dx, dy) pt. Grid search: fewest
    rows off by more than 0.1 pt, then the worst residual among the rest, then rms.
    """
    from acs_chem_svg import draw_mol

    rows = rows if rows is not None else clip_rows()
    best = None
    for beta in [x / 1000 for x in range(170, 211, 5)]:
        for dx in [x / 100 for x in range(-5, 6)]:
            for dy in [x / 100 for x in range(-5, 6)]:
                ext = [draw_mol.clip_distance_band(r["outline"], r["dir"], beta * r["size"], r["margin"], 0.0,
                                                   (dx, dy)) for r in rows]
                for c in [x / 10000 for x in range(9900, 10301, 25)]:
                    res = [e + c * r["margin"] - r["end"] for e, r in zip(ext, rows)]
                    inside = [abs(x) for x in res if abs(x) <= 0.2]  # worst of the rows the model can fit
                    key = (sum(abs(x) > 0.1 for x in res), round(max(inside), 3),
                           math.sqrt(statistics.fmean(x * x for x in res)))
                    if best is None or key < best[0]:
                        best = (key, beta, dx, dy, c, res)
    (n_bad, worst, rms), beta, dx, dy, c, res = best
    for r, x in zip(rows, res):
        r["residual"] = round(x, 3)
    return {"band_per_size": beta, "offset": (dx, dy), "per_margin": c, "n": len(rows),
            "n_over_0.1": n_bad, "worst": worst, "rms": round(rms, 3),
            "rows": [{k: v for k, v in r.items() if k not in ("outline", "dir")} for r in rows]}


# --------------------------------------------------------------------------- summary


def summarize(m: dict) -> dict:
    """Flat key-parameter table from one measure() result (medians, pt)."""
    def med(vals):
        vals = [v for v in vals if v is not None]
        return round(statistics.median(vals), 3) if vals else None

    labels = [e for e in m["labels"] if "bonds" in e]
    lb = [b for e in labels for b in e["bonds"]]
    out = {
        "bond_length": (m["bond_length"] or {}).get("median"),
        "line_width": (m["line_width"] or {}).get("median"),
        "bond_spacing_ring": (m["bond_spacing_ring"] or {}).get("median"),
        "ring_inner_inset": (m["ring_inner_line"]["inset_unlabeled_ends"] or {}).get("median"),
        "acyclic_double_spacing": med([e["spacing"] for e in m["acyclic_double"]]),
        "triple_spacing": med([(e["v"][2] - e["v"][0]) / 2 for e in m["triple"]]),
        "label_baseline_below_atom": med([e["baseline_minus_atom_y"] for e in labels]),
        "label_anchor_dx": med([abs(e["anchor_advance_center_x_minus_atom"]) for e in labels
                                if len(e["element"]) == 1]),
        "end_C=O_vertical_from_O": med([b["end_dist_from_atom"] for e in labels for b in e["bonds"]
                                        if e["label"] == "O" and b["order"] == 2 and abs(b["dir_deg"]) == 90]),
        "end_OH_30deg_from_O": med([b["end_dist_from_atom"] for e in labels for b in e["bonds"]
                                    if e["label"] == "OH" and e["glyph_run"] == "OH" and b["dir_deg"] == -30]),
        "end_horizontal_into_N_from_left": med([b["end_dist_from_atom"] for e in labels for b in e["bonds"]
                                                if e["element"] == "N" and b["dir_deg"] == 0]),
        "gap_to_glyph_outline": med([b["gap_to_elem_outline"] for b in lb]),
        "axial_gap_hull_over_margin": med([b["axial_gap_to_hull_over_margin"] for e in labels
                                           for b in e["bonds"] if e["label"] not in ("H", "I")]),
        "subscript_size_ratio": med([o["size"] / e["size"] for e in labels for o in e["other_glyphs"]
                                     if o["ch"].isdigit()]),
        "subscript_drop": med([o["dbaseline"] for e in labels for o in e["other_glyphs"] if o["ch"].isdigit()]),
    }
    for st in m["stereo"]:
        if "Hash" in st["display"]:
            out["hash_count"] = st["count"]
            out["hash_pitch"] = (st["pitch"] or {}).get("median")
            out["hash_line_thickness"] = (st["line_thickness"] or {}).get("median")
            out["hash_last_len"] = st["hashes"][-1]["len_v"] if st["hashes"] else None
        else:
            out["wedge_wide_end_width"] = st["wide_end_width"]
            out["wedge_end_u"] = st["u_range"][1]
    if m["charges"]:
        out["charge_bar_length"] = med([c["bar_length"] for c in m["charges"]])
        out["charge_bar_thickness"] = med([c["bar_thickness"] for c in m["charges"]])
    heads = [h for a in m["arrows"] for h in a["heads"].values()]
    if heads:
        out["arrowhead_length"] = med([h["length"] for h in heads])
        out["arrowhead_half_width"] = med([(h["half_width_left"] + h["half_width_right"]) / 2 for h in heads])
        out["arrowhead_inset"] = med([h["inset_from_tip"] for h in heads])
        out["arrow_shaft_width"] = med([a["shaft"]["width"] for a in m["arrows"] if "shaft" in a])
    return out



def pdf_char_origins(name: str, ref_dir: Path = REF_DIR) -> tuple[str, list[float]]:
    """(page text, x origin of each character) from the reference PDF. These are ChemDraw's own
    glyph positions, independent of our kerning (the SVG gives one origin per run only)."""
    import ctypes

    import pypdfium2.raw as raw

    tp = pypdfium2.PdfDocument(str(ref_dir / f"{name}.pdf"))[0].get_textpage()
    xs = []
    for i in range(tp.count_chars()):
        x, y = ctypes.c_double(), ctypes.c_double()
        raw.FPDFText_GetCharOrigin(tp.raw, i, x, y)
        xs.append(x.value)
    return tp.get_text_range(), xs

# --------------------------------------------------------------------------- scheme summary (§3)

SCHEME_REFS = {"1.0": ["scheme_basic_x1.0", "scheme_multi_x1.0"], "0.8": ["scheme_basic_x0.8", "scheme_multi_x0.8"],
               "0.6": ["scheme_basic_x0.6", "scheme_multi_x0.6"]}


def scheme_summary(results: dict) -> dict:
    """Scheme parameters per preset (acs-style.md §3) from measure() results of
    scheme_basic_x<k> + scheme_multi_x<k>. Medians unless noted; `n` is the sample count."""
    def med(vals):
        vals = [v for v in vals if v is not None]
        return {"median": round(statistics.median(vals), 3), "min": round(min(vals), 3),
                "max": round(max(vals), 3), "n": len(vals)} if vals else None

    out = {}
    for k, names in SCHEME_REFS.items():
        ms = [results[n]["scheme"] for n in names if n in results]
        if not ms:
            continue
        arrows = [a for m in ms for a in m["arrows"]]
        lines = [(a, ln) for a in arrows for side in ("above", "below") for ln in a[side]]
        centered = [ln for _, ln in lines if not ln["left_justified"]]
        stepped = [ln for _, ln in lines if ln["left_justified"]]
        widest = [a["widest"] for a in arrows if "widest" in a and not a["widest"]["left_justified"]]
        pitches = [p for a in arrows for p in a["pitches"]]
        size = max(ln["size"] for _, ln in lines)
        plain = [p["pitch"] - p["size"] for p in pitches if not p["sub_above"] and not p["sub_below"]]
        rows = [r for m in ms for r in m["number_rows"]]
        ring_rows = [a for a in arrows if a["structures"] and all(s["ring_centers_minus_axis"] for s in a["structures"])]
        ink_rows = [a for a in arrows if a not in ring_rows]
        out[k] = {
            "arrow_length": med([a["length"] for a in arrows]),
            "conditions_pad (widest line: tail -> origin, end of advance -> notch)":
                med([w["pad_tail"] for w in widest] + [w["pad_notch"] for w in widest]),
            "conditions_center_minus_tail_notch_mid (advance)": med([ln["center_minus_tail_notch_mid"] for ln in centered]),
            "conditions_center_minus_tail_notch_mid (ink)": med([ln["ink_center_minus_tail_notch_mid"] for ln in centered]),
            "multistep x0 - tail (glyph origin)": med([ln["x0_minus_tail"] for ln in stepped]),
            "multistep ink left - tail": med([ln["ink_left_minus_tail"] for ln in stepped]),
            "above: axis - last baseline": med([a.get("above_axis_minus_last_base") for a in arrows]),
            "above: axis - ink bottom": med([a.get("above_axis_minus_ink_bottom") for a in arrows]),
            "below: first baseline - axis": med([a.get("below_first_base_minus_axis") for a in arrows]),
            "below: ink top - axis": med([a.get("below_ink_top_minus_axis") for a in arrows]),
            "line pitch - size, no subscripts": med(plain),
            "line pitch - size - plain, sub on upper line only": med(
                [p["pitch"] - p["size"] - plain[0] for p in pitches if p["sub_above"] and not p["sub_below"]] if plain else []),
            "line pitch - size - plain, sub on both lines": med(
                [p["pitch"] - p["size"] - plain[0] for p in pitches if p["sub_above"] and p["sub_below"]] if plain else []),
            "subscript size / size": med([ln["sub"]["size_ratio"] for _, ln in lines if "sub" in ln]),
            "subscript drop / size": med([ln["sub"]["drop_over_size"] for _, ln in lines if "sub" in ln]),
            "superscript size / size": med([ln["sup"]["size_ratio"] for _, ln in lines if "sup" in ln]),
            "superscript raise / size": med([ln["sup"]["raise_over_size"] for _, ln in lines if "sup" in ln]),
            "structure ring center - axis (rows where every structure has a ring)":
                med([s["ring_centers_minus_axis"][0] for a in ring_rows for s in a["structures"]]),
            "structure ink center - axis (other rows)": med([s["ink_center_minus_axis"] for a in ink_rows for s in a["structures"]]),
            "arrow gap (structure ink -> tail, head -> structure ink)":
                med([a["gap_left_struct_to_tail"] for a in arrows] + [a["gap_head_to_right_struct"] for a in arrows]),
            "plus gap (structure ink <-> plus ink)": med([p[g] for m in ms for p in m["pluses"]
                                                          for g in ("gap_left_struct", "gap_right_struct")]),
            "plus ink center - axis": med([p["ink_center_minus_axis"] for m in ms for p in m["pluses"]]),
            "number: row's lowest structure ink -> number ink top": med([g for r in rows for g in r["gap_row_bottom_to_ink_top"]]),
            "number: baseline spread within a row (max - min)": med([max(r["bases"]) - min(r["bases"]) for r in rows]),
            "number: center - structure ink center x": med([abs(c) for r in rows for c in r["center_minus_struct_center_x"]]),
            "caption size": size,
        }
    return out


# --------------------------------------------------------------------------- mechanism summary (§4)

WRAP_REFS = {"1.0": "scheme_wrap_x1.0", "0.8": "scheme_wrap_x0.8", "0.6": "scheme_wrap_x0.6"}


def wrap_summary(name: str | Doc) -> dict:
    """Multi-row scheme parameters (acs-style.md §3.4) from one scheme_wrap_x<k> file (or a
    loaded Doc: compare.py measures our layout the same way). Every structure, arrow and caption
    joins the row of the nearest arrow axis; a row's ink is the union of its structure ink,
    arrow boxes (axis ± head half-width) and caption ink."""
    doc = load(name) if isinstance(name, str) else name
    m = measure(doc)
    structs = _structures(doc)
    arrows = sorted(m["arrows"], key=lambda a: min(a["tail"][0], a["head"][0]))
    axes = sorted({round(a["tail"][1], 2) for a in arrows})

    def nearest(y):
        return min(axes, key=lambda ax: abs(y - ax))

    items = {ax: [] for ax in axes}  # (kind, label, ink)
    for fid, s in structs.items():
        items[nearest((s["ink"][1] + s["ink"][3]) / 2)].append(("structure", fid, s["ink"]))
    for a in arrows:
        hw = max(h["half_width_left"] + h["half_width_right"] for h in a["heads"].values() if h) / 2
        x0, x1 = sorted((a["tail"][0], a["head"][0]))
        items[nearest(a["tail"][1])].append(("arrow", f"{a['length']:.2f}", (x0, a["tail"][1] - hw, x1, a["tail"][1] + hw)))
    for cap in doc.captions:
        glyphs = _caption_glyphs(doc, cap)
        if glyphs:
            ink = union(ln["ink"] for ln in text_lines(glyphs))
            # numbers sit below their row: assign by the structure above, i.e. the nearest axis above the ink top
            above = [ax for ax in axes if ax <= ink[1]]
            bold = all(ln["bold"] for ln in text_lines(glyphs))
            ax = max(above) if bold and above else nearest((ink[1] + ink[3]) / 2)
            items[ax].append(("number" if bold else "caption", cap["text"], ink))
    rows = []
    for ax in axes:
        its = sorted(items[ax], key=lambda it: it[2][0])
        ink = union(it[2] for it in its)
        sink = union(it[2] for it in its if it[0] == "structure")
        row_arrows = [it for it in its if it[0] == "arrow"]
        rows.append({"axis_y": ax, "ink": ink, "structure_ink": sink, "width": ink[2] - ink[0],
                     "sequence": [f"{k}:{lbl}" for k, lbl, _ in its if k in ("structure", "arrow")],
                     "first": its[0][0], "last": [it for it in its if it[0] in ("structure", "arrow")][-1][0],
                     "arrow_lengths": [float(lbl) for _, lbl, _ in row_arrows],
                     "first_arrow_tail_minus_row_left": row_arrows[0][2][0] - ink[0] if row_arrows else None})
    pairs = []
    for r1, r2 in zip(rows, rows[1:]):
        pairs.append({"axis_pitch": r2["axis_y"] - r1["axis_y"],
                      "ink_gap (row bottom incl. numbers -> next row top incl. conditions)": r2["ink"][1] - r1["ink"][3],
                      "structure_ink_gap": r2["structure_ink"][1] - r1["structure_ink"][3],
                      "left_edge_diff (row2 - row1)": r2["ink"][0] - r1["ink"][0],
                      "right_edge_diff (row2 - row1)": r2["ink"][2] - r1["ink"][2],
                      "center_diff (row2 - row1)": (r2["ink"][0] + r2["ink"][2] - r1["ink"][0] - r1["ink"][2]) / 2})
    return _r({"file": doc.name, "overall_ink": union(r["ink"] for r in rows), "rows": rows, "between_rows": pairs})


RESONANCE_REFS = {"1.0": "resonance_x1.0", "0.8": "resonance_x0.8", "0.6": "resonance_x0.6"}


def _charge_box(c):
    """Ink box of a filled charge (bars of bar_length x bar_thickness; '+' adds the vertical bar)."""
    (x, y), h, t = c["center"], c["bar_length"] / 2, c["bar_thickness"] / 2
    return (x - h, y - (h if c["sign"] == "+" else t), x + h, y + (h if c["sign"] == "+" else t))


def resonance_summary(name: str | Doc) -> dict:
    """Resonance parameters (acs-style.md §4.3) from one resonance_x<k> file: the double-headed
    arrows (exact heads; hand-drawn length), structure ink (bonds, labels and charges: each
    charge joins the structure of its atom) -> arrow end gaps, the bracket gaps to the content
    ink, and the contributors' ring centers against the arrow axis. `name` may be a loaded Doc."""
    doc = load(name) if isinstance(name, str) else name
    m = measure(doc)
    structs = _structures(doc)
    frag_of = {a.id: a.frag for a in doc.atoms.values()}
    for c in m["charges"]:
        atom = min(doc.atoms.values(), key=lambda a: length(sub(a.p, c["center"])))
        s = structs[frag_of[atom.id]]
        s["ink_without_charges"] = s.get("ink_without_charges", s["ink"])
        s["ink"] = union([s["ink"], _charge_box(c)])
    order = sorted(structs.values(), key=lambda s: s["ink"][0])
    arrows = sorted(m["arrows"], key=lambda a: min(a["tail"][0], a["head"][0]))
    recs = []
    for a in arrows:
        x0, x1 = sorted((a["tail"][0], a["head"][0]))
        left = max((s for s in order if s["ink"][2] <= x0), key=lambda s: s["ink"][2])
        right = min((s for s in order if s["ink"][0] >= x1), key=lambda s: s["ink"][0])
        heads = [h for h in a["heads"].values() if h]
        recs.append({"length": a["length"], "axis_y": a["tail"][1],
                     "heads": {"length": stats([h["length"] for h in heads]),
                               "half_width": stats([(h["half_width_left"] + h["half_width_right"]) / 2 for h in heads]),
                               "inset": stats([h["inset_from_tip"] for h in heads])},
                     "shaft": a.get("shaft"),
                     "gap_left_ink_to_end": x0 - left["ink"][2], "gap_end_to_right_ink": right["ink"][0] - x1,
                     "gap_end_to_right_ink_without_charges": right.get("ink_without_charges", right["ink"])[0] - x1})
    brackets = m["brackets"]
    out = {"arrows": recs,
           "arrow_length": stats([r["length"] for r in recs]),
           "arrow_gap (structure ink incl. charges <-> arrow end)":
               stats([r["gap_left_ink_to_end"] for r in recs] + [r["gap_end_to_right_ink"] for r in recs]),
           "ring_center_minus_axis": [round(c[1] - recs[0]["axis_y"], 3) for s in order for c in s["ring_centers"]],
           "brackets": [{k: b[k] for k in ("lips_point", "height", "stroke_width", "lip_length", "gap_to_content_x",
                                            "extends_above_content", "extends_below_content")} for b in brackets],
           "bracket gap_x": stats([b["gap_to_content_x"] for b in brackets]),
           "bracket gap_above": stats([b["extends_above_content"] for b in brackets]),
           "bracket gap_below": stats([b["extends_below_content"] for b in brackets])}
    return _r(out)


EQUILIBRIUM_REFS = {"1.0": "equilibrium_x1.0", "0.8": "equilibrium_x0.8", "0.6": "equilibrium_x0.6",
                    "1.0 large": "arrows_large_x1.0"}


def equilibrium_summary(name: str | Doc) -> dict:
    """Equilibrium arrows and head variants (acs-style.md §3.1, §3.6) from one reference: per
    arrow its CDXML head attributes (HeadSize / ArrowheadCenterSize / ArrowheadWidth say which
    variant), then for equilibrium arrows the kind (balanced, or unbalanced favoring products or
    reactants: the long shaft's tip on the right or the left), shaft spacing and stagger, the
    short shaft's ratio and centering, and the half heads (forward shaft: the one pointing at
    Head3D; back: the other). Full heads (forward and ↔ arrows) as _arrow_entry measures them.
    Conditions (hand-placed) are provenance only: from the outer shafts and the arrow's center."""
    doc = load(name) if isinstance(name, str) else name
    m = measure(doc)
    by_axis = {round(a["axis_y"], 1): a for a in m["scheme"]["arrows"]}
    recs = []
    for a in sorted(m["arrows"], key=lambda a: a["tail"][1]):
        c = a["cdxml"]
        rec = {"head_variant": f"HeadSize {c['HeadSize']} / CenterSize {c['ArrowheadCenterSize']} / "
                               f"Width {c['ArrowheadWidth']}"}
        e = a.get("equilibrium")
        if not e:
            hs = [h for h in a["heads"].values() if h]
            rec.update({"kind": "resonance" if len(hs) == 2 else "forward", "length": a["length"],
                        "head_length": stats([h["length"] for h in hs]),
                        "head_half_width": stats([(h["half_width_left"] + h["half_width_right"]) / 2 for h in hs]),
                        "head_inset": stats([h["inset_from_tip"] for h in hs]), "shaft": a.get("shaft")})
            recs.append(rec)
            continue
        sh = e["shafts"]
        fwd = next(r for r in sh if r["direction"] == "forward")
        back = next(r for r in sh if r["direction"] == "back")
        rec["kind"] = e["kind"]
        if e["kind"] == "unbalanced":
            long = max(sh, key=lambda r: r["length"])
            tip_x = add(a["tail"], mul(unit(sub(a["head"], a["tail"])), long["tip_u"]))[0]
            rec["kind"] += " favoring " + ("products" if tip_x > (a["tail"][0] + a["head"][0]) / 2 else "reactants")
            rec.update({k: e[k] for k in ("long_length", "short_length", "short_center_minus_long_center",
                                          "long_over_short_butt_to_notch")})
        else:
            rec.update({k: e[k] for k in ("stagger_tail", "stagger_head")})
        rec.update({
            "spacing (center to center)": e["spacing"], "shafts_center_minus_axis": e["center_v"],
            "shaft_width": [r["width"] for r in sh], "barb": [r["barb"] for r in sh],
            "half_head_length": [r["head"]["length"] for r in sh],
            "half_head_barb_from_center": [r["head"]["barb_from_center"] for r in sh],
            "half_head_tip_from_center": [r["head"]["tip_from_center"] for r in sh],
            "half_head_inset": [r["head"]["inset_from_tip"] for r in sh],
            "shaft_end_from_tip": [r["u_end_from_tip"] for r in sh],
            "forward_back_control_per_barb": fwd["head"]["back_control_per_barb"],
            "back_back_control_from_center": back["head"]["back_control_from_center"],
            "back_control_shift (forward - back)": round(fwd["head"]["back_control_from_center"]
                                                         - back["head"]["back_control_from_center"], 3),
        })
        s = by_axis.get(round(min(a["tail"][1], a["head"][1]), 1)) or min(
            m["scheme"]["arrows"], key=lambda r: abs(r["axis_y"] - a["tail"][1]))
        rec["conditions (provenance)"] = {
            "above: top shaft - last baseline": s.get("above_axis_minus_last_base"),
            "below: first baseline - bottom shaft": s.get("below_first_base_minus_axis"),
            "line centers - arrow center": [ln["center_minus_tail_notch_mid"] for ln in s["above"] + s["below"]],
            "widest pad (each end)": s.get("widest") and round((s["length"] - s["widest"]["advance"]) / 2, 3)}
        recs.append(rec)
    return _r({"arrows": recs})


MECH_REFS = {"1.0": "mech_basic_x1.0", "0.8": "mech_basic_x0.8", "0.6": "mech_basic_x0.6"}


def _box_gap(p, box):
    """Distance from p to an axis-aligned box (0 inside)."""
    dx = max(box[0] - p[0], 0.0, p[0] - box[2])
    dy = max(box[1] - p[1], 0.0, p[1] - box[3])
    return math.hypot(dx, dy)


def _dot_sets(doc):
    """Lone pairs (mutually nearest circles) and single dots (radicals), each with its atom:
    the nearest atom center to the pair midpoint / dot."""
    circles = [s for s in doc.shapes if s.kind == "circle"]

    def nearest(i):
        return min((j for j in range(len(circles)) if j != i),
                   key=lambda j: length(sub(circles[i].center, circles[j].center)), default=None)
    mate = {i: nearest(i) for i in range(len(circles))}
    spacing = min(length(sub(circles[i].center, circles[j].center)) for i, j in mate.items() if j is not None)
    pairs, singles, seen = [], [], set()
    for i, j in mate.items():
        if i in seen:
            continue
        a = circles[i]
        if j is not None and mate[j] == i and length(sub(a.center, circles[j].center)) < 1.2 * spacing:
            seen |= {i, j}
            pairs.append((a, circles[j]))
        else:
            seen.add(i)
            singles.append(a)
    return pairs, singles


def _elem_box(doc, atom):
    elem, _ = label_glyphs(doc, atom) if atom.label else ([], [])
    return union(g.box for g in elem) if elem else None


def mech_summary(name: str) -> dict:
    """Mechanism parameters from one mech_basic file (acs-style.md §4): exact dot and
    curved-arrow geometry, and the hand-placed positions (lone-pair gaps, radical distances,
    curved-arrow tails, heads and shape)."""
    doc = load(name)
    m = measure(doc)
    pairs, singles = _dot_sets(doc)
    out = {"dot_radius": stats([c.r for p in pairs for c in p] + [c.r for c in singles]),
           "dot_spacing": stats([length(sub(a.center, b.center)) for a, b in pairs])}
    lp = []
    for a, b in pairs:
        mid = mul(add(a.center, b.center), 0.5)
        atom = min(doc.atoms.values(), key=lambda at: length(sub(at.p, mid)))
        box = _elem_box(doc, atom)
        if box is None:
            continue
        gap = min(_box_gap(a.center, box), _box_gap(b.center, box)) - a.r
        d = sub(mid, _box_center(box))
        side = (("right" if d[0] > 0 else "left") if abs(d[0]) * (box[3] - box[1]) > abs(d[1]) * (box[2] - box[0])
                else ("below" if d[1] > 0 else "above"))
        lp.append({"atom": atom.label, "side": side, "gap_dot_edge_to_elem_inkbox": round(gap, 3),
                   "mid_offset_from_elem_center_along_side": round(abs(d[0]) if side in ("left", "right") else abs(d[1]), 3)})
    out["lone_pairs"] = lp
    out["lone_pair_gap"] = stats([e["gap_dot_edge_to_elem_inkbox"] for e in lp])
    rad = []
    for c in singles:
        atom = min(doc.atoms.values(), key=lambda at: length(sub(at.p, c.center)))
        box = _elem_box(doc, atom)
        d = sub(c.center, atom.p)
        rec = {"atom": atom.label or "C", "dir_deg": round(math.degrees(math.atan2(-d[1], d[0])), 1)}
        if box:
            rec["gap_dot_edge_to_elem_inkbox"] = round(_box_gap(c.center, box) - c.r, 3)
        else:
            rec["center_dist_from_atom"] = round(length(d), 3)
        rad.append(rec)
    out["radicals"] = rad
    lp_mids = [mul(add(a.center, b.center), 0.5) for a, b in pairs]
    bond_mids = [mul(add(doc.atoms[b.a].p, doc.atoms[b.b].p), 0.5) for b in doc.bonds]
    arrows = []
    for c in m["curved_arrows"]:
        pts = c["svg_path_pts"]
        p0, c1, c2, p3 = pts[:4]
        tip = add(p3, c["head"]["tip_minus_stroke_end"])
        chord = length(sub(p3, p0))
        u = unit(sub(p3, p0))
        nrm = (-u[1], u[0])
        rec = {"type": "fishhook" if c["head"]["n_vertices"] == 3 else "two_electron",
               "chord": round(chord, 3), "control_offset_over_chord": c["control_offset_max_over_chord"],
               "handle_lengths_over_chord": [round(length(sub(c1, p0)) / chord, 3), round(length(sub(c2, p3)) / chord, 3)],
               "handle_offsets_over_chord": [round(dot(sub(c1, p0), nrm) / chord, 3), round(dot(sub(c2, p0), nrm) / chord, 3)],
               "handle_along_over_chord": [round(dot(sub(c1, p0), u) / chord, 3), round(dot(sub(c2, p0), u) / chord, 3)]}
        near = min(range(len(lp_mids)), key=lambda i: length(sub(lp_mids[i], p0)))
        k = min(range(len(bond_mids)), key=lambda i: length(sub(bond_mids[i], p0)))
        if length(sub(lp_mids[near], p0)) < length(sub(bond_mids[k], p0)):  # from a lone pair, else a bond
            a, b = pairs[near]
            owner = min(doc.atoms.values(), key=lambda at: length(sub(at.p, lp_mids[near])))
            box = _elem_box(doc, owner)
            out_dir = unit(sub(lp_mids[near], _box_center(box) if box else owner.p))
            rec["tail_from_lone_pair_mid"] = round(length(sub(p0, lp_mids[near])), 3)
            rec["tail_gap_to_dot_edge"] = round(min(length(sub(p0, a.center)), length(sub(p0, b.center))) - a.r, 3)
            # along the pair's outward direction, beyond the dots' outer edge
            rec["tail_beyond_dot_edge_outward"] = round(dot(sub(p0, lp_mids[near]), out_dir) - a.r, 3)
        else:
            bd = doc.bonds[k]
            fr = frame(doc, bd.a, bd.b)
            along, off = fr.local(p0)
            rec["tail_from_bond_mid"] = round(length(sub(p0, bond_mids[k])), 3)
            rec["tail_offset_from_bond_axis"] = round(abs(off), 3)
            rec["tail_along_from_bond_mid"] = round(along - fr.L / 2, 3)
        atom = min(doc.atoms.values(), key=lambda at: length(sub(at.p, tip)))
        box = _elem_box(doc, atom)
        rec["head_target"] = atom.label or "C"
        if box:
            rec["head_tip_gap_to_elem_inkbox"] = round(_box_gap(tip, box), 3)
        else:
            rec["head_tip_from_atom"] = round(length(sub(tip, atom.p)), 3)
        arrows.append(rec)
    out["curved_arrows"] = arrows
    ca = m["curved_arrows"]
    out["curved_arrow_stroke"] = stats([c["stroke_width"] for c in ca])
    full = [c["head"] for c in ca if c["head"]["n_vertices"] == 4]
    half = [c["head"] for c in ca if c["head"]["n_vertices"] == 3]
    out["full_head"] = {"length": stats([h["length"] for h in full]),
                        "half_width": stats([(h["extent_left"] + h["extent_right"]) / 2 for h in full]),
                        "inset": stats([-max(u for u, v in h["pts_local"] if abs(v) < 0.05 and u < -0.5) for h in full])}
    out["full_head"]["stroke_end_off_axis"] = stats([h["stroke_end_local"][1] for h in full])
    out["fishhook_head"] = {"stroke_end_off_axis_toward_barb": stats(
                                [h["stroke_end_local"][1] * (1 if h["extent_left"] > h["extent_right"] else -1) for h in half]),
                            "length": stats([h["length"] for h in half]),
                            "barb": stats([max(h["extent_left"], h["extent_right"]) for h in half]),
                            "inset": stats([-max(u for u, v in h["pts_local"] if abs(v) < 0.05 and u < -0.5) for h in half])}
    return out


ABBREV_REFS = ["abbrev_x1.0", "abbrev_x0.8", "abbrev_x0.6"]


def abbrev_summary(name: str) -> dict:
    """§3.5: per abbreviation (abbrev_x*), each drawn form: displayed text, ChemDraw's
    LabelJustification, which glyph's advance center sits on the atom (its offset), baseline −
    atom y, and each bond end with its residual against the §2.1 model (preset values) on the
    anchor glyphs label_glyphs picks."""
    from acs_chem_svg import draw_mol
    from acs_chem_svg.preset import load_preset

    doc = load(name)
    asg = assign(doc)
    k = name.rsplit("_x", 1)[1]
    preset = load_preset("acs1996" if k == "1.0" else f"acs1996_x{k}")
    L, r = preset.lengths, preset.ratios
    margin = float(doc.settings["MarginWidth"])
    forms, offsets, baselines, residuals = [], [], [], []
    for a in doc.atoms.values():
        if not a.nickname:
            continue
        elem, whole = label_glyphs(doc, a)
        text = "".join(g.ch for g in whole)
        right = elem[-1] is whole[-1] and elem[0] is not whole[0]
        c = whole[-1] if right else whole[0]
        outline = [[(x - a.p[0], y - a.p[1]) for x, y in cont] for g in elem for cont in g.outline()]
        ends = []
        for b in doc.bonds:
            if a.id not in (b.a, b.b) or not asg.get(b.id, {}).get("lines"):
                continue
            fr = frame(doc, b.b if b.a == a.id else b.a, a.id)
            end = statistics.median(fr.L - _line_from_shape(ln.shape, fr).u1 for ln in asg[b.id]["lines"])
            model = draw_mol.clip_distance_band(outline, (-fr.u[0], -fr.u[1]), r.label_band_per_size * c.size, margin,
                                                r.label_clearance_per_margin,
                                                (L.label_outline_offset_x, L.label_outline_offset_y))
            ends.append({"dir_deg": round(math.degrees(math.atan2(fr.u[1], -fr.u[0])), 1), "end": _r(end),
                         "model_minus_end": _r(model - end)})
            residuals.append(model - end)
        off = c.x + c.adv / 2 - a.p[0]
        offsets.append((right, off))
        baselines.append(c.base - a.p[1])
        forms.append({"key": a.label, "text": text, "justification": a.justification, "centered": c.ch,
                      "anchor": "".join(g.ch for g in elem), "centered_dx": _r(off),
                      "baseline_minus_atom_y": _r(c.base - a.p[1]), "bonds": ends})
    return {"forms": forms,
            "centered_dx_left_justified": stats(o for rt, o in offsets if not rt),
            "centered_dx_right_justified": stats(o for rt, o in offsets if rt),
            "baseline_minus_atom_y": stats(baselines),
            "bond_end_model_minus_end": stats(residuals),
            "bond_ends_over_0.1": sum(abs(x) > 0.1 for x in residuals)}


# --------------------------------------------------------------------------- main


CDXML_SYMBOL_REFS = ["labels_x1.0", "m1_set_x1.0", "m1_set_x0.8", "m1_set_x0.6",
                     "mech_basic_x1.0", "mech_basic_x0.8", "mech_basic_x0.6",
                     "resonance_x1.0", "resonance_x0.8", "resonance_x0.6"]


def cdxml_symbols(name: str) -> dict:
    """Charge symbols (<graphic GraphicType="Symbol" SymbolType="Plus|Minus">) of one reference
    (§13): the BoundingBox's first point vs the drawn sign's center (_charges, SVG mapped to the
    CDXML frame), and the second point as (first x - dx, first y + dy). dx_per_size2 is
    dx / LabelSize^2. The scheme references draw no charges."""
    doc = load(name)
    root = etree.parse(str(REF_DIR / f"{name}.cdxml")).getroot()
    size = float(root.get("LabelSize"))
    drawn = [sub(c["center"], doc.offset) for c in measure(doc)["charges"]]
    out = []
    for g in root.iter("graphic"):
        if g.get("GraphicType") != "Symbol" or g.get("SymbolType") not in ("Plus", "Minus"):
            continue
        x0, y0, x1, y1 = (float(v) for v in g.get("BoundingBox").split())
        near = min(drawn, key=lambda c: length(sub(c, (x0, y0))), default=None)
        out.append({"sign": g.get("SymbolType"), "first": (x0, y0),
                    "first_minus_drawn_center": tuple(round(v, 3) for v in sub((x0, y0), near)) if near else None,
                    "dx": round(x0 - x1, 3), "dy": round(y1 - y0, 3), "dx_per_size2": round((x0 - x1) / size ** 2, 4)})
    res = {"label_size": size, "symbols": out}
    dots = [sub(c.center, doc.offset) for c in doc.shapes if c.kind == "circle"]
    if dots:
        res.update(_cdxml_dots(root, dots))
    return res


CHARGES_LEWIS_REFS = ["charges_lewis_x1.0", "charges_lewis_x0.8", "charges_lewis_x0.6"]


def charges_lewis(name: str) -> dict:
    """charges_lewis_x* (the user's edit of cdxml_spike.py --lewis-draft, 2026-10-28; §2.2,
    §13.4). Per symbol the user added with ChemDraw's tools (circled and plain signs, Electron
    dots): SymbolType, its attachment, BoundingBox dx, and the ink ChemDraw drew around the box's
    first point: ring radius (stroke center) and stroke width, bar length and thickness, dot
    radius. Also each node with Radical or Charge set and its label runs (face 64 = superscript),
    and the free <t> texts."""
    doc = load(name)
    root = etree.parse(str(REF_DIR / f"{name}.cdxml")).getroot()

    def center(sh):
        if sh.center:
            return sh.center
        xs, ys = [q[0] for q in sh.pts], [q[1] for q in sh.pts]
        return ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2)
    symbols = []
    for g in root.iter("graphic"):
        kind = g.get("SymbolType")
        if g.get("GraphicType") != "Symbol" or kind not in ("CirclePlus", "CircleMinus", "Plus", "Minus", "Electron"):
            continue
        x0, y0, x1, y1 = (float(v) for v in g.get("BoundingBox").split())
        c = (x0 + doc.offset[0], y0 + doc.offset[1])
        ink = [sh for sh in doc.shapes if length(sub(center(sh), c)) < 0.05]
        rec = {"symbol": kind, "first": (x0, y0), "dx": round(x0 - x1, 3),
               "attached": None if g.find("represent") is None else g.find("represent").get("attribute")}
        for sh in ink:
            if sh.kind == "stroke":
                rec["ring_radius"], rec["ring_width"] = round(sh.r, 3), round(sh.width, 3)
            elif sh.kind == "circle":
                rec["dot_radius"] = round(sh.r, 3)
            else:
                xs, ys = [q[0] for q in sh.pts], [q[1] for q in sh.pts]
                if max(xs) - min(xs) > max(ys) - min(ys):
                    rec["bar_length"], rec["bar_thickness"] = round(max(xs) - min(xs), 3), round(max(ys) - min(ys), 3)
        symbols.append(rec)
    nodes = [{"id": n.get("id"), "Element": n.get("Element"), "Charge": n.get("Charge"), "Radical": n.get("Radical"),
              "Warning": n.get("Warning"),
              "runs": [(s.text, s.get("face")) for s in n.iter("s")]}
             for n in root.iter("n") if n.get("Radical") or n.get("Charge") or n.get("Warning")]
    texts = [{"text": "".join(s.text or "" for s in t.iter("s")), "faces": [s.get("face") for s in t.iter("s")],
              "size": [s.get("size") for s in t.iter("s")]} for t in root.find("page").findall("t")]
    return {"label_size": float(root.get("LabelSize")), "symbols": symbols, "nodes": nodes, "free_texts": texts}


def _cdxml_dots(root, dots: list) -> dict:
    """Lone-pair and radical symbols (§13.2) vs the drawn dot centers (CDXML frame).
    LonePair: BoundingBox = the two dot centers, counterclockwise on screen around the atom
    (`ccw`: the cross product of (first - atom) and (second - atom) is negative, y down).
    Electron: (dot center, center x - dx, center y). Per fragment with dots and per side, how
    far its BoundingBox reaches beyond the outermost dot center, and that dot's kind: the reach
    is cdxml.lone_pair_box_reach / electron_box_reach where the dot bounds the box, more where
    other ink does."""
    def near(p):
        return min(dots, key=lambda c: length(sub(c, p)))
    nodes = {n.get("id"): _xy(n.get("p"), (0.0, 0.0)) for n in root.iter("n")}
    pairs, electrons, reach = [], [], []
    for frag in root.iter("fragment"):
        centers = []
        for g in frag.findall("graphic"):
            kind = g.get("SymbolType")
            if kind not in ("LonePair", "Electron"):
                continue
            rep = g.find("represent")
            atom = nodes[rep.get("object")]
            x0, y0, x1, y1 = (float(v) for v in g.get("BoundingBox").split())
            if kind == "LonePair":
                a, b = (x0, y0), (x1, y1)
                da, db = sub(a, atom), sub(b, atom)
                pairs.append({"attribute": rep.get("attribute"), "separation": round(length(sub(a, b)), 3),
                              "ccw": da[0] * db[1] - da[1] * db[0] < 0,
                              "first_minus_drawn": tuple(round(v, 3) for v in sub(a, near(a))),
                              "second_minus_drawn": tuple(round(v, 3) for v in sub(b, near(b)))})
                centers += [(a, "pair"), (b, "pair")]
            else:
                c = (x0, y0)
                electrons.append({"attribute": rep.get("attribute"), "dx": round(x0 - x1, 3), "dy": round(y1 - y0, 3),
                                  "first_minus_drawn": tuple(round(v, 3) for v in sub(c, near(c)))})
                centers.append((c, "electron"))
        if not centers:
            continue
        bb = [float(v) for v in frag.get("BoundingBox").split()]
        # per side: how far the box reaches beyond the outermost dot center, and that dot's kind
        rec = {"fragment": frag.get("id")}
        for side, k, sgn in (("left", 0, -1), ("top", 1, -1), ("right", 0, 1), ("bottom", 1, 1)):
            c, kind = max(centers, key=lambda ck: sgn * ck[0][k])
            rec[side] = (round(sgn * (bb[k + (2 if sgn > 0 else 0)] - c[k]), 3), kind)
        reach.append(rec)
    return {"lone_pairs": pairs, "electrons": electrons, "fragment_reach_beyond_dot_centers": reach}


CDXML_SCHEME_REFS = [f"{k}_x{s}" for k in ("scheme_basic", "scheme_multi", "scheme_wrap", "equilibrium", "resonance")
                     for s in ("1.0", "0.8", "0.6")] + ["arrows_large_x1.0"]


def _cdxml_caption(doc: Doc, cap: dict, t) -> dict:
    """One caption <t> (or a formula fragment's label) against its drawn glyphs, CDXML frame."""
    glyphs = _caption_glyphs(doc, cap)
    lines = text_lines(glyphs)
    first = lines[0]
    center = (first["x0"] + first["x1"]) / 2
    ink = union(g.box for g in glyphs)
    bb = cap["bbox"]
    rec = {"text": cap["text"], "justification": cap["justification"] or "Left", "formula_fragment": cap["formula"],
           "p_minus_first_origin": (round(cap["p"][0] - first["x0"], 3), round(cap["p"][1] - first["base"], 3)),
           "p_x_minus_first_line_center": round(cap["p"][0] - center, 3),
           "line_centers_minus_first": [round((ln["x0"] + ln["x1"]) / 2 - center, 3) for ln in lines],
           "box_minus_ink": [round(a - b, 3) for a, b in zip(bb, ink)],
           "line_pitch": [round(b["base"] - a["base"], 3) for a, b in zip(lines, lines[1:])],
           "runs": [(r[0], r[2]) for r in cap["runs"]],
           "attrs": {k: t.get(k) for k in ("Warning", "LineHeight", "LineStarts", "CaptionJustification",
                                           "InterpretChemically") if t.get(k) is not None}}
    if len(lines) > 1:
        rec["box_right_minus_line_ends"] = [round(bb[2] - ln["x1"], 3) for ln in lines]
        rec["box_right_minus_line_ink"] = [round(bb[2] - ln["ink"][2], 3) for ln in lines]
    return rec


def _cdxml_arrow(a, lw: float) -> dict:
    """One <arrow>: the axis points against Head3D / Tail3D, and its BoundingBox against the
    ends and the half width from its attributes (ArrowheadWidth x LineWidth / 100, plus half the
    ArrowShaftSpacing x LineWidth / 100)."""
    def v(k):
        return [float(x) for x in a.get(k).split()]
    h, t, c, M, m, bb = v("Head3D"), v("Tail3D"), v("Center3D"), v("MajorAxisEnd3D"), v("MinorAxisEnd3D"), v("BoundingBox")
    d = (h[0] - t[0], h[1] - t[1])
    half = int(a.get("ArrowheadWidth")) * lw / 100 + int(a.get("ArrowShaftSpacing") or 0) * lw / 200
    return {
        "attrs": {k: a.get(k) for k in ("HeadSize", "ArrowheadCenterSize", "ArrowheadWidth", "ArrowheadHead",
                                        "ArrowheadTail", "ArrowShaftSpacing", "ArrowEquilibriumRatio") if a.get(k)},
        "head3d_right_of_tail3d": h[0] > t[0], "length": round(math.hypot(*d), 3),
        "center_minus_mid": (round(c[0] - (h[0] + t[0]) / 2, 3), round(c[1] - (h[1] + t[1]) / 2, 3)),
        "major_minus_(center + head - tail)": (round(M[0] - c[0] - d[0], 3), round(M[1] - c[1] - d[1], 3)),
        "minor_minus_(center + rot90(head - tail))": (round(m[0] - c[0] + d[1], 3), round(m[1] - c[1] - d[0], 3)),
        "box_x_minus_ends": (round(bb[0] - min(h[0], t[0]), 3), round(bb[2] - max(h[0], t[0]), 3)),
        "half_width_from_attrs": round(half, 3),
        "box_pad (top, bottom)": (round(h[1] - bb[1] - half, 3), round(bb[3] - h[1] - half, 3)),
    }


def cdxml_scheme(name: str) -> dict:
    """CDXML conventions of a scheme reference (§13.1): captions (<t> on the page and formula
    fragments, _cdxml_caption), arrows (_cdxml_arrow, with the legacy <graphic> each supersedes),
    and the <scheme><step> structure: what kind of object each step cites, and the page order."""
    doc = load(name)
    off = doc.offset
    root = etree.parse(str(REF_DIR / f"{name}.cdxml")).getroot()
    page = root.find("page")
    lw = float(root.get("LineWidth"))
    texts = [t for t in page.iter("t") if t.getparent().tag != "n" or t.getparent().get("NodeType") in ("Fragment", "Nickname")]
    caps = []
    for cap in doc.captions:
        if not _caption_glyphs(doc, cap):
            continue
        p = sub(cap["p"], off)
        t = min(texts, key=lambda t: length(sub(tuple(float(x) for x in t.get("p").split()), p)))
        caps.append(_cdxml_caption(doc, cap, t))
    by_id = {el.get("id"): el for el in page.iter() if el.get("id")}
    arrows = []
    for a in page.iter("arrow"):
        rec = _cdxml_arrow(a, lw)
        g = next((g for g in page.iter("graphic") if g.get("SupersededBy") == a.get("id")), None)
        rec["graphic"] = None if g is None else {
            "ArrowType": g.get("ArrowType"), "HeadSize": g.get("HeadSize"), "same_Z": g.get("Z") == a.get("Z"),
            "box_is_head3d_tail3d": [float(x) for x in g.get("BoundingBox").split()]
            == [float(x) for x in a.get("Head3D").split()[:2] + a.get("Tail3D").split()[:2]]}
        arrows.append(rec)

    def kind(i):
        el = by_id.get(i)
        if el is None:
            return "missing"
        if el.tag == "fragment" and el.find("n") is not None and el.find("n").get("NodeType") in ("Fragment", "Nickname"):
            return "fragment (formula)"
        return el.tag
    schemes = []
    for sch in page.iter("scheme"):
        schemes.append([{k.replace("ReactionStep", ""): [kind(i) for i in st.get(k).split()]
                         for k in ("ReactionStepReactants", "ReactionStepProducts", "ReactionStepArrows",
                                   "ReactionStepObjectsAboveArrow", "ReactionStepObjectsBelowArrow") if st.get(k)}
                        for st in sch.findall("step")])
    order = [el.tag for el in page]
    return {"line_width": lw, "caption_size": float(root.get("CaptionSize")), "captions": caps, "arrows": arrows,
            "brackets": _cdxml_brackets(doc, page, lw), "schemes": schemes,
            "page_order": [k for k, _ in itertools.groupby(order)]}


def _cdxml_brackets(doc, page, lw: float) -> list[dict]:
    """Each <graphic GraphicType="Bracket"> (resonance_x*, §13.3) against the drawn bracket at
    its x (_brackets, CDXML frame): its BoundingBox as (x, y, x, y') against the bar's ends,
    which end it lists first, LipSize against the drawn lip length / LineWidth, and any
    objecttag."""
    drawn = _brackets(copy.deepcopy(doc))
    out = []
    for g in page.iter("graphic"):
        if g.get("GraphicType") != "Bracket":
            continue
        bb = [float(x) for x in g.get("BoundingBox").split()]
        d = min(drawn, key=lambda b: abs(b["x"] - doc.offset[0] - bb[0]))
        top, bottom = d["y_top"] - doc.offset[1], d["y_bottom"] - doc.offset[1]
        lip = d["lip_length"]["median"]
        out.append({
            "attrs": {k: g.get(k) for k in ("BracketType", "LipSize")},
            "lips_point": d["lips_point"], "first_y_is": "top" if bb[1] < bb[3] else "bottom",
            "x_minus_bar_x": round(bb[0] - (d["x"] - doc.offset[0]), 3),
            "box_ys_minus_bar_ends": (round(min(bb[1], bb[3]) - top, 3), round(max(bb[1], bb[3]) - bottom, 3)),
            "lip_length": lip, "lip_length_per_line_width_x100": round(lip / lw * 100, 1),
            "objecttags": [t.get("Name") for t in g.findall("objecttag")]})
    return out


CHAIR_REFS = ["chair_x1.0", "chair_x0.8", "chair_x0.6"]


def _bond_z(name: str) -> tuple[dict, set]:
    """CDXML bond id -> Z, and the CrossingBonds pairs ChemDraw stores (frozensets of bond ids), for
    drawn bonds (a nickname's internal bonds cross each other: scheme_multi_x*)."""
    root = etree.parse(str(REF_DIR / f"{name}.cdxml")).getroot()
    drawn = [b for b in root.iter("b") if b.getparent().getparent().tag != "n"]  # not nickname internals (OBn)
    z = {b.get("id"): int(b.get("Z")) for b in drawn}
    pairs = {frozenset((b.get("id"), c)) for b in drawn for c in (b.get("CrossingBonds") or "").split() if c in z}
    return z, pairs


def _poly_dist(p, pts) -> float:
    return min(seg_point_dist(p, a, b) for a, b in zip(pts, pts[1:] + pts[:1]))


def crossing_summary(name: str, doc: Doc | None = None) -> dict:
    """Bond crossings (§14.2). Every CrossingBonds pair: the bond with the higher Z is drawn whole
    (front), the other (rear) is cut. Per side of the front bond's axis, the rear polygon corners
    nearest the axis (the cut ends; unclaimed slivers in the rear's band count, so the stub left
    between a label and a gap shows): their distance to the front axis (`beyond_axis`, the same for
    both corners of an end when the cut is parallel to the front bond) and to the front bond's ink.
    `half_extent` is the front bond's drawn half-width at the crossing: lw/2 (single), the outer
    line's center offset (double), the wedge outline's half-width there (wedge). `doc`: measure
    this drawing instead (compare.py's render of the templates, same bond ids), with front and
    rear as in the reference."""
    if doc:
        doc = copy.deepcopy(doc)  # assign() claims shapes; measure() may have claimed them already
        for sh in doc.shapes:
            sh.owner = None
    else:
        doc = load(name)
    asg = assign(doc)
    z, pairs = _bond_z(name)
    bonds = {b.id: b for b in doc.bonds}
    lw, m = float(doc.settings["LineWidth"]), float(doc.settings["MarginWidth"])
    spacing = float(doc.settings["BondLength"]) * float(doc.settings["BondSpacing"]) / 100
    rows = []
    for pair in sorted(pairs, key=lambda p: sorted(p)):
        b1, b2 = (bonds[i] for i in sorted(pair))
        front, rear = (b1, b2) if z[b1.id] > z[b2.id] else (b2, b1)
        ff, rf = frame(doc, front.a, front.b), frame(doc, rear.a, rear.b)
        angle = math.degrees(math.acos(min(1.0, abs(dot(ff.u, rf.u)))))
        fa = asg.get(front.id, {})
        front_ink = [ln.shape.pts for ln in fa.get("lines", [])] + [sh.pts for sh, _ in fa.get("shapes", [])]
        # the crossing point, in the front frame
        den = ff.u[0] * rf.u[1] - ff.u[1] * rf.u[0]
        d = sub(rf.o, ff.o)
        t = (d[0] * rf.u[1] - d[1] * rf.u[0]) / den
        if fa.get("lines"):
            half = max(max(abs(ln.v) for ln in fa["lines"]), lw / 2)
        else:  # wedge: outline half-width at the crossing, between its two ends
            loc = [q for _, lc in fa["shapes"] for q in lc]
            near = [q for q in loc if q[0] < ff.L / 2]
            far = [q for q in loc if q[0] >= ff.L / 2]
            if near and far:
                (u0, h0), (u1, h1) = (max(near, key=lambda q: abs(q[1])), max(far, key=lambda q: abs(q[1])))
                half = abs(h0) + (abs(h1) - abs(h0)) * (t - u0) / (u1 - u0)
            else:
                half = None
        rear_polys = [ln.shape for ln in asg.get(rear.id, {}).get("lines", [])]
        for sh in doc.shapes:
            if sh.kind == "poly" and not sh.owner:
                loc = [rf.local(q) for q in sh.pts]
                if all(-2 <= u <= rf.L + 2 and abs(v) <= spacing + lw for u, v in loc):
                    rear_polys.append(sh)
        sides = {}
        for sh in rear_polys:
            for q in sh.pts:
                v = ff.local(q)[1]
                sides.setdefault(1 if v > 0 else -1, []).append((abs(v), q))
        ends = {}
        for sd, pts in sides.items():
            near = min(v for v, _ in pts)
            corners = [(v, q) for v, q in pts if v < near + 0.05]
            ends["+" if sd > 0 else "-"] = {
                "beyond_axis": stats([v for v, _ in corners]),
                "beyond_ink": stats([_poly_dist(q, ink) for _, q in corners for ink in [min(
                    front_ink, key=lambda k: _poly_dist(q, k))]]) if front_ink else None,
                "beyond_half_extent_per_margin": round((statistics.median(v for v, _ in corners) - half) / m, 3)
                if half is not None else None,
            }
        rows.append({"front": front.id, "rear": rear.id, "front_z": z[front.id], "rear_z": z[rear.id],
                     "front_order": front.order, "front_display": front.display or "plain",
                     "rear_order": rear.order, "angle_deg": round(angle, 2),
                     "front_half_extent": round(half, 3) if half is not None else None,
                     "rear_pieces": len(rear_polys), "cut": ends})
    return {"margin_width": m, "line_width": lw, "crossings": rows}


def _ring6(doc, frag: str) -> list[str] | None:
    """The six-ring atom ids of a fragment, in ring order starting at the first CDXML atom."""
    ring = [a for a, at in doc.atoms.items() if at.frag == frag
            and sum(1 for b in doc.bonds if b.ring and a in (b.a, b.b)) == 2]
    if len(ring) != 6:
        return None
    order = [ring[0]]
    while len(order) < 6:
        nxt = [n for n in neighbors(doc, order[-1]) if n in ring and n not in order]
        order.append(nxt[0] if len(order) == 1 else nxt[-1])
    return order


def chair_summary(name: str) -> dict:
    """Chairs (§14.1). Fragments of chair_x* in CDXML order: two bare chairs (ChemDraw's
    Conformers templates: `chair` and `flipped`), axial methyls on all six atoms, equatorial methyls
    on all six, β-D-glucopyranose, an axial OAc, then the crossing cases (crossing_summary).
    Per template: vertex offsets from vertex 1 (pt, y down) and bond lengths. Per vertex of a
    substituted chair (matched to a template by its offsets): the substituent's direction (deg,
    counterclockwise from +x, y up) and length, and for equatorials the angle to the ring bond one
    removed (the bonds parallel to it in a textbook chair). Crossing bonds in the chairs: which is
    in front (higher Z) and the mean screen y of each (larger y = lower on the page)."""
    doc = load(name)
    z, pairs = _bond_z(name)
    bl = float(doc.settings["BondLength"])
    chairs = [(fid, r) for fid in doc.fragments for r in [_ring6(doc, fid)] if r]
    templates = {}
    for key, (fid, ring) in zip(("chair", "flipped"), chairs[:2]):
        o = doc.atoms[ring[0]].p
        templates[key] = {"fragment": fid, "atoms": ring,
                          "vertices": [tuple(round(v, 3) for v in sub(doc.atoms[a].p, o)) for a in ring],
                          "bond_lengths": [round(frame(doc, a, b).L, 3) for a, b in zip(ring, ring[1:] + ring[:1])]}

    def match(ring):
        """(template, vertex offset of ring[0] in it, direction): the template whose vertices the
        ring's positions repeat (within 0.03 pt: two 0.01 pt CDXML roundings) under a translation."""
        pts = [doc.atoms[a].p for a in ring]
        for key, t in templates.items():
            tv = t["vertices"]
            for start in range(6):
                for dr in (1, -1):
                    o = sub(pts[0], tv[start])
                    if all(length(sub(sub(pts[k], o), tv[(start + dr * k) % 6])) < 0.03 for k in range(6)):
                        return key, start, dr
        return None

    def ang(v):
        return round(math.degrees(math.atan2(-v[1], v[0])), 2)

    subs = []
    for fid, ring in chairs[2:]:
        mt = match(ring)
        if not mt:
            subs.append({"fragment": fid, "error": "no template match"})
            continue
        key, start, dr = mt
        verts = []
        for k, a in enumerate(ring):
            vtx = (start + dr * k) % 6
            for s in neighbors(doc, a):
                if s in ring:
                    continue
                v = sub(doc.atoms[s].p, doc.atoms[a].p)
                far = [(ring[(k + 1) % 6], ring[(k + 2) % 6]), (ring[(k - 1) % 6], ring[(k - 2) % 6])]
                par = min(angle_between(v, sub(doc.atoms[q].p, doc.atoms[p].p)) for p, q in far)
                verts.append({"vertex": vtx + 1, "atom": a, "substituent": s, "label": doc.atoms[s].label,
                              "dir_deg": ang(v), "length": round(length(v), 3),
                              "vertical": abs(v[0]) < 0.02, "up": v[1] < 0,
                              "deg_from_ring_bond_one_removed": round(par, 2),
                              "element": doc.atoms[a].element})
        hetero = [{"vertex": (start + dr * k) % 6 + 1, "element": doc.atoms[a].element}
                  for k, a in enumerate(ring) if doc.atoms[a].element != "C"]
        cross = []
        for p in pairs:
            ids = sorted(p)
            bs = [next(b for b in doc.bonds if b.id == i) for i in ids]
            if not any(b.frag == fid for b in bs):
                continue
            front = max(bs, key=lambda b: z[b.id])
            cross.append({"front": front.id, "rear": next(b.id for b in bs if b is not front),
                          "mean_y": {b.id: round((doc.atoms[b.a].p[1] + doc.atoms[b.b].p[1]) / 2, 2) for b in bs},
                          "ring_bond": {b.id: b.ring for b in bs}})
        subs.append({"fragment": fid, "template": key, "ring_atom_1_at_vertex": start + 1, "direction": dr,
                     "substituents": verts, "heteroatoms": hetero, "crossings": cross})
    return {"bond_length": bl, "templates": templates,
            "templates_per_bond_length": {k: [tuple(round(c / bl, 4) for c in v) for v in t["vertices"]]
                                          for k, t in templates.items()},
            "substituted": subs}


def angle_between(a, b) -> float:
    """Angle (deg, 0-90) between the lines along a and b."""
    c = abs(dot(unit(a), unit(b)))
    return math.degrees(math.acos(min(1.0, c)))


def all_names():
    return sorted(p.stem for p in REF_DIR.glob("*.cdxml") if (REF_DIR / f"{p.stem}.svg").exists())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*")
    ap.add_argument("--json", type=Path)
    ap.add_argument("--summary", action="store_true", help="print the key-parameter table only")
    ap.add_argument("--fit-clip", action="store_true", help="fit the bond-to-label clip model (§2.1)")
    ap.add_argument("--scheme", action="store_true", help="scheme parameters per preset (§3)")
    ap.add_argument("--mech", action="store_true", help="mechanism parameters per preset (§4)")
    ap.add_argument("--wrap", action="store_true", help="multi-row scheme parameters per preset (§3.4)")
    ap.add_argument("--resonance", action="store_true", help="resonance arrow and bracket parameters per preset (§4.3)")
    ap.add_argument("--abbrev", action="store_true", help="abbreviation label parameters per preset (§3.5)")
    ap.add_argument("--equilibrium", action="store_true",
                    help="equilibrium arrows and arrowhead variants per preset (§3.1, §3.6)")
    ap.add_argument("--cdxml-symbols", action="store_true", help="charge-symbol <graphic> boxes per reference (§13)")
    ap.add_argument("--cdxml-scheme", action="store_true", help="CDXML captions, arrows and steps per reference (§13.1)")
    ap.add_argument("--charges-lewis", action="store_true",
                    help="circled signs, symbols and two-electron radicals in charges_lewis_x* (§2.2, §13.4)")
    ap.add_argument("--chair", action="store_true", help="chair templates, axial and equatorial bonds (§14.1)")
    ap.add_argument("--crossing", action="store_true", help="bond crossings: front by Z, rear-bond cuts (§14.2)")
    args = ap.parse_args(argv)
    if args.chair:
        print(json.dumps({n: chair_summary(n) for n in args.names or CHAIR_REFS}, indent=1, ensure_ascii=False))
        return 0
    if args.crossing:
        print(json.dumps({n: crossing_summary(n) for n in args.names or CHAIR_REFS}, indent=1, ensure_ascii=False))
        return 0
    if args.charges_lewis:
        names = args.names or CHARGES_LEWIS_REFS
        print(json.dumps({n: charges_lewis(n) for n in names}, indent=1, ensure_ascii=False))
        return 0
    if args.cdxml_scheme:
        names = args.names or CDXML_SCHEME_REFS
        print(json.dumps({n: cdxml_scheme(n) for n in names}, indent=1, ensure_ascii=False))
        return 0
    if args.cdxml_symbols:
        names = args.names or CDXML_SYMBOL_REFS
        print(json.dumps({n: cdxml_symbols(n) for n in names}, indent=1, ensure_ascii=False))
        return 0
    if args.equilibrium:
        refs = {k: n for k, n in EQUILIBRIUM_REFS.items() if not args.names or n in args.names}
        print(json.dumps({k: equilibrium_summary(n) for k, n in refs.items()}, indent=1, ensure_ascii=False))
        return 0
    if args.abbrev:
        names = args.names or ABBREV_REFS
        print(json.dumps({n: abbrev_summary(n) for n in names}, indent=1, ensure_ascii=False))
        return 0
    if args.resonance:
        refs = {k: n for k, n in RESONANCE_REFS.items() if not args.names or n in args.names}
        print(json.dumps({k: resonance_summary(n) for k, n in refs.items()}, indent=1, ensure_ascii=False))
        return 0
    if args.wrap:
        refs = {k: n for k, n in WRAP_REFS.items() if not args.names or n in args.names}
        print(json.dumps({k: wrap_summary(n) for k, n in refs.items()}, indent=1, ensure_ascii=False))
        return 0
    if args.mech:
        refs = {k: n for k, n in MECH_REFS.items() if not args.names or n in args.names}
        print(json.dumps({k: mech_summary(n) for k, n in refs.items()}, indent=1, ensure_ascii=False))
        return 0
    if args.scheme:
        names = args.names or [n for ns in SCHEME_REFS.values() for n in ns]
        print(json.dumps(scheme_summary({n: measure(load(n)) for n in names}), indent=1, ensure_ascii=False))
        return 0
    if args.fit_clip:
        fit = fit_label_clip(clip_rows(args.names or None))
        for r in fit.pop("rows"):
            print(f"{r['file']:19} {r['kind']:5} {r['element']:3} dir {r['dir_deg']:7} end {r['end']:6.3f} "
                  f"residual {r['residual']:+.3f}{'  !' if abs(r['residual']) > 0.1 else ''}")
        print(json.dumps(fit, indent=1))
        return 0
    results = {n: measure(load(n)) for n in (args.names or all_names())}
    if args.summary:
        rows = {n: summarize(m) for n, m in results.items()}
        keys = list(dict.fromkeys(k for r in rows.values() for k in r))
        names = list(rows)
        print("parameter".ljust(34) + "".join(n[:18].rjust(19) for n in names))
        for k in keys:
            print(k.ljust(34) + "".join(str(rows[n].get(k, "")).rjust(19) for n in names))
        return 0
    text = json.dumps(results, indent=1, default=str)
    if args.json:
        args.json.write_text(text)
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
