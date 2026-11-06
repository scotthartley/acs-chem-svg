"""The reference model and loaders: CDXML atoms and bonds, SVG ink, bond frames, label glyphs."""

from __future__ import annotations

import math
import re
import statistics
from dataclasses import dataclass, field
from pathlib import Path

import pypdfium2
from lxml import etree
from rdkit import Chem as _Chem

from acs_chem_svg.fontmetrics import load_face
from acs_chem_svg.geometry import add, dot, length, mul, sub, unit


REPO = Path(__file__).resolve().parents[2]
REF_DIR = REPO / "reference_renders" / "chemdraw"
# The ACS presets and the scale in their reference file names (<family>_x<scale>): one table, for
# the scripts and tests that go from a preset to its references or back. Presets made from a style
# probe (scripts/new_preset.py) have no references and are not listed.
ACS_PRESETS = {"acs1996": "1.0", "acs1996_x0.8": "0.8", "acs1996_x0.6": "0.6"}
ELEMENTS = {n: _Chem.GetPeriodicTable().GetElementSymbol(n) for n in range(1, 104)}  # CDXML Element -> symbol


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
    alignment: str = ""  # CDXML LabelAlignment (Above / Below: a label's tail on the line above / below, ring NTs)
    node_type: str = ""  # CDXML NodeType (Nickname, Fragment, Unspecified; "" for an element)


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
            # a label typed into a structure that ChemDraw expanded (metals_x* PPh3) is Unspecified
            # with the expansion inside: drawn and justified as a nickname
            # (or ChemDraw could not interpret it: L, Ln, dppf, XPhos, Bpin, p-Tol in abbrev_2_x* and
            # ligands_x* are Unspecified with no inner fragment; they are placed the same way)
            nick = n.get("NodeType") in ("Fragment", "Nickname") or (
                n.get("NodeType") == "Unspecified" and t is not None)
            if nick:  # the attachment atom of the nickname sits at the node (OMe: O)
                el = next((ELEMENTS.get(int(m.get("Element", "6")), "?") for m in n.iter("n")
                           if m is not n and m.get("p") == n.get("p")), "X")
            doc.atoms[n.get("id")] = Atom(
                n.get("id"), _xy(n.get("p"), off), el, label,
                int(n.get("Charge", "0")), int(n.get("NumHydrogens")) if n.get("NumHydrogens") else None,
                fid, nick, 1 if n.get("Radical") == "Doublet" else 0,
                t.get("LabelJustification", "") if t is not None else "",
                t.get("LabelAlignment", "") if t is not None else "", n.get("NodeType", ""))
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


def _kern(left: Glyph, right: Glyph) -> float:
    """Pair kerning between two neighboring glyphs of one run style (T|s in Ts: fontmetrics,
    as resvg and ChemDraw apply it); 0 across a size or face change."""
    if left.size != right.size or left.bold != right.bold or left.italic != right.italic:
        return 0.0
    face = load_face("Bold" if left.bold else "Oblique" if left.italic else "Regular")
    return face.kerning(left.ch, right.ch, left.size)


def _glyph_chain(doc: Doc, seed: list) -> list:
    """`seed` extended both ways by glyphs set contiguously with it on (about) its baseline,
    pair kerning included (Ts, OTs, NHTs: abbrev_2_x*)."""
    whole = list(seed)
    tol = 0.05
    changed = True
    while changed:
        changed = False
        first, last = whole[0], whole[-1]
        for g in doc.glyphs:
            if g in whole or abs(g.base - first.base) > 0.35 * first.size:
                continue
            if abs(g.x - (last.x + last.adv + _kern(last, g))) < tol:
                whole.append(g)
                changed = True
                break
            if abs(g.x + g.adv + _kern(g, first) - first.x) < tol:
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
    if atom.alignment == "Right":  # bond leaving to the right: the last glyph is on the atom (P(o-tol)₃'s ₃)
        side = "left"
    elif atom.alignment in ("Left", "Above", "Below"):
        side = "right"
    else:
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


def angle_between(a, b) -> float:
    """Angle (deg, 0-90) between the lines along a and b."""
    c = abs(dot(unit(a), unit(b)))
    return math.degrees(math.acos(min(1.0, c)))


def all_names():
    return sorted(p.stem for p in REF_DIR.glob("*.cdxml") if (REF_DIR / f"{p.stem}.svg").exists())
