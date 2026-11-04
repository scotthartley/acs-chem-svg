"""CDXML conventions: charge-symbol boxes, Lewis references, captions, arrows, steps (acs-style.md §13)."""

from __future__ import annotations

import copy
import itertools
import math

from lxml import etree

from acs_chem_svg.geometry import length, sub

from .core import Doc, REF_DIR, _xy, load, union
from .measurements import _brackets, _caption_glyphs, measure, text_lines


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
    """charges_lewis_x* (the user's edit of cdxml_gate.py --lewis-draft, 2026-10-28; §2.2,
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
