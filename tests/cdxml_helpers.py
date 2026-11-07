"""Shared by the CDXML export tests (test_cdxml_*.py; style-spec.md §13–§13.4): the DTD check,
cached exports, and the node-pairing helpers the reference diffs use.

- DTD: reference_renders/cdxml/CDXML.dtd (Revvity's file, fetched once from the DOCTYPE URL;
  ChemDraw.app doesn't ship it). Local only, like the references: those tests skip without it. It
  dates from 2003, so libxml2 needs one enumeration loosened (AS lists +, - and ?, which aren't
  name tokens), and ChemDraw 26 writes two attributes it doesn't declare (DTD_GAPS).
"""

import io
import re
from functools import lru_cache

import pytest
from compare import cited_ids, preset_for, reference_cdx_scheme, render_reference
from lxml import etree
from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions

import measure_refs as mr
from chem_svg_cdxml.cdxml_writer import to_cdxml
from chem_svg_cdxml.figure import spec_cdxml
from chem_svg_cdxml.chem import ABBREVIATIONS
from chem_svg_cdxml.draw_mol import Drawer
from chem_svg_cdxml.geometry import sub
from chem_svg_cdxml.spec import ResonanceStep, load_spec, resolve_molfile


DTD_PATH = mr.REPO / "reference_renders" / "cdxml" / "CDXML.dtd"  # Revvity's: local only, never committed
DTD_GAPS = {"No declaration for attribute ShowPrimeLabels of element CDXML",
            "No declaration for attribute IsNickname of element t"}
TOL = 0.1
CHARGE_TOL = 1.0  # hand-placed sign centers (preset charge.labeled_offset_*: round defaults)
PAIR_TOL = 0.01
CAPTION_PAIR_TOL = 0.02  # a caption's p from its SVG glyph origin: scheme_wrap_x1.0's LiAlH4 is 97.305, stored 97.30
SPECS = [f"examples/m1/m1_{k}_x{s}.yaml" for k in ("smiles", "charges") for s in ("1.0", "0.8", "0.6")]
ABBREV_REFS = [f"abbrev_x{s}" for s in ("1.0", "0.8", "0.6")]
STYLE_ATTRS = ("BondLength", "LineWidth", "BoldWidth", "BondSpacing", "HashSpacing", "MarginWidth", "LabelSize",
               "CaptionSize", "LabelFont", "CaptionFont", "LabelFace", "ChainAngle")
PARSER = etree.XMLParser(load_dtd=False, no_network=True)
needs_dtd = pytest.mark.skipif(not DTD_PATH.exists(), reason="CDXML.dtd (reference_renders/cdxml/) is local only")
needs_refs = pytest.mark.skipif(not mr.REF_DIR.exists(), reason="ChemDraw references (reference_renders/) are local only")


@lru_cache(maxsize=None)
def dtd() -> etree.DTD:
    text = DTD_PATH.read_text()

    def loosen(m):  # an enumeration with a value that isn't a name token -> CDATA
        return m[0] if all(re.fullmatch(r"[\w.:-]+", v.strip()) for v in m[1].split("|")) else "CDATA"
    return etree.DTD(io.StringIO(re.sub(r"\(([^()]*\|[^()]*)\)(?=\s+(?:\"|#))", loosen, text)))


def dtd_errors(root) -> list[str]:
    d = dtd()
    d.validate(root)
    return [e.message for e in d.error_log if e.message not in DTD_GAPS]


@lru_cache(maxsize=None)
def exported(spec: str) -> str:
    return spec_cdxml(mr.REPO / spec, warn=lambda m: None)


@lru_cache(maxsize=None)
def reference_export(ref: str) -> str:
    """The reference's templates exported at ChemDraw's CDXML coordinates (the templates hold
    its SVG frame, svg = cdxml + Doc.offset)."""
    r = render_reference(ref)
    off = (-r.ref.offset[0], -r.ref.offset[1])
    return to_cdxml([(m, Drawer(m, r.preset), off) for m in r.molecules], r.preset.name, f"{ref}.cdxml",
                    warn=lambda m: None)


STEP_KINDS = {"ReactionStepReactants": ("fragment",), "ReactionStepProducts": ("fragment",),
              "ReactionStepArrows": ("graphic",),  # above / below: captions, and species a mechanism step adds
              "ReactionStepObjectsAboveArrow": ("t", "fragment"), "ReactionStepObjectsBelowArrow": ("t", "fragment")}


def _check_structure(text: str) -> None:
    root = etree.fromstring(text.encode(), PARSER)
    assert dtd_errors(root) == []
    # object ids; font ids are their own namespace (clearance_x1.0 has font 21 and bond 21)
    els = {el.get("id"): el for el in root.find("page").iter() if el.get("id")}
    ids = [el.get("id") for el in root.find("page").iter() if el.get("id")]
    assert len(ids) == len(set(ids)), "duplicate ids"
    nodes = {n.get("id") for n in root.iter("n")}
    for b in root.iter("b"):
        assert {b.get("B"), b.get("E")} <= nodes, f"bond {b.get('id')} points to a missing node"
    for r in root.iter("represent"):
        assert r.get("object") in nodes, f"represent points to missing node {r.get('object')}"
    for st in root.iter("step"):  # every id a step cites exists and is of the right kind
        for attr, tags in STEP_KINDS.items():
            for i in (st.get(attr) or "").split():
                assert i in els and els[i].tag in tags, f"step {st.get('id')} {attr}: {i} is not a <{tags}>"
    for g in root.iter("graphic"):
        if g.get("SupersededBy"):
            a = els.get(g.get("SupersededBy"))
            assert a is not None and a.tag == "arrow" and a.get("Z") == g.get("Z")
    page = root.find("page")  # the references' page order
    order = ["fragment", "t", "graphic", "curve", "scheme", "arrow"]
    tags = [order.index(el.tag) for el in page]
    assert tags == sorted(tags)


def _canonical(smiles: str) -> str:
    m = Chem.MolFromSmiles(smiles)
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(m)


def _xy(s):
    return tuple(float(v) for v in s.split()[:2])


def _xy4(s):
    return tuple(float(v) for v in s.split())


def _nodes(root):
    """Fragment nodes (not nickname internals) by id; mech_basic's tBu-Cl is in a <group>."""
    return {n.get("id"): n for f in root.iter("fragment") if f.getparent().tag != "n" for n in f.findall("n")}


def _bonds(root, nodes):
    """frozenset(node ids) -> (Order, display, narrow-end node id); WedgeEnd -> Begin from E."""
    out = {}
    for b in root.iter("b"):
        if b.get("B") not in nodes:
            continue
        disp = b.get("Display", "")
        narrow = None
        if disp:
            narrow = b.get("B") if disp.endswith("Begin") else b.get("E")
            disp = disp.replace("End", "Begin")
        out[frozenset((b.get("B"), b.get("E")))] = (b.get("Order", "1"), disp, narrow)
    return out


def _charges(root):
    return {g.find("represent").get("object"): (g.get("SymbolType"), [float(v) for v in g.get("BoundingBox").split()])
            for g in root.iter("graphic") if g.get("SymbolType") in ("Plus", "Minus")}


def _text(t):
    return "".join(s.text or "" for s in t.iter("s"))


INNER_ATTRS = ("Element", "Charge", "NodeType")  # and NumHydrogens on heteroatoms: a carbon's is optional


def _inner_pairs(theirs, ours, key: str) -> tuple[list[str], list[tuple[str, str, float, float]]]:
    """A nickname's inner fragment (ours None: not written): nodes paired by p (INNER_ATTRS
    equal, p as "nickname atoms"), bonds by their paired ends. scheme_multi_x1.0's inner nodes
    also carry Z, NeedsClean and ExternalConnectionNum, and its methyl carbons no
    NumHydrogens="3", unlike abbrev_x*'s: ChemDraw's bookkeeping, not chemistry."""
    if ours is None:
        return [f"{key}: no inner fragment"], []
    bad, nums, pair = [], [], {}
    for o in ours.findall("n"):
        m = [t for t in theirs.findall("n") if _near(_xy(t.get("p")), _xy(o.get("p")), TOL)]
        if len(m) != 1:
            bad.append(f"{key} inner node at {o.get('p')}: {len(m)} ChemDraw nodes")
            continue
        pair[o.get("id")] = m[0].get("id")
        attrs = INNER_ATTRS + (("NumHydrogens",) if m[0].get("Element") else ())
        if {k: o.get(k) for k in attrs} != {k: m[0].get(k) for k in attrs}:
            bad.append(f"{key} inner node {m[0].get('id')}: ChemDraw {dict(m[0].attrib)} ours {dict(o.attrib)}")
        for k, (a, b) in enumerate(zip(_xy(m[0].get("p")), _xy(o.get("p")))):
            nums.append(("nickname atoms", f"{key} inner node {m[0].get('id')} {'xy'[k]}", a, b))
    if len(pair) != len(theirs.findall("n")):
        bad.append(f"{key}: {len(pair)} inner nodes paired, ChemDraw has {len(theirs.findall('n'))}")

    def bonds(f, ids):
        return {frozenset(ids.get(b.get(e), b.get(e)) for e in "BE"): b.get("Order", "1") for b in f.findall("b")}
    if bonds(ours, pair) != bonds(theirs, {}):
        bad.append(f"{key}: inner bonds differ")
    return bad, nums


def reference_pairs(ref: str, text: str | None = None) -> tuple[list[str], list[tuple[str, str, float, float]]]:
    """(mismatches, [(category, key, ChemDraw, ours)] for numbers); `text` defaults to
    reference_export(ref)."""
    theirs = etree.parse(str(mr.REF_DIR / f"{ref}.cdxml"), PARSER).getroot()
    ours = etree.fromstring((text or reference_export(ref)).encode(), PARSER)
    tn, on = _nodes(theirs), _nodes(ours)
    to_ref = {}
    bad, nums = [], []
    for oid, o in on.items():
        p = _xy(o.get("p"))
        match = [tid for tid, t in tn.items() if max(abs(a - b) for a, b in zip(_xy(t.get("p")), p)) < PAIR_TOL]
        if len(match) != 1:
            bad.append(f"node at {p}: {len(match)} ChemDraw nodes")
            continue
        to_ref[oid] = match[0]
    for oid, tid in to_ref.items():
        o, t = on[oid], tn[tid]
        key = f"node {tid}"
        for attr, default in (("Element", "6"), ("Charge", "0")):
            if o.get(attr, default) != t.get(attr, default):
                bad.append(f"{key} {attr}: ChemDraw {t.get(attr, default)} ours {o.get(attr, default)}")
        ot, tt = o.find("t"), t.find("t")
        if (ot is None) != (tt is None):
            bad.append(f"{key}: label {'missing' if ot is None else 'extra'}")
            continue
        if tt is None:
            continue
        if o.get("NumHydrogens") != t.get("NumHydrogens"):
            bad.append(f"{key} NumHydrogens: ChemDraw {t.get('NumHydrogens')} ours {o.get('NumHydrogens')}")
        for attr in ("LabelJustification", "Justification"):
            if ot.get(attr) != tt.get(attr):
                bad.append(f"{key} {attr}: ChemDraw {tt.get(attr)} ours {ot.get(attr)}")
        if _text(ot) != _text(tt):
            bad.append(f"{key} text: ChemDraw {_text(tt)} ours {_text(ot)}")
        # LabelAlignment on multi-character labels only: 28 of ChemDraw's 282 one-glyph labels have it
        for el, ref_el, attr in ((o, t, "NodeType"), (ot, tt, "IsNickname"), (ot, tt, "LabelAlignment")):
            if el.get(attr) != ref_el.get(attr) and (attr != "LabelAlignment" or len(_text(tt)) > 1):
                bad.append(f"{key} {attr}: ChemDraw {ref_el.get(attr)} ours {el.get(attr)}")
        if t.find("fragment") is not None:
            b2, n2 = _inner_pairs(t.find("fragment"), o.find("fragment"), f"{key} {_text(tt)}")
            bad += b2
            nums += n2
        for k, (a, b) in enumerate(zip(_xy(tt.get("p")), _xy(ot.get("p")))):
            nums.append(("label p", f"{key} {'xy'[k]} {_text(tt)}", a, b))
        for k, (a, b) in enumerate(zip(_xy4(tt.get("BoundingBox")), _xy4(ot.get("BoundingBox")))):
            nums.append(("label box", f"{key} [{k}] {_text(tt)}", a, b))
    # every node of a fragment we drew is paired (templates skip only nickname fragments)
    drawn = {tn[t].getparent() for t in to_ref.values()}
    missing = [tid for tid, t in tn.items() if t.getparent() in drawn and tid not in to_ref.values()]
    if missing:
        bad.append(f"ChemDraw nodes not exported: {missing}")
    tb = _bonds(theirs, tn)
    for pair, (order, disp, narrow) in _bonds(ours, on).items():
        tpair = frozenset(to_ref.get(i) for i in pair)
        if tpair not in tb:
            bad.append(f"bond {sorted(tpair)}: not in ChemDraw's file")
            continue
        t_order, t_disp, t_narrow = tb[tpair]
        if (order, disp) != (t_order, t_disp) or (narrow and to_ref[narrow] != t_narrow):
            bad.append(f"bond {sorted(tpair)}: ChemDraw {t_order} {t_disp} from {t_narrow}, "
                       f"ours {order} {disp} from {to_ref.get(narrow)}")
    for f in ours.iter("fragment"):  # the drawn ink box; a hand-placed charge sign or dot moves it
        if f.getparent().tag == "n":  # a nickname's inner fragment: _inner_pairs
            continue
        t = tn[to_ref[f.find("n").get("id")]].getparent()
        syms = {g.get("SymbolType") for g in f.findall("graphic")}
        cat = ("fragment box (with dots)" if syms & {"LonePair", "Electron"}
               else "fragment box (with a charge)" if syms else "fragment box")
        for k, (a, b) in enumerate(zip(_xy4(t.get("BoundingBox")), _xy4(f.get("BoundingBox")))):
            nums.append((cat, f"fragment {t.get('id')} [{k}]", a, b))
    tc, oc = _charges(theirs), _charges(ours)
    if sorted(tc) != sorted(to_ref[i] for i in oc):
        bad.append(f"charge symbols on ChemDraw nodes {sorted(tc)}, ours {sorted(to_ref[i] for i in oc)}")
    for oid, (sym, box) in oc.items():
        tsym, tbox = tc.get(to_ref[oid], (None, None))
        if tsym != sym:
            bad.append(f"charge on node {to_ref[oid]}: ChemDraw {tsym} ours {sym}")
            continue
        for k in range(2):
            nums.append(("charge center", f"node {to_ref[oid]} {'xy'[k]}", tbox[k], box[k]))
        nums.append(("charge box dx", f"node {to_ref[oid]}", tbox[0] - tbox[2], box[0] - box[2]))
    return bad, nums


def _floats(s):
    return [float(v) for v in s.split()]


@lru_cache(maxsize=None)
def reference_scheme_export(ref: str) -> str:
    """compare.reference_cdx_scheme written out."""
    return to_cdxml(reference_cdx_scheme(ref), preset_for(ref), f"{ref}.cdxml", warn=lambda m: None)


def _near(p, q, tol=PAIR_TOL):
    return max(abs(a - b) for a, b in zip(p, q)) < tol


def scheme_pairs(ref: str) -> tuple[list[str], list[tuple[str, str, float, float]]]:
    """(mismatches, [(category, key, ChemDraw, ours)]) for the scheme objects."""
    theirs = etree.parse(str(mr.REF_DIR / f"{ref}.cdxml"), PARSER).getroot().find("page")
    ours = etree.fromstring(reference_scheme_export(ref).encode(), PARSER).find("page")
    bad, nums = [], []
    # arrows, paired by Head3D
    o2t = {}
    t_arrows = list(theirs.iter("arrow"))
    for oa in ours.iter("arrow"):
        m = [ta for ta in t_arrows if _near(_floats(ta.get("Head3D")), _floats(oa.get("Head3D")))]
        if len(m) != 1:
            bad.append(f"arrow at {oa.get('Head3D')}: {len(m)} ChemDraw arrows")
            continue
        ta = m[0]
        og = next(g for g in ours.iter("graphic") if g.get("SupersededBy") == oa.get("id"))
        tg = next(g for g in theirs.iter("graphic") if g.get("SupersededBy") == ta.get("id"))
        o2t[og.get("id")] = tg.get("id")
        key = f"arrow {ta.get('id')}"
        for attr in ("FillType", "ArrowheadHead", "ArrowheadTail", "ArrowheadType", "HeadSize", "ArrowheadCenterSize",
                     "ArrowheadWidth", "ArrowShaftSpacing", "ArrowEquilibriumRatio"):
            if oa.get(attr) != ta.get(attr):
                bad.append(f"{key} {attr}: ChemDraw {ta.get(attr)} ours {oa.get(attr)}")
        for attr in ("ArrowType", "HeadSize"):
            if og.get(attr) != tg.get(attr):
                bad.append(f"{key} graphic {attr}: ChemDraw {tg.get(attr)} ours {og.get(attr)}")
        for attr in ("Head3D", "Tail3D", "Center3D", "MajorAxisEnd3D", "MinorAxisEnd3D", "BoundingBox"):
            cat = "arrow box" if attr == "BoundingBox" else "arrow points"
            for k, (a, b) in enumerate(zip(_floats(ta.get(attr)), _floats(oa.get(attr)))):
                nums.append((cat, (k, f"{key} {attr}[{k}]"), a, b))
        for k, (a, b) in enumerate(zip(_floats(tg.get("BoundingBox")), _floats(og.get("BoundingBox")))):
            nums.append(("arrow points", (k, f"{key} graphic box[{k}]"), a, b))
    # captions, paired by p: their <t> or formula fragment (text frags by drawn origin: the label's p)
    doc = mr.load(ref)
    tobj = {t.get("id"): t for t in theirs.iter("t") if t.getparent().tag != "n"}  # mech_basic: in <group>s
    for f in theirs.iter("fragment"):  # scheme_wrap_x*: some in <group>s
        if f.get("id") in doc.text_frags:
            tobj[f.get("id")] = f.find("n").find("t")
    # the lines of their stacked captions as drawn: we write one <t> per line where a line has a
    # sub- or superscript (style-spec.md §13.1), each paired with its line: p on its drawn
    # baseline, box its ink
    their_lines = {}
    for cap in doc.captions:
        if "\n" not in cap["text"]:
            continue
        p = sub(cap["p"], doc.offset)
        tid = next(i for i, t in tobj.items() if _near(_floats(t.get("p")), p))
        for k, ln in enumerate(mr.text_lines(mr._caption_glyphs(doc, cap))):
            ink = (*sub(ln["ink"][:2], doc.offset), *sub(ln["ink"][2:], doc.offset))
            their_lines[(tid, k)] = ((p[0], ln["base"] - doc.offset[1]), ink, cap["text"].split("\n")[k])
    split = {}
    for ot in ours.findall("t"):
        m = [i for i, t in tobj.items() if _near(_floats(t.get("p")), _floats(ot.get("p")), CAPTION_PAIR_TOL)]
        if len(m) == 1 and "\n" in _text(tobj[m[0]]) and "\n" not in _text(ot):
            m = []  # the first line of a split caption: paired with its line below
        if not m:
            lm = [k for k, (q, _, _) in their_lines.items() if _near(q, _floats(ot.get("p")), CAPTION_PAIR_TOL)]
            if len(lm) == 1:
                (tid, k), (q, ink, text) = lm[0], their_lines[lm[0]]
                o2t[ot.get("id")] = tid
                split.setdefault(tid, []).append(k)
                key = f"caption {tid} line {k + 1} {text!r}"
                if _text(ot) != text:
                    bad.append(f"{key} text: ours {_text(ot)!r}")
                if ot.get("Justification") != tobj[tid].get("Justification"):
                    bad.append(f"{key} Justification: ChemDraw {tobj[tid].get('Justification')} ours {ot.get('Justification')}")
                for j, (a, b) in enumerate(zip(q, _floats(ot.get("p")))):
                    nums.append(("caption p", (j, f"{key} p[{j}]"), a, b))
                for j, (a, b) in enumerate(zip(ink, _floats(ot.get("BoundingBox")))):
                    nums.append(("caption box", (j, f"{key} box[{j}]"), a, b))
                continue
        if len(m) != 1:
            bad.append(f"caption {_text(ot)!r} at {ot.get('p')}: {len(m)} ChemDraw captions")
            continue
        tt = tobj[m[0]]
        o2t[ot.get("id")] = m[0]
        key = f"caption {m[0]} {_text(tt)!r}"
        if _text(ot) != _text(tt):
            bad.append(f"{key} text: ours {_text(ot)!r}")
        for attr in ("Justification", "CaptionJustification", "LineStarts"):
            if ot.get(attr) != tt.get(attr):
                bad.append(f"{key} {attr}: ChemDraw {tt.get(attr)} ours {ot.get(attr)}")
        for k, (a, b) in enumerate(zip(_floats(tt.get("p")), _floats(ot.get("p")))):
            nums.append(("caption p", (k, f"{key} p[{k}]"), a, b))
        stacked = "\n" in _text(tt)
        formula_sup = any(s.get("face") == "96" and re.search(r"[A-Za-z0-9)][+\u2212-]$", s.text or "")
                          for s in tt.iter("s"))
        cat = "stacked caption box" if stacked else "formula superscript box" if formula_sup else "caption box"
        for k, (a, b) in enumerate(zip(_floats(tt.get("BoundingBox")), _floats(ot.get("BoundingBox")))):
            nums.append((cat, (k, f"{key} box[{k}]"), a, b))
    for tid, ks in split.items():
        if sorted(ks) != list(range(len(_text(tobj[tid]).split("\n")))):
            bad.append(f"caption {tid}: lines {sorted(ks)} written, ChemDraw has {_text(tobj[tid]).count(chr(10)) + 1}")
    # fragments, paired by their first node
    tnodes = {tuple(_floats(n.get("p"))[:2]): f.get("id") for f in theirs.iter("fragment") if f.getparent().tag != "n"
              for n in f.findall("n")}
    for f in ours.findall("fragment"):
        p = _floats(f.find("n").get("p"))
        m = [i for q, i in tnodes.items() if _near(q, p)]
        if len(m) == 1:
            o2t[f.get("id")] = m[0]
        else:
            bad.append(f"fragment {f.get('id')}: {len(m)} ChemDraw fragments")
    written = set(o2t.values())

    def their_ids(st, attr):  # a cited <group>: what is in it (compare.cited_ids)
        return cited_ids(theirs, st, attr)
    # steps, paired by arrow; their ids we didn't write are left out
    tsteps = {st.get("ReactionStepArrows"): st for st in theirs.iter("step")}
    osteps = list(ours.iter("step"))
    if len(osteps) != len(tsteps):
        bad.append(f"{len(osteps)} steps, ChemDraw {len(tsteps)}")
    for st in osteps:
        tst = tsteps.get(o2t.get(st.get("ReactionStepArrows")))
        if tst is None:
            bad.append(f"step {st.get('id')}: no ChemDraw step with its arrow")
            continue
        for attr in STEP_KINDS:
            got = sorted({o2t.get(i, f"?{i}") for i in (st.get(attr) or "").split()})  # split lines: one id
            want = sorted(i for i in their_ids(tst, attr) if i in written)
            if got != want:
                bad.append(f"step {tst.get('id')} {attr}: ChemDraw {want} ours {got}")
    # brackets, paired by x
    tb = [g for g in theirs.iter("graphic") if g.get("GraphicType") == "Bracket"]
    ob = [g for g in ours.iter("graphic") if g.get("GraphicType") == "Bracket"]
    if len(ob) != len(tb):
        bad.append(f"{len(ob)} brackets, ChemDraw {len(tb)}")
    for g in ob:
        m = [t for t in tb if abs(_floats(t.get("BoundingBox"))[0] - _floats(g.get("BoundingBox"))[0]) < 1.0]
        if len(m) != 1:
            bad.append(f"bracket at {g.get('BoundingBox')}: {len(m)} ChemDraw brackets")
            continue
        key = f"bracket {m[0].get('id')}"
        for attr in ("BracketType", "LipSize"):
            if g.get(attr) != m[0].get(attr):
                bad.append(f"{key} {attr}: ChemDraw {m[0].get(attr)} ours {g.get(attr)}")
        for k, (a, b) in enumerate(zip(_floats(m[0].get("BoundingBox")), _floats(g.get("BoundingBox")))):
            nums.append(("bracket box", (k, f"{key} box[{k}]"), a, b))
    # the steps sit in the same <scheme>s
    group = [sorted(o2t.get(st.get("ReactionStepArrows")) for st in sch.findall("step")) for sch in ours.iter("scheme")]
    tgroup = [sorted(st.get("ReactionStepArrows") for st in sch.findall("step")) for sch in theirs.iter("scheme")]
    if sorted(group) != sorted(tgroup):
        bad.append(f"schemes: ChemDraw {tgroup} ours {group}")
    return bad, nums


def _expand_aliases(m: Chem.Mol) -> Chem.Mol:
    """A molfile's abbreviation aliases (OMe on the O, as cdxml_to_molfile writes ChemDraw's
    nicknames) replaced by their groups (chem.ABBREVIATIONS), so it compares with what RDKit
    reads from a nickname node."""
    for a in m.GetAtoms():
        alias = a.GetProp("molFileAlias") if a.HasProp("molFileAlias") else None
        if alias not in ABBREVIATIONS:
            continue
        group = Chem.MolFromSmiles(ABBREVIATIONS[alias])
        star = next(x for x in group.GetAtoms() if x.GetAtomicNum() == 0)
        head = star.GetNeighbors()[0].GetIdx() + m.GetNumAtoms()
        (nb_bond,) = a.GetBonds()
        rw = Chem.RWMol(Chem.CombineMols(m, group))
        rw.AddBond(nb_bond.GetOtherAtomIdx(a.GetIdx()), head, nb_bond.GetBondType())
        for i in sorted((a.GetIdx(), star.GetIdx() + m.GetNumAtoms()), reverse=True):
            rw.RemoveAtom(i)
        out = rw.GetMol()
        Chem.SanitizeMol(out)
        return _expand_aliases(out)
    return m


def _spec_smiles(spec_path) -> dict[str, str]:
    """name -> canonical SMILES (no maps) of each species as the spec gives it, abbreviations
    expanded. A chair's stereo is in its geometry, which ChemDraw doesn't read either (chair_x*:
    no wedges in ChemDraw's own chairs), so a chair species is compared without stereo."""
    spec = load_spec(spec_path)
    out = {}
    for name, sp in spec.species.items():
        m = (Chem.MolFromSmiles(sp.smiles) if sp.smiles
             else _expand_aliases(Chem.MolFromMolFile(str(resolve_molfile(sp.molfile, spec_path.parent)))))
        for a in m.GetAtoms():
            a.SetAtomMapNum(0)
        if sp.chair is not None:
            Chem.RemoveStereochemistry(m)
        out[name] = Chem.MolToSmiles(m)
    return out


def check_scheme_round_trip(spec: str) -> None:
    """Every species reads back as its spec's (charges and radicals too); each step's reactants
    and products, by fragment id, are the spec step's (reversed for an arrow favoring reactants,
    as ChemDraw stores it), and the species drawn above its arrow are cited there."""
    path = mr.REPO / spec
    text = exported(spec)
    smiles = _spec_smiles(path)
    root = etree.fromstring(text.encode(), PARSER)
    frags = [f.get("id") for f in root.find("page").findall("fragment")]
    got = dict(zip(frags, (Chem.MolToSmiles(m) for m in Chem.MolsFromCDXML(text))))
    s = load_spec(path)
    drawn = [n for st in s.steps for n in (st.resonance if isinstance(st, ResonanceStep) else st.sources + st.drawn_products)]
    assert sorted(got.values()) == sorted(smiles[n] for n in dict.fromkeys(drawn))
    want_steps = []  # (reactants, products, above) per <step>: a resonance step has one per ↔
    for i, st in enumerate(s.steps, 1):
        if isinstance(st, ResonanceStep):
            want_steps += [([smiles[a]], [smiles[b]], []) for a, b in zip(st.resonance, st.resonance[1:])]
            continue
        over = s.over_names(i)
        want = ([smiles[n] for n in st.sources if n not in over], [smiles[n] for n in st.drawn_products])
        if st.arrow_type == "equilibrium_reactants":
            want = want[::-1]
        want_steps.append((*want, [smiles[n] for n in over]))
    steps = list(root.iter("step"))
    assert len(steps) == len(want_steps)
    for el, (reactants, products, over) in zip(steps, want_steps):
        have = tuple([got[i] for i in el.get(k).split()] for k in ("ReactionStepReactants", "ReactionStepProducts"))
        assert have == (reactants, products)
        above = [got[j] for j in (el.get("ReactionStepObjectsAboveArrow") or "").split() if j in got]
        assert above == over
    if len(steps) == 1:
        RDLogger.DisableLog("rdApp.*")  # RDKit looks for the <t> conditions among its fragments
        try:
            rxns = rdChemReactions.ReactionsFromCDXMLBlock(text)
        finally:
            RDLogger.EnableLog("rdApp.*")
        assert len(rxns) == 1
        # a reaction's SMILES keeps drawn H atoms ([H]C(C)(C)C): read back as plain molecules
        r, p = ([Chem.MolToSmiles(Chem.MolFromSmiles(x)) for x in side.split(".")]
                for side in rdChemReactions.ReactionToSmiles(rxns[0]).split(">>"))
        assert sorted(r) == sorted(want_steps[0][0]) and sorted(p) == sorted(want_steps[0][1])
