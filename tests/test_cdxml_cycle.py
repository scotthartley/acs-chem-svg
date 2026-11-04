"""CDXML export of catalytic cycles (acs-style.md §15.1, phase 5).

- Calibration: every arcs_x* and arcs_eq_x* arc arrow (forward, balanced and unbalanced
  equilibrium, headless), rewritten on ChemDraw's own circle, gives back its <arrow> (Head3D,
  Tail3D, Center3D, axis ends, AngularSize, head and shaft attributes) and its legacy <graphic>
  (BoundingBox head then center, ArrowType, the arc tool's AngularSize) within 0.1 pt.
- The cycle examples: each arc <arrow> is a consistent ChemDraw arc (measure_refs.arc_circle,
  from Head3D, Tail3D and AngularSize alone, gives its Center3D and radius); each, redrawn from
  its attributes alone, is a drawn arrow or tail of the layout; one <scheme> with a <step> per
  arrow, the captions, and an RDKit round trip.
"""

import math

import pytest
from compare import preset_for
from lxml import etree
from rdkit import Chem

import measure_refs as mr
from acs_chem_svg import chem, richtext
from acs_chem_svg.cdxml_writer import CdxArc, CdxScheme, load_nicknames, spec_cdxml, to_cdxml
from acs_chem_svg.cycle import build_cycle, tail_items
from acs_chem_svg.arrows import arc_step_arrow
from acs_chem_svg.geometry import arc_segments, length, sub
from acs_chem_svg.scene import Group, Path
from acs_chem_svg.spec import load_spec
from acs_chem_svg.svg_writer import to_svg

PARSER = etree.XMLParser(load_dtd=False, no_network=True, remove_blank_text=True)
TOL = 0.1
# arcs_eq_x*'s headless arcs are stored slightly elliptical (axes 0.10-0.23 pt apart), so a circle
# meets their points within this, and their AngularSize within 0.25°
HEADLESS_TOL = 0.25
# the circle from Head3D / Tail3D (0.01 pt) and AngularSize (0.1°, as ChemDraw writes it) vs
# Center3D: the coordinates' share, plus what half a 0.1° step in the sweep moves the center
# (_sweep_shift; 0.12 pt on a 101 pt arc of 41.6°); ChemDraw's own arcs, up to 0.07 (§15.1)
CENTER_TOL = 0.02
ARC_REFS = ["arcs_x1.0", "arcs_x0.8", "arcs_x0.6", "arcs_eq_x1.0", "arcs_eq_x0.8", "arcs_eq_x0.6"]
VARIANTS = {"1000": "standard", "1500": "large"}
CYCLES = sorted(str(p.relative_to(mr.REPO)) for p in (mr.REPO / "examples" / "m6").glob("m6_*_cycle_*.yaml"))
needs_refs = pytest.mark.skipif(not mr.REF_DIR.exists(), reason="ChemDraw references (reference_renders/) are local only")


def _xy(el, key):
    return tuple(float(v) for v in el.get(key).split()[:2])


def _floats(s):
    return [float(v) for v in s.split()]


def _angle(c, p):
    return math.atan2(-(p[1] - c[1]), p[0] - c[0])


def _export(spec: str):
    return etree.fromstring(spec_cdxml(mr.REPO / spec, warn=lambda m: None).encode(), PARSER)


def _kind(a) -> str:
    """The arrow_type (or "none") of a CDXML arc <arrow>, by its attributes."""
    if a.get("ArrowheadHead") is None:
        return "none"
    if a.get("ArrowheadHead") == "Full":
        return "forward"
    return "equilibrium_products" if a.get("ArrowEquilibriumRatio") else "equilibrium"


def _elliptical(a) -> bool:
    c = _xy(a, "Center3D")
    return abs((_xy(a, "MajorAxisEnd3D")[0] - c[0]) - (_xy(a, "MinorAxisEnd3D")[1] - c[1])) > 1.0


def _arcs(root):
    """[(<graphic>, <arrow>)] of the arc arrows."""
    page = root.find("page")
    arrows = {a.get("id"): a for a in page.findall("arrow")}
    return [(g, arrows[g.get("SupersededBy")]) for g in page.findall("graphic") if g.get("GraphicType") == "Arc"]


# --------------------------------------------------------------------------- calibration


@needs_refs
@pytest.mark.parametrize("ref", ARC_REFS)
def test_reference_arcs_rewritten(ref):
    root = etree.parse(str(mr.REF_DIR / f"{ref}.cdxml"), PARSER).getroot()
    # ChemDraw glitches (the user, 2026-10-30): arcs_eq_x1.0's large balanced arrows are stored as
    # ellipses (minor axis 60.16 vs major 80.84 at 60°) and drawn short; not circular arcs, left out
    theirs = [(g, a) for g, a in _arcs(root) if not _elliptical(a)]
    assert len(_arcs(root)) - len(theirs) == (6 if ref == "arcs_eq_x1.0" else 0)
    records = []
    for _, a in theirs:
        c = _xy(a, "Center3D")
        r = _xy(a, "MajorAxisEnd3D")[0] - c[0]
        a0, a1 = _angle(c, _xy(a, "Tail3D")), _angle(c, _xy(a, "Head3D"))
        sweep = float(a.get("AngularSize"))
        while math.copysign(1, sweep) * (a1 - a0) <= 0:  # unwrap the way the arc runs
            a1 += math.copysign(2 * math.pi, sweep)
        records.append(CdxArc(c, r, a0, a1, _kind(a), VARIANTS[a.get("HeadSize")]))
    ours = _arcs(etree.fromstring(to_cdxml(CdxScheme([], arcs=records), preset_for(ref)).encode(), PARSER))
    bad = []
    for k, ((g0, a0), (g1, a1)) in enumerate(zip(theirs, ours)):
        tol = HEADLESS_TOL if _kind(a0) == "none" else TOL
        for key in ("FillType", "ArrowheadHead", "ArrowheadTail", "ArrowheadType", "HeadSize", "ArrowheadCenterSize",
                    "ArrowheadWidth", "ArrowShaftSpacing", "ArrowEquilibriumRatio"):
            if a0.get(key) != a1.get(key):
                bad.append(f"arc {k} {key}: ChemDraw {a0.get(key)} ours {a1.get(key)}")
        for key in ("GraphicType", "ArrowType", "HeadSize", "AngularSize"):
            if g0.get(key) != g1.get(key):
                bad.append(f"arc {k} graphic {key}: ChemDraw {g0.get(key)} ours {g1.get(key)}")
        if abs(float(a0.get("AngularSize")) - float(a1.get("AngularSize"))) > (0.15 if tol == TOL else 0.25):
            bad.append(f"arc {k} AngularSize: ChemDraw {a0.get('AngularSize')} ours {a1.get('AngularSize')}")
        for key in ("Head3D", "Tail3D", "Center3D", "MajorAxisEnd3D", "MinorAxisEnd3D"):
            if length(sub(_xy(a0, key), _xy(a1, key))) > tol:
                bad.append(f"arc {k} {key}: ChemDraw {a0.get(key)} ours {a1.get(key)}")
        if max(abs(u - v) for u, v in zip(_floats(g0.get("BoundingBox")), _floats(g1.get("BoundingBox")))) > tol:
            bad.append(f"arc {k} graphic BoundingBox: ChemDraw {g0.get('BoundingBox')} ours {g1.get('BoundingBox')}")
    assert len(ours) == len(theirs) == {"arcs_x1.0": 12, "arcs_eq_x1.0": 24, "arcs_eq_x0.8": 18, "arcs_eq_x0.6": 18}.get(ref, 6)
    assert not bad, "\n".join(bad)


# --------------------------------------------------------------------------- the examples


@pytest.mark.parametrize("spec", CYCLES)
def test_cycle_arcs_are_chemdraw_arcs(spec):
    """ChemDraw rebuilds an arc from Head3D, Tail3D and AngularSize (§15.1): that circle is our
    Center3D and radius within CENTER_TOL, and the axis ends are center + (r, 0) and + (0, r)."""
    for _, a in _arcs(_export(spec)):
        c, R, _ = mr.arc_circle(a)
        center = _xy(a, "Center3D")
        tol = CENTER_TOL + _sweep_shift(R, float(a.get("AngularSize")))
        assert length(sub(c, center)) < tol
        assert abs(_xy(a, "MajorAxisEnd3D")[0] - center[0] - R) < tol
        assert length(sub(_xy(a, "MinorAxisEnd3D"), (center[0], center[1] + R))) < tol


def _sweep_shift(r: float, sweep_deg: float) -> float:
    """How far the center of a circle through two fixed points moves when its sweep is off by
    0.05° (AngularSize's rounding): the center's distance from the chord is (c/2)/tan(θ/2), so
    its derivative is (c/4)/sin²(θ/2) with c = 2r sin(θ/2)."""
    half = math.radians(abs(sweep_deg)) / 2
    return r / (2 * math.sin(half)) * math.radians(0.05)


@pytest.mark.parametrize("spec", CYCLES)
def test_cycle_arcs_are_the_drawing(spec):
    """Each arc <arrow>, redrawn from its attributes alone (circle: Center3D and the axis end;
    angles: Tail3D, Head3D; kind: its head attributes, an equilibrium favoring reactants read
    backwards), covers the ink of one of the layout's arrows or tails within 0.02 pt; forward
    arrows and leaving tails are Full, equilibria HalfLeft at both ends, entering tails headless."""
    path = mr.REPO / spec
    s = load_spec(path)
    cyc = build_cycle(s, path.parent, path.stem)
    _, off = to_svg(cyc.root)
    drawn = [tuple(v + o for v, o in zip(g.ink_bbox(), off + off)) for g in
             [a.group for a in cyc.arcs] + [Group("tail", [t.group.children[0]]) for t in cyc.tails]]
    arcs = _arcs(_export(spec))
    assert len(arcs) == len(drawn)
    for g, a in arcs:
        c = _xy(a, "Center3D")
        r = _xy(a, "MajorAxisEnd3D")[0] - c[0]
        a0, a1 = _angle(c, _xy(a, "Tail3D")), _angle(c, _xy(a, "Head3D"))
        sweep = float(a.get("AngularSize"))
        while math.copysign(1, sweep) * (a1 - a0) <= 0:
            a1 += math.copysign(2 * math.pi, sweep)
        kind = _kind(a)
        if kind == "none":
            box = Path(arc_segments(c, r, a0, a1), fill=False, stroke=cyc.preset.arrow.shaft_width, closed=False).ink_bbox()
        else:
            box = arc_step_arrow(kind, c, r, a0, a1, cyc.preset, "x", s.arrow_head).ink_bbox()
        assert min(max(abs(u - v) for u, v in zip(box, d)) for d in drawn) < 0.02, (kind, a.get("Tail3D"))  # 2-decimal CDXML
        assert (g.get("ArrowType") is None) == (kind == "none")
    heads = sorted(_kind(a) for _, a in arcs)
    steps = s.cycle.ordered_steps()
    n_in = sum(len(tail_items(st.in_)) for st in steps)
    n_out = sum(len(tail_items(st.out)) for st in steps)
    want = [st.arrow_type if st.arrow_type != "equilibrium_reactants" else "equilibrium_products" for st in steps]
    assert heads == sorted(want + ["forward"] * n_out + ["none"] * n_in)


@pytest.mark.parametrize("spec", CYCLES)
def test_cycle_steps_and_captions(spec):
    """One <scheme>, a <step> per arrow citing its arc graphic, from
    the step's species and those entering to its next species and those leaving; captions for
    the labels, tail texts, center text and numbers."""
    s = load_spec(mr.REPO / spec)
    root = _export(spec)
    page = root.find("page")
    frags = {f.get("id"): k for k, f in enumerate(page.findall("fragment"))}
    graphics = {g.get("id"): g for g in page.findall("graphic")}
    steps = s.cycle.ordered_steps()
    (scheme,) = root.iter("scheme")
    cdx_steps = scheme.findall("step")
    assert len(cdx_steps) == len(steps)
    names = list(build_cycle(s, (mr.REPO / spec).parent).species)
    for st, el in zip(steps, cdx_steps):
        arcs = el.get("ReactionStepArrows").split()
        assert len(arcs) == 1
        assert all(graphics[i].get("GraphicType") == "Arc" for i in arcs)
        ins = [t.species for t in st.in_ if t.species]
        outs = [t.species for t in st.out if t.species]
        assert [names[frags[i]] for i in el.get("ReactionStepReactants").split()] == [st.from_] + ins
        assert [names[frags[i]] for i in el.get("ReactionStepProducts").split()] == [st.to] + outs
    texts = sorted("".join(x.text for x in t.iter("s")) for t in page.findall("t"))
    want = [st.label for st in steps if st.label]
    want += [txt for st in steps for side in (st.in_, st.out) for txt, _ in tail_items(side) if txt]
    want += [s.cycle.center] if s.cycle.center else []
    want += [str(v) for v in (s.numbers or {}).values()]
    assert texts == sorted("".join(r.text for r in richtext.parse(w)) for w in want)  # Br– is set as Br−


def _comparable(m: Chem.Mol, plain_bonds: bool) -> str:
    """SMILES to compare: double-bond geometry dropped (the drawing fixes an unspecified one), and
    for a species with dative bonds its skeleton (CDXML has single bonds there, and RDKit reads
    the four-bonded P as PH)."""
    m = Chem.RWMol(m)
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE:
            b.SetStereo(Chem.BondStereo.STEREONONE)
            b.SetBondDir(Chem.BondDir.NONE)
        if plain_bonds and b.GetBondType() == Chem.BondType.DATIVE:
            b.SetBondType(Chem.BondType.SINGLE)
    for a in m.GetAtoms():  # @SP is in the drawing's geometry; RDKit reads no square-planar stereo
        if a.GetChiralTag() == Chem.ChiralType.CHI_SQUAREPLANAR:
            a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    if plain_bonds:
        for a in m.GetAtoms():
            a.SetNumExplicitHs(0)
            a.SetNoImplicit(True)
        m.UpdatePropertyCache(strict=False)
    return Chem.MolToSmiles(m)


def _as_rdkit_reads(m: Chem.Mol, abbreviate, preset: str) -> Chem.Mol:
    """`m` with each abbreviation written as an Unspecified node (PMe3, PPh3) a bare carbon:
    RDKit reads those as C, ChemDraw's own (metals_x* PPh3) too; Fragment and Nickname nodes it
    expands, so they stay in full."""
    kinds = {k: t.get("NodeType") for k, t in load_nicknames(preset).items()}
    keep = [k for k in abbreviate or [] if kinds.get(k) == "Unspecified"]
    if not keep:
        return m
    rw = Chem.RWMol(chem.condense(m, keep))
    for a in rw.GetAtoms():
        if a.GetAtomicNum() == 0 and a.HasProp("atomLabel"):
            a.SetAtomicNum(6)
            a.ClearProp("atomLabel")
    out = rw.GetMol()
    Chem.SetAromaticity(out)  # condense keeps the Kekulé bonds
    return out


@pytest.mark.parametrize("spec", CYCLES)
def test_cycle_round_trip(spec):
    """Every species reads back with RDKit as its spec's (_comparable; skeletons for every
    species of a spec with dative bonds)."""
    s = load_spec(mr.REPO / spec)
    mols = [_as_rdkit_reads(Chem.MolFromSmiles(sp.smiles), sp.abbreviate, s.preset) for sp in s.species.values()]
    plain = any(b.GetBondType() == Chem.BondType.DATIVE for m in mols for b in m.GetBonds())
    got = Chem.MolsFromCDXML(etree.tostring(_export(spec)).decode())
    assert sorted(_comparable(m, plain) for m in got) == sorted(_comparable(m, plain) for m in mols)
