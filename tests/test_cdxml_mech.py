"""CDXML export, phase 3 (acs-style.md §13.2): M3 mechanisms (mech_basic_x* round trip: curves,
dot symbols; M3 specs).
"""

import math

import pytest
from lxml import etree
from rdkit import Chem

import measure_refs as mr
from acs_chem_svg.cdxml_writer import spec_cdxml, to_cdxml
from acs_chem_svg.geometry import sub
from acs_chem_svg.preset import load_preset

from cdxml_helpers import (
    PARSER, TOL, _canonical, _check_structure, _floats, _near, _nodes, check_scheme_round_trip, exported,
    needs_dtd, needs_refs, reference_pairs, reference_scheme_export, scheme_pairs,
)


MECH_REFS = [f"mech_basic_x{s}" for s in ("1.0", "0.8", "0.6")]
M3_SPECS = sorted(str(p.relative_to(mr.REPO)) for p in (mr.REPO / "examples" / "m3").glob("m3_*.yaml"))
CURVE_TOL = 0.01  # anchors and inner handles: ChemDraw's numbers copied
MIRROR_TOL = 0.025  # end handles: we mirror its rounded points (2p0 - c1), it mirrors unrounded ones


def _dot_symbols(root, nodes: dict) -> dict:
    """node id (mapped by `nodes`) -> [(SymbolType, BoundingBox)] of its LonePair / Electron symbols."""
    out = {}
    for g in root.iter("graphic"):
        if g.get("SymbolType") in ("LonePair", "Electron"):
            r = g.find("represent")
            assert r.get("attribute") == "Radical"
            out.setdefault(nodes.get(r.get("object"), r.get("object")), []).append((g.get("SymbolType"), _floats(g.get("BoundingBox"))))
    return out


def _ccw(atom, box) -> bool:
    u, v = sub(box[:2], atom), sub(box[2:], atom)
    return u[0] * v[1] - u[1] * v[0] < 0


def mech_pairs(ref: str) -> tuple[list[str], list[tuple[str, str, float, float]]]:
    """(mismatches, [(category, key, ChemDraw, ours)]) for the curves and dot symbols."""
    theirs = etree.parse(str(mr.REF_DIR / f"{ref}.cdxml"), PARSER).getroot()
    ours = etree.fromstring(reference_scheme_export(ref).encode(), PARSER)
    bad, nums = [], []
    # curves, paired by first anchor
    tcurves = list(theirs.iter("curve"))
    ocurves = list(ours.iter("curve"))
    if len(ocurves) != len(tcurves):
        bad.append(f"{len(ocurves)} curves, ChemDraw {len(tcurves)}")
    for oc in ocurves:
        op = _floats(oc.get("CurvePoints"))
        m = [c for c in tcurves if _near(_floats(c.get("CurvePoints"))[2:4], op[2:4])]
        if len(m) != 1:
            bad.append(f"curve at {op[2:4]}: {len(m)} ChemDraw curves")
            continue
        tc = m[0]
        key = f"curve {tc.get('id')}"
        for attr in ("CurveType", "ArrowheadHead", "ArrowheadType"):
            if oc.get(attr) != tc.get(attr):
                bad.append(f"{key} {attr}: ChemDraw {tc.get(attr)} ours {oc.get(attr)}")
        tp = _floats(tc.get("CurvePoints"))
        if tc.get("CurveType") == "32":  # its fishhooks repeat the last anchor; we don't (§13.2)
            assert len(tp) == 18 and tp[12:] == tp[6:12], key
            tp = tp[:12]
        if len(tp) != len(op):
            bad.append(f"{key}: {len(op) // 2} points, ChemDraw {len(tp) // 2}")
            continue
        for k, (a, b) in enumerate(zip(tp, op)):
            cat = "curve end handle" if k // 2 in (0, 5) else "curve points"
            nums.append((cat, f"{key} [{k}]", a, b))
    # dot symbols, per atom (nodes paired by p)
    tn, on = _nodes(theirs), _nodes(ours)
    to_ref = {oid: tid for oid, o in on.items() for tid, t in tn.items() if _near(_floats(t.get("p")), _floats(o.get("p")))}
    tsym, osym = _dot_symbols(theirs, {}), _dot_symbols(ours, to_ref)
    if sorted(tsym) != sorted(osym):
        bad.append(f"dots on ChemDraw nodes {sorted(tsym)}, ours {sorted(osym)}")
    for nid, syms in osym.items():
        theirs_here = tsym.get(nid, [])
        atom = _floats(tn[nid].get("p"))
        for kind in ("LonePair", "Electron"):
            o = [b for k, b in syms if k == kind]
            t = [b for k, b in theirs_here if k == kind]
            if len(o) != len(t):
                bad.append(f"node {nid}: {len(o)} {kind}, ChemDraw {len(t)}")
                continue
            for ob in o:  # paired with ChemDraw's nearest of the same kind
                tb = min(t, key=lambda b: abs(b[0] - ob[0]) + abs(b[1] - ob[1]))
                key = f"node {nid} {kind}"
                for k in range(2):
                    nums.append(("dot center", f"{key} {'xy'[k]}", tb[k], ob[k]))
                if kind == "LonePair":
                    nums.append(("lone pair separation", key, math.dist(tb[:2], tb[2:]), math.dist(ob[:2], ob[2:])))
                    if not (_ccw(atom, ob) and _ccw(atom, tb)):
                        bad.append(f"{key}: not counterclockwise (ChemDraw {_ccw(atom, tb)}, ours {_ccw(atom, ob)})")
                else:
                    nums.append(("electron box", f"{key} dx", tb[0] - tb[2], ob[0] - ob[2]))
                    nums.append(("electron box", f"{key} dy", tb[3] - tb[1], ob[3] - ob[1]))
    return bad, nums


# Hand-placed: dot centers (lone_pair.gap is a round default), the charge signs the dots push
# out (Drawer.charge_center), and the fragment boxes they extend. Largest seen: 1.14 (O- sign, x1.0).
PLACE_TOL = 1.2
PLACED = ("dot center", "charge center", "fragment box (with dots)", "fragment box (with a charge)")
# what the rules don't reproduce: (ref or "*", category, key prefix) -> its measured size
MECH_KNOWN = {
    # tBu radical on an unlabeled C: ChemDraw's dot is 2.8-5.0 from ours (our largest-gap rule)
    ("*", "dot center", "node 70 Electron"): 5.1,
    # ChemDraw pads the arrow box 0.18 above the axis, not 0.5 (as scheme_basic_x0.6)
    ("mech_basic_x0.8", "arrow box", "arrow"): 0.35, ("mech_basic_x0.6", "arrow box", "arrow"): 0.35,
    # Cl. typed as text (no face 96): ChemDraw stores its box 0.1 higher than a drawn label's
    ("mech_basic_x1.0", "label box", "node 79"): 0.11,
}


def _mech_tol(ref, cat, key, default):
    return next((t for (r, c, k), t in MECH_KNOWN.items() if r in (ref, "*") and c == cat and key.startswith(k)),
                default)


@needs_refs
@pytest.mark.parametrize("ref", MECH_REFS)
def test_mech_reference_diff(ref):
    """mech_basic_x* rebuilt (compare.reference_cdx_scheme: templates, ChemDraw's CurvePoints,
    lone pairs on ChemDraw's atoms, its steps) and compared with its file: curves to 0.01 (end
    handles to MIRROR_TOL), CurveType and ArrowheadHead; per atom, dot symbol kinds and counts,
    LonePair separation and Electron dx exact, both counterclockwise; dot centers, the signs
    they push out and the boxes they extend within PLACE_TOL (hand-placed); the scheme objects
    and nodes as for the scheme references (steps citing a <group> cite its fragment)."""
    bad, nums = mech_pairs(ref)
    sbad, snums = scheme_pairs(ref)
    nbad, nnums = reference_pairs(ref, reference_scheme_export(ref))
    bad += sbad + nbad
    tols = {"curve points": CURVE_TOL, "curve end handle": MIRROR_TOL, "lone pair separation": 0.01,
            "electron box": 0.01, **dict.fromkeys(PLACED, PLACE_TOL)}
    for cat, key, a, b in nums + [(c, k, a, b) for c, (_, k), a, b in snums] + nnums:
        tol = _mech_tol(ref, cat, key, tols.get(cat, TOL))
        if abs(b - a) > tol:
            bad.append(f"{cat} {key}: ChemDraw {a} ours {b} (diff {b - a:+.3f})")
    assert not bad, "\n".join(bad)
    cats = {c for c, *_ in nums}
    assert {"curve points", "dot center", "lone pair separation", "electron box"} <= cats


@needs_refs
def test_curve_end_is_stroke_end():
    """In all 12 reference curves the last anchor is the SVG stroke end (the head's notch, our
    MechArrow.p3), not the head tip, so the writer writes p3 as it is."""
    n = 0
    for ref in MECH_REFS:
        for c in mr.measure(mr.load(ref))["curved_arrows"]:
            assert max(map(abs, c["last_anchor_minus_stroke_end"])) < 0.01, ref
            assert max(map(abs, c["head"]["tip_minus_stroke_end"])) > 0.5  # the tip is elsewhere
            n += 1
    assert n == 12


@needs_refs
def test_fishhook_side_matches_reference():
    """mechanism.fishhook_side gives +1 (cdxml_writer: HalfRight) on every reference fishhook,
    and the writer's ArrowheadHead equals ChemDraw's on every curve."""
    from compare import curve_points
    from acs_chem_svg.mechanism import fishhook_side
    for ref in MECH_REFS:
        hooks = [c for c in mr.load(ref).curves if c["type"] == "32"]
        assert len(hooks) == 2 and all(c["head"] == "HalfRight" for c in hooks)
        p = [curve_points(c) for c in hooks]
        assert all(fishhook_side(c1, c2, p3) == 1.0 for _, c1, c2, p3 in p)


def test_fishhook_left_is_half_left():
    """The references' fishhooks turn left on screen (counterclockwise, y down) and are
    HalfRight; one turning right is written HalfLeft, by symmetry (no reference has one; checked
    at the phase 3 gate)."""
    from acs_chem_svg.cdxml_writer import CdxCurve, CdxScheme
    turns_left = CdxCurve((0.0, 0.0), (10.0, 10.0), (20.0, 10.0), (30.0, 0.0), "fishhook")  # down, then up
    turns_right = CdxCurve((0.0, 0.0), (10.0, -10.0), (20.0, -10.0), (30.0, 0.0), "fishhook")
    heads = []
    for c in (turns_left, turns_right):
        root = etree.fromstring(to_cdxml(CdxScheme([], curves=[c]), "acs1996", warn=lambda m: None).encode(), PARSER)
        heads.append(root.find(".//curve").get("ArrowheadHead"))
    assert heads == ["HalfRight", "HalfLeft"]


@needs_refs
def test_dot_box_keys_match_references():
    """cdxml.electron_box_dx, lone_pair_box_reach and electron_box_reach against
    measure_refs.py --cdxml-symbols: every Electron box's dx; per fragment side bounded by a dot,
    the reach (more where other ink bounds it), each reach seen at least once."""
    for ref, preset in zip(MECH_REFS, ("acs1996", "acs1996_x0.8", "acs1996_x0.6")):
        cx = load_preset(preset).cdxml
        m = mr.cdxml_symbols(ref)
        assert all(e["dx"] == pytest.approx(cx.electron_box_dx, abs=0.01) and e["dy"] == 0 for e in m["electrons"])
        assert all(p["separation"] == pytest.approx(load_preset(preset).lone_pair.dot_spacing, abs=0.01)
                   for p in m["lone_pairs"])
        reach = {"pair": cx.lone_pair_box_reach, "electron": cx.electron_box_reach}
        seen = set()
        for r in m["fragment_reach_beyond_dot_centers"]:
            for side in ("left", "top", "right", "bottom"):
                v, kind = r[side]
                assert v >= reach[kind] - 0.01, (ref, r)
                if abs(v - reach[kind]) <= 0.01:
                    seen.add(kind)
        assert seen == {"pair", "electron"}, ref


@needs_dtd
@pytest.mark.parametrize("ref", MECH_REFS)
def test_mech_reference_exports_validate(ref):
    _check_structure(reference_scheme_export(ref))


@needs_dtd
@pytest.mark.parametrize("spec", M3_SPECS)
def test_mech_spec_exports_validate(spec):
    _check_structure(exported(spec))


@pytest.mark.parametrize("spec", M3_SPECS)
def test_mech_rdkit_round_trip(spec):
    check_scheme_round_trip(spec)


def test_mech_objects():
    """m3_sn2_x1.0: one <curve> per curved arrow, free on the page after the graphics; LonePair
    symbols on the nucleophile and leaving group, bound by represent attribute="Radical"."""
    root = etree.fromstring(exported("examples/m3/m3_sn2_x1.0.yaml").encode(), PARSER)
    page = root.find("page")
    curves = page.findall("curve")
    assert len(curves) == 2 and all(c.get("CurveType") == "8" and c.get("ArrowheadHead") == "Full" for c in curves)
    assert all(len(_floats(c.get("CurvePoints"))) == 12 for c in curves)
    assert {g.get("SymbolType") for g in root.iter("graphic") if g.get("GraphicType") == "Symbol"} >= {"LonePair"}
    root = etree.fromstring(exported("examples/m3/m3_radical_initiation_x1.0.yaml").encode(), PARSER)
    hooks = root.findall(".//curve")
    # six points, as a 2e- arrow: no repeated last anchor (a stray node in ChemDraw; gate 2026-10-28)
    assert hooks and all(c.get("CurveType") == "32" and len(_floats(c.get("CurvePoints"))) == 12 for c in hooks)
    assert all(n.get("Radical") == "Doublet" for g in root.iter("graphic") if g.get("SymbolType") == "Electron"
               for n in root.iter("n") if n.get("id") == g.find("represent").get("object"))


def test_mechanism_maps_do_not_warn():
    """A mechanism's species are all mapped (its arrows need the maps), so the map warning is
    silenced there (user decision 2026-10-28); M1 species still warn."""
    msgs = []
    spec_cdxml(mr.REPO / "examples/m3/m3_carbonyl_addition_x1.0.yaml", warn=msgs.append)
    assert msgs == []
    spec_cdxml(mr.REPO / "examples/m1/m1_charges_x1.0.yaml", warn=msgs.append)
    assert any("atom maps not written" in m for m in msgs)


def test_added_reactant_is_above_arrow():
    """m3_carbonyl_addition_x1.0 step 2: water, added by the mechanism step and drawn above its
    arrow, is cited in ReactionStepObjectsAboveArrow, not as a reactant (user decision 2026-10-28)."""
    text = exported("examples/m3/m3_carbonyl_addition_x1.0.yaml")
    root = etree.fromstring(text.encode(), PARSER)
    frags = [f.get("id") for f in root.find("page").findall("fragment")]
    smiles = dict(zip(frags, (Chem.MolToSmiles(m) for m in Chem.MolsFromCDXML(text))))
    st = list(root.iter("step"))[1]
    assert [smiles[i] for i in st.get("ReactionStepObjectsAboveArrow").split()] == ["O"]
    assert [smiles[i] for i in st.get("ReactionStepReactants").split()] == [_canonical("N#CC(C)([O-])C")]
