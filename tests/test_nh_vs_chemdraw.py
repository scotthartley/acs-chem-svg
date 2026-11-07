"""Where an NH's H goes vs ChemDraw (nh_x*: pyrrolidine turned in 15° steps, an acyclic R-NH-R
at a peak, a valley and turned 15° / 30°, aziridine up and down; style-spec.md §2). Each NH is
redrawn at ChemDraw's atom positions (its two neighbors as plain carbons): the H side (beside,
above, below), the N and H glyph origins and each bond end within 0.1 pt, and the CDXML <t>
(LabelJustification, LabelAlignment, LineStarts, p) as ChemDraw stores it."""

import math

import pytest
from lxml import etree

import measure_refs as mr
from chem_svg_cdxml import draw_mol
from chem_svg_cdxml.cdxml_writer import to_cdxml
from chem_svg_cdxml.chem import Atom, Bond, Molecule
from chem_svg_cdxml.draw_mol import Drawer
from chem_svg_cdxml.geometry import sub
from chem_svg_cdxml.preset import load_preset

pytestmark = pytest.mark.skipif(not (mr.REF_DIR / "nh_x1.0.cdxml").exists(),
                                reason="ChemDraw references (reference_renders/) are local only")

TOL = 0.1
REFS = {"nh_x1.0": "acs1996", "nh_x0.8": "acs1996_x0.8", "nh_x0.6": "acs1996_x0.6"}
# bonds 6° and 9° off horizontal or vertical: ChemDraw stops them 0.17-0.61 pt nearer N than our
# §2.1 clip; a ChemDraw error, not emulated (user decision 2026-10-30, style-spec.md Open 28)
NEAR_AXIS_KNOWN = 0.65
# x0.6 only (Open 28, not emulated either): two bonds 21-24° off horizontal 0.10 / 0.12 pt short, and the acyclic NH
# turned 30° (atom 312), whose bond at -60° ChemDraw starts 0.46 pt farther out (x1.0 / x0.8: 0.04)
KNOWN = {("nh_x0.6", "241"): 0.15, ("nh_x0.6", "43"): 0.15, ("nh_x0.6", "312"): 0.5}
SIDE = {"Left": "right", "Right": "left", "Above": "top", "Below": "bottom"}  # LabelAlignment -> h_side


def _cases(ref):
    doc = mr.load(ref)
    root = etree.parse(str(mr.REF_DIR / f"{ref}.cdxml")).getroot()
    nodes = {n.get("id"): n for n in root.iter("n")}
    for a in doc.atoms.values():
        if a.label != "NH":
            continue
        nbs = [doc.atoms[b.b if b.a == a.id else b.a] for b in doc.bonds if a.id in (b.a, b.b)]
        m = Molecule("m", [Atom(0, "1", "N", *a.p, num_h=1)] + [Atom(k + 1, str(k + 2), "C", *o.p) for k, o in enumerate(nbs)],
                     [Bond(k, 0, k + 1) for k in range(len(nbs))])
        yield doc, a, nbs, nodes[a.id], m


@pytest.mark.parametrize("ref", REFS)
def test_h_side_glyphs_and_bond_ends(ref):
    preset = load_preset(REFS[ref])
    bad, sides = [], set()
    for doc, a, nbs, node, m in _cases(ref):
        d = Drawer(m, preset)
        lab = d.labels[0]
        want = SIDE[node.find("t").get("LabelAlignment")]
        sides.add(want)
        if lab.h_side != want:
            bad.append(f"{a.id}: H {lab.h_side}, ChemDraw {want}")
            continue
        for r in lab.runs:
            g = min((g for g in doc.glyphs if g.ch == r.text[0] and g.box), key=lambda g: math.dist((g.x, g.base), (r.x, r.y)))
            if max(abs(r.x - g.x), abs(r.y - g.base)) > TOL:
                bad.append(f"{a.id} {want} glyph {r.text}: ours ({r.x:.3f}, {r.y:.3f}) ChemDraw ({g.x:.3f}, {g.base:.3f})")
        for b, o in zip(m.bonds, nbs):
            u = mr.unit(sub(o.p, a.p))
            n = math.dist(a.p, o.p)
            ends = [min((q[0] - a.p[0]) * u[0] + (q[1] - a.p[1]) * u[1] for q in s.pts) for s in doc.shapes
                    if s.kind == "poly" and not s.curved
                    and max(abs((q[0] - a.p[0]) * u[1] - (q[1] - a.p[1]) * u[0]) for q in s.pts) < 1.0
                    and all(-0.5 < (q[0] - a.p[0]) * u[0] + (q[1] - a.p[1]) * u[1] < n + 1.0 for q in s.pts)]
            off_axis = abs((math.degrees(math.atan2(u[1], u[0])) + 45) % 90 - 45)
            known = max(NEAR_AXIS_KNOWN if 3 < off_axis < 12 else TOL, KNOWN.get((ref, a.id), TOL))
            if abs(d.end_clip(b, 0) - min(ends)) > known:
                bad.append(f"{a.id} {want} bond: end ours {d.end_clip(b, 0):.3f} ChemDraw {min(ends):.3f}")
    assert sides == {"left", "right", "top", "bottom"} and not bad, "\n".join(bad)


@pytest.mark.parametrize("ref", REFS)
def test_cdxml_label(ref):
    preset = load_preset(REFS[ref])
    bad = []
    for _, a, _, node, m in _cases(ref):
        d = Drawer(m, preset)
        root = etree.fromstring(to_cdxml([(m, d, (0.0, 0.0))], REFS[ref], warn=lambda msg: None).encode())
        n = next(n for n in root.iter("n") if n.get("Element") == "7")
        ours, theirs = n.find("t"), node.find("t")
        for k in ("LabelJustification", "Justification", "LabelAlignment", "LineStarts"):
            if ours.get(k) != theirs.get(k):
                bad.append(f"{a.id}: {k} ours {ours.get(k)} ChemDraw {theirs.get(k)}")
        dp = [float(v) - float(w) for v, w in zip(ours.get("p").split(), n.get("p").split())]
        dt = [float(v) - float(w) for v, w in zip(theirs.get("p").split(), node.get("p").split())]
        if max(abs(x - y) for x, y in zip(dp, dt)) > 0.02:
            bad.append(f"{a.id} {theirs.get('LabelAlignment')}: <t> p - atom ours {dp} ChemDraw {dt}")
    assert not bad, "\n".join(bad)


@pytest.mark.parametrize("dirs, side", [
    ([-60, -120], "top"), ([60, 120], "bottom"), ([-54, -126], "top"),  # bonds both down / up
    ([-18, -126], "top"), ([-12, -120], "left"),  # 18° off vertical: above; 24°: beside (H_VERTICAL_DEG)
    ([0, 180], "right"), ([-30, 90], "left")])
def test_vertical_h_rule(dirs, side):
    """draw_mol.vertical_h: the open side within H_VERTICAL_DEG of vertical (ChemDraw: 15° above,
    30° beside); bonds opposite have no open side."""
    atoms = [Atom(0, "1", "N", 0.0, 0.0, num_h=1)] + [
        Atom(k + 1, str(k + 2), "C", 14.4 * math.cos(math.radians(t)), -14.4 * math.sin(math.radians(t)))
        for k, t in enumerate(dirs)]
    m = Molecule("m", atoms, [Bond(k, 0, k + 1) for k in range(len(dirs))])
    assert draw_mol.h_side(m, 0) == side


@pytest.mark.parametrize("element, num_h, charge", [("O", 1, 1), ("N", 2, 1), ("N", 1, 1), ("P", 1, 0)])
def test_vertical_h_only_where_measured(element, num_h, charge):
    """Only a neutral NH: an oxonium, NH2+, an iminium NH+ or a PH keeps H beside (unmeasured)."""
    atoms = [Atom(0, "1", element, 0.0, 0.0, charge=charge, num_h=num_h),
             Atom(1, "2", "C", -12.47, 7.2), Atom(2, "3", "C", 12.47, 7.2)]
    m = Molecule("m", atoms, [Bond(0, 0, 1), Bond(1, 0, 2)])
    assert draw_mol.h_side(m, 0) == "right"
