"""Arc arrows vs ChemDraw (acs-style.md §15.1), within 0.1 pt, at all three presets.

Each arcs_x* arrow is redrawn on ChemDraw's circle (measure_refs.arc_circle: from Head3D,
Tail3D and the signed AngularSize, which Center3D matches), with the head
variant from its HeadSize. Counterclockwise arcs: the shaft's start and end (the notch), the tip
and every vertex and control point of the head are compared. ChemDraw's own shaft cubics stray
up to 0.33 pt from its circle (§15.1), so the shaft between its ends is held to that instead.
Clockwise arcs: ChemDraw's heads are the measured shape scaled by up to CW_HEAD_SCALE (no
pattern found), so their tail, tip and head orientation are compared, and the head's size to that.
"""

import math
from functools import lru_cache

import pytest
from lxml import etree

import measure_refs as mr
from acs_chem_svg.arrows import arc_arrow
from acs_chem_svg.geometry import length, sub
from acs_chem_svg.preset import load_preset

pytestmark = pytest.mark.skipif(not mr.REF_DIR.exists(),
                                reason="ChemDraw references (reference_renders/) are local only")

TOL = 0.1
SHAFT_KNOWN = 0.34  # ChemDraw's cubic approximation of its circle (measure_refs --arc: 0.333)
CW_HEAD_SCALE = 0.035  # clockwise heads: length / measured length 0.972 to 1.031 (measure_refs --arc)
REFS = {"arcs_x1.0": "acs1996", "arcs_x0.8": "acs1996_x0.8", "arcs_x0.6": "acs1996_x0.6"}
VARIANTS = {"1000": "standard", "1500": "large"}


@lru_cache(maxsize=None)
def pairs(ref):
    """(ChemDraw shaft points from the tail, ChemDraw head, our shaft Path, our head Path, center,
    R, sweep in degrees) per arc arrow."""
    doc = mr.load(ref)
    preset = load_preset(REFS[ref])
    out = []
    for a in etree.parse(str(mr.REF_DIR / f"{ref}.cdxml")).getroot().iter("arrow"):
        c, R, sweep = mr.arc_circle(a, doc.offset)
        h = tuple(float(x) + o for x, o in zip(a.get("Head3D").split()[:2], doc.offset))
        t = tuple(float(x) + o for x, o in zip(a.get("Tail3D").split()[:2], doc.offset))
        shaft, spts = mr.arc_shaft(doc, t)
        head = next(s for s in doc.shapes if s.kind == "poly" and s.curved and any(length(sub(q, h)) < 0.05 for q in s.pts))
        a0 = math.atan2(-(t[1] - c[1]), t[0] - c[0])
        ours = arc_arrow(c, R, a0, a0 + math.radians(sweep), preset, "arc", VARIANTS[a.get("HeadSize")]).children[0]
        out.append((spts, shaft.width, head, ours, ours.head_leaves()[0], c, R, sweep))
    return out


@pytest.mark.parametrize("ref", REFS)
def test_arc_shaft_ends_and_circle(ref):
    for spts, width, _, ours, _, c, R, sweep in pairs(ref):
        assert length(sub(ours.points[0], spts[0])) < TOL
        if sweep > 0:
            assert length(sub(ours.points[-1], spts[-1])) < TOL  # the notch
        assert abs(length(sub(spts[-1], c)) - R) < TOL  # ChemDraw's notch is on the circle either way
        assert max(abs(length(sub(q, c)) - R) for q in mr._cubic_samples(spts)) < SHAFT_KNOWN
        assert ours.stroke == pytest.approx(width, abs=0.01)
    sweeps = sorted(p[-1] for p in pairs(ref))  # 60°, 120° and 180° each way, drawn by hand (59.6–60.1)
    assert sweeps == pytest.approx(sorted([s * k for s in (60.0, 120.0, 180.0) for k in (1, -1)]
                                          * (2 if ref == "arcs_x1.0" else 1)), abs=0.5)


@pytest.mark.parametrize("ref", REFS)
def test_arc_head_on_the_chord(ref):
    measured = mr.arc_summary(ref)["arcs"]  # same order as pairs()
    for (spts, _, head, shaft, ours, c, R, sweep), m in zip(pairs(ref), measured):
        assert abs(m["head"]["axis_minus_chord_deg"]) < 0.05  # ChemDraw's head on its own chord, both ways
        mine = ours.points
        if sweep > 0:
            for q in head.pts:  # every ChemDraw vertex and control has one of ours on it, and back
                assert min(length(sub(q, m)) for m in mine) < TOL, q
            for m in mine:
                assert min(length(sub(m, q)) for q in head.pts) < TOL, m
            continue
        tip = max(head.pts, key=lambda q: length(sub(q, spts[-1])))  # ChemDraw's tip: farthest from its notch
        my_tip = max(mine, key=lambda q: length(sub(q, shaft.points[-1])))
        assert length(sub(tip, my_tip)) < TOL
        own = load_preset(REFS[ref]).head(VARIANTS[m["head_attrs"].split(" / ")[0]])
        assert m["head"]["length"] / own.head_length == pytest.approx(1.0, abs=CW_HEAD_SCALE)
        assert m["head"]["inset (chord notch -> tip)"] / m["head"]["length"] == pytest.approx(
            own.head_inset / own.head_length, abs=0.005)  # the measured shape, scaled


# --------------------------------------------------------------------------- equilibrium arcs (arcs_eq_x*)

EQ_REFS = {"arcs_eq_x1.0": "acs1996", "arcs_eq_x0.8": "acs1996_x0.8", "arcs_eq_x0.6": "acs1996_x0.6"}
# approximate (measure_refs --arc-eq): ChemDraw's arc equilibria differ between clockwise and
# counterclockwise by up to 0.3 pt, and its short shaft sits 0.1-0.4 pt off ±spacing / 2 and
# along the arc; heads within 0.65 of ChemDraw's ink. x0.6's 60.6° unbalanced short shafts end
# 1.1 pt off (a ChemDraw outlier).
EQ_POS_TOL, EQ_HEAD_TOL = 0.45, 0.65
EQ_KNOWN = {("arcs_eq_x0.6", -60.6, False): 1.45, ("arcs_eq_x0.6", 60.6, False): 1.45}


@pytest.mark.parametrize("ref", EQ_REFS)
def test_arc_equilibrium_vs_chemdraw(ref):
    """Each arcs_eq_x* equilibrium arrow (not ChemDraw's elliptical glitches) redrawn on its
    circle (arrows.arc_step_arrow): per shaft, its offset from the circle, butt and notch along
    it, and every head point's distance to ChemDraw's head."""
    from acs_chem_svg.arrows import arc_step_arrow
    preset = load_preset(EQ_REFS[ref])
    bad, n = [], 0
    for x in mr.arc_eq_arrows(ref):
        if x["kind"] == "none" or x["glitch"]:
            continue
        c, R, a0, s = x["center"], x["R"], x["a0"], math.copysign(1.0, x["sweep"])
        g = arc_step_arrow(x["kind"], c, R, a0, a0 + math.radians(x["sweep"]), preset, "eq", x["variant"])

        def pos(q):
            d = ((math.atan2(-(q[1] - c[1]), q[0] - c[0]) - a0) * s + math.pi) % (2 * math.pi) - math.pi
            return (d + 2 * math.pi if d < -1.0 else d) * R
        for ours, forward in zip(g.children, (True, False)):
            theirs = next(sh for sh in x["shafts"] if sh["forward"] == forward)
            tol = EQ_KNOWN.get((ref, x["sweep"], forward), EQ_POS_TOL)
            start, end = ours.segments[0][1], ours.segments[-1][-1]
            got = {"offset": math.dist(start, c) - R, "butt": pos(start), "notch": pos(end)}
            for k, v in got.items():
                if abs(v - theirs[k]) > tol:
                    bad.append(f"{x['kind']} {x['sweep']} {'forward' if forward else 'back'} {k}: "
                               f"ours {v:.3f} ChemDraw {theirs[k]:.3f}")
            head = mr.path_samples(theirs["head"])
            worst = max(min(length(sub(q, p)) for p in head) for leaf in ours.head_leaves() for q in leaf.points)
            if worst > max(tol, EQ_HEAD_TOL):
                bad.append(f"{x['kind']} {x['sweep']} {'forward' if forward else 'back'} head: {worst:.3f} off")
            n += 1
    assert n == (36 if ref == "arcs_eq_x1.0" else 24)
    assert not bad, "\n".join(bad)
