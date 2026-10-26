"""M1 acceptance (design.md §10): our drawing of each reference's templates vs ChemDraw, within 0.1 pt.

Both drawings are measured by the same code (scripts/measure_refs.py). Items are paired by
CDXML atom and bond ids. The one known exception to the clip model (acs-style.md §2.1: F
with the bond leaving rightward, 0.31-0.41 pt) is checked separately as xfail.
"""

import copy
from functools import lru_cache

import pytest
from compare import render_reference

import measure_refs as mr

pytestmark = pytest.mark.skipif(not mr.REF_DIR.exists(),
                                reason="ChemDraw references (reference_renders/) are local only")

TOL = 0.1
X1_BAR_LENGTH = 3.935  # m1_set_x1.0 / labels_x1.0 charge bars (acs-style.md §2.2)
REFS = ["m1_set_x1.0", "m1_set_x0.8", "m1_set_x0.6", "double_bonds_x1.0", "labels_x1.0", "double_bonds_2",
        "stereo_x1.0", "stereo_x0.8", "stereo_x0.6", "clearance_x1.0", "clearance_x0.8", "clearance_x0.6"]


@lru_cache(maxsize=None)
def measured(ref):
    r = render_reference(ref)
    return r, mr.measure(r.ref), mr.measure(r.ours)


def _unclaimed(doc):
    """A copy with no shape claimed: measure() already set Shape.owner on the originals."""
    d = copy.deepcopy(doc)
    for sh in d.shapes:
        sh.owner = None
    return d


def _hausdorff(a, b):
    """Largest distance from a vertex of one wedge outline to the other outline's boundary."""
    pa = [p for sh, _ in a["shapes"] for p in sh.pts]
    pb = [p for sh, _ in b["shapes"] for p in sh.pts]

    def edges(pts):
        return list(zip(pts, pts[1:] + pts[:1]))

    def far(src, dst):
        return max(min(mr.seg_point_dist(p, q0, q1) for q0, q1 in edges(dst)) for p in src)
    return max(far(pa, pb), far(pb, pa))


def _known_exception(r, atom_id, other_id):
    """F with the bond leaving rightward: ChemDraw ends 0.31-0.41 pt further out (§2.1)."""
    a, o = r.ref.atoms[atom_id], r.ref.atoms[other_id]
    d = mr.unit(mr.sub(o.p, a.p))
    return a.element == "F" and d[0] > 0.96


def pairs(ref):
    """[(category, key, chemdraw, ours, known_exception)] for every compared quantity."""
    r, mc, mo = measured(ref)
    out = []

    def add(cat, key, a, b, prov=False):
        out.append((cat, key, a, b, prov))

    asg_c, asg_o = mr.assign(_unclaimed(r.ref)), mr.assign(_unclaimed(r.ours))
    for b in r.ref.bonds:
        if b.a not in r.ours.atoms:
            continue
        add("bond length", b.id, mr.frame(r.ref, b.a, b.b).L, mr.frame(r.ours, b.a, b.b).L)
        lc, lo = asg_c.get(b.id, {}).get("lines", []), asg_o.get(b.id, {}).get("lines", [])
        assert len(lc) == len(lo), f"bond {b.id}: {len(lc)} ChemDraw lines vs {len(lo)} ours"
        for k, (x, y) in enumerate(zip(lc, lo)):
            add("line width", f"{b.id}/{k}", x.width, y.width)
            add("line offset v", f"{b.id}/{k}", x.v, y.v)
        if len(lc) >= 2:
            add("double/triple spacing", b.id, lc[-1].v - lc[0].v, lo[-1].v - lo[0].v)
    ring = {e["bond"]: e for e in mo["ring_inner_line"]["entries"]}
    for e in mc["ring_inner_line"]["entries"]:
        o = ring[e["bond"]]
        if not e["a_labeled"]:
            add("ring inset", e["bond"] + "a", e["inset_a"], o["inset_a"])
        if not e["b_labeled"]:
            add("ring inset", e["bond"] + "b", e["inset_b"], o["inset_b"])
    acy = {e["bond"]: e for e in mo["acyclic_double"]}
    for e in mc["acyclic_double"]:
        o = acy[e["bond"]]
        for k, (x, y) in enumerate(zip(e["v"], o["v"])):
            add("acyclic line offset", f"{e['bond']}/{k}", x, y)
        for k, (x, y) in enumerate(zip(e["ends"], o["ends"])):
            add("acyclic line ends", f"{e['bond']}/{k}a", x[0], y[0])
            add("acyclic line ends", f"{e['bond']}/{k}b", x[1], y[1])
    tri = {e["bond"]: e for e in mo["triple"]}
    for e in mc["triple"]:
        for k, (x, y) in enumerate(zip(e["u"], tri[e["bond"]]["u"])):
            add("triple line ends", f"{e['bond']}/{k}a", x[0], y[0])
            add("triple line ends", f"{e['bond']}/{k}b", x[1], y[1])
    lab_o = {e["atom"]: e for e in mo["labels"]}
    for e in mc["labels"]:
        if e.get("nickname") or "error" in e or e["atom"] not in lab_o:
            continue
        o = lab_o[e["atom"]]
        add("label baseline", e["atom"], e["baseline_minus_atom_y"], o["baseline_minus_atom_y"])
        add("label anchor", e["atom"], e["anchor_advance_center_x_minus_atom"], o["anchor_advance_center_x_minus_atom"])
        assert e["glyph_run"] == o["glyph_run"], (e["glyph_run"], o["glyph_run"])
        for x, y in zip(e["other_glyphs"], o["other_glyphs"]):
            if x["ch"].isdigit():
                add("subscript size", f"{e['atom']}{x['ch']}", x["size"], y["size"])
                add("subscript drop", f"{e['atom']}{x['ch']}", x["dbaseline"], y["dbaseline"])
            add("label glyph x", f"{e['atom']}{x['ch']}", x["dx_origin_from_elem_origin"], y["dx_origin_from_elem_origin"])
        ob = {(b["bond"], round(b["line_v"], 1)): b for b in o["bonds"]}
        for b in e["bonds"]:
            y = ob[(b["bond"], round(b["line_v"], 1))]
            bond = next(x for x in r.ref.bonds if x.id == b["bond"])
            prov = _known_exception(r, e["atom"], bond.b if bond.a == e["atom"] else bond.a)
            add("bond end at label", f"{e['glyph_run']} atom {e['atom']} bond {b['bond']} v{b['line_v']:+.1f}",
                b["end_dist_from_atom"], y["end_dist_from_atom"], prov)
    st_o = {s["bond"]: s for s in mo["stereo"]}
    for s in mc["stereo"]:
        o = st_o[s["bond"]]
        if "Hash" in s["display"]:
            add("hash last length", s["bond"], s["hashes"][-1]["len_v"], o["hashes"][-1]["len_v"])
            add("hash line thickness", s["bond"], s["line_thickness"]["median"], o["line_thickness"]["median"])
            add("hash count", s["bond"], s["count"], o["count"])
            if s["count"] > 1:
                add("hash pitch", s["bond"], s["pitch"]["median"], o["pitch"]["median"])
            add("hash end u", s["bond"], s["hashes"][-1]["u"], o["hashes"][-1]["u"])
        else:
            bond = next(x for x in r.ref.bonds if x.id == s["bond"])
            wide = bond.b if bond.display.endswith("Begin") else bond.a
            if r.ref.atoms[wide].label or len(mr.neighbors(r.ref, wide)) == 1:
                add("wedge wide end", s["bond"], s["wide_end_width"], o["wide_end_width"])  # flat end
            add("wedge end u", s["bond"], s["u_range"][1], o["u_range"][1])
            add("wedge outline (Hausdorff)", s["bond"], 0.0, _hausdorff(asg_c[s["bond"]], asg_o[s["bond"]]))
    key = lambda c: c["center"]  # noqa: E731
    for x, y in zip(sorted(mc["charges"], key=key), sorted(mo["charges"], key=key)):
        # deliberate departure (user, 2026-10-24): ChemDraw's bars shrink faster than k; ours
        # stay proportional to the label size, i.e. the x1.0 ChemDraw length x k
        add("charge bar length (x1.0 ChemDraw x k)", x["atom_label"], X1_BAR_LENGTH * r.preset.scale, y["bar_length"])
        add("charge bar thickness", x["atom_label"], x["bar_thickness"], y["bar_thickness"])
    assert len(mc["charges"]) == len(mo["charges"])
    return out


def _fails(rows):
    return [f"{cat} [{key}]: ChemDraw {a} ours {b} (diff {b - a:+.3f})"
            for cat, key, a, b, _ in rows if a is None or b is None or abs(b - a) > TOL]


@pytest.mark.parametrize("ref", REFS)
def test_exact_geometry_within_tolerance(ref):
    rows = [p for p in pairs(ref) if not p[4]]
    cats = {p[0] for p in rows}
    assert {"bond length", "line width", "label baseline"} <= cats or ref == "double_bonds_x1.0"
    fails = _fails(rows)
    assert not fails, "\n".join(fails)


@pytest.mark.parametrize("ref", REFS)
def test_every_bond_and_label_compared(ref):
    """No silent skips: every non-nickname CDXML bond and label is paired."""
    r, mc, mo = measured(ref)
    rows = pairs(ref)
    drawn = set(r.ours.atoms)  # templates skip nickname fragments
    bonds = {k for c, k, *_ in rows if c == "bond length"}
    assert bonds == {b.id for b in r.ref.bonds if b.a in drawn}
    labels = {k for c, k, *_ in rows if c == "label baseline"}
    assert labels == {a.id for a in r.ref.atoms.values() if a.label and a.id in drawn}


@pytest.mark.xfail(reason="known exception to the clip model (acs-style.md §2.1); see the diffs", strict=True,
                   raises=AssertionError)
@pytest.mark.parametrize("ref", ["clearance_x1.0", "clearance_x0.8", "clearance_x0.6"])
def test_known_exceptions_within_tolerance(ref):
    rows = [p for p in pairs(ref) if p[4]]
    assert rows, "expected F bonds leaving rightward"
    report = "\n".join(f"{cat} [{key}]: ChemDraw {a} ours {b} (diff {b - a:+.3f})" for cat, key, a, b, _ in rows)
    print(report)
    assert not _fails(rows), report


def _bond_ink(svg_text, w, h, px, blanks=(), radius=4.0, boxes=()):
    """Alpha mask of bonds, wedges and hashes: text removed, a square of ±radius pt blanked
    around each point in `blanks` (charge signs, labeled atoms), and each box in `boxes`."""
    from lxml import etree

    from acs_chem_svg.preview import alpha
    root = etree.fromstring(svg_text.encode())
    for el in list(root.iter("{http://www.w3.org/2000/svg}text")):
        el.getparent().remove(el)
    a = alpha(etree.tostring(root, encoding="unicode"), w, h, px)
    for cx, cy in blanks:
        a[max(0, int((cy - radius) * px)):int((cy + radius) * px),
          max(0, int((cx - radius) * px)):int((cx + radius) * px)] = 0
    for x0, y0, x1, y1 in boxes:
        a[max(0, int((y0 - 2) * px)):int((y1 + 2) * px), max(0, int((x0 - 2) * px)):int((x1 + 2) * px)] = 0
    return a


@pytest.mark.parametrize("ref", ["double_bonds_x1.0", "m1_set_x0.8", "double_bonds_2", "stereo_x1.0", "stereo_x0.6"])
def test_miter_joins_match_chemdraw_pixels(ref):
    """Stroked, miter-joined polylines reproduce ChemDraw's filled bond polygons (joins,
    pointed degree-3 ends, wedge outlines): at 20 px/pt fewer than 0.2% of bond ink pixels
    differ, away from labels and charges."""
    r, mc, mo = measured(ref)
    w, h = r.ref.svg_size
    # charges and label ends are checked numerically above; this test is about the joins
    blanks = [c["center"] for c in mc["charges"] + mo["charges"]]
    blanks += [a.p for a in r.ref.atoms.values() if a.label]
    undrawn = [box for fid, box in r.ref.fragments.items()
               if not any(a.frag == fid and a.id in r.ours.atoms for a in r.ref.atoms.values())]
    # our hashes are stroked rectangles (editable), ChemDraw's are trapezoids; their geometry
    # is checked numerically above
    for b in r.ref.bonds:
        if "Hash" in b.display:
            (x0, y0), (x1, y1) = r.ref.atoms[b.a].p, r.ref.atoms[b.b].p
            undrawn.append((min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)))
    a = _bond_ink((mr.REF_DIR / f"{ref}.svg").read_text(), w, h, 20, blanks, 7 * r.preset.scale, undrawn)
    b = _bond_ink(r.svg, w, h, 20, blanks, 7 * r.preset.scale, undrawn)
    differing = (abs(a - b) > 0.5).sum() / (a > 0.5).sum()
    assert differing < 0.002, f"{differing:.3%} of ink pixels differ"
