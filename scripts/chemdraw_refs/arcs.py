"""Arc arrows and arc equilibrium arrows (acs-style.md §15.1)."""

from __future__ import annotations

import math
import statistics

from lxml import etree

from acs_chem_svg.geometry import add, dot, length, mul, sub, unit

from .core import REF_DIR, Shape, load, stats
from .measurements import _r


ARC_REFS = {"1.0": "arcs_x1.0", "0.8": "arcs_x0.8", "0.6": "arcs_x0.6"}


def _cubic_samples(pts, n: int = 10) -> list:
    """Points along a stroked path of cubics (pts as parse_path gives them: start, then c1 c2 p3
    per segment)."""
    out, cur = [pts[0]], pts[0]
    for i in range(1, len(pts) - 2, 3):
        c1, c2, p3 = pts[i], pts[i + 1], pts[i + 2]
        for k in range(1, n + 1):
            t, s = k / n, 1 - k / n
            out.append(tuple(s ** 3 * cur[j] + 3 * s * s * t * c1[j] + 3 * s * t * t * c2[j] + t ** 3 * p3[j]
                             for j in range(2)))
        cur = p3
    return out


def _signed_deg(a, b) -> float:
    """Angle (deg) from a to b, positive counterclockwise on the page (y down)."""
    return -math.degrees(math.atan2(a[0] * b[1] - a[1] * b[0], dot(a, b)))


def arc_circle(a, off=(0.0, 0.0)):
    """(center, radius, sweep in degrees, + counterclockwise on the page) of a CDXML arc <arrow>,
    from Head3D, Tail3D and the signed AngularSize alone: the center is on the chord's
    perpendicular bisector, on the left of tail -> head (as seen on the page) for a positive
    sweep. Center3D agrees within 0.07 when the arc was drawn, not flipped (acs-style.md §15.1)."""
    def p(k):
        v = [float(x) for x in a.get(k).split()[:2]]
        return (v[0] + off[0], v[1] + off[1])
    h, t = p("Head3D"), p("Tail3D")
    sweep = float(a.get("AngularSize"))
    d = sub(h, t)
    half = math.radians(abs(sweep)) / 2
    R = length(d) / 2 / math.sin(half)
    left = unit((d[1], -d[0]))  # the page's left of the travel (y down)
    c = add(mul(add(h, t), 0.5), mul(left, math.copysign(R * math.cos(half), sweep)))
    return c, R, sweep


def arc_shaft(doc, tail):
    """(stroke, its points from the tail): the stroked path with an end on the tail. ChemDraw
    strokes every arc counterclockwise on the page, so a clockwise arc's runs notch -> tail."""
    def gap(s):
        return min(length(sub(s.pts[0], tail)), length(sub(s.pts[-1], tail)))
    shaft = min((s for s in doc.shapes if s.kind == "stroke"), key=gap)
    pts = shaft.pts if length(sub(shaft.pts[0], tail)) <= length(sub(shaft.pts[-1], tail)) else shaft.pts[::-1]
    return shaft, pts


def arc_summary(name: str) -> dict:
    """Arc arrows (acs-style.md §15.1): per CDXML <arrow> with AngularSize, the circle
    (arc_circle: from Head3D, Tail3D and AngularSize), the sweep Tail3D -> Head3D on the page,
    where Center3D and the axis ends sit relative to it, and the drawn shaft (the
    stroked cubic path from Tail3D) and head (the curved fill polygon with a vertex on Head3D).
    Head frame: corners = the head's two vertices farthest from the tip on either side of the
    tip -> notch line (notch = the shaft's end); axis = tip - corner midpoint. The axis is
    compared with the chord notch -> tip and the circle's tangents at the tip and the notch."""
    doc = load(name)
    off = doc.offset
    root = etree.parse(str(REF_DIR / f"{name}.cdxml")).getroot()
    lw = float(root.get("LineWidth"))
    graphics = {g.get("SupersededBy"): g for g in root.iter("graphic") if g.get("GraphicType") == "Arc"}
    recs = []
    for a in root.iter("arrow"):
        if a.get("AngularSize") is None:
            continue

        def p(k):
            v = [float(x) for x in a.get(k).split()[:2]]
            return (v[0] + off[0], v[1] + off[1])
        h, t, c3, major, minor = p("Head3D"), p("Tail3D"), p("Center3D"), p("MajorAxisEnd3D"), p("MinorAxisEnd3D")
        c, R, _ = arc_circle(a, off)
        shaft, spts = arc_shaft(doc, t)
        head = next(s for s in doc.shapes if s.kind == "poly" and s.curved and any(length(sub(q, h)) < 0.05 for q in s.pts))
        notch = spts[-1]
        uniq = []
        for q in head.pts:
            if all(length(sub(q, u)) > 1e-3 for u in uniq):
                uniq.append(q)
        chord = unit(sub(h, notch))
        cn = (-chord[1], chord[0])
        c1 = max((q for q in uniq if dot(sub(q, h), cn) > 0.3), key=lambda q: length(sub(q, h)))
        c2 = max((q for q in uniq if dot(sub(q, h), cn) < -0.3), key=lambda q: length(sub(q, h)))
        mid = mul(add(c1, c2), 0.5)
        axis = unit(sub(h, mid))
        n = (-axis[1], axis[0])

        def tangent(q):  # along the travel at q on the circle
            r = sub(q, c)
            tg = unit((-r[1], r[0]))
            return tg if dot(tg, axis) > 0 else mul(tg, -1.0)
        ang = [math.atan2(-(q[1] - c[1]), q[0] - c[0]) for q in (t, h)]
        samples = _cubic_samples(spts)
        turn = math.copysign(1.0, _signed_deg(sub(t, c), sub(samples[len(samples) // 2], c)))  # +1: ccw on the page
        sweep = turn * (math.degrees(turn * (ang[1] - ang[0])) % 360)
        notch_arc = abs(math.atan2(-(notch[1] - c[1]), notch[0] - c[0]) - ang[1])
        notch_arc = min(notch_arc, 2 * math.pi - notch_arc) * R
        bb = [float(x) for x in a.get("BoundingBox").split()]
        ink = (min(q[0] for q in spts + head.pts) - off[0], min(q[1] for q in spts + head.pts) - off[1],
               max(q[0] for q in spts + head.pts) - off[0], max(q[1] for q in spts + head.pts) - off[1])
        g = graphics.get(a.get("id"))
        recs.append({
            "id": a.get("id"),
            "head_attrs": f"{a.get('HeadSize')} / {a.get('ArrowheadCenterSize')} / {a.get('ArrowheadWidth')}",
            "AngularSize": float(a.get("AngularSize")),
            "graphic_AngularSize": g.get("AngularSize") if g is not None else None,
            "graphic_box_minus_(head, center)": g is not None and tuple(
                v - w for v, w in zip([float(x) for x in g.get("BoundingBox").split()],
                                      [h[0] - off[0], h[1] - off[1], c[0] - off[0], c[1] - off[1]])),
            "radius": R, "head_radius_minus_R": length(sub(h, c)) - R, "tail_radius_minus_R": length(sub(t, c)) - R,
            "sweep_tail_to_head (ccw on page +)": sweep,
            "Center3D_minus_center": sub(c3, c),
            "major_minus_Center3D": sub(major, c3), "minor_minus_Center3D": sub(minor, c3),
            "shaft": {"width": shaft.width, "start_minus_tail": length(sub(spts[0], t)),
                      "stroked_from": "tail" if spts is shaft.pts else "notch",
                      "cubics": (len(spts) - 1) // 3,
                      "max_radial_deviation": max(abs(length(sub(q, c)) - R) for q in _cubic_samples(spts)),
                      "end_radius_minus_R": length(sub(notch, c)) - R, "end_arc_from_tip": notch_arc},
            "head": {"length": dot(sub(h, mid), axis), "half_width": (dot(sub(c1, mid), n), -dot(sub(c2, mid), n)),
                     "inset (chord notch -> tip)": length(sub(h, notch)),
                     "notch_off_axis": dot(sub(notch, h), n),
                     "axis_minus_chord_deg": _signed_deg(chord, axis),
                     "axis_minus_tip_tangent_deg": _signed_deg(tangent(h), axis),
                     "axis_minus_notch_tangent_deg": _signed_deg(tangent(notch), axis),
                     "corner_radius_minus_R": (length(sub(c1, c)) - R, length(sub(c2, c)) - R)},
            "outline_stroke": stats([s.width for s in doc.shapes if s.kind == "outline"]),
            "box_minus_ink (x0, y0, x1, y1)": tuple(bb[k] - ink[k] for k in range(4)),
            "half_width_from_attrs": int(a.get("ArrowheadWidth")) * lw / 100,
        })
    return _r({"line_width": lw, "arcs": recs})


ARC_EQ_REFS = {"1.0": "arcs_eq_x1.0", "0.8": "arcs_eq_x0.8", "0.6": "arcs_eq_x0.6"}


def path_anchors(s: Shape) -> list:
    """The on-curve points of a parsed path: its pts are p0, (c1, c2, p3)* when curved."""
    return [q for k, q in enumerate(s.pts) if k % 3 == 0] if s.curved else list(s.pts)


def path_samples(s: Shape, n: int = 8) -> list:
    """Points along a parsed path, its cubics sampled n times each."""
    if not s.curved:
        return list(s.pts)
    out, pts = [s.pts[0]], s.pts
    for k in range(1, len(pts) - 2, 3):
        p0, c1, c2, p3 = pts[k - 1], pts[k], pts[k + 1], pts[k + 2]
        for j in range(1, n + 1):
            t = j / n
            u = 1 - t
            out.append(tuple(u ** 3 * p0[i] + 3 * u * u * t * c1[i] + 3 * u * t * t * c2[i] + t ** 3 * p3[i]
                             for i in range(2)))
    return out


def arc_eq_arrows(name: str) -> list[dict]:
    """The arc arrows of an arcs_eq_x* reference (acs-style.md §15.1): per CDXML <arrow>, its
    circle (Center3D and the major axis end), kind (equilibrium, equilibrium_products, none),
    Tail3D angle a0 and sweep, glitch (minor and major axis radii more than 1 pt apart: an
    ellipse, which ChemDraw draws short), and per drawn shaft (the strokes and half heads whose
    every point is within 4.5 of the circle and the sweep): its offset from the circle, butt and
    notch (its stroke's ends, the notch the one nearer a head), the head's tip (its sampled point
    farthest along the shaft's travel) and the head polygon. Positions are arc lengths along the
    circle from Tail3D, in the direction of travel."""
    doc = load(name)
    root = etree.parse(str(REF_DIR / f"{name}.cdxml")).getroot()
    out = []
    for a in root.iter("arrow"):
        if not a.get("AngularSize"):
            continue

        def p(k, a=a):
            v = [float(x) for x in a.get(k).split()[:2]]
            return (v[0] + doc.offset[0], v[1] + doc.offset[1])
        c = p("Center3D")
        R = p("MajorAxisEnd3D")[0] - c[0]
        sweep = float(a.get("AngularSize"))
        s = 1.0 if sweep > 0 else -1.0
        a0 = math.atan2(-(p("Tail3D")[1] - c[1]), p("Tail3D")[0] - c[0])
        L = math.radians(abs(sweep)) * R

        def pos(q, c=c, a0=a0, s=s, R=R):
            d = ((math.atan2(-(q[1] - c[1]), q[0] - c[0]) - a0) * s + math.pi) % (2 * math.pi) - math.pi
            return (d + 2 * math.pi if d < -1.0 else d) * R
        mine = [sh for sh in doc.shapes if sh.kind in ("stroke", "poly")
                and all(abs(math.dist(q, c) - R) < 4.5 and -12 < pos(q) < L + 12 for q in path_samples(sh))]
        heads = [sh for sh in mine if sh.kind == "poly"]
        shafts = []
        for sh in (sh for sh in mine if sh.kind == "stroke"):
            an = path_anchors(sh)
            ends = [an[0], an[-1]]
            k = min((0, 1), key=lambda k: min((min(math.dist(ends[k], q) for q in path_samples(h)) for h in heads),
                                               default=0.0)) if heads else 1
            notch, butt = ends[k], ends[1 - k]
            sgn = 1.0 if pos(notch) > pos(butt) else -1.0
            head = min(heads, key=lambda h: min(math.dist(notch, q) for q in path_samples(h))) if heads else None
            tip = max(path_samples(head), key=lambda q: sgn * pos(q)) if head else None
            shafts.append({"offset": statistics.fmean(math.dist(q, c) for q in an) - R, "butt": pos(butt),
                           "notch": pos(notch), "forward": sgn > 0, "tip": pos(tip) if tip else None,
                           "tip_offset": math.dist(tip, c) - R if tip else None, "head": head})
        kind = ("none" if not a.get("ArrowheadHead") else
                "equilibrium_products" if a.get("ArrowEquilibriumRatio") else "equilibrium")
        minor = p("MinorAxisEnd3D")[1] - c[1]
        out.append({"center": c, "R": R, "a0": a0, "sweep": sweep, "L": L, "kind": kind,
                    "variant": {"1000": "standard", "1500": "large"}[a.get("HeadSize")],
                    "glitch": abs(minor - R) > 1.0, "shafts": shafts})
    return out


def arc_eq_summary(name: str) -> list[dict]:
    """arc_eq_arrows without the shapes, rounded: --arc-eq."""
    return [{"kind": x["kind"], "variant": x["variant"], "sweep": x["sweep"], "R": round(x["R"], 3),
             "L": round(x["L"], 3), "glitch": x["glitch"],
             "shafts": [{k: (round(v, 3) if isinstance(v, float) else v) for k, v in sh.items() if k != "head"}
                        for sh in x["shafts"]]} for x in arc_eq_arrows(name)]
