"""M0 spike: RDKit MolDraw2DSVG (ACS1996 mode) vs ChemDraw on the M1 set (design.md §4.1).

Each M1 molecule is rebuilt in RDKit from the ChemDraw CDXML: same atoms, coordinates, bond
orders, charges, H counts, and wedge/hash flags. Any difference is then a drawing convention,
not a layout difference. RDKit's SVG is mapped into the ChemDraw frame with a least-squares
fit on atom draw coordinates, and the bond-level measurements from measure_refs.py are run on
both. Writes:

    spike/m1_x{k}_rdkit.svg               RDKit output (as drawn)
    spike/m1_x{k}_compare.png             ChemDraw | RDKit | overlay, equal pt scale
    spike/metrics.json                    numeric table used by references/renderer-spike.md

Usage: uv run python scripts/spike_rdkit.py
"""

from __future__ import annotations

import json
import re
import statistics
import sys
from pathlib import Path

import numpy as np
from lxml import etree
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import measure_refs as mr  # noqa: E402

from acs_chem_svg.preview import compare_png, zoom_crops  # noqa: E402

REPO = mr.REPO
OUT = REPO / "spike"
PX_PER_PT = 4
ANUM = {v: k for k, v in mr.ELEMENTS.items()}
DIRS = {"WedgeBegin": Chem.BondDir.BEGINWEDGE, "WedgedHashBegin": Chem.BondDir.BEGINDASH}


# --------------------------------------------------------------------------- build + draw


def build_mol(doc: mr.Doc, pt_per_unit: float = 9.6, atom_ids=None, map_ids: bool = False,
              nicknames: bool = False):
    """RDKit mol with the CDXML atoms/bonds/coords; returns (mol, cdxml atom ids by index).

    Coordinates are pt / pt_per_unit with y up (9.6 turns 14.4 pt bonds into 1.5 A). Pass
    `atom_ids` to keep one fragment, and `map_ids` to store each CDXML id as the atom map.
    `nicknames` keeps nickname nodes (OMe) as their attachment atom; the caller sets the alias."""
    rw = Chem.RWMol()
    ids = []
    for aid, a in doc.atoms.items():
        if (a.nickname and not nicknames) or (atom_ids is not None and aid not in atom_ids):
            continue
        at = Chem.Atom(ANUM[a.element])
        at.SetFormalCharge(a.charge)
        if a.label and a.num_h is not None and not a.nickname:
            at.SetNumExplicitHs(a.num_h)
            at.SetNoImplicit(True)
        if a.radicals:  # tBu·, Cl· (mech_basic): no implicit H on the radical center
            at.SetNumRadicalElectrons(a.radicals)
            at.SetNoImplicit(True)
            if not a.label:
                at.SetNumExplicitHs(4 - a.radicals - sum(1 for b in doc.bonds if aid in (b.a, b.b)))
        if map_ids:
            at.SetAtomMapNum(int(aid))
        rw.AddAtom(at)
        ids.append(aid)
    idx = {aid: i for i, aid in enumerate(ids)}
    order = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    for b in doc.bonds:
        if b.a not in idx or b.b not in idx:
            continue
        a, e = (b.a, b.b) if not b.display.endswith("End") else (b.b, b.a)
        rw.AddBond(idx[a], idx[e], order[b.order])
        display = b.display.replace("End", "Begin")  # *End: the bond was reversed above
        if display in DIRS:
            rw.GetBondBetweenAtoms(idx[a], idx[e]).SetBondDir(DIRS[display])
    mol = rw.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for aid, i in idx.items():
        x, y = doc.atoms[aid].p
        conf.SetAtomPosition(i, Point3D(x / pt_per_unit, -y / pt_per_unit, 0.0))
    mol.AddConformer(conf, assignId=True)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                     ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    Chem.AssignChiralTypesFromBondDirs(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol, ids


def draw(mol, k: float) -> tuple[str, list]:
    d2d = rdMolDraw2D.MolDraw2DSVG(-1, -1)
    opts = d2d.drawOptions()
    Draw.SetACS1996Mode(opts, Draw.MeanBondLength(mol))
    opts.fixedBondLength *= k
    opts.fixedFontSize = int(round(10 * k))
    opts.minFontSize = 4
    opts.clearBackground = False
    opts.prepareMolsBeforeDrawing = False  # keep CDXML wedges and Kekulé bonds
    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()
    coords = [d2d.GetDrawCoords(i) for i in range(mol.GetNumAtoms())]
    return d2d.GetDrawingText(), [(p.x, p.y) for p in coords]


# --------------------------------------------------------------------------- parse RDKit SVG

_TOK = re.compile(r"[MLQCZ]|-?\d*\.?\d+(?:e-?\d+)?")


def _segments(d):
    """Polylines from a path (Q/C flattened), each a list of points; closed flag per run."""
    toks = _TOK.findall(d)
    runs, cur, cmd, i, last = [], [], None, 0, None
    while i < len(toks):
        t = toks[i]
        if t in "MLQCZ":
            cmd = t
            i += 1
            if t == "Z" and cur:
                runs.append((cur, True))
                cur = []
            continue
        if cmd == "M":
            if cur:
                runs.append((cur, False))
            last = (float(toks[i]), float(toks[i + 1]))
            cur = [last]
            i += 2
            cmd = "L"
        elif cmd == "L":
            last = (float(toks[i]), float(toks[i + 1]))
            cur.append(last)
            i += 2
        elif cmd in ("Q", "C"):
            n = 2 if cmd == "Q" else 3
            pts = [last] + [(float(toks[i + 2 * j]), float(toks[i + 2 * j + 1])) for j in range(n)]
            for s in range(1, 7):
                t_ = s / 6
                q = pts
                while len(q) > 1:
                    q = [((1 - t_) * a[0] + t_ * b[0], (1 - t_) * a[1] + t_ * b[1]) for a, b in zip(q, q[1:])]
                cur.append(q[0])
            last = pts[-1]
            i += 2 * n
        else:
            i += 1
    if cur:
        runs.append((cur, False))
    return runs


def _style(el):
    st = dict(kv.split(":", 1) for kv in el.get("style", "").split(";") if ":" in kv)
    for k in ("fill", "stroke", "stroke-width"):
        if el.get(k):
            st[k] = el.get(k)
    return st


def parse_rdkit(svg: str, transform, scale: float = 1.0):
    """Shapes and per-atom glyph boxes, mapped into the ChemDraw frame. `scale` is the
    transform's scale, applied to stroke widths."""
    root = etree.fromstring(svg.encode())
    shapes, glyphs = [], {}
    for el in root.iter("{http://www.w3.org/2000/svg}path"):
        cls = el.get("class", "")
        st = _style(el)
        sw = float(st.get("stroke-width", "0").rstrip("px")) * scale if st.get("stroke", "none") != "none" else 0.0
        filled = st.get("fill", "none") not in ("none", "")
        runs = _segments(el.get("d", ""))
        if cls.startswith("atom-"):
            a = int(cls.split()[0].split("-")[1])
            for pts, _ in runs:
                glyphs.setdefault(a, []).append([transform(p) for p in pts])
            continue
        if not cls.startswith("bond-"):
            continue
        for pts, closed in runs:
            pts = [transform(p) for p in pts]
            if filled and closed:
                sh = mr.Shape("poly", pts)
                sh.width = sw  # stroked outline widens a filled wedge by ~sw
                shapes.append(sh)
            else:
                shapes += mr.stroke_rects(pts, sw)  # stroked line, butt caps -> rectangles
    return shapes, glyphs


def fit_transform(src, dst):
    """Least-squares x' = s*x + tx, y' = s*y + ty."""
    a = np.array(src)
    b = np.array(dst)
    ac, bc = a.mean(0), b.mean(0)
    s = float(((a - ac) * (b - bc)).sum() / ((a - ac) ** 2).sum())
    t = bc - s * ac
    resid = float(np.abs(a * s + t - b).max())
    return (lambda p: (s * p[0] + t[0], s * p[1] + t[1])), s, t, resid


# --------------------------------------------------------------------------- metrics


def label_metrics(doc, asg, boxes_by_atom):
    """Per labeled atom: element glyph box (nearest glyph to the atom), bond end distance,
    and gap to that box along the bond axis. Same definition for both renderers."""
    rows = []
    for aid, boxes in boxes_by_atom.items():
        a = doc.atoms[aid]
        if not boxes:
            continue
        elem = min(boxes, key=lambda b: mr.length(mr.sub(mr._box_center(b), a.p)))
        others = [b for b in boxes if b is not elem]
        side = None
        if others:
            ox = statistics.fmean(mr._box_center(b)[0] for b in others)
            side = "right" if ox > mr._box_center(elem)[0] else "left"
        for b in doc.bonds:
            if aid not in (b.a, b.b) or b.id not in asg or "lines" not in asg[b.id]:
                continue
            other = b.b if b.a == aid else b.a
            fr = mr.frame(doc, other, aid)
            ln = mr._line_from_shape(asg[b.id]["lines"][0].shape, fr)
            end = fr.world(ln.u1, 0.0)
            gap = mr.ray_box_gap(end, fr.u, elem)
            rows.append({"key": f"{aid}-{b.id}", "atom": a.label, "bond_order": b.order,
                         "dir_deg": round(float(np.degrees(np.arctan2(-fr.u[1], fr.u[0]))), 1),
                         "end_dist": round(fr.L - ln.u1, 3),
                         "gap_elem_box": round(gap, 3) if gap is not None else None,
                         "elem_box_h": round(elem[3] - elem[1], 3),
                         "elem_center_minus_atom": tuple(round(v, 3) for v in mr.sub(mr._box_center(elem), a.p)),
                         "h_side": side})
    return rows


def bond_metrics(doc, asg):
    m = {"line_width": [], "ring_spacing": [], "ring_inset": [], "acyclic": [], "triple": [],
         "wedge": None, "hash": None}
    for b in doc.bonds:
        r = asg.get(b.id)
        if not r:
            continue
        if "stereo" in r:
            fr = r["frame"]
            locs = [loc for _, loc in r["shapes"]]
            if "Hash" in b.display:
                us = sorted(statistics.fmean(p[0] for p in loc) for loc in locs)
                m["hash"] = {"count": len(locs),
                             "pitch": round(statistics.median([q - p for p, q in zip(us, us[1:])]), 3) if len(us) > 1 else None,
                             "last_len": round(max(max(p[1] for p in loc) - min(p[1] for p in loc) for loc in locs), 3) if locs else None}
            elif locs:
                pts = [p for loc in locs for p in loc]
                umax = max(p[0] for p in pts)
                wide = [p[1] for p in pts if p[0] > umax - 0.3]
                extra = max((sh.width for sh, _ in r["shapes"]), default=0.0)
                m["wedge"] = {"wide_end_width": round(max(wide) - min(wide) + extra, 3), "end_u": round(umax, 3),
                              "outline_stroke": extra}
            continue
        lines = r.get("lines", [])
        m["line_width"] += [ln.width for ln in lines]
        if b.order == 2 and len(lines) == 2:
            sp = lines[1].v - lines[0].v
            if b.ring:
                inner = max(lines, key=lambda ln: abs(ln.v))
                m["ring_spacing"].append(sp)
                A, B = doc.atoms[b.a], doc.atoms[b.b]
                if not A.label:
                    m["ring_inset"].append(inner.u0)
                if not B.label:
                    m["ring_inset"].append(r["frame"].L - inner.u1)
            else:
                v = sorted(ln.v for ln in lines)
                m["acyclic"].append({"ends": f"{doc.atoms[b.a].label or 'C'}={doc.atoms[b.b].label or 'C'}",
                                     "v": [round(x, 3) for x in v], "spacing": round(sp, 3),
                                     "centered": abs(v[0] + v[1]) < 0.3})
        if b.order == 3:
            m["triple"].append({"n_lines": len(lines), "v": [round(ln.v, 3) for ln in lines]})
    med = lambda xs: round(statistics.median(xs), 3) if xs else None  # noqa: E731
    return {"line_width": med(m["line_width"]), "ring_spacing": med(m["ring_spacing"]),
            "ring_inset": med(m["ring_inset"]), "acyclic": m["acyclic"], "triple": m["triple"],
            "wedge": m["wedge"], "hash": m["hash"]}


def native(m, s):
    """Divide every length in a metrics dict by the RDKit->ChemDraw fit scale, giving the
    values RDKit actually drew (pt). Directions and booleans are untouched."""
    skip = {"key", "dir_deg", "count", "n_lines", "bond_order", "centered", "h_side", "atom", "ends"}

    def walk(x, key=None):
        if isinstance(x, dict):
            return {k: walk(v, k) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [walk(v, key) for v in x]
        if isinstance(x, bool) or key in skip or x is None or isinstance(x, str):
            return x
        if isinstance(x, (int, float, np.floating)):
            return round(float(x) / s, 3)
        return x
    return walk(m)


# --------------------------------------------------------------------------- main


def run(k: str) -> dict:
    name = f"m1_set_x{k}"
    doc = mr.load(name)
    mol, ids = build_mol(doc)
    svg, coords = draw(mol, float(k))
    (OUT / f"m1_x{k}_rdkit.svg").write_text(svg)
    transform, s, t, resid = fit_transform(coords, [doc.atoms[a].p for a in ids])
    shapes, glyph_runs = parse_rdkit(svg, transform, s)

    rd = mr.Doc(name + "_rdkit", atoms=doc.atoms, bonds=doc.bonds, shapes=shapes, settings=doc.settings)
    cd_asg = mr.assign(doc)
    rd_asg = mr.assign(rd)

    # glyph boxes per labeled atom
    cd_boxes = {}
    for aid, a in doc.atoms.items():
        if a.label and not a.nickname:
            _, whole = mr.label_glyphs(doc, a)
            cd_boxes[aid] = [g.box for g in whole]
    rd_boxes = {}
    for i, polys in glyph_runs.items():
        aid = ids[i]
        # RDKit emits one path per glyph (possibly several contours); group contours by overlap
        boxes = []
        for pts in polys:
            bb = (min(p[0] for p in pts), min(p[1] for p in pts), max(p[0] for p in pts), max(p[1] for p in pts))
            for j, ob in enumerate(boxes):
                if bb[0] >= ob[0] - 0.05 and bb[2] <= ob[2] + 0.05 and bb[1] >= ob[1] - 0.05 and bb[3] <= ob[3] + 0.05:
                    break
                if ob[0] >= bb[0] - 0.05 and ob[2] <= bb[2] + 0.05 and ob[1] >= bb[1] - 0.05 and ob[3] <= bb[3] + 0.05:
                    boxes[j] = bb
                    break
            else:
                boxes.append(bb)
        rd_boxes[aid] = boxes
    unclaimed = []
    for sh in shapes:
        if sh.owner:
            continue
        c = sh.centroid
        b = min(doc.bonds, key=lambda b: mr.seg_point_dist(c, doc.atoms[b.a].p, doc.atoms[b.b].p))
        fr = mr.frame(doc, b.a, b.b)
        us = [fr.local(p)[0] for p in sh.pts]
        unclaimed.append({"near_bond": f"{doc.atoms[b.a].label or 'C'}{'=' if b.order == 2 else '-'}"
                                       f"{doc.atoms[b.b].label or 'C'}",
                          "drawn_length_native": round((max(us) - min(us)) / s, 3)})

    w_pt, h_pt = doc.svg_size
    cd_svg = (mr.REF_DIR / f"{name}.svg").read_text()
    to_cd = f"matrix({s} 0 0 {s} {t[0]} {t[1]})"
    compare_png(cd_svg, svg, w_pt, h_pt, OUT / f"m1_x{k}_compare.png", px_per_pt=PX_PER_PT, labels=(
        f"ChemDraw 26.1 (x{k})", f"RDKit ACS1996 mode (x{k})"), transform=to_cd)
    if k == "1.0":
        zoom_crops(cd_svg, svg, w_pt, h_pt, list(doc.fragments.values()), OUT / "m1_x1.0_zoom.png",
                   transform=to_cd)
    rdkit_mapped = {"bonds": bond_metrics(rd, rd_asg), "labels": label_metrics(rd, rd_asg, rd_boxes)}
    return {
        "scale": k,
        "fit": {"scale": round(s, 4), "max_atom_residual_pt": round(resid, 4)},
        "rdkit_native_bond_length": round(float(doc.settings["BondLength"]) / s, 3),
        "rdkit": native(rdkit_mapped, s),
        "rdkit_unassigned_bond_shapes": unclaimed,
        "chemdraw": {"bonds": bond_metrics(doc, cd_asg), "labels": label_metrics(doc, cd_asg, cd_boxes)},
        "rdkit_font": "BuiltinRobotoRegular (ACS1996 mode default)",
    }


def main() -> int:
    OUT.mkdir(exist_ok=True)
    res = {k: run(k) for k in ("1.0", "0.8", "0.6")}
    (OUT / "metrics.json").write_text(json.dumps(res, indent=1, default=str))
    for k, r in res.items():
        print(f"x{k}: fit scale {r['fit']['scale']}, atom residual {r['fit']['max_atom_residual_pt']} pt, "
              f"RDKit native bond {r['rdkit_native_bond_length']}, unassigned RDKit shapes {r['rdkit_unassigned_bond_shapes']}")
        for who in ("chemdraw", "rdkit"):
            b = r[who]["bonds"]
            print(f"  {who:8} lw {b['line_width']} ring_sp {b['ring_spacing']} inset {b['ring_inset']} "
                  f"wedge {b['wedge']} hash {b['hash']}")
            print(f"           acyclic {[(a['ends'], a['centered'], a['spacing']) for a in b['acyclic']]} triple {b['triple']}")
        rd_by_key = {row["key"]: row for row in r["rdkit"]["labels"]}
        for lc in r["chemdraw"]["labels"]:
            lr = rd_by_key.get(lc["key"])
            if lr is None:
                print(f"  label {lc['atom']:4} dir {lc['dir_deg']:6}  (no RDKit bond line found)")
                continue
            print(f"  label {lc['atom']:4} dir {lc['dir_deg']:6}  end CD {lc['end_dist']:6} RD {lr['end_dist']:6} | "
                  f"gap CD {lc['gap_elem_box']} RD {lr['gap_elem_box']} | box_h CD {lc['elem_box_h']} RD {lr['elem_box_h']} "
                  f"| H side CD {lc['h_side']} RD {lr['h_side']} | anchor RD {lr['elem_center_minus_atom']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
