"""CDXML export spike gate (acs-style.md §13): what does ChemDraw recompute when it opens our file?

Step 1 (default): write build/cdxml_spike/<case>.cdxml for the gate cases: the M1 references'
templates at ChemDraw's own CDXML coordinates (so a re-save can be compared with the original
reference), and render.py --cdxml's output for examples/m1_smiles_x1.0.yaml.

Step 2 (by hand, in ChemDraw): open each file, note any warning, and save it as CDXML and as SVG
(File > Save As, same name) into build/cdxml_spike/resaved/. Then drag a labeled atom (does its
charge follow?) and run Clean Structure on a copy; neither is saved.

Step 3: --diff compares each re-saved file with our export, nodes paired by position:
- the attributes ChemDraw added, removed or rewrote, counted per element and attribute, with the
  largest numeric change (label p, BoundingBox, justification, NeedsClean, charge boxes,
  DoublePosition on bonds, ...);
- for the reference cases with a re-saved SVG, its ink vs the original reference's, both drawn
  by ChemDraw at the same positions (measure_refs): label baselines and anchors, bond ends at
  labels, charge centers. Equal ink means ChemDraw drew our file exactly as its own.

Usage: uv run python scripts/cdxml_spike.py          # write the gate files
       uv run python scripts/cdxml_spike.py --diff   # after the ChemDraw re-saves
"""

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

from lxml import etree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import measure_refs as mr  # noqa: E402
from compare import render_reference  # noqa: E402

from acs_chem_svg.cdxml_writer import spec_cdxml, to_cdxml  # noqa: E402
from acs_chem_svg.draw_mol import Drawer  # noqa: E402
from acs_chem_svg.geometry import length, sub  # noqa: E402

OUT = mr.REPO / "build" / "cdxml_spike"
RESAVED = OUT / "resaved"
REFS = ["labels_x1.0", "m1_set_x1.0", "m1_set_x0.6", "stereo_x1.0", "double_bonds_x1.0", "clearance_x1.0"]
SPECS = ["examples/m1_smiles_x1.0.yaml"]
PAIR_TOL = 0.05  # pt: node positions should survive a re-save unchanged
NUM_TOL = 0.005
PARSER = etree.XMLParser(load_dtd=False, no_network=True, remove_blank_text=True)


def write() -> list[Path]:
    OUT.mkdir(parents=True, exist_ok=True)
    RESAVED.mkdir(exist_ok=True)
    paths = []
    for ref in REFS:
        r = render_reference(ref)
        off = (-r.ref.offset[0], -r.ref.offset[1])  # templates hold its SVG frame: back to CDXML's
        text = to_cdxml([(m, Drawer(m, r.preset), off) for m in r.molecules], r.preset.name, f"{ref}.cdxml",
                        warn=lambda m: None)
        paths.append(OUT / f"{ref}.cdxml")
        paths[-1].write_text(text)
    for s in SPECS:
        paths.append(OUT / f"{Path(s).stem}.cdxml")
        paths[-1].write_text(spec_cdxml(mr.REPO / s))
    return paths


# --------------------------------------------------------------------------- attribute diff


def _floats(v: str | None):
    try:
        return [float(x) for x in v.split()] if v else None
    except ValueError:
        return None


class Changes:
    """(where, attribute) -> counts of added / removed / changed, and the largest numeric change."""

    def __init__(self):
        self.counts = collections.defaultdict(collections.Counter)
        self.delta = collections.defaultdict(float)
        self.examples = {}

    def compare(self, where: str, ours, theirs, skip=()) -> None:
        for k in sorted(set(ours.attrib) | set(theirs.attrib)):
            if k in skip:
                continue
            a, b = ours.get(k), theirs.get(k)
            if a == b:
                continue
            fa, fb = _floats(a), _floats(b)
            if fa and fb and len(fa) == len(fb):
                d = max(abs(x - y) for x, y in zip(fa, fb))
                if d <= NUM_TOL:
                    continue
                self.delta[(where, k)] = max(self.delta[(where, k)], d)
            kind = "added" if a is None else "removed" if b is None else "changed"
            self.counts[(where, k)][kind] += 1
            self.examples.setdefault((where, k), f"ours {a!r} -> ChemDraw {b!r}")

    def report(self) -> list[str]:
        out = []
        for key in sorted(self.counts):
            c = self.counts[key]
            d = f", max |delta| {self.delta[key]:.3f}" if key in self.delta else ""
            out.append(f"  {key[0]:9} {key[1]:20} " + ", ".join(f"{n} {k}" for k, n in sorted(c.items()))
                       + f"{d}  (e.g. {self.examples[key]})")
        return out or ["  no attribute changes"]


def _top_nodes(root):
    return [n for f in root.iter("fragment") if f.getparent().tag == "page" for n in f.findall("n")]


def _p(el):
    return tuple(_floats(el.get("p"))[:2])


def attribute_diff(ours_root, theirs_root) -> list[str]:
    ch, lines = Changes(), []
    ch.compare("CDXML", ours_root, theirs_root, skip=("Name",))
    on, tn = _top_nodes(ours_root), _top_nodes(theirs_root)
    pair = {}
    for o in on:
        m = [t for t in tn if length(sub(_p(t), _p(o))) < PAIR_TOL]
        if len(m) == 1:
            pair[o.get("id")] = m[0]
        else:
            lines.append(f"  node {o.get('id')} at {_p(o)}: {len(m)} nodes in the re-save (moved?)")
    lines.append(f"  nodes: {len(on)} ours, {len(tn)} re-saved, {len(pair)} paired")
    for o in on:
        t = pair.get(o.get("id"))
        if t is None:
            continue
        el = o.get("Element", "6")
        ch.compare(f"n[{el}]", o, t, skip=("id", "Z"))
        ot, tt = o.find("t"), t.find("t")
        if (ot is None) != (tt is None):
            ch.counts[(f"n[{el}]", "<t>")]["added" if ot is None else "removed"] += 1
        elif ot is not None:
            ch.compare("t", ot, tt)
            os_, ts_ = ot.findall("s"), tt.findall("s")
            if [s.text for s in os_] != [s.text for s in ts_]:
                ch.counts[("s", "text")]["changed"] += 1
                ch.examples.setdefault(("s", "text"), f"ours {[s.text for s in os_]} -> {[s.text for s in ts_]}")
            for a, b in zip(os_, ts_):
                ch.compare("s", a, b)
    tid = {t.get("id"): o_id for o_id, t in pair.items()}
    tb = {frozenset((tid.get(b.get("B")), tid.get(b.get("E")))): b for b in theirs_root.iter("b")}
    for b in ours_root.iter("b"):
        t = tb.get(frozenset((b.get("B"), b.get("E"))))
        if t is None:
            # a wedge whose ends ChemDraw swapped (WedgeEnd) shows up here
            ch.counts[("b", "<b> unpaired")]["removed"] += 1
            continue
        ch.compare("b", b, t, skip=("id", "Z", "B", "E"))
        if t.get("B") != pair[b.get("B")].get("id"):
            ch.counts[("b", "B/E swapped")]["changed"] += 1
    tg = {}
    for g in theirs_root.iter("graphic"):
        r = g.find("represent")
        tg[tid.get(r.get("object")) if r is not None else None] = g
    for g in ours_root.iter("graphic"):
        t = tg.get(g.find("represent").get("object"))
        if t is None:
            ch.counts[("graphic", "<graphic>")]["removed"] += 1
            continue
        ch.compare("graphic", g, t, skip=("id", "Z"))
    of = [f for f in ours_root.iter("fragment") if f.getparent().tag == "page"]
    tf = {f: f for f in theirs_root.iter("fragment") if f.getparent().tag == "page"}
    for f in of:
        n = pair.get(f.find("n").get("id"))
        if n is not None and n.getparent() in tf:
            ch.compare("fragment", f, n.getparent(), skip=("id", "Z"))
    extra = collections.Counter(el.tag for el in theirs_root.find("page").iter()
                                if el.tag not in ("page", "fragment", "n", "t", "s", "b", "graphic", "represent"))
    if extra:
        lines.append(f"  new elements in the re-save: {dict(extra)}")
    return lines + ch.report()


# --------------------------------------------------------------------------- ink diff


def _cdxml_xy(doc, p):
    return (round(p[0] - doc.offset[0], 2), round(p[1] - doc.offset[1], 2))


def ink_diff(ref: str) -> list[str]:
    """The re-saved SVG vs the original reference, both drawn by ChemDraw at the same positions."""
    a, b = mr.load(ref), mr.load(ref, RESAVED)
    ma, mb = mr.measure(a), mr.measure(b)
    pos_b = {_cdxml_xy(b, at.p): at.id for at in b.atoms.values()}
    a_to_b = {}
    for at in a.atoms.values():
        q = _cdxml_xy(a, at.p)
        m = [i for p, i in pos_b.items() if abs(p[0] - q[0]) < PAIR_TOL and abs(p[1] - q[1]) < PAIR_TOL]
        if m:
            a_to_b[at.id] = m[0]
    worst = collections.defaultdict(float)
    lb = {e["atom"]: e for e in mb["labels"]}
    for e in ma["labels"]:
        o = lb.get(a_to_b.get(e["atom"]))
        if o is None or "error" in e or "error" in o:
            continue
        if e["glyph_run"] != o["glyph_run"]:
            worst[f"label text {e['glyph_run']} -> {o['glyph_run']}"] += 1
        for k in ("baseline_minus_atom_y", "anchor_advance_center_x_minus_atom"):
            worst[f"label {k}"] = max(worst[f"label {k}"], abs(e[k] - o[k]))
        ends_a = sorted(x["end_dist_from_atom"] for x in e.get("bonds", []))
        ends_b = sorted(x["end_dist_from_atom"] for x in o.get("bonds", []))
        if len(ends_a) == len(ends_b):
            for x, y in zip(ends_a, ends_b):
                worst["bond end at label"] = max(worst["bond end at label"], abs(x - y))
        else:
            worst["bond ends at label: count differs"] += 1
    ca = [_cdxml_xy(a, c["center"]) for c in ma["charges"]]
    cb = [_cdxml_xy(b, c["center"]) for c in mb["charges"]]
    if len(ca) != len(cb):
        worst[f"charge count {len(ca)} -> {len(cb)}"] += 1
    for c in ca:
        if cb:
            worst["charge center"] = max(worst["charge center"], min(length(sub(c, d)) for d in cb))
    sa, sb = mr.summarize(ma), mr.summarize(mb)
    for k in sa:
        if isinstance(sa[k], (int, float)) and isinstance(sb.get(k), (int, float)) and abs(sa[k] - sb[k]) > 0.01:
            worst[f"summary {k} ({sa[k]} -> {sb[k]})"] = abs(sa[k] - sb[k])
    return [f"  {k:45} {v:.3f}" for k, v in worst.items()] or ["  ink identical"]


def diff() -> int:
    found = False
    for path in sorted(OUT.glob("*.cdxml")):
        saved = RESAVED / path.name
        if not saved.exists():
            print(f"{path.stem}: no re-save in {RESAVED.relative_to(mr.REPO)}")
            continue
        found = True
        print(f"{path.stem}: attributes (ours -> ChemDraw's re-save)")
        for line in attribute_diff(etree.parse(str(path), PARSER).getroot(), etree.parse(str(saved), PARSER).getroot()):
            print(line)
        if path.stem in REFS and (RESAVED / f"{path.stem}.svg").exists():
            print(f"{path.stem}: ink (original reference -> re-save)")
            for line in ink_diff(path.stem):
                print(line)
    return 0 if found else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--diff", action="store_true", help="compare ChemDraw's re-saves with our exports")
    args = ap.parse_args(argv)
    if args.diff:
        return diff()
    for p in write():
        print(p.relative_to(mr.REPO))
    print(f"Open each in ChemDraw; save as CDXML and SVG into {RESAVED.relative_to(mr.REPO)}/, same names.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
