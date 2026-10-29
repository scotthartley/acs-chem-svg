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

Phase 2 (--phase 2, acs-style.md §13.1; files in build/cdxml_spike/phase2/): the scheme
references' round-trip exports (compare.reference_cdx_scheme) and render.py --cdxml's output for
M2 specs. By hand in ChemDraw, besides the
re-saves: select each arrow (do the reaction tools see its reactants and products?), drag a
reactant, and retype a condition (does it stay in the step?). --diff then reports, per file:
1. whether every <scheme> and <step> survived with its ids and citations;
2. the arrow and <graphic> attributes ChemDraw rewrote;
3. caption p, BoundingBox, LineStarts and <s> runs;
4. for the reference cases, the caption glyphs ChemDraw drew in the re-saved SVG against the
   original reference's: line pitch of each stacked caption, and the largest glyph origin shift
   (sub- and superscripts included).

Phase 3 (--phase 3, acs-style.md §13.2; files in build/cdxml_spike/phase3/): mech_basic_x1.0 /
x0.6 round trips (compare.reference_cdx_scheme) and render.py --cdxml's output for M3 specs
(two-electron arrows, fishhooks with HalfRight and HalfLeft barbs, an added reactant above its
arrow). By hand in ChemDraw, besides the re-saves (CDXML and SVG): drag an atom with dots (do they
follow?), move a curve, and run Clean Structure on a copy. --diff reports, besides phase 2's items:
5. the curves: CurveType, ArrowheadHead and CurvePoints kept or rewritten;
6. the LonePair / Electron symbols per atom: boxes and <represent> kept;
7. the re-saved SVG's ink: each curve's stroke end, head tip and barb side against our MechArrow
   (spec cases: tip and head polygon from the layout) or the original reference (reference
   cases), and the drawn dot centers against the centers our symbols declare.

Phase 4a (--phase 4a, acs-style.md §13.3; files in build/cdxml_spike/phase4a/): round trips of
abbrev_x1.0 (templates, as phase 1), resonance_x1.0 and scheme_wrap_x1.0 (as phase 2), and
render.py --cdxml's output for m1_abbrev, m4_resonance, m4_wrap, m4_nitration and m2_sequence.
By hand in ChemDraw, besides the re-saves (CDXML and SVG): expand and contract a nickname (Edit >
Expand Label), drag a contributor (do the brackets and ↔ steps hold?), and move a bracket.
--diff reports, besides phase 2's items:
8. the nickname nodes: NodeType, IsNickname, inner atoms kept or moved, label <t> rewrites;
9. the brackets: BoundingBox and LipSize (phase 2's <graphic> diff), and any <bracketedgroup>
   or objecttag ChemDraw adds (new elements in the re-save).

Phase 4b (--phase 4b, acs-style.md §13.4; files in build/cdxml_spike/phase4b/): render.py
--cdxml's output for m1_charges (circled signs, 2+ / 2-, +3, triplet CH2), m1_lewis,
m1_lewis_lines and m1_lewis_assigned (oxidation states). By hand in ChemDraw, besides the
re-saves (CDXML and SVG): note any warning, drag a Lewis atom and a charged atom (what follows?).
--diff reports, besides phase 2's items:
10. the symbols per type, bound or free, with the largest move, and the free digit texts; node
   Charge / Radical rewrites show in the species attributes.

Draft (--lewis-draft, phase 4b; files in build/cdxml_spike/lewis_draft/): a first draft of the
charges_lewis_x* reference, for the user to edit in ChemDraw and save as the reference, at each
preset. Row 1: the charge structures with their bonds, without charges or unpaired electrons
(ChemDraw's tools place those, so they are measured, not ours read back). Rows 2-4: there is no
ChemDraw Lewis reference to copy (Open 17), so they are our Lewis drawing in ChemDraw's objects:
every atom labeled; shared pairs as dots (bonds kept but Visible="no") and as bond lines (lewis:
{bonds: lines}); formal-charge and oxidation-state views. Each column of two dots is a LonePair
symbol (shared pairs free, assigned and lone pairs bound to their atom), an odd electron an
Electron.

Phase 5 (--phase 5, acs-style.md §14; files in build/cdxml_spike/phase5/): chair_x1.0's templates
(as phase 1, with ChemDraw's front bonds) and render.py --cdxml's output for m1_chair_x1.0. We
write each crossing's rear bond with the lower Z (cdxml_writer.z_order) and no CrossingBonds. By
hand in ChemDraw, besides the re-saves (CDXML and SVG): does each crossing show its gap on the
same bond as our SVG, and does moving an atom keep it? --diff reports, besides phase 2's items:
11. every crossing ChemDraw records (CrossingBonds) with the bond it draws in front (higher Z),
   against ours, by node positions.

Usage: uv run python scripts/cdxml_spike.py [--phase N]          # write the gate files
       uv run python scripts/cdxml_spike.py --lewis-draft        # the Lewis reference draft
       uv run python scripts/cdxml_spike.py [--phase N] --diff   # after the ChemDraw re-saves
"""

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

from lxml import etree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import measure_refs as mr  # noqa: E402
from compare import preset_for, reference_cdx_scheme, render_reference  # noqa: E402

from acs_chem_svg.cdxml_writer import spec_cdxml, to_cdxml  # noqa: E402
from acs_chem_svg.draw_mol import Drawer  # noqa: E402
from acs_chem_svg.geometry import add, dot, length, sub, unit  # noqa: E402
from acs_chem_svg.layout import build_scheme  # noqa: E402
from acs_chem_svg.spec import load_spec  # noqa: E402
from acs_chem_svg.svg_writer import to_svg  # noqa: E402

OUT = mr.REPO / "build" / "cdxml_spike"
RESAVED = OUT / "resaved"
REFS = ["labels_x1.0", "m1_set_x1.0", "m1_set_x0.6", "stereo_x1.0", "double_bonds_x1.0", "clearance_x1.0"]
SPECS = ["examples/m1_smiles_x1.0.yaml"]
PHASE2_REFS = ["scheme_basic_x1.0", "scheme_multi_x1.0", "equilibrium_x1.0"]
PHASE2_SPECS = ["examples/m2_scheme_basic_x1.0.yaml", "examples/m2_stages3_x1.0.yaml",
                "examples/m2_equilibrium_x1.0.yaml", "examples/m2_suzuki_x1.0.yaml"]
PHASE3_REFS = ["mech_basic_x1.0", "mech_basic_x0.6"]
PHASE3_SPECS = ["examples/m3_sn2_x1.0.yaml", "examples/m3_radical_initiation_x1.0.yaml",
                "examples/m3_radical_propagation_x1.0.yaml", "examples/m3_carbonyl_addition_x1.0.yaml"]
PHASE4A_REFS = ["abbrev_x1.0", "resonance_x1.0", "scheme_wrap_x1.0"]
PHASE4A_SPECS = ["examples/m1_abbrev_x1.0.yaml", "examples/m4_resonance_x1.0.yaml", "examples/m4_wrap_x1.0.yaml",
                 "examples/m4_nitration_x1.0.yaml", "examples/m2_sequence_x1.0.yaml"]
PHASE4B_SPECS = ["examples/m1_charges_x1.0.yaml", "examples/m1_lewis_x1.0.yaml", "examples/m1_lewis_lines_x1.0.yaml",
                 "examples/m1_lewis_assigned_x1.0.yaml"]
PHASE5_REFS = ["chair_x1.0"]
PHASE5_SPECS = ["examples/m1_chair_x1.0.yaml"]
PHASES = {"2": (PHASE2_REFS, PHASE2_SPECS), "3": (PHASE3_REFS, PHASE3_SPECS), "4a": (PHASE4A_REFS, PHASE4A_SPECS),
          "4b": ([], PHASE4B_SPECS), "5": (PHASE5_REFS, PHASE5_SPECS)}
PAIR_TOL = 0.05  # pt: node positions should survive a re-save unchanged
NUM_TOL = 0.005
PARSER = etree.XMLParser(load_dtd=False, no_network=True, remove_blank_text=True)


def out_dir(phase: str) -> Path:
    return OUT if phase == "1" else OUT / f"phase{phase}"


def species_reference(ref: str) -> str:
    """An M1 reference's templates exported at ChemDraw's CDXML coordinates."""
    r = render_reference(ref)
    off = (-r.ref.offset[0], -r.ref.offset[1])  # templates hold its SVG frame: back to CDXML's
    return to_cdxml([(m, Drawer(m, r.preset), off) for m in r.molecules], r.preset.name, f"{ref}.cdxml",
                    warn=lambda m: None)


def write(phase: str = "1") -> list[Path]:
    out = out_dir(phase)
    (out / "resaved").mkdir(parents=True, exist_ok=True)
    if phase in PHASES:
        refs, specs = PHASES[phase]
        paths = []
        for ref in refs:
            paths.append(out / f"{ref}.cdxml")
            paths[-1].write_text(species_reference(ref) if ref.startswith(("abbrev_", "chair_")) else
                                 to_cdxml(reference_cdx_scheme(ref), preset_for(ref), f"{ref}.cdxml", warn=lambda m: None))
        for s in specs:
            paths.append(out / f"{Path(s).stem}.cdxml")
            paths[-1].write_text(spec_cdxml(mr.REPO / s, warn=lambda m: None))
        return paths
    paths = []
    for ref in REFS:
        paths.append(OUT / f"{ref}.cdxml")
        paths[-1].write_text(species_reference(ref))
    for s in SPECS:
        paths.append(OUT / f"{Path(s).stem}.cdxml")
        paths[-1].write_text(spec_cdxml(mr.REPO / s))
    return paths


LEWIS_SHARED = {"methanol": "CO", "co2": "O=C=O", "nitrate": "[O-][N+](=O)[O-]", "hcn": "C#N",
                "ccl4": "ClC(Cl)(Cl)Cl"}
LEWIS_ASSIGNED = {"co2_formal": ("O=C=O", "formal_charge"), "hcn_formal": ("C#N", "formal_charge"),
                  "methanol_ox": ("CO", "oxidation_state"), "water_ox": ("O", "oxidation_state")}
# the charge items, each written without its charge and radicals: the user places them with
# ChemDraw's tools (circled signs, 2+ / 2-, a "+3" text, two unpaired electrons)
DRAFT_CHARGES = {"nitromethane": "C[N+](=O)[O-]", "isopropyl": "C[CH+]C",
                 "enolate": "C=C[O-:1]", "iron_ii": "[Fe+2]", "oxide": "[O-2]", "iron_iii": "[Fe]",
                 "methylene": "[CH2]"}
LEWIS_ROW_GAP = 20.0  # draft layout only, not a style value
LEWIS_PRESETS = {"x1.0": "acs1996", "x0.8": "acs1996_x0.8", "x0.6": "acs1996_x0.6"}


def lewis_draft(out: Path = OUT / "lewis_draft") -> list[Path]:
    """The charges_lewis draft files (module docstring), one per preset. Never overwrites: the
    user saved the reference over the first draft (2026-10-28)."""
    from acs_chem_svg.layout import m1_placed
    (out / "specs").mkdir(parents=True, exist_ok=True)
    rows_species = (
        ("charges", {k: f'smiles: "{v}"' + (", show_lone_pairs: [1]" if ":1]" in v else "")
                     for k, v in DRAFT_CHARGES.items()}),
        ("shared", {k: f'smiles: "{v}", lewis: true' for k, v in LEWIS_SHARED.items()}),
        ("lines", {k: f'smiles: "{v}", lewis: {{bonds: lines}}' for k, v in LEWIS_SHARED.items()}),
        ("assigned", {k: f'smiles: "{v}", lewis: {{electrons: {m}}}' for k, (v, m) in LEWIS_ASSIGNED.items()}))
    paths = []
    for scale, preset in LEWIS_PRESETS.items():
        species, charges, y = [], [], 0.0
        for row, items in rows_species:
            spec = out / "specs" / f"{row}_{scale}.yaml"
            spec.write_text(f"preset: {preset}\nspecies:\n" + "".join(f"  {k}: {{{v}}}\n" for k, v in items.items()))
            placed = m1_placed(load_spec(spec), spec.parent)
            boxes = [p.group.bbox() for p in placed]
            top, bottom = min(b[1] for b in boxes), max(b[3] for b in boxes)
            off = (LEWIS_ROW_GAP - min(b[0] for b in boxes), y + LEWIS_ROW_GAP - top)
            for p in placed:
                if row == "charges":
                    charges.append(len(species))
                species.append((p.mol, p.drawer, add(p.offset, off)))
            y += LEWIS_ROW_GAP + bottom - top
        root = etree.fromstring(to_cdxml(species, preset, f"charges_lewis_draft_{scale}.cdxml", warn=lambda m: None).encode())
        frags = root.find("page").findall("fragment")
        for i in charges:  # left for the user: charges and unpaired electrons
            for n in frags[i].findall("n"):
                n.attrib.pop("Charge", None)
                n.attrib.pop("Radical", None)
            for g in frags[i].findall("graphic"):
                if g.get("SymbolType") in ("Plus", "Minus", "Electron"):
                    frags[i].remove(g)
        paths.append(out / f"charges_lewis_draft_{scale}.cdxml")
        if paths[-1].exists():
            raise SystemExit(f"{paths[-1].relative_to(mr.REPO)} exists; not overwritten")
        paths[-1].write_bytes(etree.tostring(root.getroottree(), xml_declaration=True, encoding="UTF-8",
                                             doctype=root.getroottree().docinfo.doctype, pretty_print=True))
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
        if b.getparent().getparent().tag == "n":
            continue  # a nickname's inner bond: nickname_diff
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
        if g.get("SymbolType") in DOTS:
            continue
        r = g.find("represent")
        tg[tid.get(r.get("object")) if r is not None else None] = g
    for g in ours_root.iter("graphic"):
        if g.find("represent") is None or g.get("SymbolType") in DOTS:
            continue  # an arrow's legacy <graphic>, ovals: scheme_diff / symbol_diff; dots: mech_diff
        t = tg.get(g.find("represent").get("object"))
        if t is None:
            ch.counts[("graphic", "<graphic>")]["removed"] += 1
            ch.examples.setdefault(("graphic", "<graphic>"), f"{g.get('SymbolType')} on node {g.find('represent').get('object')}")
            continue
        ch.compare("graphic", g, t, skip=("id", "Z"))
    of = [f for f in ours_root.iter("fragment") if f.getparent().tag == "page"]
    tf = {f: f for f in theirs_root.iter("fragment") if f.getparent().tag == "page"}
    for f in of:
        n = pair.get(f.find("n").get("id"))
        if n is not None and n.getparent() in tf:
            ch.compare("fragment", f, n.getparent(), skip=("id", "Z"))
    extra = collections.Counter(el.tag for el in theirs_root.find("page").iter()
                                if el.tag not in ("page", "fragment", "n", "t", "s", "b", "graphic", "represent")
                                and el.tag not in [x.tag for x in ours_root.find("page").iter()])
    if extra:
        lines.append(f"  new elements in the re-save: {dict(extra)}")
    return lines + ch.report()


# --------------------------------------------------------------------------- phase 4a: nicknames


def nickname_diff(ours_root, theirs_root) -> list[str]:
    """What ChemDraw's re-save did to the nickname nodes (module docstring, item 8): each of ours
    paired by p; NodeType and IsNickname kept, its inner fragment's atoms (by Element and
    position) and bonds, and the largest inner atom move."""
    lines, ch = [], Changes()
    theirs = [n for n in _top_nodes(theirs_root)]
    worst, kept = 0.0, 0
    ours = [n for n in _top_nodes(ours_root) if n.find("fragment") is not None]
    for o in ours:
        m = [t for t in theirs if length(sub(_p(t), _p(o))) < PAIR_TOL]
        if len(m) != 1:
            lines.append(f"  nickname {_text(o.find('t'))} at {_p(o)}: {len(m)} nodes in the re-save")
            continue
        t = m[0]
        key = f"nickname[{o.get('NodeType')}]"
        for attr in ("NodeType", "Element", "NeedsClean"):
            if o.get(attr) != t.get(attr):
                ch.counts[(key, attr)]["changed"] += 1
                ch.examples.setdefault((key, attr), f"{_text(o.find('t'))}: {o.get(attr)} -> {t.get(attr)}")
        if (o.find("t") is None) == (t.find("t") is None) and o.find("t").get("IsNickname") != t.find("t").get("IsNickname"):
            ch.counts[(key, "IsNickname")]["changed"] += 1
        fi, ti = o.find("fragment"), t.find("fragment")
        if ti is None:
            ch.counts[(key, "inner <fragment>")]["removed"] += 1
            continue
        on_, tn_ = fi.findall("n"), ti.findall("n")
        if sorted(n.get("Element", "6") for n in on_) != sorted(n.get("Element", "6") for n in tn_) or \
                len(fi.findall("b")) != len(ti.findall("b")):
            ch.counts[(key, "inner atoms or bonds")]["changed"] += 1
            ch.examples.setdefault((key, "inner atoms or bonds"), _text(o.find("t")))
            continue
        kept += 1
        for n in on_:
            worst = max(worst, min(length(sub(_p(n), _p(x))) for x in tn_ if x.get("Element", "6") == n.get("Element", "6")))
    lines.append(f"  nicknames: {len(ours)} ours, {kept} with the same inner atoms and bonds; "
                 f"largest inner atom move {worst:.3f}")
    return lines + ch.report()


# --------------------------------------------------------------------------- phase 2: schemes


def _text(t) -> str:
    return "".join(x.text or "" for x in t.iter("s"))


def _anchors(page) -> dict[str, tuple]:
    """id -> a position that survives ChemDraw's renumbering: a fragment's first node, a <t>'s
    p (a formula fragment ChemDraw made from a caption: its label's p), a <graphic>'s Head3D.
    Stand-in for ids, which ChemDraw renumbers on every re-save."""
    out = {}
    species = {i for st in page.iter("step") for k in ("ReactionStepReactants", "ReactionStepProducts")
               for i in (st.get(k) or "").split()}
    for f in page.findall("fragment"):
        n = f.find("n")
        t = n.find("t") if n is not None else None
        # one labeled node that no step has as a species: a formula ChemDraw made from a caption
        # (K2CO3; H3O+ became hydronium), keyed by its label's p, as our <t>
        caption_like = n is not None and len(f.findall("n")) == 1 and t is not None and f.get("id") not in species
        out[f.get("id")] = ("t", *_floats(t.get("p"))[:2]) if caption_like else ("f", *_floats(n.get("p"))[:2])
    for t in page.findall("t"):
        out[t.get("id")] = ("t", *_floats(t.get("p"))[:2])
    for g in page.iter("graphic"):
        if g.get("SupersededBy"):
            out[g.get("id")] = ("a", *_floats(g.get("BoundingBox"))[:2])
    return out


def _step_sets(page) -> list[tuple]:
    """Each <scheme>'s steps as (kind, rounded position) sets per citation attribute."""
    anchors = _anchors(page)

    def key(i):
        a = anchors.get(i)
        return (a[0], round(a[1], 1), round(a[2], 1)) if a else ("?", i)
    attrs = ("ReactionStepReactants", "ReactionStepProducts", "ReactionStepArrows", "ReactionStepObjectsAboveArrow",
             "ReactionStepObjectsBelowArrow")
    return [frozenset(tuple(tuple(sorted(key(i) for i in (st.get(k) or "").split())) for k in attrs)
                      for st in sch.findall("step")) for sch in page.iter("scheme")]


def scheme_diff(ours_root, theirs_root) -> list[str]:
    """What ChemDraw's re-save did to the scheme objects (module docstring, items 1-3). Ids are
    renumbered on re-save, so steps are compared by what they cite, each object by position."""
    lines, ch = [], Changes()
    op, tp = ours_root.find("page"), theirs_root.find("page")
    ours_steps, their_steps = _step_sets(op), _step_sets(tp)
    same = sorted(map(sorted, ours_steps)) == sorted(map(sorted, their_steps))
    lines.append("  steps and schemes: " + ("the same objects in each step and scheme" if same else
                 f"CHANGED: ours {[len(x) for x in ours_steps]} steps per scheme, re-save {[len(x) for x in their_steps]}"))
    if not same:
        theirs = {st for sch in their_steps for st in sch}
        for sch in ours_steps:
            for st in sch:
                if st not in theirs:
                    lines.append(f"  step not in the re-save: {st}")
    for tag, key in (("arrow", "Head3D"), ("graphic", "BoundingBox")):
        for el in op.iter(tag):
            if el.get("GraphicType") == "Oval":
                continue  # ChemDraw rewrites an oval's BoundingBox: symbol_diff pairs them by Center3D
            t = next((x for x in tp.iter(tag) if x.get(key) and max(
                abs(u - v) for u, v in zip(_floats(x.get(key)), _floats(el.get(key)))) < PAIR_TOL), None)
            if t is None:
                ch.counts[(tag, f"<{tag}>")]["removed"] += 1
                ch.examples.setdefault((tag, f"<{tag}>"), f"ours at {el.get(key)}")
                continue
            ch.compare(tag, el, t, skip=("id", "Z", "SupersededBy"))
    their_anchors = _anchors(tp)
    for el in op.findall("t"):
        p = _floats(el.get("p"))
        tid = next((i for i, a in their_anchors.items() if a[0] == "t" and max(abs(a[1] - p[0]), abs(a[2] - p[1])) < 0.02),
                   None)
        t = next((x for x in tp.iter() if x.get("id") == tid), None) if tid else None
        if t is None:
            ch.counts[("caption", "<t>")]["removed"] += 1
            ch.examples.setdefault(("caption", "<t>"), f"ours {_text(el)!r}")
            continue
        if t.tag == "fragment":  # ChemDraw read the caption as a formula (K2CO3)
            ch.counts[("caption", "became a formula fragment")]["changed"] += 1
            ch.examples.setdefault(("caption", "became a formula fragment"), f"ours {_text(el)!r}")
            continue
        ch.compare("caption", el, t, skip=("id", "Z"))
        os_, ts_ = el.findall("s"), t.findall("s")
        if [(x.text, x.get("face"), x.get("font")) for x in os_] != [(x.text, x.get("face"), x.get("font")) for x in ts_]:
            ch.counts[("caption s", "runs")]["changed"] += 1
            ch.examples.setdefault(("caption s", "runs"),
                                   f"ours {[(x.text, x.get('face')) for x in os_]} -> {[(x.text, x.get('face')) for x in ts_]}")
    return lines + ch.report()


def caption_ink_diff(ref: str, resaved: Path) -> list[str]:
    """Caption glyphs in the re-saved SVG against the original reference's (item 4)."""
    a, b = mr.load(ref), mr.load(ref, resaved)
    out = []
    worst = 0.0
    for cap in a.captions:
        ga = mr._caption_glyphs(a, cap)
        if not ga:
            continue
        p = (cap["p"][0] - a.offset[0], cap["p"][1] - a.offset[1])
        cb = min(b.captions, key=lambda c: length(sub((c["p"][0] - b.offset[0], c["p"][1] - b.offset[1]), p)))
        gb = mr._caption_glyphs(b, cb)
        if len(ga) != len(gb):
            out.append(f"  {cap['text']!r}: {len(ga)} glyphs, re-save {len(gb)}")
            continue
        shift = max(length(sub((x.x - a.offset[0], x.base - a.offset[1]), (y.x - b.offset[0], y.base - b.offset[1])))
                    for x, y in zip(sorted(ga, key=lambda g: (round(g.base), g.x)), sorted(gb, key=lambda g: (round(g.base), g.x))))
        sizes = sorted({round(g.size, 2) for g in ga}), sorted({round(g.size, 2) for g in gb})
        worst = max(worst, shift)
        la, lb = mr.text_lines(ga), mr.text_lines(gb)
        pitch = ""
        if len(la) > 1:
            pa = [round(y["base"] - x["base"], 3) for x, y in zip(la, la[1:])]
            pb = [round(y["base"] - x["base"], 3) for x, y in zip(lb, lb[1:])]
            pitch = f", line pitch {pa} -> {pb}"
        if shift > NUM_TOL or pitch or sizes[0] != sizes[1]:
            out.append(f"  {cap['text']!r}: largest glyph shift {shift:.3f}{pitch}"
                       + (f", sizes {sizes[0]} -> {sizes[1]}" if sizes[0] != sizes[1] else ""))
    return out + [f"  largest caption glyph shift {worst:.3f}"]


# --------------------------------------------------------------------------- phase 3: mechanisms

DOTS = ("LonePair", "Electron")


def _curves(root) -> list:
    return list(root.iter("curve"))


def _dot_boxes(root) -> list[tuple]:
    """(SymbolType, node p, BoundingBox) per dot symbol; nodes by position (ids are renumbered)."""
    nodes = {n.get("id"): _p(n) for n in root.iter("n")}
    return [(g.get("SymbolType"), nodes.get(g.find("represent").get("object")) if g.find("represent") is not None else None,
             _floats(g.get("BoundingBox"))) for g in root.iter("graphic") if g.get("SymbolType") in DOTS]


def mech_diff(ours_root, theirs_root) -> list[str]:
    """Items 5 and 6: curves paired by first anchor, dot symbols by atom and first point."""
    lines, ch = [], Changes()
    tc = _curves(theirs_root)
    oc = _curves(ours_root)
    lines.append(f"  curves: {len(oc)} ours, {len(tc)} re-saved")
    for c in oc:
        p = _floats(c.get("CurvePoints"))
        t = next((x for x in tc if max(abs(u - v) for u, v in zip(_floats(x.get("CurvePoints"))[2:4], p[2:4])) < PAIR_TOL),
                 None)
        if t is None:
            ch.counts[("curve", "<curve>")]["removed"] += 1
            ch.examples.setdefault(("curve", "<curve>"), f"ours from {p[2:4]}")
            continue
        ch.compare("curve", c, t, skip=("id", "Z"))
    od, td = _dot_boxes(ours_root), _dot_boxes(theirs_root)
    lines.append(f"  dot symbols: {len(od)} ours, {len(td)} re-saved")
    for kind, node, box in od:
        m = [b for k, n, b in td if k == kind and n and node and length(sub(n, node)) < PAIR_TOL
             and max(abs(u - v) for u, v in zip(b[:2], box[:2])) < PAIR_TOL]
        if not m:
            ch.counts[(kind, "moved or removed")]["changed"] += 1
            ch.examples.setdefault((kind, "moved or removed"), f"ours {box} on the atom at {node}")
            continue
        d = max(abs(u - v) for u, v in zip(m[0], box))
        if d > NUM_TOL:
            ch.counts[(kind, "BoundingBox")]["changed"] += 1
            ch.delta[(kind, "BoundingBox")] = max(ch.delta[(kind, "BoundingBox")], d)
            ch.examples.setdefault((kind, "BoundingBox"), f"ours {box} -> ChemDraw {m[0]}")
    reps = collections.Counter(r.get("attribute") for g in theirs_root.iter("graphic") if g.get("SymbolType") in DOTS
                               for r in g.findall("represent"))
    lines.append(f"  re-saved dot <represent> attributes: {dict(reps)}")
    return lines + ch.report()


def _head_record(tip, p_end, tangent, pts) -> dict:
    n = (-tangent[1], tangent[0])
    loc = [(dot(sub(q, tip), tangent), dot(sub(q, tip), n)) for q in pts]
    left, right = max(v for _, v in loc), -min(v for _, v in loc)
    return {"tip": tip, "end": p_end, "side": "-" if abs(left - right) < 0.05 else "left" if left > right else "right",
            "local": sorted((round(u, 2), round(v, 2)) for u, v in loc)}


def _drawn_heads(doc) -> list[dict]:
    """ChemDraw's drawn curves (measure_refs), CDXML frame: stroke end, tip, barb side, head shape."""
    out = []
    off = doc.offset
    for c in mr.measure(doc)["curved_arrows"]:
        if "head" not in c:
            continue
        end = c["svg_path_pts"][3]
        tip = add(end, c["head"]["tip_minus_stroke_end"])
        tangent = unit(sub(tip, end)) if c["head"]["n_vertices"] == 4 else None
        h = c["head"]
        side = "-" if abs(h["extent_left"] - h["extent_right"]) < 0.05 else (
            "left" if h["extent_left"] > h["extent_right"] else "right")
        out.append({"tip": sub(tip, off), "end": sub(end, off), "side": side,
                    "local": sorted((round(u, 2), round(v, 2)) for u, v in h["pts_local"]), "tangent": tangent})
    return out


def _our_heads(spec_path: Path) -> list[dict]:
    """Our MechArrows in the exported file's frame (to_svg's offset, as spec_cdxml writes it)."""
    spec = load_spec(spec_path)
    lay = build_scheme(spec, spec_path.parent, spec_path.stem)
    _, off = to_svg(lay.root)
    out = []
    for m in lay.curved:
        t = unit(sub(m.p3, m.c2 if length(sub(m.p3, m.c2)) > 1e-6 else m.c1))
        out.append(_head_record(add(m.tip, off), add(m.p3, off), t, [add(q, off) for q in m.head]))
    return out


def mech_ink_diff(name: str, resaved: Path, ours_root) -> list[str]:
    """Item 7."""
    b = mr.load(name, resaved)
    drawn = _drawn_heads(b)
    if name in PHASE3_REFS:
        want, what = _drawn_heads(mr.load(name)), "original reference"
    else:
        want, what = _our_heads(mr.REPO / next(s for s in PHASE3_SPECS if Path(s).stem == name)), "our layout"
    lines = [f"  {len(want)} heads ({what}), {len(drawn)} drawn in the re-save"]
    worst = collections.defaultdict(float)
    for w in want:
        d = min(drawn, key=lambda x: length(sub(x["end"], w["end"])), default=None)
        if d is None:
            continue
        worst["stroke end"] = max(worst["stroke end"], length(sub(d["end"], w["end"])))
        worst["head tip"] = max(worst["head tip"], length(sub(d["tip"], w["tip"])))
        if d["side"] != w["side"]:
            lines.append(f"  barb side at {tuple(round(v, 2) for v in w['end'])}: {what} {w['side']}, re-save {d['side']}")
        if len(d["local"]) == len(w["local"]):
            worst["head shape (local)"] = max(worst["head shape (local)"], max(
                length(sub(p, q)) for p, q in zip(d["local"], w["local"])))
        else:
            lines.append(f"  head at {tuple(round(v, 2) for v in w['end'])}: {len(w['local'])} vertices, re-save {len(d['local'])}")
    declared = [tuple(bb[:2]) for _, _, bb in _dot_boxes(ours_root)]
    declared += [tuple(bb[2:]) for k, _, bb in _dot_boxes(ours_root) if k == "LonePair"]
    circles = [sub(c.center, b.offset) for c in b.shapes if c.kind == "circle"]
    lines.append(f"  dots: {len(declared)} declared, {len(circles)} drawn")
    for c in declared:
        if circles:
            worst["dot center (declared -> drawn)"] = max(worst["dot center (declared -> drawn)"],
                                                          min(length(sub(c, q)) for q in circles))
    return lines + [f"  {k:45} {v:.3f}" for k, v in worst.items()]


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


def symbol_diff(ours_root, theirs_root) -> list[str]:
    """Phase 4b (module docstring, item 10): the symbols per (SymbolType, bound or free) in ours
    and the re-save, the largest move of a box's first point (paired by type), and the free
    <t> texts (charge digits)."""
    def syms(root):
        out = collections.defaultdict(list)
        for g in root.iter("graphic"):
            if g.get("GraphicType") == "Symbol":
                out[(g.get("SymbolType"), "bound" if g.find("represent") is not None else "free")].append(
                    _floats(g.get("BoundingBox"))[:2])
        return out
    a, b = syms(ours_root), syms(theirs_root)
    lines = []
    for k in sorted(set(a) | set(b)):
        move = max((min(length(sub(p, q)) for q in b[k]) for p in a[k]), default=0.0) if b[k] else float("nan")
        lines.append(f"  {k[0]:12} {k[1]:5} ours {len(a[k]):3}  re-save {len(b[k]):3}  largest move {move:.3f}")

    def texts(root):
        return sorted("".join(x.text or "" for x in t.iter("s")) for t in root.find("page").findall("t"))
    def ovals(root):  # ChemDraw rewrites an oval's BoundingBox and OvalType: pair by Center3D
        return [_floats(g.get("Center3D"))[:2] for g in root.iter("graphic") if g.get("GraphicType") == "Oval"]
    oa, ob = ovals(ours_root), ovals(theirs_root)
    if oa or ob:
        move = max((min(length(sub(p, q)) for q in ob) for p in oa), default=0.0) if ob else float("nan")
        lines.append(f"  {'Oval':12} {'':5} ours {len(oa):3}  re-save {len(ob):3}  largest move {move:.3f}")
    ta, tb = texts(ours_root), texts(theirs_root)
    lines.append(f"  free texts: {'the same' if ta == tb else f'{ta} -> {tb}'}")
    return lines


def _crossings(root) -> dict[frozenset, tuple]:
    """{the two bonds' node-position pairs: front bond's node positions} for every crossing: the
    CrossingBonds ChemDraw records, else (our files: none recorded) every pair of bonds whose axes
    cross, the front one by Z."""
    pos = {n.get("id"): tuple(round(v, 1) for v in _floats(n.get("p"))) for n in root.iter("n") if n.get("p")}
    bonds = {b.get("id"): b for b in root.iter("b") if b.getparent().getparent().tag != "n"}

    def key(b):
        return frozenset((pos[b.get("B")], pos[b.get("E")]))
    pairs = {frozenset((b.get("id"), c)) for b in bonds.values() for c in (b.get("CrossingBonds") or "").split()
             if c in bonds}
    if not pairs:
        from acs_chem_svg.geometry import segments_cross
        items = list(bonds.values())
        for i, b in enumerate(items):
            for c in items[i + 1:]:
                ends = [pos[b.get("B")], pos[b.get("E")], pos[c.get("B")], pos[c.get("E")]]
                if len(set(ends)) == 4 and segments_cross(*ends):
                    pairs.add(frozenset((b.get("id"), c.get("id"))))
    out = {}
    for p in pairs:
        b, c = (bonds[i] for i in p)
        front = b if int(b.get("Z")) > int(c.get("Z")) else c
        out[frozenset((key(b), key(c)))] = key(front)
    return out


def crossing_diff(ours_root, theirs_root) -> list[str]:
    ours, theirs = _crossings(ours_root), _crossings(theirs_root)
    lines = [f"  {len(ours)} crossings ours, {len(theirs)} recorded by ChemDraw"]
    for k in sorted(set(ours) | set(theirs), key=lambda k: sorted(sorted(x) for x in k)):
        o, t = ours.get(k), theirs.get(k)
        if o != t:
            lines.append(f"  crossing {sorted(sorted(x) for x in k)}: front ours {sorted(o) if o else None}, "
                         f"ChemDraw {sorted(t) if t else None}")
    if len(lines) == 1:
        lines.append("  every crossing has the same front bond")
    return lines


def diff(phase: str = "1") -> int:
    found = False
    out = out_dir(phase)
    resaved = out / "resaved"
    for path in sorted(out.glob("*.cdxml")):
        saved = resaved / path.name
        if not saved.exists():
            print(f"{path.stem}: no re-save in {resaved.relative_to(mr.REPO)}")
            continue
        found = True
        if phase in PHASES:
            ours, theirs = (etree.parse(str(f), PARSER).getroot() for f in (path, saved))
            print(f"{path.stem}: scheme objects (ours -> ChemDraw's re-save)")
            for line in scheme_diff(ours, theirs):
                print(line)
            print(f"{path.stem}: species attributes")
            for line in attribute_diff(ours, theirs):
                print(line)
            if phase == "4b":
                print(f"{path.stem}: symbols and charge digits")
                for line in symbol_diff(ours, theirs):
                    print(line)
            if phase == "4a":
                print(f"{path.stem}: nicknames")
                for line in nickname_diff(ours, theirs):
                    print(line)
            if phase == "5":
                print(f"{path.stem}: crossings (front bond by Z; ours -> ChemDraw's re-save)")
                for line in crossing_diff(ours, theirs):
                    print(line)
            if phase == "3":
                print(f"{path.stem}: curves and dots")
                for line in mech_diff(ours, theirs):
                    print(line)
                if (resaved / f"{path.stem}.svg").exists():
                    print(f"{path.stem}: drawn heads and dots (re-saved SVG)")
                    for line in mech_ink_diff(path.stem, resaved, ours):
                        print(line)
            if path.stem in PHASE2_REFS and (resaved / f"{path.stem}.svg").exists():
                print(f"{path.stem}: caption ink (original reference -> re-save)")
                for line in caption_ink_diff(path.stem, resaved):
                    print(line)
            continue
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
    ap.add_argument("--lewis-draft", action="store_true", help="write the Lewis reference draft (phase 4b)")
    ap.add_argument("--phase", choices=("1", "2", "3", "4a", "4b", "5"), default="1",
                    help="1: M1 species; 2: M2 schemes; 3: M3 mechanisms; 4a: abbreviations, resonance, wrap")
    args = ap.parse_args(argv)
    if args.lewis_draft:
        for p in lewis_draft():
            print(p.relative_to(mr.REPO))
        return 0
    if args.diff:
        return diff(args.phase)
    for p in write(args.phase):
        print(p.relative_to(mr.REPO))
    resaved = (out_dir(args.phase) / "resaved").relative_to(mr.REPO)
    print(f"Open each in ChemDraw; save as CDXML and SVG into {resaved}/, same names.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
