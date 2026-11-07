"""Measure ChemDraw reference renders into numbers for references/style-spec.md.

Atom positions, bond orders, stereo flags, and object types come from the
.cdxml saved next to each export (exact, no clustering heuristics). Drawn
geometry comes from the .svg, mapped into the CDXML frame by the path
transform (svg = cdxml + (e, f)). SVG user units are checked against the PDF
page size, so every number is in pt.

ChemDraw draws every bond line, wedge, hash, charge sign, and arrow shaft as a
filled polygon; lone pairs as filled circles; curved arrows as stroked cubics.

Usage:
    uv run python scripts/measure_refs.py                 # all references
    uv run python scripts/measure_refs.py m1_set_x1.0     # one file
    uv run python scripts/measure_refs.py --json out.json
    uv run python scripts/measure_refs.py --scheme        # scheme parameters per preset (§3)
    uv run python scripts/measure_refs.py --mech          # mechanism parameters per preset (§4)
    uv run python scripts/measure_refs.py --resonance     # resonance arrows and brackets per preset (§4.3)
    uv run python scripts/measure_refs.py --wrap          # multi-row scheme parameters per preset (§3.4)
    uv run python scripts/measure_refs.py --abbrev        # abbreviation labels per preset (§3.5)
    uv run python scripts/measure_refs.py --equilibrium   # equilibrium arrows, head variants (§3.1, §3.6)
    uv run python scripts/measure_refs.py --cdxml-symbols # charge-symbol <graphic> boxes (§13, CDXML export)
    uv run python scripts/measure_refs.py --cdxml-scheme  # CDXML captions, arrows and steps (§13.1)
    uv run python scripts/measure_refs.py --charges-lewis # circled signs, symbols, two-electron radicals (§2.2, §13.4)
    uv run python scripts/measure_refs.py --chair         # chair templates, axial / equatorial bonds (§14.1)
    uv run python scripts/measure_refs.py --crossing      # bond crossings: front by Z, rear cuts (§14.2)
    uv run python scripts/measure_refs.py --arc           # arc arrows: circle, sweep, head on the chord (§15.1)
    uv run python scripts/measure_refs.py --nh            # where an NH's H goes: beside, above, below (§2)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from chemdraw_refs import arcs, cdxml, chair, clip, core, labels, measurements, mech, scheme
from chemdraw_refs.core import all_names, load
from chemdraw_refs.measurements import measure
from chemdraw_refs.clip import clip_rows, fit_label_clip, summarize
from chemdraw_refs.scheme import (
    EQUILIBRIUM_REFS, RESONANCE_REFS, SCHEME_REFS, WRAP_REFS, equilibrium_summary, resonance_summary,
    scheme_summary, wrap_summary,
)
from chemdraw_refs.mech import MECH_REFS, mech_summary
from chemdraw_refs.labels import ABBREV_REFS, NH_REFS, abbrev_summary, nh_summary
from chemdraw_refs.cdxml import (
    CDXML_SCHEME_REFS, CDXML_SYMBOL_REFS, CHARGES_LEWIS_REFS, cdxml_scheme, cdxml_symbols, charges_lewis,
)
from chemdraw_refs.chair import CHAIR_REFS, chair_summary, crossing_summary
from chemdraw_refs.arcs import ARC_EQ_REFS, ARC_REFS, arc_eq_summary, arc_summary

# Re-export every module's names, helpers included: compare.py, cdxml_gate.py, cdxml_to_molfile.py
# and the tests use `mr.<name>` (mr._caption_glyphs, mr._dot_sets, ...).
for _m in (core, measurements, clip, scheme, mech, labels, cdxml, chair, arcs):
    globals().update({k: v for k, v in vars(_m).items() if not k.startswith('__') and k not in globals()})



# References --fit-clip leaves out unless named (style-spec.md §2.1): bond ends the model is not
# fitted to. nh_x*: bonds 6-9° off the axes end up to 0.6 pt short, a ChemDraw error (Open 28);
# metals_x*: the bond straight down from a centered Pd (Open 25; a refit with them still leaves it
# 0.16-0.25 off and moves the other parameters); chair_x*: a rear bond cut at a crossing (§14.2).
CLIP_FIT_EXCLUDE = ("nh_x", "metals_x", "chair_x")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="*")
    ap.add_argument("--json", type=Path)
    ap.add_argument("--summary", action="store_true", help="print the key-parameter table only")
    ap.add_argument("--fit-clip", action="store_true", help="fit the bond-to-label clip model (§2.1; CLIP_FIT_EXCLUDE unless named)")
    ap.add_argument("--scheme", action="store_true", help="scheme parameters per preset (§3)")
    ap.add_argument("--mech", action="store_true", help="mechanism parameters per preset (§4)")
    ap.add_argument("--wrap", action="store_true", help="multi-row scheme parameters per preset (§3.4)")
    ap.add_argument("--resonance", action="store_true", help="resonance arrow and bracket parameters per preset (§4.3)")
    ap.add_argument("--abbrev", action="store_true", help="abbreviation label parameters per preset (§3.5)")
    ap.add_argument("--equilibrium", action="store_true",
                    help="equilibrium arrows and arrowhead variants per preset (§3.1, §3.6)")
    ap.add_argument("--cdxml-symbols", action="store_true", help="charge-symbol <graphic> boxes per reference (§13)")
    ap.add_argument("--cdxml-scheme", action="store_true", help="CDXML captions, arrows and steps per reference (§13.1)")
    ap.add_argument("--charges-lewis", action="store_true",
                    help="circled signs, symbols and two-electron radicals in charges_lewis_x* (§2.2, §13.4)")
    ap.add_argument("--chair", action="store_true", help="chair templates, axial and equatorial bonds (§14.1)")
    ap.add_argument("--crossing", action="store_true", help="bond crossings: front by Z, rear-bond cuts (§14.2)")
    ap.add_argument("--arc", action="store_true", help="arc arrows: circle, sweep, head on the chord (§15.1)")
    ap.add_argument("--arc-eq", action="store_true", help="arc equilibrium arrows: shafts, butts, notches, tips (§15.1)")
    ap.add_argument("--nh", action="store_true", help="where an NH's H goes: beside, above, below (§2)")
    args = ap.parse_args(argv)
    if args.nh:
        print(json.dumps({n: nh_summary(n) for n in args.names or NH_REFS}, indent=1, ensure_ascii=False))
        return 0
    if args.arc_eq:
        refs = {k: n for k, n in ARC_EQ_REFS.items() if not args.names or n in args.names}
        print(json.dumps({k: arc_eq_summary(n) for k, n in refs.items()}, indent=1, ensure_ascii=False))
        return 0
    if args.arc:
        refs = {k: n for k, n in ARC_REFS.items() if not args.names or n in args.names}
        print(json.dumps({k: arc_summary(n) for k, n in refs.items()}, indent=1, ensure_ascii=False))
        return 0
    if args.chair:
        print(json.dumps({n: chair_summary(n) for n in args.names or CHAIR_REFS}, indent=1, ensure_ascii=False))
        return 0
    if args.crossing:
        print(json.dumps({n: crossing_summary(n) for n in args.names or CHAIR_REFS}, indent=1, ensure_ascii=False))
        return 0
    if args.charges_lewis:
        names = args.names or CHARGES_LEWIS_REFS
        print(json.dumps({n: charges_lewis(n) for n in names}, indent=1, ensure_ascii=False))
        return 0
    if args.cdxml_scheme:
        names = args.names or CDXML_SCHEME_REFS
        print(json.dumps({n: cdxml_scheme(n) for n in names}, indent=1, ensure_ascii=False))
        return 0
    if args.cdxml_symbols:
        names = args.names or CDXML_SYMBOL_REFS
        print(json.dumps({n: cdxml_symbols(n) for n in names}, indent=1, ensure_ascii=False))
        return 0
    if args.equilibrium:
        refs = {k: n for k, n in EQUILIBRIUM_REFS.items() if not args.names or n in args.names}
        print(json.dumps({k: equilibrium_summary(n) for k, n in refs.items()}, indent=1, ensure_ascii=False))
        return 0
    if args.abbrev:
        names = args.names or ABBREV_REFS
        print(json.dumps({n: abbrev_summary(n) for n in names}, indent=1, ensure_ascii=False))
        return 0
    if args.resonance:
        refs = {k: n for k, n in RESONANCE_REFS.items() if not args.names or n in args.names}
        print(json.dumps({k: resonance_summary(n) for k, n in refs.items()}, indent=1, ensure_ascii=False))
        return 0
    if args.wrap:
        refs = {k: n for k, n in WRAP_REFS.items() if not args.names or n in args.names}
        print(json.dumps({k: wrap_summary(n) for k, n in refs.items()}, indent=1, ensure_ascii=False))
        return 0
    if args.mech:
        refs = {k: n for k, n in MECH_REFS.items() if not args.names or n in args.names}
        print(json.dumps({k: mech_summary(n) for k, n in refs.items()}, indent=1, ensure_ascii=False))
        return 0
    if args.scheme:
        names = args.names or [n for ns in SCHEME_REFS.values() for n in ns]
        print(json.dumps(scheme_summary({n: measure(load(n)) for n in names}), indent=1, ensure_ascii=False))
        return 0
    if args.fit_clip:
        fit = fit_label_clip(clip_rows(args.names or [n for n in all_names() if not n.startswith(CLIP_FIT_EXCLUDE)]))
        for r in fit.pop("rows"):
            print(f"{r['file']:19} {r['kind']:5} {r['element']:3} dir {r['dir_deg']:7} end {r['end']:6.3f} "
                  f"residual {r['residual']:+.3f}{'  !' if abs(r['residual']) > 0.1 else ''}")
        print(json.dumps(fit, indent=1))
        return 0
    results = {n: measure(load(n)) for n in (args.names or all_names())}
    if args.summary:
        rows = {n: summarize(m) for n, m in results.items()}
        keys = list(dict.fromkeys(k for r in rows.values() for k in r))
        names = list(rows)
        print("parameter".ljust(34) + "".join(n[:18].rjust(19) for n in names))
        for k in keys:
            print(k.ljust(34) + "".join(str(rows[n].get(k, "")).rjust(19) for n in names))
        return 0
    text = json.dumps(results, indent=1, default=str)
    if args.json:
        args.json.write_text(text)
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
