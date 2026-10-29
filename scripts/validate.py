"""Validate a spec (design.md §6): schema, chemistry (RDKit sanitize, abbreviations, alignment),
electron bookkeeping for curved arrows (electrons.check_spec), then scheme and mechanism
geometry (geometry_checks). Prints one line per problem; exit 1 if any.

When the arrows of a step don't give its declared products, it also prints what they do give,
one `suggested: <name>: {smiles: "..."}` line per molecule, ready to paste into `species:`
(name "?": a molecule no declared product shares atoms with). --suggest prints these for every
step with arrows, passing or not.

In a mechanism of two or more steps with curved arrows, it also prints a `note:` (not a problem)
for each extra product of a mechanism step that is drawn and never used again: a chemist omits it
(omit: [...]) and, if the loss matters, writes it on the arrow ("−H2O"). It also notes a loss
("−X") written below an equilibrium arrow, where it reads as the reverse direction.

Usage: uv run python scripts/validate.py [--suggest] spec.yaml [...]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pydantic import ValidationError

from acs_chem_svg.electrons import StepCheck, check_spec
from acs_chem_svg.geometry_checks import (check_labels, check_lewis, check_line_end_dots, check_mechanism,
                                          check_scheme)
from acs_chem_svg.layout import build_scheme, m1_drawers
from acs_chem_svg.spec import ResonanceStep, Spec, load_spec


def validate(path: Path) -> list[str]:
    return checked(path)[0]


def checked(path: Path) -> tuple[list[str], list[StepCheck]]:
    """(problems, the electron checks of the steps with arrows)."""
    try:
        spec = load_spec(path)
    except ValidationError as e:
        return [f"schema: {err['loc']}: {err['msg']}" for err in e.errors()], []
    except (OSError, ValueError) as e:
        return [f"schema: {e}"], []
    try:
        if not spec.steps:
            return [f"geometry: {p}" for d in m1_drawers(spec, path.parent)
                    for p in check_labels(d.mol.name, d.draw(d.mol.name), d.mol) + check_lewis(d.mol.name, d)
                    + check_line_end_dots(d.mol.name, d)], []
        checks = check_spec(spec, path.parent)
        electrons = [f"electrons: {p}" for c in checks for p in c.problems]
    except (ValueError, NotImplementedError, RuntimeError, FileNotFoundError) as e:
        return [f"chemistry: {e}"], []
    try:
        lay = build_scheme(spec, path.parent, path.stem)
    except (ValueError, NotImplementedError, RuntimeError, FileNotFoundError) as e:
        return electrons + [f"chemistry/layout: {e}"], checks
    return electrons + [f"geometry: {p}" for p in check_scheme(lay) + check_mechanism(lay)], checks


def notes(spec: Spec) -> list[str]:
    """A loss ("−X") written below an equilibrium arrow; drawn products after the first of a
    mechanism step that no later step takes, in a scheme with two or more mechanism steps (a
    one-step mechanism's leaving group stays)."""
    steps = spec.steps or []
    out = []
    for i, st in enumerate(steps, 1):
        if isinstance(st, ResonanceStep) or st.arrow_type == "forward":
            continue
        for ln in st.lines()[1]:
            if ln.startswith(("−", "-")):
                out.append(f"note: step {i}: {ln!r} below an equilibrium arrow reads as lost in the reverse "
                           "direction; write it above, or leave it out if it is the conjugate of a reagent "
                           "drawn over the arrow")
    mech = [i for i, st in enumerate(steps, 1) if not isinstance(st, ResonanceStep) and st.arrows]
    if len(mech) < 2:
        return out
    for i in mech:
        st = steps[i - 1]
        later = {n for s in steps[i:] for n in (s.resonance if isinstance(s, ResonanceStep) else s.sources)}
        side = "below" if spec.over_names(i) and st.arrow_type == "forward" else "above"
        for n in st.drawn_products[1:]:
            if n not in later:
                out.append(f"note: step {i}: {n} is drawn and never used again; omit it (omit: [{n}]) and, if the "
                           f"loss matters, write it on the arrow ({side}: \"−<formula>\", e.g. \"−H2O\")")
    return out


def suggestions(checks: list[StepCheck], every: bool = False) -> list[str]:
    """Paste-ready product lines for the steps whose arrows don't give the declared products
    (every step with arrows if `every`)."""
    out = []
    for c in checks:
        if not c.suggested or not (every or c.result != c.expected):
            continue
        head = f"step {c.step} ({', '.join(c.sources)} ->)"
        out.append(f"{head}: the arrows give{'' if c.result == c.expected else ' (paste them if they are the products you meant; if not, fix the arrows)'}:")
        out += [f"  suggested: {n}: {{smiles: \"{s}\"}}" for n, s in c.suggested]
        if any(n == "?" for n, _ in c.suggested):
            out.append("  (?: no declared product shares atoms with it: add it as a species and to the step's `to`)")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("specs", nargs="+", type=Path)
    ap.add_argument("--suggest", action="store_true",
                    help="print the products the arrows give for every step with arrows, as species lines")
    args = ap.parse_args(argv)
    bad = 0
    for p in args.specs:
        problems, checks = checked(p)
        bad += bool(problems)
        print(f"{p}: {'OK' if not problems else f'{len(problems)} problem(s)'}")
        for msg in problems:
            print(f"  {msg}")
        if not problems:
            for msg in notes(load_spec(p)):
                print(f"  {msg}")
        for line in suggestions(checks, args.suggest):
            print(f"  {line}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
