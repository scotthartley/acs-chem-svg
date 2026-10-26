---
name: acs-chem-svg
description: Draw chemical structures, reaction schemes, and arrow-pushing mechanisms in ACS (ChemDraw ACS Document 1996) style as editable SVG, with a PNG preview Claude inspects before delivering. Use when the user asks to draw a molecule, a reaction scheme, a mechanism with curved arrows, or resonance structures.
---

# acs-chem-svg

> **Status: M5. Single structures, schemes, arrow-pushing mechanisms and resonance work** (all
> three presets): rings, bonds, wedges, labels with H counts and charges, abbreviations (OMe,
> OBn, CO₂Me, NO₂, Ph, …), reaction arrows with conditions, plus signs, bold compound numbers,
> products aligned to their reactants, curved arrows (2e⁻ and fishhook), lone pairs and
> radicals with electron bookkeeping, resonance brackets with ↔ arrows, several sequences, and
> rows wrapped to a width (504 pt by default). Product suggestions, species reused by name, and
> a map-number overlay help turn a one-sentence request into a spec.

`<repo>` is this skill's directory: the target of the project's `.claude/skills/acs-chem-svg` symlink.
Run every command as `uv run --project <repo> python <repo>/scripts/<script> ...`.

## Workflow

1. **Write a spec** (full schema: `references/scheme-schema.md`):
   ```yaml
   preset: acs1996          # or acs1996_x0.8, acs1996_x0.6
   species:
     ester:   {smiles: "COC(=O)c1ccc(O)cc1", abbreviate: [OMe], orientation: {rotate: 120}}
     bnbr:    {smiles: "BrCc1ccccc1"}
     ether:   {smiles: "COC(=O)c1ccc(OCc2ccccc2)cc1", abbreviate: [OMe, OBn]}
     alcohol: {smiles: "OCc1ccc(OCc2ccccc2)cc1", abbreviate: [OBn]}
   steps:
     - {from: [ester, bnbr], to: ether, above: K2CO3, below: "acetone, reflux"}
     - {to: alcohol, above: LiAlH4, below: "THF, 0 °C"}
   numbers: {ester: 1, bnbr: 2, ether: 3, alcohol: 4}
   ```
   - Conditions text: digits after a letter or `)` become subscripts (K2CO3, Pd(PPh3)4);
     a trailing `+`/`−` is a superscript (H3O+); `*n*-BuLi` is italic; `\` escapes.
   - Use only the conditions the user gave. A reagent named in the request (BnBr) goes over
     the arrow unless the user wants it drawn and numbered.
   - Several lines: `above: [line1, line2]`. Two-stage arrows: `stages: ["*n*-BuLi, THF, −78 °C", "H3O+"]`.
   - Products are aligned to the species they come from. Orient only the first species (and
     separate reactants) with `orientation: {rotate: …, flip_h: …}`; `align: false` opts out.
   - Without `steps`, the species are drawn side by side (single structures).
   - `kekule: ["1=2", ...]` (atom maps) chooses the Kekulé form; `show_carbons: terminal`
     labels CH₃; a molfile from `templates/` keeps its layout.
2. **Validate:** `scripts/validate.py spec.yaml`. Fix the spec until it prints OK (schema,
   RDKit chemistry, abbreviations, electron bookkeeping, overlaps, curved arrows crossing
   bonds or labels, dots on bonds, conditions wider than their arrow).
3. **Render:** `scripts/render.py spec.yaml --out DIR [--crops]` writes `spec.svg` and
   `spec.png` (4 px/pt; `--px-per-pt 10` for close inspection) and prints the drawing size.
   `--crops` adds 16 px/pt crops of every curved arrow, charged or dotted label, ↔ arrow and
   bracket corner. `--affinity` writes for Affinity (separate arrowhead shapes, bold/italic face names).
4. **Inspect** the PNG (and crops): orientation of the first species, label sides, double-bond
   sides, wedges, charges, conditions centered on their arrows, numbers on one line; for
   mechanisms, that each arrow leaves its lone pair or bond and points at its target.
5. **Adjust** the spec only (orientation, abbreviations, line breaks, `align`, `overrides`),
   never the SVG. Stop after about 5 visual iterations and report what is left.
6. **Deliver** `.svg` (editable text, named groups; verified in Inkscape; `--affinity` for
   Affinity, which drops marker arrowheads), `.yaml`, `.png`.
   Mention it if stderr warned that the fallback font was used for metrics.

## From a request (mechanisms)

Worked example and common mistakes: `references/workflow.md`.
1. Write the steps in words (who attacks what, what leaves).
2. Assign map numbers once, in a table, to every atom that moves or is an arrow end; keep them
   for the whole scheme. Map the H atoms that move as atoms (`[H:17]`).
3. Write the reactants and each step's arrows.
4. Write rough products, run `validate.py --suggest`, and paste the `suggested:` lines.
5. Render with `--maps` from the first render (map numbers beside the atoms, preview only):
   it shows which end of each species faces its partner.
6. Orient, render, inspect. Reuse a species by name when it appears again (`nitronium` in
   both rows): each later occurrence is drawn as its own copy. Keep the user's preset
   (default acs1996): a long mechanism wrapping to three rows is expected.

## Mechanisms (curved arrows)

Map every atom an arrow touches (`[O-:1]`) and keep the map numbers across the scheme. H atoms
that move must be written as atoms (`[H:4][CH2:3]…`). Arrows are drawn on the step's reactants:
```yaml
species:
  hydroxide: {smiles: "[H][O-:1]"}
  mebr:      {smiles: "[CH3:2][Br:3]", show_carbons: terminal}
  methanol:  {smiles: "[OH:1][CH3:2]", show_carbons: terminal}
  bromide:   {smiles: "[Br-:3]"}
steps:
  - from: [hydroxide, mebr]
    to: [methanol, bromide]
    arrows:
      - {from: "lp(1)", to: "between(1,2)"}            # type: two_electron (default) | fishhook
      - {from: "bond(2-3)", to: "atom(3)"}
overrides:
  arrows: [{step: 1, index: 0, bow: left, curvature: 0.6, tail_side: top}]
```
- Endpoints: from `lp(n)`, `bond(a-b)` / `bond(a=b)` / `bond(a#b)`, `electron(n)` (a radical);
  to `atom(n)` or `between(a,b)` (a new or stronger bond). Full vocabulary and drawing rules:
  `references/scheme-schema.md` (Mechanisms).
- Every step's arrows must turn its reactants into its products exactly (charges, radicals,
  bond orders, H moved). Read the validator's message; it names the arrow and the atom, and
  on a product mismatch it prints paste-ready `suggested:` products.
- Generic atoms (X·, R): a stand-in element in the SMILES plus `labels: {9: X}` (map → symbol).
- Lone pairs are drawn on arrow-tail and arrow-target atoms and on atoms that gain a pair; add
  others with `show_lone_pairs: [n]` per species. Radicals are always drawn.
- `render.py` prints each arrow's default `bow`. If an arrow crosses something, try
  `overrides` (the other `bow`, `curvature`, `tail_side`), or turn a
  species (`orientation`) so the atoms face each other. Order the reactants so partners sit
  side by side (`from: [mebr, hydroxide]` swaps them).
- The bookkeeping proves internal consistency, **not** mechanistic correctness: say so, and ask
  the user to check research-level mechanisms.

## Resonance, sequences, rows

```yaml
layout: {width: 504}                                  # optional; pt
steps:
  - {from: [benzene, nitronium], to: sigma_a, arrows: [...]}
  - resonance: [sigma_a, sigma_b, sigma_c]            # first = the previous step's first product
    arrows: [[{from: "bond(13=14)", to: "between(12,13)"}], [...]]   # list k: contributor k -> k+1
  - {from: [sigma_c, bisulfate], to: [nitrobenzene, sulfuric], arrows: [...]}
  - {new_row: true, from: [a, b], to: c}             # another sequence, on a new row
```
- Contributors keep their atom maps and connectivity; give them the same orientation by
  aligning (the default) and orient only the first.
- `from` may keep just the first products of the previous step; the others stay as by-products.
  In a step with curved arrows, the reactants it adds (a base) are drawn above its arrow.
- Rows break only between steps; a continuation row starts with an arrow, or, after a
  mechanism step's arrow, with its products. If a step plus a
  bracket can't fit, raise `layout.width`. Full rules: `references/scheme-schema.md` (M4).

## Known limits (raise errors rather than guess)

- An abbreviation whose bond leaves to the right (MeO–, BnO–) is unmeasured: turn the
  structure so the bond comes from the left or below, or don't abbreviate that group.
- At most 2 stages per arrow.
- Explicit CH₂/CH labels on internal carbons, charges beyond ±1, carbenes.

## Rules

- Never edit `src/`, `scripts/`, `presets/`, or `references/` from a use session.
- Problems and workarounds go to `<repo>/issues/inbox.md` (create it if missing), appended, with:
  - the date,
  - a spec excerpt,
  - the symptom,
  - the workaround.
- Validation proves internal consistency, not mechanistic correctness. Say so for
  research-level mechanisms.
