---
name: acs-chem-svg
description: Draw chemical structures, Lewis structures, reaction schemes, and arrow-pushing mechanisms in ACS (ChemDraw ACS Document 1996) style as editable SVG, with a PNG preview Claude inspects before delivering. Use when the user asks to draw a molecule, a Lewis structure, a reaction scheme, a mechanism with curved arrows, or resonance structures.
---

# acs-chem-svg

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
   - Several lines: `above: [line1, line2]`. Multistep arrows: `stages: ["*n*-BuLi, THF, −78 °C", "H3O+"]`.
   - Reversible steps: `arrow_type: equilibrium` (⇌), `equilibrium_products` / `equilibrium_reactants`
     (unbalanced). `arrow_head: standard` (top level) gives ChemDraw's smaller default head at ×1.0.
   - Products are aligned to the species they come from. Orient only the first species of each
     sequence (and separate reactants) with `orientation: {rotate: …, flip_h: …}`; `align: false` opts out.
   - Without `steps`, the species are drawn side by side (single structures).
   - `kekule: ["1=2", ...]` (atom maps) chooses the Kekulé form; `show_carbons: terminal`
     labels CH₃; a molfile from `templates/` keeps its layout.
   - Lewis structures (every atom labeled, bonds as dot pairs): `lewis: true`, `lewis: {bonds: lines}`
     (bond lines), or `lewis: {electrons: formal_charge | oxidation_state}`; see `references/scheme-schema.md`.
   - Chairs: `chair: true` draws a single six-membered sp³ ring as a chair, axial/equatorial from the SMILES
     stereo (`chair: {flip: true}`, `{axial: [maps]}`, `{start: map}`); see `references/scheme-schema.md`.
   - Charges: `charges: {circled: true}` (spec) or `circled_charges` (species) circles ±1 signs;
     `charge_labels: {1: "+3"}` (species, by atom map) draws a text, display only.
2. **Validate:** `scripts/validate.py spec.yaml`. Fix the spec until it prints OK (schema,
   RDKit chemistry, abbreviations, electron bookkeeping, overlaps, curved arrows crossing
   bonds or labels, dots on bonds, conditions wider than their arrow).
3. **Render:** `scripts/render.py spec.yaml --out DIR [--crops] [--maps]` writes `spec.svg` and
   `spec.png` (4 px/pt; `--px-per-pt 10` for close inspection) and prints the drawing size.
   `--crops` adds 16 px/pt crops of curved arrows, charged or dotted labels, ↔ arrows,
   equilibrium arrow ends and bracket corners.
4. **Inspect** the PNG (and crops): orientation of the first species, label sides, double-bond
   sides, wedges, charges, conditions centered on their arrows, numbers on one line; for
   mechanisms, that each arrow leaves its lone pair or bond and points at its target.
5. **Adjust** the spec only (orientation, abbreviations, line breaks, `align`, `overrides`),
   never the SVG. At most 3 tries per curved arrow and 5 render rounds overall; then deliver
   and report what is left.
6. **Deliver** `.svg` (editable text, named groups; `--affinity` writes for Affinity: separate
   arrowhead shapes, bold/italic face names), `.yaml`, `.png`. Mention a fallback-font warning.
   For ChemDraw, `--cdxml` also writes `.cdxml` for every spec kind (M1–M4, Lewis structures,
   abbreviations, charges). Report its stderr warnings (atom maps are not written).

## Mechanisms (curved arrows)

Worked example and common mistakes: `references/workflow.md`.
1. Write the steps in words (who attacks what, what leaves).
2. Assign map numbers once, in a table, to every atom an arrow touches (`[O-:1]`); keep them
   for the whole scheme. Map each H that moves as an atom (`[H:4][CH2:3]…`). Writing one of a
   carbon's H's as an atom draws all of them as bonds, so write only the H that moves.
3. Write the reactants and each step's arrows (drawn on the step's reactants).
4. Write rough products, run `validate.py --suggest`, and paste the `suggested:` lines. Every
   step's arrows must turn its reactants into its products exactly (charges, radicals, bond
   orders, H moved); the validator's message names the arrow and the atom.
5. **Account for every product:**
   - An acid/base that only gives or takes a proton goes `over:` the arrow; its conjugate
     (AcOH from AcO−) and the acid/base regenerated last are `omit:`ted, without text.
   - Omit other products the next step doesn't use; note a loss that matters on the arrow
     (`above: "−H2O"`, or `below:` a → with a structure over it).
   - Never write text below ⇌ (it reads as the reverse direction).
6. Render with `--maps` (map numbers beside the atoms, preview only): it shows which end of
   each species faces its partner.
7. Orient, render, inspect. `render.py` prints each arrow's default `bow`; if an arrow crosses
   something, try `overrides` (the other `bow`, then `tail_side`, `curvature`), or turn a species
   (`orientation`) or swap the reactants (`from: [mebr, hydroxide]`) so the atoms face each other.
   Keep the user's preset (default acs1996); three rows is fine.

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
- Generic atoms (X·, R): a stand-in element in the SMILES plus `labels: {9: X}` (map → symbol).
- Lone pairs the arrows don't involve: `show_lone_pairs: [n]` per species.

## Resonance, sequences, rows

```yaml
layout: {width: 504}                                  # optional; pt
steps:
  - {from: [benzene, nitronium], to: sigma_a, arrows: [...]}
  - resonance: [sigma_a, sigma_b, sigma_c]            # first = the previous step's first product
    arrows: [[{from: "bond(13=14)", to: "between(12,13)"}], [...]]   # list k: contributor k -> k+1
  - {from: [sigma_c, bisulfate], to: [nitrobenzene, sulfuric], omit: [sulfuric], arrows: [...]}
  - {new_row: true, from: [a, b], to: c}             # another sequence, on a new row
```
- Contributors keep their atom maps and connectivity; they align to the first (orient only it).
- Reuse a species by name when it appears again: each later occurrence is its own copy.
- `omit:` products stay in the bookkeeping, out of the drawing; a later step that names one
  adds it. A mechanism step's added reactants (a base) and its `over:` sit above its arrow.
- Rows break only between steps; if a step plus a bracket can't fit, raise `layout.width`.
  Full rules: `references/scheme-schema.md` (Resonance, sequences, rows).

## Known limits (raise errors rather than guess)

- Explicit CH₂/CH labels on internal carbons, carbenes.
- Lewis structures: rings, more than 4 neighbors (PCl₅), Lewis species in schemes.
- Chairs: fused or spiro rings (decalins, steroids), sp² ring atoms, boats.
- Condensed labels other than the listed abbreviations (CH₃Cl as text), aromatic circles.

## Rules

- Never edit `src/`, `scripts/`, `presets/`, or `references/` from a use session.
- Problems and workarounds go to `<repo>/issues/inbox.md` (create it if missing), appended, with:
  - the date,
  - a spec excerpt,
  - the symptom,
  - the workaround.
- Validation proves internal consistency, not mechanistic correctness. Say so for
  research-level mechanisms.
