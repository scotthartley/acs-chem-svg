# Scheme spec schema

Validated by `src/acs_chem_svg/spec.py` (pydantic). Unknown keys are errors. M2 covers one row of
species, plus signs, reaction arrows with conditions, and compound numbers. A spec without `steps`
is an M1 spec: its species are drawn side by side. M3 adds curved arrows, lone pairs and
radicals; M4 adds resonance brackets, several sequences and wrapping into rows (below).

## M2: linear scheme

```yaml
preset: acs1996            # acs1996 (default) | acs1996_x0.8 | acs1996_x0.6
species:
  ester:    {smiles: "COC(=O)c1ccc(O)cc1", abbreviate: [OMe], orientation: {rotate: 120}}
  bnbr:     {smiles: "BrCc1ccccc1", orientation: {rotate: 0}}
  ether:    {smiles: "COC(=O)c1ccc(OCc2ccccc2)cc1", abbreviate: [OMe, OBn]}
  alcohol:  {smiles: "OCc1ccc(OCc2ccccc2)cc1", abbreviate: [OBn]}
  benzene:  {smiles: "c1ccccc1", align: false}
steps:
  - {from: [ester, bnbr], to: ether, above: K2CO3, below: "acetone, reflux"}
  - {to: alcohol, above: LiAlH4, below: ["THF, 0 °C", "2 h"]}      # from: the previous `to`
  - {to: benzene, stages: ["*n*-BuLi, THF, −78 °C", "H3O+"]}       # "1. …" above, "2. …" below
numbers: {ester: 1, bnbr: 2, ether: 3, alcohol: "4a"}
align: auto                # or false: every species on its own layout
```

| Key | Type | Meaning |
|---|---|---|
| `steps[].from` | name or list | Reactants. Required in the first step; later steps default to the previous step's `to`. It starts with the previous step's first drawn product (then any more of its drawn products, in order). Further names are added reactants (`[alcohol, X]` draws "+ X" before the arrow; above the arrow in a mechanism step). A product the next step doesn't take should be omitted (`omit`), not left drawn. |
| `steps[].to` | name or list | Products (a list gets plus signs). |
| `steps[].omit` | list of names | Products of `to` that the electron bookkeeping uses but the drawing leaves out (water lost, the conjugate acid/base of a reagent drawn over the arrow, the acid/base regenerated in the last step). At least one product stays drawn; no `numbers:` on them. Write a loss on the arrow yourself if it matters (`above: "−H2O"`; `below: "−AcO−"` under a → arrow with a reagent drawn above it). Never below an equilibrium arrow: text there reads as the reverse direction (the validator notes it). The conjugate of an `over` reagent needs no text: reagents over the arrow aren't balanced. A later step that names an omitted product adds it as a reactant (drawn there for the first time). |
| `steps[].above` / `below` | string or list | Conditions, one string per line, centered on the arrow. |
| `steps[].stages` | list of 1–2 strings | Numbered stages, left-justified: "1. …" above, "2. …" below. More than 2 is unmeasured. Exclusive with above/below. |
| `steps[].arrow_type` | `forward` (default) \| `equilibrium` \| `equilibrium_products` \| `equilibrium_reactants` | The step's arrow: → , balanced ⇌, or unbalanced ⇌ with the long shaft pointing to the favored side (products: right; reactants: left). Equilibrium conditions sit above and below the outer shafts, centered on the arrow (acs-style.md §3.6). Not on resonance steps. Works with `stages` and curved `arrows` (a reversible proton transfer). |
| `arrow_head` | `standard` \| `large` | Head of every reaction, equilibrium and ↔ arrow. Default: the preset's (`large` at ×1.0, `standard` at ×0.8 and ×0.6; acs-style.md §3.1). `standard` at ×1.0 matches ChemDraw's default head. Curved arrows are unaffected. |
| `numbers.<name>` | int or string | Bold compound number under the species, on one shared baseline. |
| `species.<name>.abbreviate` | list | Collapse these terminal groups to labels: OMe, OEt, OBn, OPh, OAc, SMe, NHAc, CO2Me, CO2Et, CO2H, CF3, CN, NO2, Ph, Bn, Et, Ac, TMS, TES, TBS, TIPS. A group whose bond leaves to the right is set as ChemDraw does: right-justified, in reversed text (MeO, MeO₂C, O₂N, F₃C, AcHN; Ph, Bn, Et, Ac and the silyl groups keep their text). |
| `species.<name>.align` | bool (default true) | A product takes the orientation of the species it is drawn from: atoms shared by atom map (≥ 2 common maps), else by maximum common substructure, keep their coordinates. `false` = own layout. Molfile species keep their layout unless `align_to` is set. |
| `species.<name>.align_to` | name | Align to this species instead (it must be drawn earlier). |
| `align` | `auto` \| `false` | `false` turns alignment off for the whole scheme. |
| `labels` | map number → symbol | Draw a generic atom (X, Y, R, Nu, R′ as `R'`): `labels: {9: X}` draws atom 9 as X in every species. Write a real element with the same bonds, H count and electrons as its stand-in in the SMILES (`[Cl:9]` for a radical X·, `[H:6][Cl:9]` for HX); the electron bookkeeping uses the stand-in. Each map number must be on some atom. |

An aligned species can't also set `orientation` (set `align: false` to orient it yourself). Each
species is drawn once; give a repeated compound a second name.

### Conditions text

| Input | Drawn |
|---|---|
| `K2CO3`, `Pd(PPh3)4`, `CH2Cl2`, `C12H25` | digits after a letter or `)` are subscripts |
| `H3O+`, `Cl−` (U+2212) | a trailing `+` / `−` after a formula is a superscript |
| `*n*-BuLi` | `*…*` is italic |
| `'SN\2'`, `'\*'` (YAML single quotes) | `\` escapes the next character (set plain) |
| `0 °C`, `−78 °C`, `(4:1)`, `1,4-dioxane`, `12 h`, `n-BuLi` | unchanged (ASCII `-` is a hyphen) |

### Layout (automatic, acs-style.md §3)

Structures sit on the arrow axis by their ring centers (or ink centers if a species has no
ring); arrows grow to fit their widest line; numbers share a baseline. SVG groups (editor
layers): `species-<name>`, `plus-<i>`, `arrow-step-<n>`, `conditions-step-<n>` (with `-above`,
`-below`), `numbers` (with `number-<name>`). `scripts/validate.py` checks the schema, the
chemistry, and overlaps.

## M3: mechanisms

```yaml
preset: acs1996
species:
  hydroxide: {smiles: "[H][O-:1]"}
  mebr:      {smiles: "[CH3:2][Br:3]", show_carbons: terminal}
  methanol:  {smiles: "[OH:1][CH3:2]", show_carbons: terminal}
  bromide:   {smiles: "[Br-:3]", show_lone_pairs: [3]}
steps:
  - from: [hydroxide, mebr]
    to: [methanol, bromide]
    arrows:
      - {type: two_electron, from: "lp(1)", to: "between(1,2)"}
      - {from: "bond(2-3)", to: "atom(3)"}
overrides:
  arrows:
    - {step: 1, index: 0, bow: left, curvature: 0.6, tail_side: top}
```

| Key | Type | Meaning |
|---|---|---|
| `steps[].arrows[].type` | `two_electron` (default) \| `fishhook` | Moves 2 or 1 electrons. |
| `steps[].arrows[].from` | endpoint | `lp(n)`: a lone pair of atom n (not for fishhooks). `bond(a-b)`, `bond(a=b)`, `bond(a#b)`: the bond, whose order must be as written. `electron(n)`: n's unpaired electron. |
| `steps[].arrows[].to` | endpoint | `atom(n)`: the electrons become nonbonding on n. `between(a,b)`: a new bond, or one order more. |
| `species.<name>.show_lone_pairs` | list of map numbers | Also draw these atoms' lone pairs (all of them). |
| `overrides.arrows[]` | `{step, index, …}` | Adjust arrow `index` (0-based) of step `step` (1-based): `bow: left\|right` (the side of the tail → head direction it bulges to), `curvature` (control offset ÷ chord; default 0.5, short loops 1.3), `tail_side: top\|right\|bottom\|left` (which lone pair or radical the tail leaves). No coordinates. |

Map numbers are atom identities. They must be unique among a step's reactants (and among its
products) and keep the same atom across steps. H atoms that move must be atoms in the SMILES
(`[H:4]`). A SMILES carbon with one H written as an atom gets its other H's as unmapped atoms
(`[H:4][CH2:3]` draws three C–H bonds), so write only the H that moves. An arrow can't start or end inside an abbreviation.

### Electron bookkeeping (`electrons.py`)

All of a step's arrows apply at once to the reactants taken together. The result must equal
the declared products, compared as canonical SMILES with atom maps, charges and radicals, with
H atoms made implicit and stereo dropped. These are errors: a source without the electrons
(no lone pair, no such bond or the wrong order written, no radical); a fishhook from a lone
pair; a bond left with an odd number of electrons; more than 8 electrons around a second-row
atom or more than 2 around H; a duplicate map number; an endpoint that doesn't resolve.
Passing proves internal consistency, not that the mechanism is right.

On a mismatch the message names map numbers found on one side only, and `validate.py` prints
what the arrows do give, one paste-ready line per molecule, named after the declared product
it shares the most atoms with (`?`: none does, so add a species):

```
step 1 (hydroxide, mebr ->): the arrows give (paste to replace the products):
  suggested: methanol: {smiles: "[OH:1][CH3:2]"}
  suggested: bromine: {smiles: "[Br-:3]"}
```

`validate.py --suggest` prints these for every step with arrows (and each resonance list:
the next contributor), passing or not. The SMILES keep map numbers, charges, radicals, mapped
H atoms and the Kekulé form the arrows made; unmapped H are implicit. So the usual way to write
products is: write the reactants and arrows, put rough products, and paste the suggestions.

`render.py --maps` writes `<stem>_maps.png`: the drawing with each mapped atom's number beside
it in red, to check which atom is which before choosing arrows, `orientation` and overrides. It
never changes the SVG. `render.py` also prints the `bow` each curved arrow got, so an override
starts from the default (`bow: left` when it already bows left changes nothing).

### Drawing (automatic, acs-style.md §4)

- **Dots:** lone pairs are drawn on lp tails, on arrow targets (the atoms in `atom(n)` and
  `between(a,b)`), on atoms that gain a pair (Br⁻ after SN2), and on `show_lone_pairs`. On a
  labeled atom each pair takes a free side of the element (top, right, bottom, left; sides with
  a bond within 45°, the H run or an NH₃⁺ charge are taken). The pair an arrow leaves faces the
  arrow's target. Arrows never move the other pairs; they work around them. On an unlabeled atom
  pairs and radicals go on largest-gap bisectors. Radicals are always drawn.
- **Curved arrows:** a cubic stroke and a filled head (half head for fishhooks). The tail starts
  beyond its lone pair or radical dot, or beside the bond's midpoint. The head stops short of the
  target's label and its dots (above a lone pair, never between its dots), atom center, or (for `between` on an existing bond)
  the bond. A two-electron `between` across species points at the atom that isn't the source.
  A fishhook `between` points at the midpoint, where the bond forms; when a plus sign sits
  there, just past it, and the fishhooks that meet there share a side. Arrows bow away from
  their species, or upward over the row. Chords shorter than a bond curl into loops. Each arrow
  tries both bow sides, every lone pair, each solution for its end, a flatter curve and a tip
  stopped farther out, and the step keeps the combination that crosses the least. Two arrows from one bond go to opposite
  sides.
- **Groups:** `species-<name>-lone-pairs`, `species-<name>-radicals` (SVG `<circle>`s),
  `curved-arrows-step-<n>` with `curved-arrow-step-<n>-<index>` (a stroked path and a head).
- **Checks** (`validate.py`): curved arrows crossing bonds, labels, charges, dots, other curved
  arrows, reaction arrows, conditions or numbers (not plus signs); heads and tails closer than
  `margin_width` to a label; dots overlapping bonds, labels, charges or each other, or a pair
  straddling a bond.

## M4: resonance, several sequences, wrapping

```yaml
layout: {width: 504}                       # pt; default: the preset's scheme.wrap_width (504 = 7 in)
steps:
  - {from: [nitric, sulfuric], over: [sulfuric], to: [nitric_h, bisulfate], omit: [bisulfate], arrows: [...]}
  - {from: nitric_h, to: [nitronium, water], omit: [water], above: "−H2O", arrows: [...]}
  - {new_row: true, from: [benzene, nitronium], to: sigma_a, arrows: [...]}  # nitronium again
  - resonance: [sigma_a, sigma_b, sigma_c]                        # first = the previous step's first product
    arrows:                                                       # one list per contributor but the last
      - [{from: "bond(13=14)", to: "between(12,13)"}]
      - [{from: "bond(15=16)", to: "between(14,15)"}]
  - {from: [sigma_c, bisulfate], to: [nitrobenzene, sulfuric], omit: [sulfuric], arrows: [...]}
```

| Key | Type | Meaning |
|---|---|---|
| `steps[].new_row` | bool | Start a new sequence on a new row: `from` is required and the chain starts again. |
| `steps[].over` | list of names | Reactants drawn above the arrow of a step with curved `arrows` (their arrows reach them there). Default: a later mechanism step's added reactants. Use it on the first step or a `new_row` step for an acid or base that only gives or takes a proton (`from: [nitric, sulfuric], over: [sulfuric]`). Names must be `from` names the step adds (not kept products); at least one reactant stays inline. |
| repeated names | — | A species named again later in the scheme (in any step) is drawn again: occurrence k is a copy with the id `name~k`, aligned to its own predecessor like any product. `numbers:` apply to the first occurrence. To give a copy its own options, declare it: `species: {"nitronium~2": {smiles: ..., align: false}}`. A copy that is a product drops the original's `orientation` (it aligns). |
| `steps[].resonance` | list of ≥ 2 names | Contributors in square brackets, joined by double-headed arrows. The first must be the previous step's first product (or the step is the first); that step's other products are drawn after the closing bracket ("] + HX"). The next step continues from the last contributor. No conditions. |
| `steps[].arrows` (resonance) | list of lists | List k is drawn on contributor k and must give contributor k + 1. Omit `arrows` for none. `overrides.arrows[].index` counts across the lists in order. |
| `layout.width` | pt | Wrap width. Rows break only between steps; a unit that can't fit alone is an error. |

Bookkeeping: contributors keep the first one's atom maps, connectivity (bonds ignoring order,
H counts) and total charge; each arrow list must turn contributor k into k + 1.

Drawing (acs-style.md §3.4, §4.3):
- **Brackets** enclose every contributor's ink, including its dots and curved arrows; the ↔
  arrows are `resonance.arrow_length` long, `resonance.arrow_gap` from the contributors (their
  curved arrows count). The arrow into the bracket stops `arrow_gap` before it. The next
  step's curved arrows may reach into the bracket.
- **Added reactants in mechanisms:** in a step with curved arrows, the reactants it adds
  (`from` names after the kept products, e.g. the base HSO₄⁻), or its `over` names, are drawn
  above its arrow, not as "+ X", side by side with "+" between them; the arrow lengthens to
  fit them, and conditions above the arrow stack above them. Steps without curved arrows keep
  "+ X". Omitted products are never drawn (nor their "+").
- **Rows:** units (the opening reactants; each step's arrow, conditions and products, with a
  following bracket and the next step's added reactants) fill each row up to `layout.width`.
  A continuation row starts with the reaction arrow at the left margin, except after a step
  with curved arrows: its arrow stays with its reactants, and the next row starts with its
  products; rows are left-aligned,
  `scheme.row_gap` apart (full ink, conditions to numbers). Each row has its own number
  baseline and its own vertical rule (ring centers on the axis when every structure in the row
  has a ring). A step's curved arrows stay with its reactants, even when its arrow starts the
  next row.
- **Groups:** `row-<r>` (only with several rows), `bracket-step-<n>` (`-left`, `-right`),
  `resonance-arrow-step-<n>-<k>`, `numbers-row-<r>` (several rows; else `numbers`).

## M1: single-species spec

```yaml
preset: acs1996
species:
  aspirin:                 # name -> SVG group ids species-aspirin, -bonds, -labels, -charges
    smiles: "CC(=O)Oc1ccccc1C(=O)O"
  benzene:
    smiles: "[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1"
    kekule: ["1=2", "3=4", "5=6"]          # which aromatic bonds are double (atom-map numbers)
  alanine:
    molfile: templates/m1/m1_set_x1.0__3.mol
    orientation: {rotate: 30, flip_h: false, flip_v: false}
```

| Key | Type | Meaning |
|---|---|---|
| `preset` | string | Preset file in `presets/`. All lengths and sizes come from it. |
| `species.<name>.smiles` | string | SMILES. Kekulé bonds are kept as written. Aromatic input is kekulized by RDKit unless `kekule` is given. Explicit `[H]` atoms are kept and drawn. Stereo becomes wedges via RDKit. |
| `species.<name>.molfile` | path | V2000/V3000 molfile, relative to the spec file or the repo root. Its layout is kept and scaled so the median bond is the preset bond length (a 1.5 Å layout maps exactly). Its wedges are kept. An atom alias naming an abbreviation (OMe on an O) draws that label. |
| `species.<name>.kekule` | list of `"m=n"` | Atom-map pairs that are double bonds. Every other aromatic bond becomes single. |
| `species.<name>.orientation.rotate` | number \| `auto` | Degrees, counter-clockwise on screen, applied to the raw CoordGen layout. SMILES default `auto` (bonds at ±30°, rings with vertical sides, linear molecules horizontal). Molfile default 0. |
| `species.<name>.orientation.flip_h` / `flip_v` | bool | Mirror. Wedges are re-derived so the stereo stays correct. |
| `species.<name>.show_carbons` | `none` (default), `terminal`, or a list of atom ids | Explicit labels on terminal carbons (CH₃, H₃C, =CH₂). Internal CH₂/CH raise an error until measured (acs-style.md §6 Open 8). Atom ids are atom-map numbers if set, else 1-based atom indices. A molfile atom alias on a carbon (`A` line: CH3, H3C, …) also shows its label. |
| `charges.circled` | bool (default false) | Draw every ±1 sign in a ring (⊕ / ⊖). Unmeasured ring size (`charge.circle_radius`). |
| `species.<name>.circled_charges` | bool | Overrides `charges.circled` for this species. |
| `species.<name>.charge_labels` | map: atom-map number → text | Draw this charge text instead of the formal charge: one sign before or after optional digits (`+3`, `2+`, `-`; `-` draws the minus bar). Display only: electron bookkeeping keeps the SMILES charge; a neutral atom draws it too (an oxidation state). Never circled. |
| `species.<name>.lewis` | `true` \| `{electrons: shared \| formal_charge \| oxidation_state}` | Draw a Lewis structure (below). `true` = `shared`. SMILES only, specs without `steps` only; not with `kekule`, `abbreviate`, `show_carbons` or `show_lone_pairs` (the last two are implied). `orientation.rotate` must be a multiple of 90. |

Give exactly one of `smiles` or `molfile`. Without `steps`, several species are placed left to
right with `scheme.plus_gap` between their ink boxes.

### Lewis structures (`lewis`)

```yaml
species:
  ccl4:     {smiles: "ClC(Cl)(Cl)Cl", lewis: true}
  co2:      {smiles: "O=C=O", lewis: {electrons: formal_charge}}
  methanol: {smiles: "CO", lewis: {electrons: oxidation_state}}
```

- **Atoms:** every H becomes an atom and every atom is labeled with its bare element symbol
  (C included); every lone pair and radical is drawn. Charges come from the SMILES (±2 or
  more draws "2+"); `charges.circled`, `circled_charges` and `charge_labels` apply as for any
  species.
- **Layout:** a square lattice, every bond `lewis.bond_length` long whatever its order. The
  root is the heavy-atom tree center; its heavy neighbors go right, left, up, down. Each other
  atom continues its largest heavy branch straight and puts the others up, then down; H atoms
  fill what is left (CH₃OH: H–C–O–H with C's H's above and below; H₂O and CO₂ linear; CCl₄ a
  cross). The search keeps non-bonded atoms off neighboring lattice points where it can; where
  it can't (acetic acid's CH₃ beside C=O), the sides facing each other get no lone pairs.
- **Bonds as electrons**, between the two element ink boxes:
  - `shared`: one column per bond order (O : : C, N : : : N), each column two dots across the
    bond, centered in the gap and at least 2 × `lone_pair.gap` from each label, so a shared
    pair sits visibly farther out than a lone pair.
  - `formal_charge`: each atom gets `order` electrons; `oxidation_state`: the more
    electronegative atom (Pauling, `chem.ELECTRONEGATIVITY`) gets all of them, equal atoms
    split. Assigned electrons pack outward from the owner's label like lone pairs (a pair, then
    one dot on the axis if odd). The labels don't move between the three modes.
- **Limits (errors):** rings; atoms with more than 4 neighbors (PCl₅); trees that don't fit
  on the lattice (neopentane); salts written as one SMILES; Lewis species in schemes.
  `validate` reports a triple bond assigned to one atom in `oxidation_state` (its 3 columns
  cross the middle of the bond) and wide `labels:` symbols that crowd the shared pairs.
- **SVG groups:** `species-<name>-labels`, `-lone-pairs`, `-radicals`, `-charges`, plus
  `-bonding-pairs` (shared) and `-bond-electrons` (assigned). There is no `-bonds` group.

## Drawn conventions (automatic)

- **Labels:** heteroatoms and explicit H atoms are labeled; carbons are not, except isolated
  ones and those chosen with `show_carbons` or a molfile alias. The H count goes on the side away from the bonds (H₂N– vs –NH₂);
  a shown CH₃ that would run into another label or bond is set the other way (H₃C), and `validate`
  reports a label that still overlaps.
- **Charges:** charges of ±1 are drawn as bars (⊕/⊖ geometry from the preset), in a ring with
  `charges.circled`; ±2 or more as digits then the sign ("2+", "2−", digits at `charge.digit_scale` ×
  label size). On a labeled atom they sit upper right, or above if the H is on the right. On an
  unlabeled atom they go on the bisector of the largest gap. A sign beside drawn lone pairs (or a
  circled / "2+" sign) moves outward until it clears the dots and bonds by 1 pt.
- **Not yet supported:** condensed labels other than the abbreviations above (CH₃Cl as text),
  aromatic circles.
