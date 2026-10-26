# From a request to a mechanism (worked example)

The order in SKILL.md "From a request", applied to: *"Draw the full arrow-pushing mechanism for
the nitration of benzene with HNO3 and H2SO4."* Commands are `uv run --project <repo> python
<repo>/scripts/<script>`.

## 1. Steps in words

1. HNO3's OH takes H+ from H2SO4 (lp of O4 to H10; the O9–H10 bond goes to O9).
2. The protonated acid loses water (lp of O3 forms N=O; the N–O4 bond goes to O4).
3. Benzene's C11=C12 bond attacks N; one N=O goes to O1. The sigma complex is a cation.
4. Resonance: the + moves 12 → 14 → 16 (bond 13=14 to 12–13, then 15=16 to 14–15).
5. HSO4⁻ takes H17 from C11; the C11–H17 bond forms C11=C16.

## 2. Map numbers, once

| Atom | Map | Atom | Map |
|---|---|---|---|
| HNO3 O (=O) | 1 | H2SO4 O, S, O, OH | 5, 6, 7, 8 |
| HNO3 N | 2 | H2SO4 O–H oxygen | 9 |
| HNO3 O⁻ | 3 | the H that moves (step 1) | 10 |
| HNO3 OH oxygen | 4 | benzene C1…C6 | 11…16 |
| the ring H that leaves | 17 | | |

Map every atom an arrow touches and every H that moves (`[H:10]`, `[H:17]`). Map the other
atoms too when they are few: an unmapped atom in a product can't be matched, and the message is
harder to read. Keep the numbers for the whole scheme.

## 3. Reactants and arrows

```yaml
preset: acs1996
species:
  nitric:   {smiles: "[O:1]=[N+:2]([O-:3])[OH:4]"}
  sulfuric: {smiles: "[O:5]=[S:6](=[O:7])([OH:8])[O:9][H:10]"}
  benzene:  {smiles: "[CH:11]1=[CH:12][CH:13]=[CH:14][CH:15]=[CH:16]1"}
steps:
  - from: [nitric, sulfuric]
    to: [nitric_h, bisulfate]
    arrows:
      - {from: "lp(4)", to: "between(4,10)"}
      - {from: "bond(9-10)", to: "atom(9)"}
  ...
```

Benzene needs no `[H:17]`: the H appears in the sigma complex, where it is written as an atom.

## 4. Rough products, then paste the suggestions

Write each product as any SMILES with the right map numbers, then
`validate.py --suggest spec.yaml` prints what the arrows give:

```
step 1 (nitric, sulfuric ->): the arrows give (paste to replace the products):
  suggested: nitric_h: {smiles: "[O:1]=[N+:2]([O-:3])[OH+:4][H:10]"}
  suggested: bisulfate: {smiles: "[O:5]=[S:6](=[O:7])([OH:8])[O-:9]"}
```

Paste them, step by step (each step's products are the next step's reactants). For a
resonance step the suggestion is the next contributor.

## 5. Check atoms

`render.py spec.yaml --out DIR --maps` writes `spec_maps.png` with each map number beside its
atom. Use it to check that each arrow starts where you meant, and to decide orientations.

## 6. Orient, render, inspect

```yaml
  nitric:  {smiles: "...", orientation: {flip_h: true}}    # OH toward sulfuric acid's H
  benzene: {smiles: "...", orientation: {rotate: 180}}     # C11=C12 toward the nitronium
steps:
  - new_row: true                                          # row B: a new sequence
    from: [benzene, nitronium]                             # nitronium again: drawn as a copy
    ...
  - from: [sigma_c, bisulfate]                             # bisulfate again
    to: [nitrobenzene, sulfuric]                           # sulfuric again
```

Orient only the first species of each sequence and separate reactants; products follow their
reactants. The approved result is `examples/m4_nitration_x1.0.yaml` (written there with `_2`
copies; reusing the names gives the same drawing). At acs1996 it takes three rows (504 pt);
that is expected.

## Common mistakes (M5 trials)

- **Changing the preset to make rows fit.** Use the preset the user asked for (default
  acs1996). Rows wrap at 504 pt; a long mechanism takes three rows, and the products alone on
  the last row are fine. Raise `layout.width` only if the user wants a wider figure.
- **Brute-force override loops.** `render.py` prints each arrow's `bow`; try the other side
  once, then `tail_side`, then turn a species. If two or three tries fail, deliver and report
  the arrow, rather than looping.
- **Unmapped atoms in a product that an arrow changed.** Map every atom an arrow touches, on
  both sides; the mismatch message lists the map numbers found on one side only.
- **A carbanion or other carbon with lone pairs drawn without a label** (C≡N⁻): its dots sit
  on the bonds. `show_carbons: [n]` labels it; the validator says so.
- **The H that moves, drawn too early.** Write `[H:17]` only where the H is an arrow end or
  changes partner; on benzene it would be drawn as a stray H.
- **A partner on the far side.** If an arrow arcs over a whole molecule, swap the reactants
  (`from: [b, a]`) or flip the species (`orientation: {flip_h: true}`) so the reacting atoms
  face each other, and check with `--maps`.
- **Arrows that pass a by-product.** When a step continues a row, the previous step's
  by-products stay between the species and the next arrow, and a reactant the step adds is
  drawn above its arrow. If the curved arrows of that step have to reach across a by-product
  (radical propagation: HBr between t-Bu• and Br2), start the step on its own row
  (`new_row: true`, `from: [tbu, br2]`).
- **Pasting a suggestion that looks wrong.** The suggestions are what the arrows give. If
  they show odd charges (`[OH+:1]`, `[CH3-:2]`), the arrows are wrong, not the products: a
  new bond is `to: between(a,b)`, while `to: atom(b)` makes a lone pair on b.
