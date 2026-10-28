# acs-chem-svg: design brief

The original design brief for development sessions (CLAUDE.md). New users: start with
[README.md](../README.md). Section numbers (§) are cited throughout the code and
`references/acs-style.md`; where this brief and `acs-style.md` disagree, `acs-style.md` wins.

A Claude Code skill for producing chemical structures, reaction schemes, and arrow-pushing mechanisms in ACS style as **SVG**, for direct use in vector editors (Inkscape; Affinity via `--affinity`) and documents.

Target usage: in any Claude Code session, the user asks in plain language, e.g.

> Draw the mechanism for electrophilic nitration of benzene, including all curved arrows and resonance structures.

Claude writes a declarative scheme spec (YAML), runs the scripts in this repo to validate and render it, **looks at the rendered PNG**, iterates until it is right, and returns `.svg` + `.yaml` (+ `.png` preview).

**Why SVG only.** Claude can rasterize SVG and inspect the result itself, which closes the visual feedback loop. CDXML output cannot be visually verified without ChemDraw. The cost is that output is not ChemDraw-editable chemistry. The internal layout model is kept backend-agnostic so that a CDXML writer could be added later, but that is out of scope.

Items marked **[VERIFY]** are assumptions that have not been confirmed and must be checked (against reference renders or by measurement) before code depends on them.

---

## 1. Design principles

1. **Claude writes chemistry; scripts write geometry.** Claude never hand-edits SVG coordinates. It edits the YAML spec and re-renders.
2. **Our renderer owns ACS style.** No external program does the drawing, so the renderer must implement ACS conventions exactly (bond offsets, label clearance, wedge shapes), measured in pt.
3. **Look before delivering.** Every render produces a PNG, and Claude inspects it, including zoomed crops of dense regions (arrowheads, labels, lone pairs).
4. **Numbers before pictures.** Structural, electron-bookkeeping, and geometric checks run first, so visual inspection is spent on aesthetics rather than on catching errors a script can find.
5. **Editable output.** Clean, flat, grouped SVG with real text (Helvetica by default) that opens sensibly in Inkscape (Affinity: render with `--affinity`, which writes arrowheads as separate shapes and bold/italic as face names).
6. **Curated memory.** Use sessions log problems; development sessions turn them into fixes, tests, or lessons (§9).

---

## 2. Repository layout

```
acs-chem-svg/
├── SKILL.md                 # Entry point for use sessions (keep < ~150 lines)
├── CLAUDE.md                # Instructions for development sessions in this repo
├── README.md                # Introduction and installation
├── docs/design.md           # This file (design brief)
├── pyproject.toml / uv.lock # Pinned environment
├── fonts/                   # Metric font(s) for label measurement (see §4.3)
├── presets/
│   ├── acs1996.yaml
│   ├── acs1996_x0.8.yaml
│   └── acs1996_x0.6.yaml
├── references/
│   ├── acs-style.md         # Numeric style spec + drawing conventions
│   ├── scheme-schema.md     # Authoritative YAML schema with examples
│   └── lessons.md           # Curated, dated lessons (dev sessions only)
├── reference_renders/       # local only (gitignored)
│   └── chemdraw/            # User-supplied ChemDraw SVG/PDF exports (visual ground truth)
├── src/acs_chem_svg/
│   ├── spec.py              # Load + validate YAML
│   ├── chem.py              # RDKit: parse, 2D coords, alignment, stereo, ring perception
│   ├── electrons.py         # Electron-bookkeeping validator for curved arrows
│   ├── fontmetrics.py       # Text measurement for labels and captions
│   ├── draw_mol.py          # ACS molecule renderer → scene graph
│   ├── layout.py            # Scheme layout: species, arrows, text, brackets
│   ├── mechanism.py         # Lone pairs, radicals, curved-arrow geometry
│   ├── geometry_checks.py   # Collisions, crossings, clearances
│   ├── svg_writer.py        # Scene graph → SVG
│   └── preview.py           # SVG → PNG, zoomed crops, side-by-side comparisons
├── scripts/                 # CLI wrappers: render, validate, preview, compare
├── templates/               # Curated molfiles with hand-set 2D coordinates
├── examples/                # Benchmark specs (nitration, SN2, E2, ...)
├── tests/
│   ├── gallery/             # Approved renders for regression comparison
│   └── test_*.py
└── issues/                  # local only (gitignored)
    └── inbox.md             # Append-only log written by use sessions
```

Installation (per project, by symlink): README.md.

---

## 3. Style

### 3.1 ACS Document 1996 base values (ChemDraw)

| Parameter | Value |
|---|---|
| Fixed bond length | 14.4 pt |
| Line width | 0.6 pt |
| Bold width | 2.0 pt |
| Margin width | 1.6 pt |
| Hash spacing | 2.5 pt |
| Bond spacing | 18 % of bond length (2.592 pt at base) |
| Chain angle | 120° |
| Atom label font / size | Helvetica / 10 pt |
| Caption font / size | Helvetica / 10 pt |

Font family is a preset parameter (`font_family`). The default is Helvetica, one of the two fonts named in ACS guidance (Arial or Helvetica); presets can switch to Arial.

### 3.2 Scaled presets

Every length and font size is multiplied by *k*, **except** parameters listed in `exempt` (line width at minimum). The renderer computes geometry at the scaled values directly; no post-hoc SVG scaling is needed.

| Parameter | ×1.0 | ×0.8 | ×0.6 |
|---|---|---|---|
| Bond length | 14.4 | 11.52 | 8.64 pt |
| Bond spacing | 2.592 | 2.074 | 1.555 pt |
| Bold width | 2.0 | 1.6 | 1.2 pt |
| Margin width | 1.6 | 1.28 | 0.96 pt |
| Hash spacing | 2.5 | 2.0 | 1.5 pt |
| Label/caption size | 10 | 8 | 6 pt |
| Line width (exempt) | 0.6 | 0.6 | 0.6 pt |

**[VERIFY]** against the user's ChemDraw presets: exact values, rounding, and whether arrows or arrowheads are also exempt.

### 3.3 Rendering conventions the renderer must implement

These are the rules ChemDraw normally applies. Several are **[VERIFY]**: measure them from `reference_renders/chemdraw/` and record the results in `references/acs-style.md`.

**Bonds**
- Single bonds: stroked line, `line_width`, butt caps.
- Double bonds in rings: second line inside the ring, offset by `bond_spacing`, shortened at both ends. **[VERIFY]** shortening amount.
- Acyclic double bonds: centered (two lines symmetric about the bond axis) or offset toward the substituent side, following ChemDraw's choice. **[VERIFY]** the rule, especially for terminal and C=O bonds.
- Triple bonds: centered, three lines.
- Bonds to labeled atoms stop `margin_width` short of the label's glyph box.

**Stereo**
- Wedge: filled triangle from a point at the stereocenter to width `bold_width` at the far end. **[VERIFY]** wide-end width against ChemDraw.
- Hashed wedge: lines perpendicular to the bond at `hash_spacing`, widening along the bond.

**Labels**
- Heteroatoms and terminal carbons with attached H are labeled (e.g., OH, NH₂, CH₃ where shown explicitly). Implicit H are attached on the side away from bonds (H₂N vs. NH₂ orientation).
- H counts are subscripts; charges are superscripts at the upper right, or circled per user preference.
- The label's anchor is the element letter's visual center, not the text baseline.

**Schemes**
- Compound numbers: bold, centered below the structure.
- Reaction arrows: `line_width`, filled arrowhead. **[VERIFY]** head dimensions against ChemDraw.
- Conditions: reagents above the arrow, solvent/temperature/time below, with true subscripts.
- Resonance: square brackets around the set; ↔ arrows between contributors.
- Ring orientation consistent across a scheme (structures aligned to a common core).

**Mechanisms**
- Lone pairs: two dots placed in the largest angular gap around the atom, clear of the label. Radical: one dot.
- Two-electron curved arrows: cubic Bézier, full arrowhead. One-electron: fishhook (half arrowhead).
- The tail starts at the electrons (lone-pair dots or bond midpoint). The head ends at an atom (new lone pair or charge) or at the midpoint between two atoms (new bond), stopping short of labels.
- Arrows bow away from the molecule's centroid unless the spec overrides; curvature is consistent across a mechanism.

---

## 4. Renderer architecture

### 4.1 Own renderer, RDKit for chemistry

RDKit provides parsing, sanitization, 2D coordinates, ring info, stereo-bond assignment, and core alignment. **Drawing is our own code** (`draw_mol.py`), because mechanisms need:

- exact atom coordinates and label bounding boxes in the final SVG frame, so arrows, lone pairs, and clearance checks share one coordinate system;
- full control over ACS conventions (double-bond placement, label clearance, wedge geometry).

**M0 spike:** render the M1 test set with RDKit's `MolDraw2DSVG` (ACS1996 mode) as a baseline, and compare it against ChemDraw references. If RDKit output is close enough that wrapping it is simpler, document why and revise this section before building the custom renderer.

### 4.2 Scene graph

`draw_mol.py`, `layout.py`, and `mechanism.py` produce a backend-agnostic scene graph: lines, filled polygons, Bézier paths, text runs, and groups, all in pt. `svg_writer.py` serializes it. Geometric checks run on the scene graph, not on SVG text.

### 4.3 Text measurement

Label layout and clearance require real glyph metrics.

- Default font: **Helvetica**. On macOS, measure with `fontTools` from the system Helvetica (`/System/Library/Fonts/Helvetica.ttc`, a TrueType collection; select the Regular and Bold faces by name). **[VERIFY]** the path and face names on the user's macOS version.
- Fallback for environments without Helvetica (e.g., cloud sessions): a freely redistributable Helvetica clone such as TeX Gyre Heros or Nimbus Sans, bundled in `fonts/`. **[VERIFY]** that its advance widths match system Helvetica closely enough. Compare widths for the label character set (C H N O S P F Cl Br I, digits, +, −, parentheses) and record the maximum deviation in pt at 10 pt.
- `fontmetrics.py` records which font file was actually used for every render, and warns when the fallback was used.
- Emit real `<text>` elements with `font-family="Helvetica, 'TeX Gyre Heros', 'Nimbus Sans', Arial, sans-serif"`, so text stays editable.
- Previews must render with the same font that was used for measurement. Point the rasterizer at the same font file, or the PNG will not match the metrics.
- **[VERIFY]** that Affinity's SVG import positions `<text>` identically, including baseline and `dominant-baseline`. If not, use `x`/`y` baselines only (no `dominant-baseline`), or provide an optional text-to-outlines export.

### 4.4 SVG output rules

- Root `width`/`height` in pt, with `viewBox` in pt, so 1 user unit = 1 pt.
- No `transform` attributes in output (bake all coordinates), no clip paths, no `vector-effect`.
- Exception: arrowheads are `<marker orient="auto">` (`markerUnits="userSpaceOnUse"`) in one `<defs>`,
  the one implicit transform. The renderer places each head on its path's end along the end
  tangent, so a head follows when the curve is edited (Inkscape, browsers). Affinity drops markers
  on import, so `render.py --affinity` writes each head as a separate filled shape instead.
- Bold and italic text: the font family list with `font-weight` / `font-style`. `--affinity` writes
  a bare face name in `font-family` (`Helvetica Bold`, `Helvetica Oblique`) instead, the only form
  Affinity imports as bold or italic; Inkscape and Chrome can't resolve a face name.
- Strokes: `stroke-width` explicit on every element, butt caps unless a convention requires round.
- Groups with meaningful `id`s, e.g. `species-benzene`, `species-benzene-bonds`, `species-benzene-labels`, `arrows-step-1`, `conditions-step-1`, `numbers`. These appear as layer or object names in the editor.
- Black on transparent background; no white rectangle.
- Coordinates rounded to 0.01 pt for clean diffs.

---

## 5. Scheme spec (YAML) — sketch

The authoritative schema lives in `references/scheme-schema.md`. Key ideas:

- **Atom maps** (`[C:3]`) give atoms stable identities across the whole scheme, and arrows reference them. Map numbers are global within a scheme.
- Hydrogens or other atoms involved in arrows must be explicit and mapped.
- Species can come from SMILES or from a curated molfile in `templates/`.

```yaml
preset: acs1996
layout: {width: 504}         # pt (7 in); rows wrap between steps to stay within it
species:
  benzene:     {smiles: "[CH:11]1=[CH:12][CH:13]=[CH:14][CH:15]=[CH:16]1", orientation: {rotate: 180}}
  nitronium_2: {smiles: "[O:1]=[N+:2]=[O:3]"}
  sigma_a:     {smiles: "[H:17][C:11]1([N+:2](=[O:3])[O-:1])[CH+:12][CH:13]=[CH:14][CH:15]=[CH:16]1"}
  # sigma_b, sigma_c, bisulfate, nitrobenzene, sulfuric_2 ... (examples/m4_nitration_x1.0.yaml)
steps:
  # row A (nitronium formation) first; then a second sequence on a new row:
  - new_row: true
    from: [benzene, nitronium_2]
    to: sigma_a
    arrows:
      - {from: "bond(11=12)", to: "between(11,2)"}
      - {from: "bond(2=1)", to: "atom(1)"}
  - resonance: [sigma_a, sigma_b, sigma_c]      # in square brackets, joined by ↔
    arrows:                                     # list k turns contributor k into k + 1
      - [{from: "bond(13=14)", to: "between(12,13)"}]
      - [{from: "bond(15=16)", to: "between(14,15)"}]
  - from: [sigma_c, bisulfate]                  # the base sits above the arrow
    to: [nitrobenzene, sulfuric_2]
    omit: [sulfuric_2]                          # regenerated: in the bookkeeping, not drawn
    arrows:
      - {from: "lp(9)", to: "between(9,17)"}
      - {from: "bond(11-17)", to: "between(11,16)"}
numbers: {nitrobenzene: 1}
overrides:                   # escape hatch for aesthetics, applied after auto-layout
  arrows:
    - {step: 4, index: 0, bow: left, curvature: 0.6}
```

Arrow endpoint vocabulary: `lp(n)`, `atom(n)`, `bond(n-m)`, `between(n,m)`, and `electron(n)` for radicals. Arrow types: `two_electron`, `fishhook`. The `overrides` block is how visual iteration happens: Claude adjusts semantic parameters (bow side, curvature, label side, rotation), never raw coordinates.

---

## 6. Validation and iteration loop

Run in order on every render. Failures return actionable messages; Claude fixes the spec and re-runs.

1. **Spec validation.** Schema, references resolve, map numbers consistent.
2. **Chemistry validation (RDKit).** Species sanitize; valences and formal charges are sane.
3. **Electron bookkeeping** (`electrons.py`):
   - Apply each curved arrow to an explicit electron model (bond orders, lone pairs, formal charges) of the step's starting species.
   - Compare the result, as atom-mapped canonical SMILES with charges, to the declared product(s).
   - For resonance sets: contributors share atoms, connectivity, and total charge, and the arrows interconvert them.
4. **Geometric checks** (`geometry_checks.py`, numeric, on the scene graph):
   - Curved arrows crossing bonds, labels, or other arrows.
   - Arrowheads or tails within `margin_width` of label glyph boxes.
   - Lone pairs overlapping bonds or labels.
   - Text overlapping structures.
   - Species overlapping, or gaps outside tolerance.
5. **Visual inspection.**
   - Render the full PNG at ≥ 4× (≈ 300 px per inch of output at minimum).
   - Render automatic zoomed crops of every curved arrow and every labeled atom with charges or lone pairs.
   - Claude reviews them against `references/acs-style.md` and adjusts `overrides` as needed.
   - Cap at ~5 visual iterations per request; then deliver with a note of remaining issues.
6. **Reference comparison** (development only). `scripts/compare` renders our output side by side with a ChemDraw reference PNG at the same scale, plus a difference overlay.

Limitation: validation proves internal consistency, not mechanistic correctness. Research-level mechanisms need user review.

---

## 7. Environment

- Python managed with `uv`; dependencies pinned in `uv.lock`.
- Core dependencies: `rdkit`, `fonttools`, `pyyaml`, `pydantic`, `resvg` (preferred, via CLI or Python binding) or `cairosvg` for rasterization, `pillow`, `pytest`.
- **[VERIFY]** that the rasterizer honors the chosen font and renders `<text>` faithfully. Pick one rasterizer and use it everywhere.
- Use sessions call scripts only via `uv run --project <repo> ...` (`<repo>`: the skill symlink's target).

---

## 8. Reference renders needed from the user

Export from ChemDraw (ACS Document 1996 unless noted) as both **SVG and PDF** into `reference_renders/chemdraw/`:

1. `m1_set`: benzene, aspirin, (R)-alanine with a wedge and a hash, nitronium, acetic acid (C=O and OH labels), a terminal alkyne. Export in all three presets.
2. `scheme_basic`: two structures, a reaction arrow with conditions above and below, a plus sign, bold compound numbers.
3. `mech_basic_x1.0` / `_x0.8` / `_x0.6`: hydroxide + bromomethane SN2 with curved arrows and lone pairs on O; a fishhook pair. Export in all three presets.
4. `resonance_x1.0` / `_x0.8` / `_x0.6`: the nitration σ-complex contributors in brackets with ↔ arrows. Export in all three presets.
5. `scheme_wrap_x1.0` / `_x0.8` / `_x0.6`: a two-row scheme drawn as it would go in a paper (where the arrow breaks, row gap, row alignment). Not yet supplied: optionally `nitration_x1.0`, the full mechanism (HNO₃ + H₂SO₄ → NO₂⁺ + H₂O + HSO₄⁻, attack, bracketed σ-complex, deprotonation by HSO₄⁻) as the M4 style target.
6. Optionally, a figure from a recent manuscript as an overall style target.

Also note the ChemDraw version. These are measured (bond spacing, wedge widths, arrowhead size, label offsets) and used for side-by-side comparison; they are not parsed as chemistry.

---

## 9. Session roles and memory

**Use sessions** (a project folder with the skill linked in `.claude/skills/`):

- Follow `SKILL.md`: write spec → validate → render → inspect PNG and crops → adjust `overrides` → deliver `.svg`, `.yaml`, `.png`.
- **Never edit** `src/`, `scripts/`, `presets/`, or `references/`.
- Append problems and workarounds to `issues/inbox.md` with date, spec excerpt, symptom, and workaround.

**Development sessions** (opened in this repo, governed by `CLAUDE.md`):

- Triage `issues/inbox.md`: each item becomes a fix + regression test, a `lessons.md` entry, or is discarded.
- Run tests and the regression gallery before every commit. Gallery changes require user approval.
- Keep `SKILL.md` short; push detail into `references/`.

---

## 10. Milestones

Each milestone ends with renders the user approves, which are added to `tests/gallery/`.

| # | Scope | Acceptance |
|---|---|---|
| **M0** | Scaffold, `uv` env, presets, `SKILL.md`/`CLAUDE.md` stubs, font metrics, rasterizer choice. Measure the reference renders and record values in `references/acs-style.md`. RDKit-vs-custom renderer spike (§4.1). | Measured-conventions table with values in pt and the reference file each came from. Spike result documented with side-by-side images. |
| **M1** | Single structures, all three presets. | M1 set renders with measured bond lengths, line widths, double-bond spacing, and label clearance matching references within 0.1 pt. Side-by-side comparisons approved. Opens cleanly in Affinity with editable text and named groups. |
| **M2** | Linear schemes: arrows, conditions, plus signs, bold numbers, aligned orientation. | 3-step sequence in all presets, approved. |
| **M3** | Mechanisms: lone pairs, radicals, two-electron and fishhook arrows; electron bookkeeping; geometric checks; zoomed-crop inspection. | SN2, E2, carbonyl addition, and radical bromination pass validation and visual review. Deliberately broken specs are caught. |
| **M4** | Resonance brackets, ↔ arrows, multi-row wrapping. | Full nitration mechanism approved. |
| **M5** | Natural-language workflow polish; benchmark prompts in `examples/`. | From a one-sentence request, Claude produces an approved nitration mechanism without the user editing the spec. |
| Later | Catalytic cycles; foldamer templates; optional CDXML backend from the same scene graph. | — |

Status: M0–M5 are done. M5 was accepted on 2026-10-25: a use session turned the one-sentence
nitration prompt (`examples/benchmarks/prompts.yaml`) into `examples/benchmarks/nitration_accepted.yaml`,
approved without edits; fresh-agent trials of all six benchmark prompts are summarized in
`references/lessons.md`.

CDXML export (the "Later" backend, `render.py --cdxml`, `references/acs-style.md` §13): phase 1,
M1 species, and phase 2, M2 schemes (arrows, conditions, plus signs, numbers, `<scheme><step>`,
§13.1), both passed their ChemDraw gates on 2026-10-28. Abbreviations, curved arrows and
brackets come later.

---

## 11. Sources

- ChemDraw ACS Document 1996 preset values: https://support.revvitysignals.com/hc/en-us/articles/4408234173332
- ACS-style drawing parameters (journal guidelines): https://chemistry.org.tw/downFile/files/news_14_1206068950.pdf
- RDKit drawing (ACS1996 mode): https://rdkit.org/docs/source/rdkit.Chem.Draw.html
- Claude Code skills: https://code.claude.com/docs/en/skills.md
