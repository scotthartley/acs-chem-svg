# Lessons

Curated, dated lessons from development sessions.

## 2026-10-24 (M1)

- **ChemDraw's label clip is offset from the drawn glyph.** A hull-plus-margin model leaves
  a systematic left/right asymmetry of about 0.1 pt. Shifting the hull by (−0.05, +0.08) pt
  removes it (worst 0.060 pt, `measure_refs.py --fit-clip`). Identical geometries in one file
  already differ by 0.05 pt, so don't chase residuals below that.
- **Wedges are not triangles.** The sides extrapolate to ±LineWidth/2 at the stereocenter,
  like the hash outline. The narrow end is cut along bisectors with the neighboring *plain*
  bonds; stereo neighbors are ignored. Drawing a point-apex triangle leaves a visible sliver.
- **Stroked miter polylines reproduce ChemDraw's filled bond polygons.** That covers degree-3
  points and ring corners, with no extra fill wedges (< 0.2% of ink pixels differ at 20 px/pt).
- **`measure_refs.measure()` mutates the Doc** (`Shape.owner`). Run `assign()` again only on
  a copy with the owners cleared.
- **Templates are in ChemDraw's SVG frame**, not CDXML's: `load_cdxml` already adds the path
  offset. Scaling a 1.5 Å molfile by the *median* bond (11.51 vs 11.52 at ×0.8) drifted far
  atoms by 0.03 pt. `chem.to_molecule` maps 1.5 Å layouts exactly.
- **RDKit with `sanitize=False` parsing never perceives E/Z.** Call `AssignStereochemistry`
  before CoordGen, or Z-alkenes come out E.
- **`Path("acs1996_x0.8").suffix` is `".8"`.** Don't detect file paths by suffix.
- **CDXML fragment order is not visual order**: labels_x1.0 fragment 1 is –NH₂, not H₂N–.

## 2026-10-24 (M1 follow-up: stereo_x*, clearance_x*, double_bonds_2)

- **Ask for references that sweep the variable.** The first M1 references showed slanted bonds
  into H/I/N as "outliers". Bonds at 12 directions for 7 labels showed they weren't outliers
  of a hull model: the hull model was wrong. ChemDraw measures the glyph's extent along the
  bond within a band around it (§2.1).
- **ChemDraw right-justifies two-letter elements** when the bonds point right (Cl–: the l is
  on the atom). The same test sets the H side. A 4.7 pt "clip error" was a placement error.
- **ChemDraw ties follow the drawn bond direction**, not the screen. My "upper side" guess
  for E-alkenes matched one drawing and failed at 180°.
- **RDKit's molfile writer drops wedges on non-stereocenters.** ChemDraw wedges are drawing
  marks, so `cdxml_to_molfile.with_wedges` writes them into the bond block. `build_mol`
  ignored `WedgeEnd` / `WedgedHashEnd` displays until this was fixed.
- **Thresholds measured on one side only need a bracket.** The hash count step lies between
  D/HashSpacing = 1.997 and 2.05; a cut exactly at 2.0 flipped when the clip model moved by
  0.01 pt.
- **An "off-anchor" label is a rule until shown otherwise.** I excluded double_bonds_2's
  aldehyde O labels as hand-moved; they are ChemDraw's default: a terminal label on an offset
  double bond moves bond_spacing/2 toward the second line. Measure the offset in the bond
  frame (perpendicular 1.30 = s/2 at both orientations) before calling anything manual.
- **Compare unions, not shapes, at joins.** Our butane wedge matched ChemDraw's polygon
  (Hausdorff < 0.01) yet the join showed a gap: ChemDraw's *neighbor* bond starts 1.5 pt
  behind the atom to fill the notch. Look at the zoomed overlay, not just per-shape numbers.

## 2026-10-24 (M2: schemes)

- **LineHeight "auto": the pitch depends on the content.** Plain lines are size + 0.05 (hhea
  ascent + descent = 1 em), but a subscript on the upper line adds 0.775 pt and one on the lower
  line a further 0.15 (×1.0). Two data points per preset fit it; neither ink nor font extents do.
- **Check kerning against an independent source.** `measure_refs` places glyphs inside a run
  with our own kerning, so "the 11 ink width matches kerned Helvetica" was circular. ChemDraw's
  PDF glyph origins (`pdf_char_origins`) settled it: ChemDraw kerns.
- **ChemDraw stores typed formulas as fragments.** K₂CO₃, LiAlH₄ and PCC are one-node fragments
  with a label; treating them as structures broke rows and gaps. They are captions.
- **Align by ring center, not by ink, when every structure has a ring.** Row 1 of scheme_multi
  (OMe on top, OBn below) has ring centers on one line and ink centers 9 pt apart; scheme_basic,
  with ethyl iodide, is ink-centered. Both sit 0.25 above the axis: the center of the arrow's
  bounding box, which ChemDraw's asymmetric head pulls up.
- **Hand placement shows up in the scaled references.** The ×0.8/×0.6 conditions overhang
  their arrows (pad down to −1.4 pt): ChemDraw scaled the text about its anchor. Use ×1.0 × k
  for such values and say so in the preset comment.
- **rdAbbreviations matches aromatic patterns.** On our Kekulé mol it silently found OMe but not
  OBn. Find groups on an aromatic copy, then collapse them on the Kekulé mol.
- **CoordGen's coordMap holds atoms to ~0.005 Å, not exactly.** RDKit's matching depiction is
  exact but lays new atoms out at 1.5 Å; scale the reference to 1.5 Å and back.
- **The conservative stroke bbox is not ink.** `Polyline.bbox` pads miter joins by a full line
  width; scheme gaps measured 0.3 pt wide until the layout used `ink_bbox`. M1 kept `bbox`, so
  approved renders didn't move.
- **Affinity and resvg disagree on how to name a face.** Affinity ignores font-weight and
  font-style and needs "Helvetica Bold" / "Helvetica Oblique" in font-family; resvg matches only
  the family "Helvetica" (the n of n-BuLi vanished from ChemDraw's reference PNG). We write face
  names and `preview.rasterize` maps them back. The PNG review can't catch an import problem:
  test new text styles in Affinity with a one-row-per-encoding file.

## 2026-10-24 (placement defaults)

- **Hand-placed reference positions are not data.** Arrow gaps, conditions offsets, number
  gaps and charge positions in the references were set by eye (arrow gap 10.3–22.0 pt); their
  medians (13.11) looked precise but weren't. Presets now use round defaults on a 0.5 pt grid,
  each preset rounded from its own measurement, with the median kept in the comment
  (`placement:` list). Don't copy artifacts of ChemDraw's drawing either: its 0.25 pt
  structure offset came from the asymmetric arrowhead's bounding box, so with our symmetric
  head structures sit on the axis.

## 2026-10-24 (M3: mechanisms)

- **A one-atom fragment isn't always a formula.** measure_refs read every single-node fragment
  with text as a caption (K₂CO₃), which dropped Br⁻ and Cl· from mech_basic. A formula has a
  `NodeType` (Nickname, Fragment, Unspecified); a lone atom has an `Element` and none.
- **ChemDraw scales the hand-drawn curve but not its head or dots.** The curved arrows' heads
  and strokes are the same at ×1.0/×0.8/×0.6, while the curves scale with the document. Its ×0.8
  dots are the ×0.6 size (0.533), so dot radius, spacing and lone-pair gaps are not k-scaled.
- **Overlays find what tables miss.** The fishhook's stroke ends 0.25 pt beside the half
  head's axis, toward the barb. That showed up as red and blue slivers in the compare overlay
  before any measured number did.
- **seg_seg_dist is not a crossing test.** It is the minimum endpoint-to-segment distance, so
  two crossing segments report a positive distance. Mechanism checks use
  `geometry.segment_distance`, which tests the crossing first.
- **A lone pair can straddle a bond.** With the bond between its dots, no dot touches the bond,
  but the pair still reads as on top of it. The check tests the segment between the two dots too.
- **Pick the arrow's free choices by what it would hit.** Fixed rules (bow away from the
  centroid, the pair facing the bow) still landed heads on C=O and triple bonds. The planner
  builds each candidate (both bows × each lone pair) and keeps the one with the fewest clashes.
  The dot set an arrow leaves is drawn facing its target before the layout runs.
- **The linear-chain rule rules out by-products.** Each step's `from` must start with the
  previous `to`, so HBr from the H abstraction would have to be carried into the Br₂ step.
  Propagation is two specs until M4 layouts. **Closed in M4 (2026-10-25):** `from` starts with
  a prefix of the previous products; the rest stay behind as by-products.

## 2026-10-24 (M3: gallery review)

- **A fixed-point solve can settle on a curl that points backwards.** The tip, stroke end and
  end tangent were iterated together. For C=O → O and O–H → O, it converged (or stalled) on a
  short loop beside the tail whose head pointed back at the bond. `mechanism.solve_end` now
  scans approach directions round the target and bisects the end-tangent condition. Several
  roots can exist, and the planner tries each. A sign change at a jump in the tip distance (a
  dot pushing the tip out) is not a root.
- **Greedy arrow choices box each other in.** Fishhooks meeting between two species share a
  side, and a bond's second arrow takes the other side. So the first arrow's bow decided where
  the third could go, and in the H abstraction that forced it across a methyl. `plan_step` now
  tries every combination of bows and keeps the one with the fewest clashes.
- **Lone pairs stay put; the arrows are decorations.** Lone pairs on target atoms (the user's
  call) put dots exactly where heads land. Moving the pairs out of the way was rejected at
  review: pairs belong in their usual places. Instead the planner may stop the tip 1.5 or
  3.0 pt farther out or flatten the curve (×0.75, ×0.5). A head more than 10° off its atom
  counts as a clash. Pushing the tip past each dot was tried and dropped: it makes the tip
  distance jump, so the end tangent has no exact solution and the head points at the pair
  instead of the atom.
- **Clear the dots, not only the label, and measure the gap to the nearest point.** head_gap
  was measured along the head's line through the box center. A head coming down onto a pair
  stopped between its dots (E2, SN2 ×0.6), and a diagonal approach sat about head_gap / √2
  from the box corner (SN2 ×0.8, ×1.0). The user: "It should end just above."
  `mechanism.clear_along` now keeps head_gap from the box and every dot of the target atom.
  A short bond-to-own-atom arrow often can't then aim exactly at the center: the shape's end
  tangent is fixed relative to its chord. So for labeled targets the head only has to point
  into the ink box past the dots (`aims`). Off a root, the head's own direction moves its
  tip, so the end is re-solved for the fixed tip (a plain fixed-point iteration oscillated
  there). Flipping a bow now ranks behind all shape choices, so the default side wins when
  it is clean.
- **ChemDraw can put a Bézier handle on its end point.** The new C–Br → Br arrow in mech_basic
  has c2 = p3, so the end tangent comes from c1. Both measure_refs and `head_polygon` fall
  back to it.
- **Fishhooks that form a bond meet in the middle, but a plus sign is often there.** They
  point at a spot just past the plus sign on their shared side, rather than at either atom.

## 2026-10-25 (M4: resonance, sequences, rows)

- **Charges belong to the structure's ink.** The σc gap read 12.56 against 7.68 elsewhere
  because `_structures` left out charge bars, and σc's ⊕ faces the arrow. With charges joined
  to their structure every ×1.0 gap is 7.68 (`resonance_summary`).
- **Measure units without placing them twice in y.** Rows depend only on x widths, and a
  species' width doesn't depend on its vertical anchor. So `build_scheme` places once to
  measure, breaks greedily, and places again with the breaks and each row's own ring rule.
  Single-row output stayed byte-identical (all 37 gallery SVGs) through the refactor.
- **A Bézier's control hull is not its ink.** Rows were stacked with curved arrows' control
  points, which reach far past the curve: the first nitration row gap came out ~58 pt instead
  of 34.5. `Path.ink_bbox` solves the cubic's extremes.
- **Plan a step's curved arrows with its reactants' row.** When a step's arrow starts a new
  row, its reactants (and so its curved arrows) stay on the row above. Planning them with the
  arrow's row moved them with the wrong row.
- **Abbreviations must keep the attach atom's map.** `condense` dropped the map number, so
  contributors aligned by atom maps left NO₂ free, and CoordGen put it on the ring after a
  flip. The dummy atom now keeps the map.
- **Adding obstacles can't change a clean plan.** Reaction and ↔ arrows joined the curved-arrow
  planner's obstacles. The planner keeps the lowest score, and new obstacles only raise other
  candidates' scores, so every approved (clash-free) plan stays the same.
- **Ring bonds next to the tail's bond need the far side.** A head on an existing bond sits on
  the tail's side of it; for adjacent ring bonds that side is inside the ring. The planner now
  also tries the far side.
- **Added species belong above an arrow (user review).** HSO₄⁻ floating after the bracket read
  as coming from nowhere. Reactants a mechanism step adds now sit above its arrow, and the
  arrow stays on its reactants' row, so a row may end with an arrow.
- **Don't pin atoms whose neighbors change.** Nitrobenzene copied σc's N position, where C11 was
  sp3, so its C–N bond was off the ring bisector. `chem.stable_core` pins only the largest
  connected group of matched atoms with unchanged neighbor counts.
- **Fallback heads don't hold the planned stop.** Water's O–H → O loop above an arrow can only
  aim (no exact end-tangent root); its tip lands 3.77 pt out, between the 2.5 and 4.0 stops.
  The head-gap test now accepts the whole FARTHER range.
- **A rule from the references can block a fix.** Moving a crowded N⁺ charge to the largest
  gap changed m1_set_x*: ChemDraw keeps N⁺ upper right with bonds at 0° and 180° (45° off
  every corner). The fallback only applies when a bond is within 45° of every corner.


## 2026-10-25 (M5: baseline trials, round r0)

Six fresh agents ran the benchmark prompts (`examples/benchmarks/`) on a frozen copy of the
skill. Validate runs / renders: synthesis 2/2, SN2 5/5, E2 4/4, radical 4/4, carbonyl 7/6,
nitration 21/6. All delivered; nitration at x0.6 after giving up on the default preset. Inbox
triage (each item: fix + test, lesson, or discarded):

- **Unbracketed atoms lost their H in the product comparison** (3 of 6 trials: `C[C:3](C)C`,
  `S(=O)(=O)O`). `canonical` set NoImplicit before folding the implicit H into the explicit
  count. Fixed; `test_unmapped_atoms_keep_their_implicit_h`.
- **The mismatch message didn't say which atom was missing.** It now lists map numbers found
  on one side only, and validate prints paste-ready `suggested:` products. Tests in
  `test_electrons.py`, `test_broken.py`.
- **fontTools parse notes on stderr** looked like the fallback-font warning (4 trials).
  Silenced the fontTools logger; `test_validate_prints_nothing_on_stderr`.
- **Lone pair on an unlabeled carbanion** (C≡N⁻): the message now says
  `add show_carbons: [n]`; `tests/broken/lone_pair_unlabeled_carbon.yaml`. Auto-labeling was
  left for later: it changes renders.
- **Overrides were blind.** render.py prints each curved arrow's bow;
  `test_the_reported_bow_reproduces_the_default`.
- **Nitronium drawn bent** (r0 nitration, x0.8): aligned to HNO3's bent atoms because its
  stable core is 2 atoms and the fallback pinned everything. An alignment that bends an sp
  atom by more than 10° is now dropped for the species' own layout, so `align: false` is no
  longer needed there. Gallery unchanged; `test_an_alignment_that_bends_nitronium_is_dropped`.
- **Validate took 13–21 s on nitration**, which made override loops expensive. The arrow
  planner's `clashes` now skips pairs whose boxes are farther apart than the threshold
  (exact: such a pair can't be closer). 21.7 s → 4.2 s; gallery byte-identical.
- **Switching presets to make rows fit** (nitration): a doc fix. SKILL.md and
  `references/workflow.md` say to keep the user's preset; three rows is expected.
- Discarded as cosmetic defaults, for gallery review later rather than M5: the SN2 C–Br loop
  sits above the bond by default; cramped homolysis fishhooks on Br–Br; an E2 partner facing
  away (documented in workflow.md's common mistakes instead).

## 2026-10-25 (M5: re-trial, round r1)

Same six prompts on the fixed skill. Validate runs / renders (r0 → r1): synthesis 2/2 → 2/1,
SN2 5/5 → 4/4, E2 4/4 → 4/3, radical 4/4 → 3/3, carbonyl 7/6 → 7/6, nitration 21/6 → 7/4. The
nitration trial stayed at acs1996 (three rows), used `--suggest`, and with its one override
removed renders pixel-identical to the approved `m4_nitration_x1.0`. Inbox triage:

- **A suggestion can hide a wrong arrow** (SN2: `lp(1) -> atom(2)` gave `[OH+:1]`, `[CH3-:2]`).
  The mismatch now names such an arrow and says `to: "between(1,2)"`; the suggestion header
  says to fix the arrows if the products look wrong. `test_lp_to_another_atom_points_to_between`.
- **Suggested water drew as HO–H** (E2). An atom holding a mapped H atom now gets its other H as
  atoms: `[O:1]([H:4])[H]`. `test_suggested_water_draws_both_h`.
- **RDKit's valence error named an atom index** (nitration: `[H:17][CH:11]`). It now names the
  species and map number and says to drop the bracket H; RDKit's own log line is blocked.
  `tests/broken/h_written_twice.yaml`.
- Doc fixes: reagents named in the request go over the arrow and no unstated conditions
  (synthesis); `--maps` from the first render (carbonyl); a step whose arrows would cross a
  by-product goes on its own row (radical).
- **O4's lone pair touches the N–O4 bond end in nitric_h** (nitration). Also in the approved
  gallery render; the dots-on-bonds check misses it by a hair. Open for the user: fixing it
  changes an approved image.
- Discarded: hν reads as "hv" at 4 px/pt (the SVG has U+03BD; Helvetica's ν looks like that);
  "arrow crosses the O4 charge" (the crop shows it passes beside it).

## 2026-10-25 (M5: acceptance)

- **Accepted.** In a use session (`~/Claude/struct_svg_gen/`, skill symlinked) the user gave the
  one-sentence nitration prompt and approved the render without editing the spec. Saved as
  `examples/benchmarks/nitration_accepted.yaml`, gallery case `nitration_accepted`.
- **The session reused the approved example.** The accepted spec is
  `examples/m4_nitration_x1.0.yaml` line for line (`_2` names, `align: false` and its comment),
  so it shows the skill finding and adapting its own example, not building the mechanism.
  `references/workflow.md` names that file as the approved result, which invites this. The
  evidence for building from scratch is the r1 trial, which could not read `examples/` (7
  validate runs, 4 renders, same drawing). For a future acceptance of a new mechanism, use a
  prompt with no example in the repo.
- **Inbox: no silyl abbreviations** (same session, `abbreviate: [TMS]`). Added TMS, TES, TBS,
  TIPS to `chem.ABBREVIATIONS`; `test_silyl_abbreviations`. They follow the existing rule: a
  bond leaving to the right is unmeasured, so turn the structure.

## 2026-10-26 (inbox: ring carbanion lone pair vs ⊖)

- **Place dots before the charge on an unlabeled atom.** A ring [CH-] with an arrow from its
  pair had the ⊖ on the outward bisector first; the pair then split the 240° exterior gap and
  landed 60° from the ⊖, touching its bar, and `tail_side` couldn't help (one set only). Now the
  dots take the bisectors and the charge goes beside them, inside the largest bond gap. Scoring
  by clearance alone first put the ⊖ inside the ring (nothing is near there): restrict charge
  directions to the exterior gap. `test_mechanism.py::test_ring_carbanion_pair_on_bisector_charge_beside`.
- **A clean geometry check is not a readable drawing.** The fix moved unlabeled cyanide's pair
  onto the axis, so `tests/broken/lone_pair_unlabeled_carbon.yaml` stopped failing, though
  ⁻:≡N with no C reads badly. validate now flags dots on any unlabeled atom with ≤ 1 bond
  directly (user decision), rather than relying on a dot/bond overlap to catch it.

## 2026-10-26 (arrowheads as SVG markers)

- **Markers are the default; Affinity needs `--affinity`.** Every arrowhead is now a
  `scene.Marker` on its shaft, written as `<marker orient="auto" markerUnits="userSpaceOnUse">`
  in one `<defs>`, so a head follows the curve when it is edited (Inkscape, browsers).
  `Marker.place` returns the leaf the heads used to be, and `head_leaves()` places them, so
  bboxes, `head_polygon`, `clashes` and `layout._arrow_parts` read the same geometry.
- **Probe (user check):** Affinity imports the shafts with no heads at all (both marker unit
  modes, `auto` and `auto-start-reverse`). Inkscape → copy → paste into Affinity brings the heads
  over, but as loose shapes that no longer follow. resvg draws markers like the separate shapes
  (gallery: ≤ 0.26 alpha, no pixel past 0.3 at 4 px/pt).
- **`render.py --affinity`** (`to_svg(affinity=True)`) writes each head as a separate shape after
  its path: byte-identical to the pre-marker output, so the approved gallery PNGs match it
  exactly. compare.py uses it too, since `measure_refs.load_rendered` parses the head paths.
- The resonance start head is a turned-around shape with `orient="auto"` (SVG 1.1), not
  `auto-start-reverse` (SVG 2).

## 2026-10-26 (bold and italic: face names only for Affinity)

- **Face names broke bold and italic in Inkscape and Chrome.** `font-family="Helvetica Bold"`
  (written since 2026-10-24 for Affinity) is no family to Inkscape or Chrome: Inkscape drew a
  wide substitute sans, Chrome drew Times (user report; reproduced with Inkscape's CLI export and
  headless Chrome). Safari and Affinity resolve it; resvg never saw it (`rasterize` maps it
  back), so no test or PNG review caught it.
- **No single encoding serves both.** Probe (one row per encoding): the face name first in the
  family list plus `font-weight` / `font-style` renders correctly in Inkscape and Chrome, but
  Affinity imports it as regular (user check). PostScript names (`Helvetica-Bold`) fail in
  Inkscape. So the default writes the family list with the attributes, and `--affinity` keeps the
  bare face names, alongside its separate arrowhead shapes.
- Editor checks of a new encoding need all three: Inkscape, Chrome and Affinity. Inkscape and
  Chrome can be checked here (`inkscape f.svg --export-type=png`; headless Chrome needs to run
  outside the sandbox); Affinity needs the user.

## 2026-10-26 (inbox: phenol H abstraction by X·)

- **Generic atoms are a drawing label over a real stand-in.** SMILES can't put a radical on
  `*`, so `[*:9]` with `electron(9)` failed. `labels: {9: X}` (spec level, by map number)
  replaces the drawn symbol only. Valence, lone pairs, radicals and the octet checks still use
  the stand-in (`[Cl:9]` for X·, `[H:6][Cl:9]` for HX). Pick a stand-in with the bonds and H
  count the label should show: a `[CH3:9]` stand-in for R would draw RH₃.
- **A resonance bracket may follow a step with by-products.** The first product opens the
  bracket and the rest follow it with "+" (`… ] + HX`). Before, abstraction and resonance had to
  be two specs.
- **Fishhook bookkeeping stays strict, and its message now names the missing arrow.** H
  abstraction takes three fishhooks (`m3_radical_propagation_x1.0`). With two, both odd bonds
  now say `add {type: fishhook, from: "bond(1-2)", to: "between(2,9)"}`.
- **Two fishhooks from a ring double bond.** In phenoxy → o-radical, the fishhooks from
  `bond(3=4)` took opposite sides, so the one to `atom(3)` started inside the ring on the
  double bond's inner line, and no `bow` combination validated. Both fixes are allowed (user
  decision): an inside tail starts `bond_tail_gap` past the inner line (offset double bonds
  only; centered ones keep the gap from the axis), and the second arrow may share the first's
  side, its tail moved `bond_tail_gap` along the bond toward its target. Same-side scores as two
  flips. At one flip it tied with the opposite side's flip and changed
  `m3_radical_propagation_x1.0` on rank alone. Phenol itself needs `orientation: {flip_h: true}`
  so its O–H faces X· (`tests/specs/phenol_x_abstraction.yaml`).

## 2026-10-26 (inbox: acetone aldol condensation)

- **`validate` never looked inside a species.** `check_scheme` compared whole species boxes,
  so a CH₃ label running into an O⁻ in the same molecule passed. `geometry_checks.check_labels`
  now reports a label that overlaps another atom's label or a bond not on its atom, for schemes
  and M1 rows (`tests/broken/label_overlap.yaml`). It is clean on every example and on all 482
  ChemDraw templates.
- **A crowded CH₃ flips to H₃C** (user decision, unmeasured; ChemDraw keeps `h_side`'s rule).
  `Drawer.flip_crowded_labels` keeps the flip only when the other side is clear. A gem-dimethyl
  beside O⁻ at `{flip_h: true, rotate: 180}` is blocked both ways (H₃C would hit the other
  methyl's H₃C), so the fix there is a different `orientation`, and the validator says so.
- **`align` / `orientation` on a step** got pydantic's "extra inputs not permitted". It now says
  they belong on the species.
- **One explicit C–H drew only that H.** `[H:6][CH:5]([H:20])…` showed two of three H's.
  `chem.parse_smiles` now adds a carbon's remaining H's as unmapped atoms when any of its H's is
  written as an atom (SMILES only; molfiles keep their own). The electron bookkeeping already
  treats unmapped H as implicit, so no spec needed changes. `m3_e2_x1.0` gained two C–H bonds,
  and `tests/broken/arrow_crosses_bond.yaml` needed `orientation: {rotate: -60}` on etbr to keep
  its arrow crossing the chain.

## 2026-10-26 (inbox: Lewis structure of CCl₄)

- **A user asked for CCl₄ with shared electrons as dots and no bond lines.** Bonds were always
  lines and internal carbons couldn't be labeled. `lewis:` on a species now draws a Lewis
  structure on a square lattice: every atom labeled, bonds as 1–3 shared-pair columns, or as
  electrons assigned per `formal_charge` / `oxidation_state` (references/scheme-schema.md).
- **"No non-bonded atoms one lattice step apart" can't be a hard rule.** A CH₃ on a three-way
  center (acetic acid, isopropanol) always puts one of its H's next to another branch. The
  layout minimizes such pairs, and the facing sides get no lone pairs, so nothing sits where
  a shared pair would. Four branches with H's (neopentane) don't fit at all and raise an error.
- **`no:` in YAML is `false`.** A species named `no` (nitric oxide) fails with a baffling
  "Input should be a valid string". Use another name (`nitric_oxide`).


## 2026-10-26 (curved arrows: hooks, aim, smoother arc)

- **Chord length alone is the wrong switch between shapes; adjacency is.** Arrows onto a bond
  next to their source (E2's C–H → C–C, a sigma complex's C=C → C–C⁺) were loops only because
  their chords are short, so they curled. They are now hooks built from end tangents: leave
  along the source bond's normal, land along the target bond's normal, circular-arc handles.
- **A C–H → C–C hook has no room on its own side: it goes through the forming π bond and
  back (user choice 2026-10-26).** The tail and the near side of the target sit in the same
  ~110° angle between the bonds, so the head squeezed in beside the tail. A head that merely
  crossed the bond was rejected. Six options were rendered side by side (`build/arrow_review/
  cramped_options.png`, same framing for each). Chosen: through the bond, a large curl
  (handles ≥ 0.9 bond lengths), the head pointing back at the bond turned 30° toward the
  shared atom. Define such a turn geometrically (toward the shared atom), not as a screen
  rotation, or a mirrored drawing gets the other option. Ring C=C → C–C hooks start outside
  the ring and keep the near side. Keep both plain sides as fallbacks: phenoxy's C=C → C–O
  needs one, because the radical's fishhook meets it there.
- **Show options with fixed framing.** Crops centered on each arrow's own box made six
  variants look alike; centering every crop on the target bond at one scale made them
  comparable.
- **Bond → own-atom arrows stay loops.** A hook was tried (end tangents chosen, 36 landing
  directions, the shortest clean stroke first). It failed: between a bond's midpoint and its
  own atom there is barely room for the 4.8 pt head. On H–C → C the tip and tail share one
  90° gap between bonds, so the head lay across its own tail. On Br–Br and C=O the clean
  candidates swung wide around the dots. ChemDraw's loop (out and back) fits. What the loops
  needed was the tail at the *drawn* bond's midpoint (a label shortens it) and a stricter aim.
- **"Enters the ink box" is too loose an aim.** Water's O–H → O head clipped the O's box 28°
  off its center and read as pointing at the lower pair; a clean candidate with a 0° miss
  existed at curvature × 0.5 but ranked lower. Labeled targets now also need ≤ 20° from the
  box center. A miss costs 0.5 + miss/360 of a clash. With a flat penalty, ties went by rank
  (nitration ×0.6 jumped from 22° to 48°). With 1 + miss/180, an aiming head that crossed a
  charge beat a clean miss, and validation failed. A miss must cost less than a clash.
- **ChemDraw's arc can't be smoothed to strictly rising curvature with its tangent angles
  kept.** No handle pair does it. Refitting only the handle lengths gives an almost flat
  profile (2% wobble) and removes the head spike. The arc comes out 16% taller per chord.
