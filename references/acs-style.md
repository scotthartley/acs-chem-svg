# ACS style: measured conventions

Every number here was measured from the ChemDraw exports in `reference_renders/chemdraw/` by
`scripts/measure_refs.py`. To reproduce, run `uv run python scripts/measure_refs.py --summary`
for the key table, or `--json build/measurements.json` for everything. No value was filled in
from memory; anything the references don't show is listed under **Open**.

- **Source:** ChemDraw 26.1.0.6801, ACS Document 1996, exported as SVG + PDF + CDXML (2026-10-24).
- **Units:** ChemDraw's SVG says `px`, but 1 SVG unit = 1 pt. The viewBox equals the PDF
  page size in pt for all 10 files.
- **SVG text:** text is exported as real `<text>` (Helvetica). Each label is split into runs
  at size changes: `H`, `2` (7.5 pt), `N`.
- **How ChemDraw draws:** bond lines, wedges, hashes, charge signs and arrow shafts are filled
  polygons. Arrowheads are filled paths plus a 0.05 pt hairline outline. Lone pairs and radicals
  are filled circles. Curved arrows and brackets are stroked paths.
- **Method:** `measure_refs.py` takes atom positions, bond orders, stereo flags and object types
  from the `.cdxml` saved with each export. CDXML maps to SVG by the path transform
  (`svg = cdxml + (e, f)`). Each SVG polygon is assigned to its bond in that bond's frame
  (u along a→b, v perpendicular). Glyph boxes and outlines come from fontTools (system
  Helvetica). Every value is a median over all instances, with the range given when it is
  wider than 0.01 pt.
- **Classes:** **exact** means ChemDraw-generated geometry, measured to 0.01 pt and held to
  the M1 0.1 pt tolerance. **approx** means hand-placed positions (by eye in ChemDraw); they
  are recorded with their spread. Presets use round **placement defaults** for them (0.5 pt
  grid, user decision 2026-10-24), with the measured value kept as provenance.

k = preset scale (1.0 / 0.8 / 0.6). The ×0.8 and ×0.6 files were built with scaled
Document Settings, not Object › Scale. Their CDXML headers show:

| Document setting | ×1.0 | ×0.8 | ×0.6 |
|---|---|---|---|
| BondLength | 14.40 | 11.52 | 8.64 |
| LineWidth | 0.60 | 0.60 | 0.60 |
| BoldWidth | 2 | 1.60 | 1.20 |
| MarginWidth | 1.60 | 1.28 | 0.96 |
| HashSpacing | 2.50 | 2 | 1.50 |
| BondSpacing (% of length) | 18 | 18 | 18 |
| LabelSize / CaptionSize | 10 | 8 | 6 |
| ChainAngle | 120 | 120 | 120 |

## 1. Bonds

| Parameter | ×1.0 | ×0.8 | ×0.6 | Source | Method | Class |
|---|---|---|---|---|---|---|
| Bond length (atom–atom) | 14.40 | 11.52 | 8.64 | m1_set_x*, double_bonds, labels | CDXML atom distance, all bonds (n = 32–36 per file) | exact |
| Line width | 0.60 | 0.60 | 0.60 | all | perpendicular extent of each bond polygon in bond frame (n = 44/file) | exact, **exempt** |
| Double-bond spacing (centerline to centerline) | 2.592 | 2.074 | 1.556 | m1_set_x*, double_bonds, resonance_x* | Δv of the two line polygons (= 0.18 × bond length) | exact |
| Triple-bond spacing | ±2.59 | ±2.07 | ±1.56 | m1_set_x* (propyne) | v of the 3 lines (−s, 0, +s); all three run the full length | exact |
| Ring inner line: inset from each unlabeled atom | 1.497 (1.493–1.501) | 1.198 | 0.899 | m1_set_x*, resonance_x*, scheme_basic_x* | centerline u of the inner polygon's ends, n = 12/file; equals s·tan 30° (1.496 / 1.197 / 0.898) | exact |
| Ring inner line end shape | butt (perpendicular) | | | m1_set_x1.0 | edge vectors of the inner polygon | exact |
| Joins at unlabeled atoms | mitered | | | all | see below | exact |

**Joins.** Each bond end at an unlabeled atom is cut along the bisector of the angle to the
neighboring bond on each side.
- At a degree-2 atom this is a single miter cut, so adjacent bonds share an edge.
- At a degree-3 atom the two cuts meet at the atom center, which gives a pointed end.
- At a degree-1 atom the end is butt, perpendicular to the bond.
- Overall, the ChemDraw output looks like miter-joined polylines. It does not look like
  separate butt-capped strokes (design.md §3.3 discrepancy, below).

### 1.1 Double-bond placement (source: double_bonds_x1.0, double_bonds_2, m1_set_x1.0, labels_x1.0)

| Case | Ends (degree) | Placement | Second-line insets | Class |
|---|---|---|---|---|
| propene, terminal C=CH₂ | C(2)=C(1) | offset; 2nd line on the substituent side | 1.499 at the substituted end, 0 at CH₂ | exact |
| 2-methylpropene | C(3)=C(1) | centered (v = ±1.296) | both lines run into the C(3) joins | exact |
| (E)-2-butene | C(2)=C(2), substituents on opposite sides | offset; 2nd line on one side (v = −2.589) | 1.496 only at the end whose substituent is on that side; 0 at the other | exact |
| (Z)-2-butene | C(2)=C(2), substituents same side | offset toward both substituents | 1.499 / 1.494 at both ends | exact |
| acetaldehyde (H drawn), acetone, cyclohexanone C=O | C(3)=O(1) | centered | both lines end 5.52 from the O center | exact |
| methylenecyclohexane (exocyclic C=CH₂) | C(3)=C(1) | centered | — | exact |
| nitronium O=N⁺=O | O(1)=N(2), collinear | centered | clipped at both labels | exact |
| (E)-2-butene at 0°, 90°, 180°, 270° | C(2)=C(2), tie | offset; 2nd line at v < 0 of the bond's drawn direction (B→E) every time, toward E's substituent | 1.49–1.50 at E, 0 at B | exact (double_bonds_2) |
| acetaldehyde, H implicit, two orientations | C(2)=O(1) | offset toward the methyl | 1.49–1.50 at C; clipped at O | exact (double_bonds_2) |
| 1-methylcyclohexene, cyclohex-2-enone, 3-methylcyclohex-2-enone | ring C=C with exocyclic substituents | 2nd line inside the ring | 1.49–1.50 both ends, also at degree-3 atoms | exact (double_bonds_2) |
| naphthalene (Kekulé), fused C=C | shared by two 6-rings, 3 double bonds each | v < 0 side (tie) | 1.497 / 1.498 | exact (double_bonds_2) |
| indene, fused C=C | 6-ring (3 double bonds) vs 5-ring (2) | inside the 6-ring | 1.497 / 1.498 | exact (double_bonds_2) |
| indene, 5-ring C=C | ring angle 108° | inside | **1.880 / 1.886** = s·cot 54° (1.884) | exact (double_bonds_2) |

**Rule** (fits every case above):
- **Centered:** a bond whose terminal end (degree 1, or only collinear substituents) faces a
  degree-3 atom or a label is drawn centered; so is a bond with both ends open.
- **Offset otherwise:** the second line goes on the side that holds more substituents. A tie
  (E-alkenes; tetrasubstituted) goes to the v < 0 side of the bond's drawn direction, at any
  rotation, so it follows how the bond was drawn, not the screen.
- **Rings:** the second line goes inside the ring. A bond shared by two rings uses the ring
  with more double bonds (indene: the benzene ring), then the larger ring, then v < 0
  (naphthalene).
- **Insets:** the second line stops on the bisector of the angle to the substituent on its
  side: s·cot(φ/2), i.e. s·tan 30° at 120° and 1.884 at 108° (both measured). No substituent
  on that side: no inset. Labeled end: the label's clip distance.
- **Centered-line ends:** at a degree-3 atom, both lines of a centered bond run on until they
  meet the neighbor bond's centerline (u = −h/tan 60° = −0.748 for h = 1.296). At a
  degree-1 carbon they end square at the atom.

### 1.2 Stereo bonds (source: m1_set_x* (R)-alanine, stereo_x*)

| Parameter | ×1.0 | ×0.8 | ×0.6 | Method | Class |
|---|---|---|---|---|---|
| Wedge sides | linear from half-width LineWidth/2 at the stereocenter (side lines extrapolated to u = 0: ±0.30 at every scale) to 1.5 × BoldWidth / 2 at the drawn end | | | wedge polygon vertices in the bond frame | exact |
| Wedge narrow end | cut along the bisector with the neighboring plain bond on each side, meeting at the atom center (alanine: toward CH₃, not toward the hashed NH₂) | | | narrow-end vertices (0.062, 0.306), (0, 0), (0.298, −0.348) at ×1.0 | exact |
| Wedge wide-end width (wedge ends at an H label) | 3.000 | 2.400 | 1.802 | v extent of the polygon at max u; = 1.5 × BoldWidth | exact |
| Wedge drawn length (to H label) | 8.207 | 6.576 | 4.944 | u range | exact |
| Hash line thickness | 0.60 | 0.60 | 0.60 | perpendicular to each hash's long edge | exact, **exempt** |
| Hash count | 3 | 2 | 2 | polygons in the wedge corridor | exact |
| Hash pitch (center to center) | 2.56 | 3.98 | 2.855 | Δu of hash centers | exact |
| Hash first center u | 2.806 | 2.313 | 1.803 | | exact |
| Hash first near edge u | 2.505 | 2.013 | 1.504 | = HashSpacing | exact |
| Hash last far edge | flush with the drawn end (8.227 / 6.593 / 4.958) | | | | exact |
| Hash lengths | 0.854, 1.927, 3.000 | 0.836, 2.400 | 0.806, 1.801 | v extent per hash; last = 1.5 × BoldWidth | exact |
| Hash outline | each hash is a trapezoid slice of an outline whose half-width runs linearly from LineWidth/2 at u = HashSpacing to 1.5 × BoldWidth / 2 at the drawn end | | | hash vertices in the bond frame | exact |
| Hash count | floor(D / HashSpacing − 0.025) + 1, D = end − LineWidth − HashSpacing; hashes evenly spaced | | | 27 of 27 hashes (stereo_x*: 0.56–2.2 × bond length to CH₃, to OH/NH₂/Cl; m1_set_x*). The step from 2 to 3 lies between D/HashSpacing = 1.997 and 2.05; 0.025 puts the cut mid-bracket | exact |
| Wedge and hash end at a degree-1 unlabeled atom | wide end at the atom (u = L), flat, 1.5 × BoldWidth | | | stereo_x*, 6 lengths each | exact |
| Wedge narrow end at a degree-1 atom | flat, LineWidth wide | | | stereo_x* chain wedge | exact |
| Wedge wide end at an unlabeled atom with other bonds | sides keep their slope past the atom until they meet the nearer edge (±LineWidth/2) of the neighbor bond on their side; the end is cut back to the atom center (V notch). A side with no neighbor stops on the line through the other corner and the center. Runs 0.55 / 0.36 / 0.18 pt past the atom | | | stereo_x* (ring carbon, chain carbon), within 0.01 pt | exact |
| Bond beside a notched wedge end | cut along the wedge's end line: when it is the atom's only plain bond (butane) it starts 1.51 pt behind the atom, at the wedge's far corner, filling the notch; in rings the miter fills it | | | stereo_x1.0 (chain and ring carbon) | exact; we start the stroke at the far corner's projection (butt end inside the wedge) |
| Stereo bond into a label | clipped like a plain bond (§2.1) at its axis | | | stereo_x* (OH, NH₂, Cl), m1_set_x* | exact |

Hashes to unlabeled atoms with other bonds are not in the references (unmeasured): the last
hash ends at the atom.

## 2. Atom labels (source: labels_x1.0, m1_set_x*, mech_basic_x1.0, double_bonds_x1.0, clearance_x*)

| Parameter | ×1.0 | ×0.8 | ×0.6 | Method | Class |
|---|---|---|---|---|---|
| Label font / size | Helvetica 10 | 8 | 6 | SVG `<text>` font-size × matrix | exact |
| Explicit carbon labels (CH₃) | laid out and clipped like any label: C centered, H₃ after it; wedge and hash end flat at the label | | | stereo_x1.0 (wedge and hash into CH₃), within 0.03 pt | exact (terminal CH₃ only; CH₂/CH in chains unmeasured) |
| Horizontal anchor | first letter of the element symbol centered on the atom (advance center = ink center within 0.005) | | | glyph origin from SVG + fontTools advance | exact |
| Right-justified labels | when the bonds point right (the same test that puts H on the left), the **last** letter of the element is centered instead: –Cl has C on the atom, Cl– has l on the atom (Br likewise) | | | clearance_x*: Cl/Br at 12 directions | exact |
| Baseline below atom center | 3.586 (3.582–3.591) | 2.869 | 2.156 | glyph baseline − atom y; = 0.3585 × size = cap height / 2 (Helvetica capHeight 1469/2048 = 0.7173 em) | exact |
| H after element (NH₂, OH, NH₃, CH₃) | origin + advance(element) | | | dx of the H origin | exact |
| Left-side labels (H₂N, HO, H₃C) | set right-to-left, contiguous: H₂N gives H at −11.393, 2 at −4.171 from the N origin | | | labels_x1.0, mech_basic_x1.0 | exact |
| Subscript size | 0.75 × label size (7.5 / 6 / 4.5) | | | SVG font-size | exact |
| Subscript baseline drop | 2.30 | 1.84 | 1.38 | = 0.23 × size | exact |
| Terminal label on an offset double bond | moved bond_spacing/2 (1.30) perpendicular to the bond, toward the second line, so it sits centered between the lines; the lines don't move | | | double_bonds_2 acetaldehyde O at 90° and −30° (1.300, 1.377 perpendicular; ≤0.1 along) | exact (×1.0 only; assumed bond_spacing/2 at ×0.8/×0.6) |
| Subscript x | directly after the previous glyph's advance (no gap) | | | | exact |
| Kerning | none observed in labels (runs split at every size change) | | | | exact |

### 2.1 Bond-to-label clearance (source: all references; clearance_x*: H, I, Cl, Br, F, OH, NH₂ at 12 directions × 3 presets)

Every line of a multiple bond stops at the same axial distance (vertical C=O: both lines and a
single C–O end 5.519 from the O center at ×1.0).

**Rule (refit 2026-10-24, `measure_refs.py --fit-clip`).** Take the element glyph outline as
drawn, shifted by (−0.02, +0.03) pt. Keep only the parts within ±0.195 × label size of the bond
axis (edges clipped to that band). The bond stops 1.0075 × MarginWidth beyond the furthest
point of those parts along the bond.
- **Fit:** 396 of 399 bond ends are within 0.1 pt (worst 0.092, rms 0.055). These are all
  lines, wedges and hashes into labels in every reference, at every direction and preset.
- **Why a band:** this is the glyph's extent along the bond, limited to what lies near the
  bond's path.
  - **Round glyphs (O):** identical to a convex-hull ray model.
  - **Flat-sided glyphs:** a bond passing an H, N, I or F corner obliquely stops beyond the
    corner, so it gets the 1.2–1.75 × margin "outliers" of the earlier hull model.
  - **Narrow glyphs:** an I or l beside the path is ignored.
- **Known exception:** F with the bond leaving rightward. ChemDraw stops 0.31–0.41 pt further
  out, between F's middle bar and its longer top bar (just outside the band). No band shape
  fits this together with I at 30°. xfail in the acceptance test.
- **Noise floor:** nominally identical geometries in one file differ by up to 0.05 pt
  (m1_set_x0.8 OH at 30°: 4.199–4.268).
- **Superseded:** the M1 hull model (1.025 × margin, hull shift (−0.05, +0.08)), and the
  provisional box for oblique H/I/N. It mispredicted F, I and Br by up to 1.2 pt and Cl/Br
  bonds to the right by 4.7 pt (right-justification, §2).

| Measurement | ×1.0 | ×0.8 | ×0.6 | Source | Class |
|---|---|---|---|---|---|
| Bond end, vertical bond into O from below (C=O and C–O) | 5.518 (5.469 in labels_x1.0) | 4.433 | 3.300 | m1_set_x*, double_bonds, labels | exact |
| Bond end, 30° bond into O / OH | 5.27 (5.23–5.30) | 4.225 | 3.194 | m1_set_x*, labels | exact |
| Bond end, horizontal bond into N from the left (NH₂, NH₃⁺, O=N) | 4.50–4.55 | 3.659 | 2.719 | m1_set_x*, labels | exact |
| Bond end, explicit H at 30° / 60° | 5.907 / 6.099 | | | double_bonds_x1.0, resonance_x1.0, clearance_x* | exact |

### 2.2 Charges (filled bars; source: labels_x1.0, m1_set_x*, resonance_x*, mech_basic_x1.0)

| Parameter | ×1.0 | ×0.8 | ×0.6 | Method | Class |
|---|---|---|---|---|---|
| ⊕ / ⊖ bar length | 3.935 | 2.733 | 1.798 | polygon bbox | exact (not proportional to k). **Not followed:** by user decision (2026-10-24), presets scale it with the label size, 3.935 × k = 3.148 / 2.361 |
| Bar thickness | 0.50 | 0.50 | 0.50 | polygon bbox | exact, **exempt**, ≠ line width |
| On an unlabeled C (resonance σ-complex) | center 3.84 from the atom, on the bisector of the largest angular gap (33°, 147°, −91° for gaps bisected at 30°, 150°, −90°). **Default 4.0** | 3.24, at −0.4°, −179.6°, −90.5° (the side carbons are **not** on the 30°/150° bisectors). **Default 3.0** | 2.87, same angles as ×0.8. **Default 3.0** | resonance_x* | placement default (`charge.unlabeled_distance`, `nonlinear`: each preset from its own measurement, user decision 2026-10-25); the bisector is used at every preset (Open 15) |
| On a labeled atom with a bond within 60° of the upper right and none near the upper left (HO⁻ drawn O–H) | mirrored: center − element ink box top-left = (−x, y) of the row below | | | — | rule (user review 2026-10-24); unmeasured (`draw_mol.CHARGE_CORNER_DEG`) |
| On a labeled atom with free upper right (O⁻, N⁺) | center − element ink box top-right: (+1.98, −1.10) O⁻; (+1.55, −2.40) N⁺ ×1.0; (+2.05, −1.45) ×0.8; (+1.78, −1.33) ×0.6. Measured ×1.0 mean (+1.77, −1.75); × k = (+1.42, −1.40) / (+1.06, −1.05) (the measured ×0.8 / ×0.6 offsets belong to ChemDraw's smaller signs; see §7 item 6). **Default (+2.0, −2.0)** | (+1.5, −1.5) | (+1.0, −1.0) | labels_x1.0, m1_set_x* | placement default (`charge.labeled_offset_x/_y`; auto-placed, varies by context) |
| On a labeled atom with bonds near both upper corners (nitro N⁺ in a ring compound) | the first lower corner with no bond within 60° (lower right, then lower left; offsets mirrored down); if every corner has a bond within 60°: the first with none within 45° (m1_set_x*: N⁺ with bonds at 0° and 180° keeps the upper right, as ChemDraw); else on the largest gap's bisector, the corner offset's length beyond the element ink box | | | — | rule (M4, nitration); unmeasured (`draw_mol.CHARGE_CORNER_MIN_DEG`, Open 16) |
| On an unlabeled C that draws a lone pair or radical (ring carbanion with an arrow from its pair) | the dots keep the largest-gap bisectors (§4.2; the charge does not count as a direction); the charge takes the direction inside the largest bond gap, every 5°, whose ink clears the dots and the atom's drawn bond lines the most, at `unlabeled_distance` or 0.5 pt steps beyond it (up to +2.0) until some direction clears by 1.0; ties go nearest the bisector | | | — | rule (inbox 2026-10-26); unmeasured (`draw_mol.CHARGE_DOT_CLEARANCE`, `CHARGE_DOT_STEP_*`, Open 16) |
| On a labeled atom that draws lone pairs or radicals (HO⁻, Br⁻, bisulfate O⁻ in m3_sn2, m4_nitration) | the corner is chosen as above; the sign then moves out along that corner's diagonal (straight up for NH₃⁺ / single-atom HO⁻) in 0.5 pt steps, up to +2.0, until its ink clears the dots and the atom's drawn bond lines by 1.0 (`CHARGE_DOT_CLEARANCE`); a charge without dots keeps its corner | | | — | rule (user decision 2026-10-26); unmeasured (`draw_mol.Drawer._charge_clear`, Open 18) |
| Circled ⊕ / ⊖ (`charges.circled`, `circled_charges`) | the plain bars plus a ring of `circle_radius`, stroked at `bar_thickness`. **Default 3.0** | 2.5 | 2.0 | — | unmeasured placement default (user decision 2026-10-26, Open 18). The starting position moves out by the ring's extra reach beyond a plain sign; a circled sign always runs the clearance search (labeled: along the corner diagonal; unlabeled: the largest-gap search above) |
| Multiple charges ("2+", "2−") and `charge_labels` texts ("+3") | digits in the label font at `digit_scale` × label size (**0.75**, = `ratios.superscript_scale`), ink `digit_gap` (**0.5**) from the bars, in the written order, the whole centered on the sign position and vertically on the bars; never circled. Positions move out by the extra width / height and run the clearance search | | | — | unmeasured (user decision 2026-10-26, Open 18) |
| NH₃⁺ (H on the right) | ⊕ centered above the N: center − atom = (+0.03, −7.43). **Default 7.5** | 6.0 | 4.5 | labels_x1.0 | placement default (`charge.above_label`) |

## 3. Schemes (source: scheme_basic_x1.0 / x0.8 / x0.6, scheme_multi_x1.0 / x0.8 / x0.6)

Method: `uv run python scripts/measure_refs.py --scheme` (per preset, medians over both
references; `_scheme` assigns each caption to an arrow by horizontal overlap, then nearest in y,
and groups structures into rows). Exact text positions come from ChemDraw's SVG run origins and,
within runs, its PDF glyph origins (`measure_refs.pdf_char_origins`).

### 3.1 Arrows (exact)

| Parameter | ×1.0 | ×0.8 | ×0.6 | Method | Class |
|---|---|---|---|---|---|
| Reaction arrow shaft width | 0.60 | 0.60 | 0.60 | shaft polygon | exact, **exempt** |
| Arrowhead length (tip → back corners) | 6.00–6.01 | 6.00 | 5.98–6.00 | fill path in the arrow frame | exact, **exempt** |
| Arrowhead half-width | 1.58 (1.55 / 1.60, ChemDraw is asymmetric) | 1.575 | 1.569–1.575 | | exact, **exempt** |
| Arrowhead inset (tip → notch on axis) | 5.25–5.26 | 5.25 | 5.23–5.25 | | exact, **exempt** |
| Arrowhead back edge | cubic corner → notch, controls at the corner and at (notch, ±0.4375 × corner offset) | same | same | path data, all 16 heads | exact (`ratios.arrow_back_control_per_half_width`) |
| Arrowhead hairline outline | 0.05 stroke over the fill | | | | exact |
| Shaft end | 5.20–5.21 from the tip (runs into the notch) | | | | exact |
| CDXML arrow attributes | HeadSize 1000, ArrowheadCenterSize 875, ArrowheadWidth 250 (= 10, 8.75, 2.5 × LineWidth) | | | | exact |
| Arrow bounding-box center | 0.25 / 0.20 / 0.155 above the axis (the asymmetric head) | | | CDXML BoundingBox | exact |

Encoding: the head is an SVG `<marker>` on the shaft's end (`render.py --affinity` expands it into
a separate filled path, as above). The geometry is the same either way.

### 3.2 Text (exact)

| Parameter | Value | Method | Class |
|---|---|---|---|
| Conditions and numbers size | = caption size (10 / 8 / 6) | | exact |
| Runs | ChemDraw starts a new SVG run at every word, size and style change; each starts at the previous origin + advance (no kerning across a boundary) | run origins, all captions, 0.1 pt | exact |
| Kerning within a run | applied: bold "11" steps 5.562 − 0.547 (×k); every other measured pair matches Helvetica's kern table | PDF glyph origins | exact (closes Open 6) |
| Subscript (K₂CO₃, Pd(PPh₃)₄, CH₂Cl₂, H₂O) | 0.75 × size, baseline dropped 0.23 × size | | exact |
| Superscript (H₃O⁺) | 0.75 × size, baseline raised 0.3586 × size (3.586 / 2.869 / 2.152) | | exact |
| Italic (*n*-BuLi) | separate run, `font-family="Helvetica Oblique"` `font-style="italic"` (CDXML face 2); n advance 5.562, no kerning across | scheme_multi_x1.0 only: the ×0.8/×0.6 files have an upright n | exact |
| Line pitch, no subscripts | size + 0.05 (10.05 / 8.05 / 6.05) = hhea ascent + descent (1.0 em) + 0.05 | | exact |
| Line pitch, subscript on the upper line | + 0.775 / 0.575 / 0.433 | 1,4-dioxane/H₂O (4:1) → 90 °C, 12 h | approx (one pitch per preset; Open 9) |
| … and on the lower line too | + a further 0.15 / 0.15 / 0.10 | PhB(OH)₂ → Pd(PPh₃)₄ | approx (Open 9) |
| Compound number | Helvetica Bold at caption size | | exact font |

### 3.3 Placement (hand-placed in the references; round defaults)

The references' positions here were set by eye in ChemDraw, so their medians are provenance,
not rules: presets use round defaults on a 0.5 pt grid, each preset rounded from its own
measurement (user decision 2026-10-24; preset `placement:` list). Structures and '+' sit on
the arrow axis: ChemDraw's 0.25 offset is the center of its asymmetric arrowhead's bounding
box, and our head is symmetric.

| Parameter | Default ×1.0 / ×0.8 / ×0.6 | Measured ×1.0 | ×0.8 | ×0.6 | Range at ×1.0 (n) | Preset key |
|---|---|---|---|---|---|---|
| Arrow length = max(min, widest line + 2 pad + head inset) | min 40 / 31.5 / 23.5 | 39.5 | 31.6 | 23.72 | 39.5–112.5 (6) | `arrow_min_length` (shortest reference arrow) |
| Conditions pad (tail → first origin; end of advance → notch, widest line) | 2.5 / 2.0 / 1.5 | 2.33 | (1.505) | (0.752) | 0.82–6.04 (10) | `conditions_pad`; ×0.8/×0.6 from ×1.0 × k (their text overhangs: ChemDraw scaled text about its anchor) |
| Conditions horizontal | advance center on the tail → notch midpoint | median offset −0.28 (advance), −0.13 (ink) | | | −1.28 to +0.58 (13) | rule |
| Multistep stages ("1." above, "2." below): ink left of each line after the tail | 3.0 / 2.5 / 2.0 | 3.195 | 2.553 | 1.921 | 3.192 / 3.197 | `multistep_indent` |
| Above: axis − last baseline | 5.5 / 4.5 / 3.5 | 5.444 | 4.361 | 3.278 | 3.89–7.18 (6) | `conditions_above_baseline` |
| Below: first baseline − axis | 11.5 / 9.5 / 7.0 | 11.638 | 9.312 | 6.982 | 11.36–12.61 (6) | `conditions_below_baseline` |
| Structure vertical: ring centers (rows where every structure has a ring; else ink centers) above the axis | 0 / 0 / 0 | 0.25 | 0.20 | 0.155 | ring −0.25 (19); ink −0.25 to −0.30 (3) | `center_above_axis` |
| Plus sign: text "+", regular, caption size; ink center above the axis | 0 / 0 / 0 | 0.25 | | | −0.251 to −0.254 | `center_above_axis` |
| Gap, structure ink ↔ arrow ends | 13.0 / 10.5 / 8.0 | 13.11 | 10.47 | 7.82 | 10.26–22.03 (12) | `arrow_gap` |
| Gap, structure ink ↔ plus ink | 16.0 / 13.0 / 9.5 | 15.96 | 12.78 | 9.59 | 2.26–18.27 (4) | `plus_gap` |
| Compound numbers: ink-centered under the structure (≤ 0.006, one at 0.5); one baseline per row; row's lowest structure ink → tallest number's ink top | 6.5 / 5.5 / 4.0 | 6.709 | 5.314 | 3.945 | 6.39–6.87 (11) | `number_gap` |

Notes:
- Above/below placement by **baseline**, not ink: row 1 of scheme_multi has K₂CO₃, LiAlH₄ and PCC
  on one baseline although their ink bottoms differ by 2.3 pt.
- Formula text typed in ChemDraw (K₂CO₃, LiAlH₄, PCC, CH₂Cl₂) is stored as a one-node fragment;
  `measure_refs` treats it as a caption (`Doc.text_frags`).
- Abbreviations (OMe, OBn) are `NodeType="Fragment"` nicknames: one run, left-justified, first
  letter centered on the atom, on the label baseline. Their bonds fit the §2.1 band model on
  the O outline (worst 0.08 pt, all presets).

### 3.4 Multi-row schemes (source: scheme_wrap_x1.0 / x0.8 / x0.6)

Method: `uv run python scripts/measure_refs.py --wrap`. Each structure, arrow and caption joins
the row of the nearest arrow axis (compound numbers: the row above them). A row's ink is the
union of its structure ink, arrow boxes (axis ± head half-width) and caption ink. The drawing
is 1 + 2 → 3 → 4 / → 5 → 6 → 7 (scheme_multi row 1 plus a Wittig and a hydrogenation). Every
distance in the ×0.8 / ×0.6 files is the ×1.0 value × k (ratios 0.8005 / 0.6005), so this is
**one hand placement**, scaled by ChemDraw.

| Parameter | ×1.0 | ×0.8 | ×0.6 | Class |
|---|---|---|---|---|
| Break | row 1 ends on a structure (4); row 2 starts with a bare arrow (→ 5); 4 is not repeated | same | same | rule (user decision 2026-10-25): break only between steps; every later row starts with a reaction arrow, so the scheme visibly continues |
| Row alignment | left: row 2's first arrow tail is at row 1's left ink edge (−0.004) | −0.007 | −0.012 | rule (user drawing) |
| Row widths (ink) | 360.1 / 377.2 (row 2 wider by 17.0) | 288.3 / 301.9 | 216.3 / 226.5 | (follows from the content) |
| Row gap: row ink bottom (incl. numbers) → next row ink top (incl. conditions) | 34.73. **Default 34.5** | 27.80. **Default 28.0** | 20.86. **Default 21.0** | placement default (user decision 2026-10-25; `scheme.row_gap`); hand-placed, n = 1 |
| Row gap: structure ink → structure ink | 48.41 | 38.76 | 29.07 | approx (here row 2's top is structure 6, not conditions) |
| Axis pitch | 142.37 | 113.97 | 85.49 | approx (depends on the structures) |
| First arrow of row 2 | 39.5, PCC / CH₂Cl₂ (= `scheme_multi` and `arrow_min_length`) | 31.63 | 23.73 | ordinary arrow-length rule |
| Overall ink width | 377.2 | 301.9 | 226.5 | (follows from the content) |
| Width limit that triggers a wrap | **504 pt** (7 in, ACS double column) at every preset; the user can ask for another width | | | user decision 2026-10-25 (`scheme.wrap_width`, `assumed_exempt`: not k-scaled); `layout.width` in a spec overrides |

Notes:
- Arrow ↔ structure gaps in this file: row 1 10.27 (n = 4), row 2 12.48 (n = 5). Both are within
  the §3.3 range; `scheme.arrow_gap` is unchanged.
- Conditions sit as in §3.3 (above 4.50–5.75, below 11.36–12.50); compound numbers share one
  baseline per row (3 is 0.18 pt high, by hand).

## 4. Mechanisms and resonance (source: mech_basic_x1.0 / x0.8 / x0.6, resonance_x1.0 / x0.8 / x0.6)

Method: `uv run python scripts/measure_refs.py --mech` (per preset). Dots are ChemDraw's
circle paths. A lone pair is two mutually nearest dots; a single dot is a radical. Its atom is
the nearest atom center, and gaps are measured to that atom's element ink box (fontTools).
Curved arrows are matched to their CDXML `curve` (CurvePoints = in-handle, anchor, out-handle
per anchor) and measured in the tangent frame at the stroke end. mech_basic draws HO⁻ + CH₃Br
with one 2e⁻ arrow (O lone pair → C of CH₃), and tBu–Cl homolysis with two fishhooks from the
C–Cl bond.

### 4.1 Exact (ChemDraw-generated)

| Parameter | ×1.0 | ×0.8 | ×0.6 | Class |
|---|---|---|---|---|
| Lone-pair / radical dot radius | 0.833 | 0.533 | 0.533 | exact; **not k-scaled** (ChemDraw's ×0.8 dots are the ×0.6 size); copied exactly by user decision (`nonlinear`) |
| Lone-pair dot spacing (center to center) | 3.75 | 2.40 | 2.40 | exact, `nonlinear` |
| Curved-arrow stroke | 0.60, stroked cubic Bézier | 0.60 | 0.60 | exact, **exempt** |
| Full head: length / half-width / inset; straight 4-vertex polygon | 4.81 / 1.50 / 4.06 | 4.80 / 1.50 / 4.05 | 4.78 / 1.49 / 4.03 | exact, **exempt** (presets use the ×1.0 values) |
| Stroke end | at the head's notch; the head axis is the stroke's end tangent | same | same | exact |
| Fishhook head: 3 vertices, tip, barb at (−4.81, 1.755), notch at −4.06; barb on the outside of the curve's end | barb 1.755 | 1.739 | 1.734 | exact, **exempt** |
| Curve shape scales with the document | control offset / chord 0.531 (2e⁻), 1.307–1.312 (fishhooks); chord 53.25 / 42.60 / 31.30 (2e⁻), 5.04 / 4.03 / 3.02 (fishhooks) | | | the hand-drawn curve is scaled by ChemDraw, the heads are not |
| Resonance arrow | double-headed, both heads identical to the reaction arrowhead (6.02 / 1.58 / 5.26); shaft runs between the notches; tip to tip 28.80 (= 2 bond lengths) | heads 6.00 / 1.58 / 5.25; 23.04 | heads 5.98 / 1.57 / 5.23; 17.28 | heads exact, **exempt**; length hand-drawn (scaled with the document) |
| Square bracket stroke | 0.60 (stroked, 3 segments) | 0.60 | 0.60 | exact, **exempt** |
| Bracket lip length | 4.007 | 3.212 | 2.417 | exact, scaled by k (× k: 3.206 / 2.404) |
| Bracket → content ink (x gap; above; below) | 5.82; 2.11; 2.20 | 4.63 (4.60–4.66); 1.69; 1.81 | 3.43 (3.37–3.50); 1.26; 1.04 | approx (hand-placed); placement defaults in §4.3 |

Encoding: curved-arrow and resonance heads are SVG `<marker>`s on the stroke (`render.py
--affinity` expands them into separate filled shapes). The geometry is the same either way.

### 4.2 Placement (hand-placed in the references; round defaults)

Each preset is rounded to the 0.5 pt grid from its own measurement (user decision 2026-10-24).
Curvature is dimensionless and gets a round ratio, the same at every preset.

| Parameter | Default ×1.0 / ×0.8 / ×0.6 | Measured ×1.0 | ×0.8 | ×0.6 | Preset key |
|---|---|---|---|---|---|
| Lone-pair gap, dot edge → element ink box | 2.0 / 2.5 / 1.0 | 2.14 (1.60–2.60, n = 10) | 2.29 (1.65–2.51) | 1.06 (0.91–1.35) | `lone_pair.gap` (not k-scaled, like the dots) |
| Radical dot center from an unlabeled atom (tert-butyl), on a largest-gap bisector of the bonds and the sets already placed (a charge on the atom is placed after the dots, §2.2; a lone pair's midpoint sits here too, unmeasured) | 3.0 / 2.5 / 2.0 | 2.79 (31°) | 2.48 (−89°) | 1.80 (−88°) | `radical.distance` |
| Radical on a labeled atom (Cl·): a free side, like a lone pair | `lone_pair.gap` | 1.94 | 1.36 | 1.28 | (rule) |
| Tail from a lone pair: beyond the dots' outer edge, outward | 2.5 / 1.5 / 1.5 | 2.61 | 1.37 | 1.41 | `curved_arrow.tail_gap` |
| Tail from a bond: off the bond axis at the midpoint of the *drawn* bond (a label at either end shortens it; user decision 2026-10-26, unmeasured) | 2.5 / 2.0 / 1.5 | 2.12, 2.82 (tails 2.0–2.5 along the bond from the midpoint) | 1.70, 2.25 | 1.28, 1.69 | `curved_arrow.bond_tail_gap` |
| Two arrows from one bond; a tail beside a double bond's second line | Opposite sides preferred. On the inner side of an offset (ring) double bond the tail starts `bond_tail_gap` beyond the second line. On the same side (used only when the opposite side clashes, or by `bow`), the second tail also moves `bond_tail_gap` along the bond toward its target | | | | rule (user decision 2026-10-26; unmeasured) |
| Head tip → target element ink box (unlabeled atom: its center); the nearest distance, to the box or one of the atom's dots, whichever is nearer | 2.5 / 1.5 / 1.5 | 2.59 (2e⁻ → CH₃); fishhooks: Cl 1.70, C 2.89 | 1.31; 1.01, 1.58 | 1.30; 0.31, 0.56 | `curved_arrow.head_gap` |
| Arc: control offset / chord | 0.5 | 0.531 | 0.531 | 0.531 | `curved_arrow.curvature` |
| Loop (chord shorter than a bond): control offset / chord | 1.3 | 1.311 | 1.310 | 1.310 | `curved_arrow.loop_curvature` |
| Arc handles (in the chord frame, along / offset ÷ chord) | c1 (0.258, 0.987 h), c2 (0.874, 0.770 h): user decision 2026-10-26, smoothed. ChemDraw's tangent angles kept (63.8° and 72.8° off the chord at h 0.531), handle lengths refit (0.584, 0.428 of the chord; ChemDraw 0.592, 0.268) so curvature runs evenly to the head (at h 0.5: 1.3 → 1.8 → 2.2 per chord) instead of dipping mid-arc and spiking at the head (1.4 → 2.1 → 1.5 → 5.1). The arc is 16% taller at one chord (0.33 vs 0.29 of it) | (0.261, 0.531), (0.921, 0.256) | same | (0.288, 0.531), (0.930, 0.264) | `mechanism.ARC_HANDLES` (shape) |
| Loop handles | c1 (1.19, h), c2 (1.45, 0.59 h) | (1.194, 1.311), (1.446, 0.769) | same | same | `mechanism.LOOP_HANDLES` (shape) |
| Hook: an arrow onto a bond next to its source (a bond sharing one atom; lp(n) onto a bond on n; E2's C–H → C–C, a sigma complex's C=C → C–C⁺) | leaves along the source bond's normal on the bow side (lp: the pair's direction) and lands along the target bond's normal (±10°, ±20° tried); both handles the circular arc's for that turn, 4/3 tan(θ/4) R, at least 0.35 × bond_length; × 0.75, × 0.5 tried. Chosen from the endpoints, not the chord length. A `curvature` override keeps the arc or loop. When the tail lies inside the angle between the two bonds (C–H → C–C forming a π bond), there's no room on that side: the arrow passes through the target bond and curls back, its tip head_gap past the axis on the far side and its head pointing back at the bond, turned 30° toward the shared atom so the loop opens away from it; handles at least 0.9 × bond_length; that crossing isn't a clash (user choice 2026-10-26, option D of six rendered). The tail's side and the far side without the curl are the fallbacks (phenoxy's C=C → C–O, where the radical's fishhook meets it). Between two bonds of one ring (a sigma complex's C=C → C–C⁺) the hook bows into the ring: leaving one bond's inward normal and landing across the next is a 120° turn, drawn as a circular arc (no 0.35 × bond_length minimum, which made it pointed), with both ends slid 0.25 × bond_length along their bonds away from the shared atom so the fixed-size head leaves a stroke at ×0.6; an explicit `bow:` still wins, and the outside stays a fallback (user decision 2026-10-27) | | | | rule (user decisions 2026-10-26, 2026-10-27); unmeasured (`mechanism.HOOK_SCAN`, `HOOK_REACH`, `hook_controls`, `Planner.cramped`, `CROSS_REACH`, `CROSS_TILT`, `MechArrow.crosses`, `Planner.ring_center`, `RING_SLIDE`) |
| Lone-pair midpoint distance from the element ink-box center | 6.2–8.0 | | 5.8–5.9 | 3.5–4.5 | (follows from the gap) |
| Charge on O⁻ with lone pairs: center − element ink box top-right | (+3.85, −2.02), moved right past the top pair | (+1.48, −1.25) | (+1.35, −1.18) | | followed in spirit: the §2.2 corner, pushed clear of the dots (§2.2 clearance rule) |
| Lone pairs on an arrow's target atom | in their usual places; the arrow works around them. The tip keeps head_gap from the dots too, so it ends above a pair, never between its dots (user review 2026-10-24: "end just above"). The head must point into the element ink box without passing through a dot, within 20° of the box's center (`MISS_LABEL_DEG`, inbox 2026-10-26: water's O–H → O pointed at the lower pair); an exact aim at the center is preferred, and among heads that miss, the closest. The planner also tries the tip 1.5 and 3.0 pt farther out and the curvature × 0.75 and × 0.5, and flips a bow only to remove a clash | | | | rule (user review 2026-10-24: "the arrows should be decorations") |
| 2e⁻ arrow from a bond to its own atom (C–Br → Br, hand-drawn) | tail off the bond axis 5.00, tip → Br ink box 1.11; chord 8.73, control/chord 1.34, its c2 on the end point | 3.08, 1.23; 8.85, 0.91 | 4.71, 1.05; 5.87, 1.54 | mech_basic_x* | not followed: inconsistent between presets; `bond_tail_gap`, `head_gap` and the loop shape kept (a hook was tried and dropped, 2026-10-26: lessons.md). Head exact (§4.1, tested) |
| Fishhook to `between(a,b)` (bond forming) | tip head_gap from the midpoint of a and b; with a plus sign there, from a point head_gap past its ink box on the bow side; arrows meeting there share a side | | | | rule (user review 2026-10-24); unmeasured |

Notes:
- mech_basic shows all lone pairs on HO⁻, Br⁻ and Cl·, and none on the neutral molecules.
- The 2e⁻ arrow leaves the top lone pair of O and bows upward over the '+'; it ends pointing
  down onto the C of CH₃. The two fishhooks each start beside the C–Cl bond, one above and one
  below. Each loops back onto one atom (C, Cl), with the barbs on the outside of the loop.

### 4.3 Resonance placement (resonance_x*; hand-placed; round defaults)

Method: `uv run python scripts/measure_refs.py --resonance` (`resonance_summary`). Structure ink
here includes each contributor's charge bars (a charge joins the structure of its nearest
atom): without them the σc gap reads 12.56, because σc's charge sits on its left, facing the
arrow. With them every gap at ×1.0 is 7.68. Bracket gaps are to the stroke center line
(`_brackets`), from the union of all content ink (bonds, charges). resonance_x* draws the
three σ-complex contributors of benzene + E⁺; ChemDraw scaled one hand placement by k.

| Parameter | Default ×1.0 / ×0.8 / ×0.6 | Measured ×1.0 | ×0.8 | ×0.6 | Preset key |
|---|---|---|---|---|---|
| ↔ arrow length, tip to tip | 29.0 / 23.0 / 17.5 | 28.80 (= 2 bond lengths) | 23.04 | 17.28 | `resonance.arrow_length` |
| Contributor ink (incl. charges) ↔ arrow tip | 7.5 / 6.0 / 4.5 | 7.68 (7.680–7.684, n = 4) | 6.11 (5.68–6.15) | 4.55 (3.94–4.61) | `resonance.arrow_gap` |
| Contributor vertical: ring centers on the arrow axis | 0 | −0.002 (n = 3) | −0.002 | −0.002 | `scheme.center_above_axis` (the §3.3 rule) |
| Bracket stroke center ↔ content ink, x | 6.0 / 4.5 / 3.5 | 5.82 (n = 2) | 4.63 (4.60–4.66) | 3.43 (3.37–3.50) | `bracket.gap_x` |
| Bracket top above the content ink | 2.0 / 1.5 / 1.5 | 2.105 | 1.687 | 1.263 | `bracket.gap_above` |
| Bracket bottom below the content ink | 2.0 / 2.0 / 1.0 | 2.204 | 1.808 | 1.042 | `bracket.gap_below` |

Rules (M4; unmeasured, user decisions 2026-10-25, flagged at review):
- The bracket encloses the union of every contributor's ink, **including its lone-pair dots
  and curved arrows** (resonance_x* has neither), plus the gaps above. Both brackets share
  the top and bottom.
- The reaction arrow into a bracket ends `scheme.arrow_gap` before the bracket ink; the next
  step's arrow starts `scheme.arrow_gap` after it.
- Reactants a mechanism step adds (the base of a deprotonation) sit above its arrow, centered
  on tail → notch, their ink bottom `conditions_above_baseline` above the axis (or above the
  step's conditions); the arrow grows to fit them like a conditions line (user review
  2026-10-25: "structures being added to the mechanism should come in above an arrow"). A
  mechanism step's arrow stays on its reactants' row, so that row may end with the arrow and
  the next start with the products.

## 5. Fonts and rasterizer

- **System Helvetica:** `/System/Library/Fonts/Helvetica.ttc` (macOS 26.6.2) contains 6 faces:
  #0 Regular, #1 Bold, #2 Oblique, #3 Bold Oblique, #4 Light, #5 Light Oblique (family
  "Helvetica"). `fontmetrics.py` selects them by name-table family and subfamily.
- **Fallback font:** TeX Gyre Heros 2.004 (GUST Font License), `fonts/texgyreheros-{regular,bold,italic}.otf`
  (italic added in M2, from CTAN, for *n*-BuLi-style runs; Helvetica Oblique is face #2).
  - **Advances:** the max |Δ advance| vs Helvetica at 10 pt is **0.002 pt**, Regular and Bold,
    for C H N O S P F l B r I 0–9 + − ( ). All characters exist in both fonts.
  - **Kerning differs:** Helvetica's legacy `kern` table kerns "11" by −0.737 pt (Bold −0.547).
    Heros GPOS kerns Fr −0.35, FO −0.20, CH −0.10, CO −0.10 (Regular). The max |Δ kern| is
    0.737 pt at "11". This matters for multi-digit compound numbers and condensed text. It
    doesn't matter for ChemDraw-style labels, which split runs at every subscript.
  - **Verdict:** use the fallback only when Helvetica is missing; Nimbus Sans was not needed.
  - Reproduce with `uv run python scripts/compare_fonts.py`.
- **Rasterizer: resvg 0.48.1 via resvg-py 0.5.0** (`acs_chem_svg.preview.rasterize`).
  - **Font database:** the rasterizer loads only explicit font files (`skip_system_fonts`).
  - **Pt workaround:** resvg-py rejects root `width`/`height` in `pt`, so `rasterize()`
    rewrites the root size to pixels from the pt viewBox.
  - **Kerning:** resvg applies both GPOS kerning and the legacy `kern` table, and ignores
    `font-kerning="none"`. `fontmetrics` therefore applies the same pair kerning.
  - **Acceptance (`scripts/check_rasterizer.py`):** ink boxes of 11 test runs at 10 pt, 20 px/pt,
    match the fontTools predictions within **0.006 pt**. This holds for Helvetica Regular and
    Bold and for Heros Regular and Bold (tolerance 0.1 pt). M2 adds Helvetica Oblique and Heros
    Italic (worst 0.012 pt), and Oblique through the default font list: `font-style="italic"` with
    our family list picks Helvetica Oblique, not Heros Italic.
  - **Bold and italic encoding depends on the editor.** Affinity imports bold and italic only
    from a bare face name in font-family (user checks: 2026-10-24, `Helvetica Bold` and
    `Helvetica Oblique` work, a family list with `font-weight="bold"` / `font-style="italic"`
    imports as regular; 2026-10-26, the face name first in the family list also imports as
    regular). ChemDraw writes `Helvetica Oblique` too. Inkscape and Chrome can't resolve a face
    name and substitute another font (2026-10-26; Safari resolves it). So the default output
    writes the family list with the attributes, and `--affinity` writes the bare face name.
    resvg matches family names only (the TTC family is `Helvetica`), so `preview.rasterize` maps
    the two face names to the family list; the attributes pick the face. Both render alike.
  - **Negative control:** with no fonts loaded, nothing is drawn, so resvg never substitutes a
    font silently. PASS.
  - cairosvg is not usable (libcairo is not installed), and it isn't needed.

## 6. [VERIFY] items

### Resolved
- **§3.2 preset values:** ChemDraw's scaled Document Settings match the design-brief table exactly,
  and line width stays 0.60.
- **§3.2 exemptions:** the ×0.8 and ×0.6 references show these are **not** scaled: line width,
  hash line thickness, arrow shaft width, arrowhead (length, width, inset) and charge bar
  thickness. Charge bar length scales, but not linearly.
- **§3.3 ring inner-line shortening:** s·tan 30° per end (1.496 / 1.197 / 0.898), with butt ends.
- **§3.3 acyclic double-bond rule:** the cases are in §1.1 (C=O is centered when the carbon
  has two substituents; terminal C=CH₂ follows the substituent count).
- **§3.3 wedge wide end:** 1.5 × BoldWidth (3.0 pt at ×1.0) at the drawn end next to a label.
- **§3.3 reaction arrowhead dimensions:** 6.0 × 1.58 (half-width) × 5.25 (inset), exempt.
- **§4.3 Helvetica path and face names:** confirmed on macOS 26.6.2.
- **§4.3 fallback advance deviation:** 0.002 pt max; the kerning difference is noted above.
- **§7 rasterizer:** resvg honors the chosen font file and renders `<text>` within 0.006 pt.
- **§4.3 Affinity import:** confirmed by the user (2026-10-24, M1 SVGs at all three presets).
  Text imports as editable Helvetica at the same position, using `x`/`y` baselines with no
  `dominant-baseline`. The `species-…` groups appear as named layers, bonds are editable
  strokes, and nothing moves. These checks apply to `--affinity` output: Affinity drops `<marker>`
  arrowheads on import and needs bare face names for bold and italic (2026-10-26, lessons.md).

### Resolved by stereo_x*, clearance_x*, double_bonds_2 (2026-10-24)
- **Hash count and spacing** (§1.2): measured on 27 hashes.
- **Wedge and hash wide end at unlabeled atoms** (§1.2): flat at degree 1, notched at
  degree 2–3.
- **Clearance shape** (§2.1): one band model for all elements and directions. Cl/Br
  right-justification (§2).
- **Terminal O on an offset C=O** (§2): centered between the two lines (first misread as a
  hand-moved label).
- **Double-bond side** (§1.1): E-alkene ties, aldehyde C=O, ring bonds with exocyclic
  substituents, fused bonds, and the 108° inset.

### Open (not measurable from the supplied files)
1. **F with the bond leaving rightward** (§2.1): 0.31–0.41 pt further out than the band model.
2. **Hashes to unlabeled atoms with other bonds**: not in the references.
3. **Label-to-label bonds** at angles other than 0° (clearance_x* has HO–OH and H₂N–NH₂ at
   0–90°). They are covered by the band model, but only as single-bond cases.
4. **Charge placement rule** on labeled atoms (auto-placed, context-dependent). Presets use
   round default offsets, (+2.0, −2.0) from the element ink box top-right at ×1.0 (§2.2).
5. ~~**Scaling of mechanism elements.**~~ **Resolved (M3, mech_basic_x*):** curved arrows are exempt; dots are
   0.833 / 0.533 / 0.533 (not proportional to k). **Brackets (resonance_x*, 2026-10-24):** stroke
   exempt (0.60), lips scaled by k, ↔ heads exempt.
6. ~~Bold compound numbers ≥ 10.~~ **Resolved (M2, scheme_multi_x*):** ChemDraw kerns "11" by
   −0.547 × k (its PDF glyph origins), as fontmetrics and resvg do.
7. **Terminal-label shift on offset double bonds** (§2) is measured only for O at ×1.0;
   labels with H (=NH) and the ×0.8/×0.6 presets are assumed to follow bond_spacing/2.
8. **Explicit labels on internal carbons:** where ChemDraw puts the H of a CH₂/CH with bonds
   on both sides (our general rule runs a bond through the H run), and CH₃ at ×0.8/×0.6, are
   not in the references. Internal carbon labels raise NotImplementedError until measured.
9. **Line pitch with subscripts** (§3.2): one two-line pitch per case per preset. The fitted
   increments (upper line subscripted; lower line too) fit both, but a reference with a
   subscript on the lower line only, and lines with parentheses but no subscripts, would tell
   whether the rule is about subscripts or about ink.
10. **Right-justified abbreviations** (MeO–, BnO–: bond leaving the label to the right): not
    in the references; `draw_mol.layout_abbreviation` raises NotImplementedError. Abbreviations
    other than OMe/OBn (CO₂Me, Ph, …) follow the OMe rule unmeasured.
11. **Multistep arrows with more than 2 stages** (where "3." goes): not in the references; the
    spec allows at most 2.
12. **Italic at ×0.8 / ×0.6:** the scaled scheme_multi files have an upright n in n-BuLi, so
    the italic run is measured at ×1.0 only.
13. **Wrapping and width limits** (multi-row schemes, conditions wider than a page column).
    **Resolved (scheme_wrap_x*, §3.4; user decisions 2026-10-25):** break only between steps,
    every later row starts with a reaction arrow, rows left-aligned, row gap 34.5 / 28.0 / 21.0,
    wrap at 504 pt (7 in) unless the user asks for another width. Still open: conditions wider
    than the width limit.
14. ~~**Placement medians** (§3.3) pool hand placement from two references; arrow gaps range
    10.3–22.0 pt.~~ **Resolved by user decision (2026-10-24):** hand-placed positions are not
    data; presets use round defaults (0.5 pt grid), with the medians kept as provenance.
15. **Charge direction on an unlabeled C at ×0.8 / ×0.6** (§2.2): resonance_x0.8/x0.6 put the
    side-carbon ⊕ at 0° / 180° (3.24 / 2.87 pt), not on the largest-gap bisector used at ×1.0
    (33° / 147°, 3.84 pt). Whether ChemDraw re-placed them on rescaling or by its own rule at
    smaller sizes is unknown. By user decision (2026-10-25) `draw_mol` keeps the bisector at every
    preset. A σ-complex drawn fresh in a ×0.6 document would settle it.

16. **M4 rules not in any reference** (user decisions or defaults, flagged at the M4 review):
    the charge on a crowded labeled atom (nitro N⁺ next to a ring: lower corner or largest gap,
    §2.2); the charge on an unlabeled atom that draws dots (ring carbanion: the pair on the
    bisector, the charge in the clearest direction of the largest gap, §2.2); brackets enclosing the contributors' curved arrows and dots, and ↔ gaps measured
    from contributor ink including its curved arrows (§4.3); structures added above a mechanism arrow
    (reusing `conditions_above_baseline`, §4.3); a lone pair that slides along its side, at most half the side,
    off a bond (×0.6 crowding); the vertical rule per row in mixed rows (ink centers when a row
    has a structure without a ring, so σ-complex rings sit off the ↔ axis next to HSO₄⁻). A
    reference with a mechanism inside a bracket (curved arrows on contributors) would settle
    the bracket and gap rules.

17. **Lewis structures** (§12): there is no ChemDraw Lewis reference. The lattice layout,
    `lewis.bond_length` and the shared/assigned dot positions are unmeasured user decisions
    (2026-10-26). A ChemDraw Lewis reference (CCl₄, CH₃OH, CO₂, HCN at three presets) would let
    them be measured.

18. **Charge clearance, circled signs, multiple charges** (§2.2): no reference has a circled
    ⊕ / ⊖ or a "2+", and the push clear of lone pairs is a rule, not a measurement
    (`charge.circle_radius`, `digit_scale`, `digit_gap`, `draw_mol.CHARGE_DOT_*`; user decision
    2026-10-26). A ChemDraw reference at three presets with circled ⊕ / ⊖ and "2+" / "2−" on a
    labeled atom, an unlabeled carbon, and an O⁻ with its lone pairs would let them be measured.

## 7. Design-brief (docs/design.md) discrepancies (references win)

1. **Bond drawing.** design.md §3.3 says "stroked line, butt caps". ChemDraw draws filled polygons
   with **mitered joins** at unlabeled atoms (bisector cuts, pointed at degree-3 atoms). Butt
   ends occur only at degree-1 atoms, at labels, and on inner/offset double-bond lines. M1
   should reproduce the join geometry, either as polygons or as polylines with
   `stroke-linejoin: miter`.
2. **Wedge shape.** design.md §3.3 says "filled triangle … width `bold_width`". Measured: a
   trapezoid from LineWidth/2 half-width at the stereocenter (cut to a point by the neighbor
   bisectors) to **1.5 × bold_width** at the drawn end, notched at atoms with other bonds.
3. **Hash lines.** design.md §3.3 says "at `hash_spacing`". Measured: the first hash at
   HashSpacing, the last flush with the end, floor(D / HashSpacing) + 1 hashes evenly spaced
   (pitch ≥ HashSpacing). Hash line thickness = line width (exempt).
4. **Label clearance.** design.md says "stop `margin_width` short of the label's glyph box".
   Measured: 1.0075 × margin_width beyond the glyph outline's reach along the bond, counting
   only outline within ±0.195 × label size of the bond axis (§2.1).
5. **Label anchor.** design.md says "element letter's visual center". Measured: horizontally the
   first letter's advance center (the last letter's when right-justified: Cl–); vertically
   baseline = atom y + 0.3585 × size (cap-height center, ignoring O overshoot).
6. **Charge sign size (deliberate departure from ChemDraw).** ChemDraw's ⊕/⊖ bars shrink
   faster than the labels (0.69× and 0.46× at ×0.8 and ×0.6). By user decision (2026-10-24),
   we keep them proportional: 3.935 × k. Their offsets from labeled atoms are also the ×1.0
   values × k. The bar thickness stays 0.50 pt, like the line width.
7. **§3.2 exempt list.** design.md names "line width at minimum". Measured exempt: line width,
   hash thickness, arrow shaft, arrowhead dimensions, charge bar thickness. Charge bar length
   is non-linear in ChemDraw (3.935 / 2.733 / 1.798), but we scale it with k (item 6).
8. **Curved-arrow heads** are smaller than reaction arrowheads (4.81 × 1.50 × 4.06 vs
   6.0 × 1.58 × 5.25). design.md treats them as one "full arrowhead".
9. **SVG units.** ChemDraw writes `px` for what are pt. Our writer must keep using explicit `pt`.
10. ~~**Install path.**~~ README.md (Installation) gives the per-project install (fixed in M3).
11. **§4.3 text measurement.** Pair kerning must be applied to match resvg. design.md is silent
    on kerning.

## 8. Implemented as (M1)

Module `src/acs_chem_svg/draw_mol.py` unless noted.

| Rule | Section | Implementation |
|---|---|---|
| Bond length, 2D layout scaled to the preset | §1 | `chem.to_molecule` (median bond → `bond_length`; 1.5 Å molfiles map exactly) |
| Mitered joins, pointed degree-3 ends, butt ends | §1 | `Drawer.chains`: main lines joined into `stroke-linejoin: miter` polylines through unlabeled atoms, most-collinear pairs first; an odd end stops at the atom center. `tests/test_m1_vs_chemdraw.py::test_miter_joins_match_chemdraw_pixels` |
| Ring inner line, inset | §1 | `Drawer.double_placement` (ring side), `Drawer.inner_inset` |
| Centered / offset rule, ties, fused rings | §1.1 | `Drawer.double_placement`, `Drawer.ring_side` |
| Centered lines run into degree-3 joins | §1.1 | `Drawer.bond_geometry` (u = h·cot φ) |
| Triple bonds | §1 | `Drawer.bond_geometry` |
| Wedge outline | §1.2 | `Drawer.wedge_polygon`, `Drawer.wedge_wide_end` (notch), `Drawer.behind_wedge` (chain bond fills the notch) |
| Hashed wedge | §1.2 | `Drawer.hash_lines` (stroked segments, length = outline width at the far edge; ChemDraw draws trapezoids) |
| Which atoms are labeled; H side and justification | §2 | `is_labeled` (`Atom.show_label` from `chem.show_carbons` or a molfile alias), `h_side` (then `Drawer.flip_crowded_labels`, unmeasured), `element_x` |
| Label anchor, baseline, subscripts, right-to-left forms, run splitting | §2 | `layout_label` |
| Terminal label centered between offset double-bond lines | §2 | `Drawer.center_terminal_labels`, `Label.shift` |
| Bond-to-label clearance | §2.1 | `clip_distance` → `clip_distance_band` (`band_extent`); fit in `measure_refs.fit_label_clip`; F-rightward exception xfail |
| Charges | §2.2 | `Drawer.charge_center`, `Drawer.charge_bars` (stroked bars); placement offsets are round defaults |
| SVG output rules | design.md §4.4 | `svg_writer.to_svg` |

## 9. Implemented as (M2)

| Rule | Section | Implementation |
|---|---|---|
| Reaction arrow: shaft stroke, concave cubic head, hairline | §3.1 | `arrows.reaction_arrow` (`scene.Path`) |
| Rich text: auto sub/superscripts, italics, escapes; consecutive runs; kerning within runs | §3.2 | `richtext.parse`, `richtext.layout_line` |
| Line pitch | §3.2 | `richtext.line_pitch` (fitted subscript increments, Open 9) |
| Abbreviations (OMe, OBn, …) | §3.3 notes, §2.1 | `chem.condense` (rdAbbreviations on an aromatic copy, collapsed on the Kekulé mol); `draw_mol.layout_abbreviation` |
| Products aligned to their reactant | — (design.md §10) | `chem.prepare` / `chem.correspondence` (atom maps, else MCS) / `chem.align_coords` (CoordGen, else RDKit matching depiction, exact) |
| Row layout: vertical rule, gaps, plus, conditions, arrow length, numbers | §3.3 | `layout.build_scheme` (gaps on `Group.ink_bbox`, the exact stroke box) |
| Overlap and width checks | design.md §6.4 | `geometry_checks.check_scheme`, `scripts/validate.py` |
| Comparison with ChemDraw | — | `tests/test_scheme_vs_chemdraw.py` (exact, 0.1 pt), `scripts/compare.py --spec` (approx, reported) |

## 10. Implemented as (M3)

| Rule | Section | Implementation |
|---|---|---|
| Lone pairs and radicals: which atoms, sides, gaps, exact dots | §4 | `draw_mol.Drawer.dot_sets` → `label_dots` (sides of the element ink box; tail sets face their target) / `bisector_dots`; `electrons.lone_pair_atoms`; `scene.Circle` |
| Curved arrows: stroke, full head, fishhook half head with its offset spine | §4.1 | `mechanism.curved_arrow`, `mechanism.head_polygon`; `tests/test_mech_vs_chemdraw.py` (exact, 0.1 pt) |
| Tails, heads, bow, arc, loop and hook shapes; clash-scored choice | §4.2 | `mechanism.Planner` (`ARC_HANDLES`, `LOOP_HANDLES`; `hook`, `hook_controls`, `HOOK_SCAN`, `HOOK_REACH`; `ring_center`, `RING_SLIDE`: ring hooks inside the ring; `bond_mid`: bond tails at the drawn bond's midpoint), `mechanism.solve_end` (end tangent through the target), `mechanism.plan_step` (bow combinations), `layout.curved_arrows`, `layout.tail_directions` |
| Heads clear the target's label and dots and aim into it | §4.2 | `mechanism.clear_along`, `aims`, `solve_end` (roots, then aiming ends), `FARTHER`, `SOFTER`, `MISS_DEG` (unlabeled targets), `MISS_LABEL_DEG` (labeled targets) |
| Electron bookkeeping | design.md §6.3 | `electrons.check_step` / `check_spec`; `tests/test_electrons.py`, `tests/broken/` |
| Mechanism checks | design.md §6.4 | `geometry_checks.check_mechanism` (`mechanism.clashes`) |
| Zoomed crops | design.md §6.5 | `preview.crops`, `scripts/render.py --crops` |
| Comparison with ChemDraw | — | `scripts/compare.py mech_basic_x*` (our arrows on ChemDraw's CurvePoints, our dots at its centers; templates/m3/) |

## 11. Implemented as (M4)

| Rule | Section | Implementation |
|---|---|---|
| ↔ arrows: reaction heads at both ends, shaft notch to notch | §4.1 | `arrows.resonance_arrow`; `tests/test_resonance_vs_chemdraw.py` (exact, 0.1 pt) |
| Square brackets: bar and lips, stroke exempt, lips × k | §4.1 | `brackets.square_bracket`, `brackets.bracket_pair`; `tests/test_resonance_vs_chemdraw.py` |
| Resonance placement: ↔ length and gaps, bracket gaps | §4.3 | `layout._place` (`put_bracket`), `layout._finish_row`; `scripts/compare.py resonance_x* --spec examples/m4_resonance_x*.yaml` |
| Resonance bookkeeping | design.md §6.3 | `electrons.check_resonance` (`skeleton`, `check_step` per list) |
| Rows: units, greedy breaks, continuation rows start with an arrow, row gap, per-row numbers and vertical rule | §3.4 | `layout._place` (`Unit`), `layout._breaks`, `layout.build_scheme`, `layout._shift_row_records`; `scripts/compare.py scheme_wrap_x* --spec examples/m4_wrap_x*.yaml` |
| Curved arrows avoid reaction and ↔ arrows; a bond target may take its far side (ring bonds) | §4.2 | `layout._arrow_parts`, `mechanism.Planner.build(far=…)` |
| Checks: ↔ arrows and bracket segments join the overlap checks; the step after a bracket may cross it | design.md §6.4 | `geometry_checks.elements`, `check_mechanism` |

## 12. Implemented as (Lewis)

All unmeasured (Open 17): user decisions 2026-10-26, round defaults on the 0.5 pt grid.

| Rule | Value | Implementation |
|---|---|---|
| Layout: square lattice, one bond length at every order, tree center as root, largest branch straight, H last; non-bonded lattice neighbors avoided where possible | — | `chem.lewis_layout` |
| `lewis.bond_length` (atom center to center) | 26.5 / 19.5 / 15.5: the smallest 0.5 pt value where 3 columns fit between two Cl or Br labels (first letter centered on both) with the inner dot edges ≥ 2 × `lone_pair.gap` from each element ink box (Br 26.31 / 19.18 / 15.35, Cl 25.50 / 18.53 / 14.86) | `presets/*.yaml` (`placement`, `nonlinear`); `tests/test_lewis.py` rederives it |
| Shared pairs: `order` columns `dot_spacing` apart along the bond, centered in the gap between the ink boxes; each column two dots `dot_spacing` apart across it | lone_pair values | `draw_mol.Drawer.bond_dots` |
| Assigned electrons (`formal_charge`, `oxidation_state`): columns from the owner's ink box, the first `lone_pair.gap` out, then every `dot_spacing`; an odd last electron on the axis | lone_pair values | `draw_mol.Drawer.bond_dots`; Pauling values in `chem.ELECTRONEGATIVITY` |
| Lone pairs: every atom; sides facing a bond or a non-bonded lattice neighbor stay free | §4.2 rules | `draw_mol.Drawer.label_dots`, `lewis_facing` |
| Checks: shared dots ≥ 2 × gap from labels; assigned dots on the owner's half; dot sets ≥ `dot_spacing` apart | — | `geometry_checks.check_lewis` |
