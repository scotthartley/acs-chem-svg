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

## 2. Atom labels (source: labels_x1.0, m1_set_x*, mech_basic_x1.0, double_bonds_x1.0, clearance_x*, metals_x*, nh_x*)

| Parameter | ×1.0 | ×0.8 | ×0.6 | Method | Class |
|---|---|---|---|---|---|
| Label font / size | Helvetica 10 | 8 | 6 | SVG `<text>` font-size × matrix | exact |
| Explicit carbon labels (CH₃) | laid out and clipped like any label: C centered, H₃ after it; wedge and hash end flat at the label | | | stereo_x1.0 (wedge and hash into CH₃), within 0.03 pt | exact (terminal CH₃ only; CH₂/CH in chains unmeasured) |
| Horizontal anchor | first letter of the element symbol centered on the atom (advance center = ink center within 0.005) | | | glyph origin from SVG + fontTools advance | exact |
| Centered two-letter symbols | a symbol of two or more letters with no H (Pd, Si) is centered whole on the atom (advance center within 0.004) when a vertical or horizontal bond has another bond 165–195° from it (the horizontal case: the user, 2026-10-30, unmeasured); otherwise the first letter is (terminal Cl / Br, -Cl right-justified). Bond lengths unchanged (14.4 / 11.56 / 8.67). Bond ends at a centered Pd: 6.91 / 7.09 left / right, 5.25 / 5.35 up / down (x1.0); Si 5.29–5.58 | 5.54 / 5.62, 4.21 / 4.25 | 4.16 / 4.26, 3.12 / 3.25 | metals_x\* (Pd crosses, SiMe₄, a terminal Cl); the 165–195° window is the user's (2026-10-30), the references have exact verticals; `draw_mol.centered_label` | exact (our clip model: within 0.1 except straight down from Pd, Open 25) |
| Right-justified labels | when the bonds point right (the same test that puts H on the left), the **last** letter of the element is centered instead: –Cl has C on the atom, Cl– has l on the atom (Br likewise) | | | clearance_x*: Cl/Br at 12 directions | exact |
| Baseline below atom center | 3.586 (3.582–3.591) | 2.869 | 2.156 | glyph baseline − atom y; = 0.3585 × size = cap height / 2 (Helvetica capHeight 1469/2048 = 0.7173 em) | exact |
| H after element (NH₂, OH, NH₃, CH₃) | origin + advance(element) | | | dx of the H origin | exact |
| Left-side labels (H₂N, HO, H₃C) | set right-to-left, contiguous: H₂N gives H at −11.393, 2 at −4.171 from the N origin | | | labels_x1.0, mech_basic_x1.0 | exact |
| NH with H above or below | an NH with two bonds whose open side (opposite the sum of the bond directions) is within 15° of straight up has its H **above** N, within 15° of straight down **below**; from 30° off vertical the H is beside N by the rule above. The same at bond angles 60°, 108°, 120° (aziridine, pyrrolidine at 24 turns, R–NH–R). N stays where it is (first letter centered, usual baseline); the H has N's glyph origin x (±0.003), full size, one line pitch away: **8.342** (8.341–8.342) | **6.650** | **5.017** | nh_x\* (`measure_refs.py --nh`, n=11 above or below and 19 beside per preset); the cut between 15° and 30° is unmeasured, `draw_mol.H_VERTICAL_DEG` = 22.5 | exact (neutral NH only; charged N, NH₂, OH, PH unmeasured) |
| Lone atom with H (no bonds) | water is H₂O, right-justified (ChemDraw 26 re-justified our OH₂ on re-save, CDXML gate 2026-10-28); other elements unmeasured (they keep H on the right) | | | m2_equilibrium_x1.0 re-save | exact (O only) |
| Subscript size | 0.75 × label size (7.5 / 6 / 4.5) | | | SVG font-size | exact |
| Subscript baseline drop | 2.30 | 1.84 | 1.38 | = 0.23 × size | exact |
| Terminal label on an offset double bond | moved bond_spacing/2 (1.30) perpendicular to the bond, toward the second line, so it sits centered between the lines; the lines don't move | | | double_bonds_2 acetaldehyde O at 90° and −30° (1.300, 1.377 perpendicular; ≤0.1 along) | exact (×1.0 only; assumed bond_spacing/2 at ×0.8/×0.6) |
| Subscript x | directly after the previous glyph's advance (no gap) | | | | exact |
| Kerning | none observed in labels (runs split at every size change) | | | | exact |

### 2.1 Bond-to-label clearance (source: all references; clearance_x*: H, I, Cl, Br, F, OH, NH₂ at 12 directions × 3 presets)

Every line of a multiple bond stops at the same axial distance (vertical C=O: both lines and a
single C–O end 5.519 from the O center at ×1.0).

**Rule (refit 2026-10-27, `measure_refs.py --fit-clip`).** Take the element glyph outline as
drawn, shifted by (−0.01, +0.03) pt. Keep only the parts within ±0.195 × label size of the bond
axis (edges clipped to that band); a **horizontal** bond keeps the whole outline. The bond stops
1.0075 × MarginWidth beyond the furthest point of those parts along the bond.
- **Fit:** all 783 bond ends are within 0.1 pt (worst 0.094, rms 0.044; 669 ends, rms 0.045
  when fitted on 2026-10-27; the re-run of 2026-11-01 gives the same parameters). These are all
  lines, wedges and hashes into labels in every reference, at every direction and preset, except
  the references `measure_refs.CLIP_FIT_EXCLUDE` leaves out: nh_x\* (ends 6–9° off the axes,
  Open 28), metals_x\* (straight down from Pd, Open 25) and chair_x\* (a rear bond cut at a
  crossing, §14.2). With all of them the fit drifts (n=1051, 56 over 0.1).
  In the 2026-10-27 set: horizontal, 220 ends, worst 0.072 with the whole outline; angled, 449, worst 0.098 with the band.
- **Why horizontal bonds differ (abbrev_x*, 2026-10-27):** a horizontal bond into T (TMS) or
  A (Ac) stops beyond T's bar and A's foot, both outside the band: 2.44 / 0.65 pt further out
  than the band model at ×1.0, −0.02 with the whole outline. The same rule fits –F (bond
  leaving rightward, the former known exception). At 30° and steeper the band is needed
  (I, Cl, Br corners in clearance_x*). The references have no bond between 0° and 30°; the
  switch at 15° (`draw_mol.HORIZONTAL_CLIP_DEG`) is unmeasured.
- **Why a band:** this is the glyph's extent along the bond, limited to what lies near the
  bond's path.
  - **Round glyphs (O):** identical to a convex-hull ray model.
  - **Flat-sided glyphs:** a bond passing an H, N, I or F corner obliquely stops beyond the
    corner, so it gets the 1.2–1.75 × margin "outliers" of the earlier hull model.
  - **Narrow glyphs:** an I or l beside the path is ignored.
- **Former exception (resolved 2026-10-27):** F with the bond leaving rightward stopped
  0.31–0.41 pt beyond the band model (F's top bar). The horizontal rule fits it.
- **Noise floor:** nominally identical geometries in one file differ by up to 0.05 pt
  (m1_set_x0.8 OH at 30°: 4.199–4.268).
- **Superseded:** the 2026-10-24 fit (offset (−0.02, +0.03), band for every direction; 396 of
  399 within 0.1 pt, F rightward excepted); the M1 hull model (1.025 × margin, hull shift (−0.05, +0.08)), and the
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
| Every sign | nearer its own atom (center to atom position) than to any other atom: a corner that fails moves on; with no corner left, the largest bond gap's bisector (equal gaps within 1°: the one nearest the upper right), the center as far from the element ink box as a corner's (user review 2026-10-29: measured along the ray instead, the ring came within 0.07–0.81 pt of the N; now 1.64 / 1.20 / 0.78, as the corner signs; the user's by eye 1.42 / 1.46 / 1.58) | | | — | rule (phase 4b gate, 2026-10-28: ChemDraw binds each sign to the nearest atom on open and resets its charge; nitromethane's corner ⊕ sat 9.66 pt from =O, 10.04 from N). The user placed that ⊕ in the upper-right gap (charges_lewis_x\*). `draw_mol.Drawer._nearest_atom`, `_charge_gap_direction`; `tests/test_cdxml_resonance.py::test_signs_sit_nearest_their_atom` |
| Circled ⊕ / ⊖ (`charges.circled`, `circled_charges`) | the plain bars plus a ring of `circle_radius` (stroke center), stroked at `bar_thickness`: **3.59** | 2.89 | 2.24 | charges_lewis_x\* (ChemDraw's CirclePlus / CircleMinus, all four per file; stroke 0.5; `measure_refs.py --charges-lewis`) | exact ring, copied (user decision 2026-10-28). Placement stays our rule (user decision 2026-10-28; the user's by-eye positions vary): the starting position moves out by the ring's extra reach beyond a plain sign; a circled sign always runs the clearance search (labeled: along the corner diagonal; unlabeled: the largest-gap search above). ChemDraw's bars inside the ring here: 3.93 / 3.27 / 2.60 (it sizes symbols from their CDXML box; §13.4) |
| Multiple charges ("2+", "2−") and `charge_labels` texts ("+3") | digits in the label font at `digit_scale` × label size (**0.75**, = `ratios.superscript_scale`), ink `digit_gap` (**0.5**) from the bars, in the written order, the whole centered on the sign position and vertically on the bars; never circled. Positions move out by the extra width / height and run the clearance search | | | — | unmeasured (user decision 2026-10-26, Open 18). ChemDraw has no multiple-charge symbol: charges_lewis_x\*' Fe²⁺ is typed into the label as a superscript; the user prefers our digit + drawn sign (2026-10-28) |
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

**Arrowhead variants** (source: arrows_large_x1.0; `measure_refs.py --equilibrium`). The table
above is the **standard** head (HeadSize 1000 / ArrowheadCenterSize 875 / ArrowheadWidth 250).
The **large** head is HeadSize 1500 / 1313 / 375, set in ChemDraw's arrow settings: the same
shape, with the same shaft, hairline and back-edge control ratio (0.4375, all large heads):

| Parameter | Large, ×1.0 | Method | Class |
|---|---|---|---|
| Head length | 9.000 (forward, ↔, half heads) → **9.00** | fill path | exact (`arrow.heads.large`) |
| Half-width | 2.30 / 2.35 (asymmetric, as the standard head), mean 2.325 → **2.325**; half heads 2.30 / 2.35 | | exact |
| Inset | 7.875–7.878 → **7.88** | | exact |
| Shaft end from tip | 7.825–7.828 → **7.83** | | exact |

Default (user decision 2026-10-27): reaction, equilibrium and ↔ arrows take the **large** head at
×1.0 and the standard one at ×0.8 and ×0.6 (`arrow.default_head`); a spec's `arrow_head`
overrides it. No reference has a large head at ×0.8 or ×0.6: those presets copy the ×1.0 values
(`assumed_exempt`; head geometry is exempt at every measured size; Open 19). Curved-arrow heads
(§4.1) are unchanged.

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
| Multistep stages ("1." above, "2." below; n stages: first ⌈n/2⌉ above, rest below, user decision 2026-10-28): ink left of each line after the tail | 3.0 / 2.5 / 2.0 | 3.195 | 2.553 | 1.921 | 3.192 / 3.197 | `multistep_indent` |
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
- Abbreviations (OMe, OBn) are `NodeType="Fragment"` nicknames (Ph, Bn, Et, Ac, silyl:
  `Nickname`): one run of rich text on the label baseline. See §3.5.

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

### 3.5 Abbreviations (source: abbrev_x\*, abbrev_2_x\*, ligands_x\*, metals_x\* at ×1.0 / 0.8 / 0.6)

Method: `uv run python scripts/measure_refs.py --abbrev` (2026-11-03 for abbrev_2_x\* and
ligands_x\*). abbrev_x\*: 21 labels on a benzene ring, bond leaving the label to the left (ring
on the right) and to the right, plus OMe, CO₂Me and NO₂ at ±30°. abbrev_2_x\*: about 30 more
groups and formulas the same way (OTf, OTs, OMs, OTBS, OTIPS, OBz, OPiv, OMOM, OTHP, OPMB, NHBoc,
NHTs, NMe₂, N₃, SPh, SO₂Ph, CHO, CO₂*t*-Bu, B(OH)₂, Bpin / pinB, SnBu₃, Boc, Ts, Cbz, Fmoc, Me,
*t*-Bu, *i*-Pr, *n*-Bu, Cy, Mes, *p*-Tol), ring NTs (an N-tosyl pyrrolidine turned in 60° steps:
Ts beside, above and below N), OTBS / NHBoc on vertical bonds, OTf / NHBoc at ±30° and ±150°.
ligands_x\*: PMe₃, PCy₃, P(*t*-Bu)₃, P(*o*-tol)₃, L, Lₙ, OAc, Cl, Br, dppf, XPhos and IPr round a Pd
(both justifications, above and below). metals_x\*: PPh₃. ChemDraw justified every label
itself (`LabelJustification`); the CDXML keeps the key text in both forms and the SVG shows the
reversed one.

The labels are rules, not a list (`abbrev.py`; registry `src/acs_chem_svg/data/abbreviations.yaml`,
plus a spec's `abbreviations:`). Groups (Me, Boc, TBS) are named substituents; formulas (OMe,
CO₂Me, PPh₃) are labels on their attachment atom; text labels (L, Lₙ, dppf, XPhos, IPr) have no
chemistry.

| Rule | ×1.0 | ×0.8 | ×0.6 | Class |
|---|---|---|---|---|
| Bond leaving to the left or vertically | key text left-justified, **first** glyph's advance center on the atom: −0.005 to +0.004 (all four references) | −0.005 to +0.005 | −0.005 to +0.005 | exact |
| Bond leaving to the right | reversed text right-justified, **last** glyph's advance center on the atom: −0.003 to +0.003 | −0.005 to +0.005 | −0.004 to +0.005 | exact |
| Baseline below atom center | 3.582–3.590 | 2.864–2.874 | 2.147–2.157 | exact (= §2) |
| Anchor (`abbrev.anchor`) | the attachment element's glyphs when the text starts (reversed: ends) with its token (OMe, MeO: O; MeO₂C: C; Bu₃Sn: Sn; NTs: N), else the end glyph (Ph: P / h; TMS: T / S; N₃ reversed: the ₃) | | | exact |
| Bond clip | §2.1 on the anchor glyphs at ChemDraw's own place (below). Model − end: −0.09 to +0.06 (abbrev_x\*), −0.09 to +0.12 (abbrev_2_x\*, not counting ring NTs's ends 3° and 9° off the axes), −0.08 to +0.08 (ligands_x\*, not counting dppf above Pd) | −0.09 to +0.08 | −0.09 to +0.09 | exact |
| Count subscript at the centered end | N₃ and P(*o*-tol)₃ right-justified: the ₃ drawn as if full size (text ends where the full-size 3 centered on the atom would: advance center +0.698 / +0.697); ChemDraw's own place (bond clip, CDXML p) has the ₃ centered | +0.552 / +0.558 | +0.419 / +0.421 | exact |
| Letter subscript at the centered end | Lₙ right-justified (typed italic subscript): ₙ centered as drawn, +0.002 | +0.005 | +0.005 | exact |
| Kerned first glyph | Ts left-justified: T drawn with its plain advance centered (−0.004); ChemDraw's own place centers its kerned advance (T\|s −1.11 at 10 pt: +0.55), which the bond clip (model − end −0.08 with it, +0.47 without) and the CDXML label (p, box +0.55) use | −0.02 with it (+0.46) | +0.03 (+0.30) | exact |
| Ring NTs (two bonds) | the NH rule (§2): Ts beside N (NTs / TsN by the bonds' x sum), or on the line below / above N when the open side is within 15° of vertical, its origin x on N's (0.00); pitch below 8.342 (= `label_h_line_pitch`), above 8.542 (`label_tail_above_pitch`) | 6.650 / 6.800 | 5.017 / 5.117 | exact, n = 1 each |
| Runs | split at style and script changes, as in conditions: MeO \| ₂ \| C, *t* \| -Bu | | | exact |

**Reversal** (`abbrev.reverse`): the key's tokens in reverse order, also inside parentheses.
Tokens are matched longest first against the group names, then element symbols (case-sensitive),
then italic-prefixed names (*s*-Bu); digits after a token or `)` are its count. A one-token
label stays as it is. An entry's `reversed:` overrides it (Bpin → pinB: pin bridges both O, so
it is no group). `tests/test_abbrev_grammar.py` holds every right-justified text ChemDraw shows
in the four references:

| Key | ChemDraw | Key | ChemDraw | Key | ChemDraw |
|---|---|---|---|---|---|
| OMe, OEt, OBn, OPh, OAc | MeO, EtO, BnO, PhO, AcO | OTf, OTs, OTBS, OTIPS | TfO, TsO, TBSO, TIPSO | NHAc, NHBoc, NHTs | AcHN, BocHN, TsHN |
| OBz, OPiv, OMOM, OTHP, OPMB | BzO, PivO, MOMO, THPO, PMBO | SMe, SPh, SO₂Ph | MeS, PhS, PhO₂S | NMe₂, NO₂, NTs | Me₂N, O₂N, TsN |
| CO₂Me, CO₂Et, CO₂H, CO₂*t*-Bu | MeO₂C, EtO₂C, HO₂C, *t*-BuO₂C | CHO, CF₃, CN | OHC, F₃C, NC | B(OH)₂, Bpin, SnBu₃ | (HO)₂B, pinB, Bu₃Sn |
| PPh₃, PMe₃, PCy₃, P(*t*-Bu)₃ | Ph₃P, Me₃P, Cy₃P, (*t*-Bu)₃P | Ph, Bn, Et, Ac, Me, Cy, Mes, *p*-Tol | unchanged | TMS, TES, TBS, TIPS, Boc, Ts, Cbz, Fmoc | unchanged |
| *t*-Bu, *i*-Pr, *n*-Bu, N₃ | unchanged | Lₙ, dppf, XPhos | unchanged | IPr, P(*o*-tol)₃ | PrI, P(*o*-tol)₃: **ours IPr, (*o*-tol)₃P** |

**User decisions (2026-11-03):** one reversal rule rather than ChemDraw's reading of labels it
can't interpret: IPr stays IPr (ChemDraw: PrI, iodine + propyl) and P(*o*-tol)₃ reverses to
(*o*-tol)₃P (ChemDraw: unchanged). Scope: the user's own definitions plus **composition** below.

**Composition** (`chem.compose`): after groups are collapsed, a heteroatom (not C, not a d-block
metal; uncharged) bonded to collapsed groups and to one or two other atoms becomes one label:
element, H count, then the groups in SMILES order, identical ones counted, a counted group with
a hyphen in parentheses: NHBoc, OTBS, SnBu₃, P(*t*-Bu)₃, ring NTs. A formula entry takes
precedence (CO₂Me stays one entry). Placement as above; two bonds take the NH rule.

**Text labels** (L, Lₙ, dppf, XPhos, IPr, Ar, or the spec's own without SMILES): drawn on a
stand-in atom (`[*:1]` with `labels: {1: dppf}`), placed by the same rules, no chemistry.

Not matched, by reference (tests list them as known):
- ring NTs's bonds 3° and 9° off horizontal stop 0.12–0.60 pt nearer N than our model (Open 28's
  cut, not emulated);
- dppf above Pd: ChemDraw stops the bond below the p descenders, 1.9 pt farther (x0.8 / x0.6:
  no bond at all); ours clips at d (Open 31);
- *p*-Tol left-justified at ×1.0: 0.12 nearer than the oblique p's outline;
- abbrev_2_x0.8 / x0.6 set *t*-Bu, *i*-Pr, *n*-Bu, *p*-Tol and CO₂*t*-Bu upright (CDXML face 96,
  not 98), and ligands_x0.8 / x0.6 Lₙ's n upright (face 32): ours stay italic (bond ends up to
  0.74 pt apart); several of their fragments also sit 1.3–2.3 pt apart in the CDXML and the SVG,
  so the CDXML diff runs on the ×1.0 files (Open 32).

### 3.6 Equilibrium arrows (source: equilibrium_x1.0 / x0.8 / x0.6, arrows_large_x1.0)

Method: `uv run python scripts/measure_refs.py --equilibrium` (`_equilibrium_entry`: the two
straight shaft polygons along each CDXML arrow and the curved half-head polygon at one end of
each, in the frame of Tail3D → Head3D). Each reference has three rows, A ⇌ B drawn with
ChemDraw's equilibrium tools: balanced, unbalanced favoring products, and the unbalanced arrow
turned around (favoring reactants), each with one line of conditions above and below.
CDXML: `ArrowheadHead`/`ArrowheadTail` `HalfLeft`, `ArrowShaftSpacing` 300, unbalanced
`ArrowEquilibriumRatio` 300.

| Parameter | ×1.0 | ×0.8 | ×0.6 | Class |
|---|---|---|---|---|
| Shafts | two, shaft width 0.60, centered on the axis (center offset ≤ 0.006) | same | same | exact |
| Shaft spacing, center to center | 1.805 | 1.800 | 1.793 | exact, **exempt** (`equilibrium.shaft_spacing` 1.80) |
| Directions | the forward shaft (left of travel: above for a left-to-right arrow) points at the head, the back shaft at the tail | same | same | exact |
| Half head | the outer half of the full head (§3.1, same variant: length, half-width, inset, shaft end from tip all within 0.02 of the full head's), with its tip and notch on the shaft's **inner** edge (0.300 = shaft width / 2 from the center line) and the barb **outward** (`HalfLeft`) | same | same | exact |
| Back-edge control | forward shaft: at 0.4375 × barb from the center line, as the full head (0.437–0.438); back shaft: 0.60 nearer the center line (0.599–0.601; large 0.600) | same | same | exact, **exempt** (`equilibrium.back_control_shift`) |
| Balanced: stagger | each shaft's butt end 0.902 behind the other shaft's tip (ChemDraw's Tail3D / Head3D are the forward shaft's butt and tip; the back tip reaches 0.9 beyond Tail3D) | 0.900 | 0.896 | exact, **exempt** (`equilibrium.stagger` 0.90) |
| Unbalanced: long shaft | Tail3D to Head3D | | | exact |
| Unbalanced: short shaft | centered on the long one (≤ 0.001); butt → notch 1/3 of the long one's: long / short = 3.000 (74.58 / 28.37; large heads 74.50 / 30.08) | 3.000 (59.67 / 23.39) | 3.000 (44.75 / 18.40) | exact (`ratios.equilibrium_long_per_short`) |
| Favoring reactants | the favoring-products arrow turned 180°: the long shaft below, pointing back | same | same | exact |
| Large head | same rules (spacing 1.80, stagger 0.90, shift 0.600–0.601, ratio 3.000) with the large half head (§3.1) | — | — | exact |

Placement (hand-placed; provenance only, user decision 2026-10-27: the existing `scheme`
defaults apply). Our rules, unmeasured: the conditions' baselines are measured from the **outer**
shafts (`conditions_above_baseline` above the forward shaft, `conditions_below_baseline` below the
back one); the lines are centered on the arrow's middle (it is symmetric; ChemDraw −0.4 to −1.6,
by eye); the length is max(`arrow_min_length`, widest line + 2 `conditions_pad` + 2 head insets),
one notch at each end.

| Hand-placed (provenance) | ×1.0 | ×0.8 | ×0.6 | Default used |
|---|---|---|---|---|
| Above: forward shaft − last baseline | 5.00 / 5.25 / 4.75 (large 5.69 / 6.25 / 5.26) | 3.82 / 4.02 / 3.62 | 2.64 / 2.79 / 2.49 | `conditions_above_baseline` 5.5 / 4.5 / 3.5 |
| Below: first baseline − back shaft | 12.36 / 12.11 / 12.61 (large 11.69 / 11.12 / 12.12) | 9.71 / 9.51 / 9.91 | 7.06 / 6.91 / 7.21 | `conditions_below_baseline` 11.5 / 9.5 / 7.0 |

Exact check: `tests/test_equilibrium_vs_chemdraw.py` (our arrows on ChemDraw's, every shaft and
half head within 0.1 pt at all three presets, and arrows_large_x1.0's forward and ↔ heads).

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
| Tangent arrow: an arc or loop of at least 1.2 bond lengths whose exit is over 30° off square, or which clashes, gives way to a cubic set by its end directions: it leaves square (the pair's direction, a bond's normal; when square has no clean curve, 10°, 20° or 30° off it, user review 2026-11-06) and lands exactly on its aim, handles 1/3 × chord (then 0.45, 0.25 and 0.6 if no clean curve at 1/3), a C or an S by geometry; the smoothest clean one (bending energy). A clean nearly flat arc through a gap (`FLATTER`) stays. Ruled out: a stroke over 1.6 × chord long, radius under 0.3 × bond_length, under 35° of turning, a self-crossing, a stroke within margin_width of a bond or plus sign, or 2 × margin_width of a label (the ends' head_length excepted; user review 2026-11-06: a stroke 0.33 pt from an H read as touching), a fishhook whose end hardly bends. Fishhooks that meet land from their shared side, and a C among them bows to it (an S is allowed, user request 2026-11-06). `shape: tangent` forces one, `shape: arc` or `curvature` an arc | as stated | — | — | — | user decision 2026-11-05 (`mechanism.TANGENT_*`, `Planner.tangent_pass`); unmeasured, no reference contains an S |
| Loop (chord shorter than a bond): control offset / chord | 1.3 | 1.311 | 1.310 | 1.310 | `curved_arrow.loop_curvature` |
| Arc handles (in the chord frame, along / offset ÷ chord) | c1 (0.258, 0.987 h), c2 (0.874, 0.770 h): user decision 2026-10-26, smoothed. ChemDraw's tangent angles kept (63.8° and 72.8° off the chord at h 0.531), handle lengths refit (0.584, 0.428 of the chord; ChemDraw 0.592, 0.268) so curvature runs evenly to the head (at h 0.5: 1.3 → 1.8 → 2.2 per chord) instead of dipping mid-arc and spiking at the head (1.4 → 2.1 → 1.5 → 5.1). The arc is 16% taller at one chord (0.33 vs 0.29 of it) | (0.261, 0.531), (0.921, 0.256) | same | (0.288, 0.531), (0.930, 0.264) | `mechanism.ARC_HANDLES` (shape) |
| Loop handles | c1 (1.19, h), c2 (1.45, 0.59 h) | (1.194, 1.311), (1.446, 0.769) | same | same | `mechanism.LOOP_HANDLES` (shape) |
| Hook: an arrow onto a bond next to its source (a bond sharing one atom; lp(n) onto a bond on n; E2's C–H → C–C, a sigma complex's C=C → C–C⁺) | leaves along the source bond's normal on the bow side (lp: the pair's direction) and lands along the target bond's normal (±10°, ±20° tried); both handles the circular arc's for that turn, 4/3 tan(θ/4) R, at least 0.35 × bond_length; × 0.75, × 0.5 tried. Chosen from the endpoints, not the chord length. A `curvature` override keeps the arc or loop. When the tail lies inside the angle between the two bonds (C–H → C–C forming a π bond), there's no room on that side: the arrow passes through the target bond and curls back, its tip head_gap past the axis on the far side and its head pointing back at the bond, turned 30° toward the shared atom so the loop opens away from it; handles at least 0.9 × bond_length; that crossing isn't a clash (user choice 2026-10-26, option D of six rendered). On a chair a ring bond lies in that angle (E2's C3–C4): the hook first tries the near side, one arc of about 120° from the C–H bond's normal across that ring bond onto the target bond, its tip slid 0.25 × bond_length along the target bond away from the shared atom where the angle opens; the stroke may cross the shared atom's bonds between the source and target bonds, its head may not, and every other bond still clashes. The pass-through, which then also crosses those bonds, is the fallback (user decisions 2026-10-29). The tail's side and the far side without the curl are the fallbacks (phenoxy's C=C → C–O, where the radical's fishhook meets it). Between two bonds of one ring (a sigma complex's C=C → C–C⁺) the hook bows into the ring: leaving one bond's inward normal and landing across the next is a 120° turn, drawn as a circular arc (no 0.35 × bond_length minimum, which made it pointed), with both ends slid 0.25 × bond_length along their bonds away from the shared atom so the fixed-size head leaves a stroke at ×0.6; an explicit `bow:` still wins, and the outside stays a fallback (user decision 2026-10-27) | | | | rule (user decisions 2026-10-26, 2026-10-27); unmeasured (`mechanism.HOOK_SCAN`, `HOOK_REACH`, `hook_controls`, `Planner.cramped`, `CROSS_REACH`, `CROSS_TILT`, `MechArrow.crosses`, `Planner.swept`, `Planner.ring_center`, `RING_SLIDE`) |
| Short hooks: minimum stroke and excursions (Meisenheimer lp(C⁻) → C–N, where the chord is a few pt and the head hid the stroke) | Every curved arrow whose stroke shows less than 2.0 × head_length outside its head costs like a miss (below a clash, above any rank); ring hooks are exempt (their 0.25 × bond_length slide is their rule). When a hook's best candidate still costs something, the planner tries, after every default candidate, flipped or not: first slid arcs, the hook with its tip slid 0.25 × bond_length along the target bond away from the shared atom and handles at least 0.45 × bond_length (the Meisenheimer's lp(C⁻) → C–N: a small arc straight off the pair; user review 2026-10-29); then excursions: the tip as for the hook (the tail's side or the far side), the exit the natural one (bond normal on the bow side, or the pair's direction, which a lone pair's excursion keeps) turned 0, ±15° … ±75°, the landing the target bond's normal turned up to ±40°, the exit and reversed landing opening at least 60° (else it is a needle), and both handles equal, scaled by bisection so the curve's farthest point stands 0.65 × bond_length off the tail → stroke-end chord (then 0.8, 0.95). Ring hooks keep their arc | | | | rule (user decision 2026-10-29, rendered 0.5 / 0.65 / 0.8 and 1.75 / 2.0 / 2.5); unmeasured (`mechanism.MIN_STROKE`, `stroke_shown`, `SLID_REACH`, `SLID_RANK`, `EXCURSION_REACH`, `EXCURSION_LAND`, `EXCURSION_LP_SCAN`, `EXCURSION_GROW`, `EXCURSION_SCAN`, `EXCURSION_OPEN`, `EXCURSION_RANK`, `excursion_controls`, `apex`) |
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
1. ~~**F with the bond leaving rightward** (§2.1).~~ **Resolved (abbrev_x*, 2026-10-27):** horizontal
   bonds clip against the whole glyph outline. Still open: where between 0° and 30° the band takes
   over (`draw_mol.HORIZONTAL_CLIP_DEG` = 15°, unmeasured); a reference with T or F bonds at 10° and
   20° would settle it. The cut makes a bond end step by about 0.18 pt as a bond turns through
   15° (into N: 5.25 → 5.07 at ×1.0). **Accepted (user decision, 2026-10-30):** labels don't
   turn with the bonds, so bond ends vary with direction anyway, and ChemDraw has a cut of its
   own (between 0° and 6°, Open 28); don't smooth it.
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
10. ~~**Right-justified abbreviations** (MeO–, BnO–).~~ **Resolved (abbrev_x*, 2026-10-27, §3.5):**
    all 21 abbreviations measured in both forms at every preset.
11. ~~**Multistep arrows with more than 2 stages** (where "3." goes).~~ **Resolved by user
    decision (2026-10-28):** the first ⌈n/2⌉ stages go above the arrow, the rest below (3 → 1, 2 /
    3; 4 → 1, 2 / 3, 4), each at the same indent and line pitch; no upper limit. Not in the
    references.
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

17. **Lewis structures** (§12): ChemDraw has no Lewis tool, so no reference can measure them.
    The lattice layout, `lewis.bond_length` (also the bond-line form's) and the shared/assigned
    dot positions are user decisions (2026-10-26), approved unchanged in charges_lewis_x\*
    (2026-10-28).

18. **Charge clearance, circled signs, multiple charges** (§2.2): ~~no reference has a circled
    ⊕ / ⊖~~ **ring measured in charges_lewis_x\* (2026-10-28)**: `charge.circle_radius` is
    ChemDraw's. Still rules, not measurements: where a circled sign sits (the user placed them
    by eye, inconsistently), the push clear of lone pairs (`draw_mol.CHARGE_DOT_*`), and the
    digits of "2+" (`digit_scale`, `digit_gap`): ChemDraw has no multiple-charge symbol, and the
    user prefers our digit + sign to its superscript label text.

20. **CDXML** (§13): ~~what ChemDraw recomputes on open~~ **resolved by the gate,
    2026-10-28:** it keeps our labels and boxes. ~~Still unknown: whether `cdxml.charge_box_dx`
    (the second BoundingBox point) matters~~ **it does (charges_lewis_x\*, 2026-10-28):
    ChemDraw sizes a symbol from it** (Electron dot radius = dx / 9; the user's CH₂ dots at
    ×0.8 / ×0.6 came out larger for that reason). Circled signs, "2+" and Lewis: §13.4.

21. **Chair equatorial directions** (§14.1): the reference's equatorial bonds at four vertices
    sit on ChemDraw's 15° snap, 0.85–0.98° off the textbook "parallel to the ring bond one
    removed". We copy the reference. Open: whether the user wants the textbook parallel instead.

22. **Chair crossings not in the reference** (§14.2): a rear wedge or hashed bond, an offset
    double or a triple in front, and two crossings on one bond are drawn by rule, unmeasured.
    The label stub (a sliver between a label clip and a gap) is copied from ChemDraw; it is 0.52
    pt at x1.0 and has no x0.8 / x0.6 measurement (those files lack the label).

23. **Carbenes** (§10): carbene pair and triplet dot placement on skeletal carbons and the singlet CDXML
    (`LonePair` + `Radical="Singlet"`) are unmeasured; no reference (user decision 2026-10-29).

24. **Arc arrows** (§15.1): clockwise heads in arcs_x\* are the measured shape rescaled by
    −2.8 to +3.1 % with no pattern found; ours keeps the measured head both ways. Arc equilibria
    (arcs_eq_x\*) are approximate: ChemDraw's clockwise and counterclockwise ones differ by up to
    0.3 pt, its short shaft sits off ±spacing / 2, and x1.0's large balanced ones are glitched
    ellipses (§15.1). CDXML (§15.3): gate passed 2026-10-30.
25. **Bond straight down from a centered Pd** (§2, metals_x\*): our clip model (§2.1) stops it
    0.24 / 0.16 / 0.17 pt nearer the atom than ChemDraw (5.11 vs 5.35 at x1.0); the other
    directions at Pd and Si are within 0.1. A refit with metals_x\* (2026-11-01) does not close it:
    the Pd ends stay 0.16–0.25 off and the other parameters move (offset −0.03, 1.015 × margin,
    band 0.19), so metals_x\* is left out of the fit (`CLIP_FIT_EXCLUDE`).
26. ~~**Ligand abbreviations other than PPh₃** (§3.5).~~ **Resolved (ligands_x\*, 2026-11-03):**
    PMe₃, PCy₃, P(*t*-Bu)₃, P(*o*-tol)₃, L, Lₙ, dppf, XPhos and IPr on Pd, both justifications, three
    presets; PMe₃'s CDXML fragment is ChemDraw's own now.
27. ~~**H above or below an NH / OH** (§2).~~ **Resolved (nh_x\*, 2026-10-30):** an NH's H goes
    above / below N when its open side is within 15° of vertical, beside from 30° (§2). Still
    open: where between 15° and 30° it switches (22.5° assumed), and charged N (iminium NH⁺,
    NH₂⁺), PH and other elements, which keep the H beside.
28. **Bonds just off an axis** (§2.1, nh_x\*): ChemDraw's bonds 6° and 9° off horizontal or
    vertical stop 0.17–0.61 pt nearer N than our model (exactly horizontal matches), and x0.6 has
    three more ends 0.10–0.46 off. **Not emulated (user decision, 2026-10-30):** a ChemDraw
    error; our rings are regular and their bond ends should follow rotation smoothly.
    `tests/test_nh_vs_chemdraw.py` allows these at their measured size (`NEAR_AXIS_KNOWN`,
    `KNOWN`); the other nh_x\* ends stay within 0.1 pt.
29. ~~**Further group abbreviations** (§3.5).~~ **Resolved (abbrev_2_x\*, 2026-11-03):** about 30
    more groups and formulas in both justifications; labels are now rules (registry, reversal,
    anchor, composition), so a label without a reference is drawn and exported too.
30. **Repeat-unit brackets** ([ ]ₙ for an oligomer or polymer; inbox 2026-10-30): not drawn. Needs
    a reference (polymer_x0.6 / 0.8 / 1.0): brackets across the two bonds at several angles, the
    subscript n's size and place, and the CDXML objects ChemDraw writes for them.

31. **dppf above Pd** (§3.5, ligands_x1.0): ChemDraw stops the bond below dppf's p descenders,
    1.9 pt beyond our anchor-glyph clip, and draws no bond at ×0.8 / ×0.6; Lₙ above Pd (its
    subscript n below the baseline) has no bond at any preset, where ours draws a short one. Clipping against the
    whole label fits dppf but not P(*o*-tol)₃ or XPhos above Pd (2.5 / 0.28 pt off), so the rule
    is unknown: a reference with descender labels (dppf, dppp, Tol) above and below a metal at
    several offsets would settle it.
32. **abbrev_2_x0.8 / x0.6 and ligands_x0.8 / x0.6 re-export** (§3.5): their italic prefixes and
    Lₙ's n are upright, unlike ×1.0's, and several fragments differ by 1.3–2.3 pt between the
    CDXML and the SVG. Re-exporting them with italics, from one save, would let the tests hold
    those presets to 0.1 pt and run the CDXML diff on them.

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

19. **Arrowhead variants** (§3.1, §3.6): no reference has a large head at ×0.8 or ×0.6 (the presets
    copy the ×1.0 values, `assumed_exempt`). ~~arrows_large_x1.0's unbalanced arrows had
    ArrowheadCenterSize 1500 (flat-backed).~~ **Resolved (re-export, 2026-10-27):** every large
    head is 1500 / 1313 / 375.

## 9. Implemented as (M2)

| Rule | Section | Implementation |
|---|---|---|
| Reaction arrow: shaft stroke, concave cubic head, hairline | §3.1 | `arrows.reaction_arrow` (`scene.Path`) |
| Head variants (standard, large) and the per-preset default; spec `arrow_head` | §3.1 | `preset.Preset.head`, `arrows.head_marker(variant=…)`; `tests/test_equilibrium_vs_chemdraw.py` |
| Equilibrium arrows: two shafts, half heads, stagger, unbalanced short shaft; spec `arrow_type` | §3.6 | `arrows.equilibrium_arrow`, `arrows.head_marker(half=…)`, `arrows.step_arrow`; `layout._place_conditions` (outer shafts, `arrows.conditions_mid`); `tests/test_equilibrium_vs_chemdraw.py` (exact, 0.1 pt) |
| Rich text: auto sub/superscripts, italics, escapes; consecutive runs; kerning within runs | §3.2 | `richtext.parse`, `richtext.layout_line` |
| Line pitch | §3.2 | `richtext.line_pitch` (fitted subscript increments, Open 9) |
| Abbreviations (OMe, OBn, …) | §3.3 notes, §2.1 | `chem.condense` (rdAbbreviations on an aromatic copy, collapsed on the Kekulé mol); `draw_mol.layout_abbreviation` |
| Products aligned to their reactant | — (design.md §10) | `chem.prepare` / `chem.correspondence` (MCS seeded by atom maps, else MCS) / pin: `chem.stable_core` + `chem.align_coords` (CoordGen, else RDKit matching depiction, exact) / orient: `chem.orient_to` (own layout, fitted rotation snapped to the 30° grid) |
| Row layout: vertical rule, gaps, plus, conditions, arrow length, numbers | §3.3 | `layout.build_scheme` (gaps on `Group.ink_bbox`, the exact stroke box) |
| Overlap and width checks | design.md §6.4 | `geometry_checks.check_scheme`, `scripts/validate.py` |
| Comparison with ChemDraw | — | `tests/test_scheme_vs_chemdraw.py` (exact, 0.1 pt), `scripts/compare.py --spec` (approx, reported) |

## 10. Implemented as (M3)

| Rule | Section | Implementation |
|---|---|---|
| Lone pairs and radicals: which atoms, sides, gaps, exact dots | §4 | `draw_mol.Drawer.dot_sets` → `label_dots` (sides of the element ink box; tail sets face their target) / `bisector_dots`; `electrons.lone_pair_atoms`; `scene.Circle` |
| Curved arrows: stroke, full head, fishhook half head with its offset spine | §4.1 | `mechanism.curved_arrow`, `mechanism.head_polygon`; `tests/test_mech_vs_chemdraw.py` (exact, 0.1 pt) |
| Tails, heads, bow, arc, loop and hook shapes; clash-scored choice | §4.2 | `mechanism.Planner` (`ARC_HANDLES`, `LOOP_HANDLES`; `tangent_pass`, `tangent_controls`, `curve_shape`, `TANGENT_*`: curves set by their end directions; `hook`, `hook_controls`, `HOOK_SCAN`, `HOOK_REACH`; `SLID_REACH`, `excursion_controls`, `EXCURSION_REACH`, `stroke_shown`, `MIN_STROKE`: short hooks; `swept`: the ring bonds a chair hook may cross; `ring_center`, `RING_SLIDE`: ring hooks inside the ring; `bond_mid`: bond tails at the drawn bond's midpoint), `mechanism.solve_end` (end tangent through the target), `mechanism.plan_step` (bow combinations), `layout.curved_arrows`, `layout.tail_directions` |
| Heads clear the target's label and dots and aim into it | §4.2 | `mechanism.clear_along`, `aims`, `solve_end` (roots, then aiming ends), `FARTHER`, `SOFTER`, `FLATTER`, `MISS_DEG` (unlabeled targets), `MISS_LABEL_DEG` (labeled targets) |
| Electron bookkeeping | design.md §6.3 | `electrons.check_step` / `check_spec`; `tests/test_electrons.py`, `tests/broken/` |
| Carbenes / nitrenes: two RDKit radicals are a singlet's lone pair, always drawn and an `lp(n)` source, unless the species lists the atom in `triplet:` (two radicals); arrows that form one leave a singlet; placement reuses the lone-pair rules (`label_dots`, `bisector_dots` at `radical.distance`); a skeletal triplet's second dot goes in the bond gap opposite the first (Ph₂C: one below the bonds, one between them); unmeasured (Open 23) | §4 | `chem.nonbonding`, `Drawer._opposite_gap`, `Atom.singlet`, `spec.mark_triplets` / `with_spin`, `electrons.spin_problems`; `tests/test_carbenes.py` |
| Mechanism checks | design.md §6.4 | `geometry_checks.check_mechanism` (`mechanism.clashes`) |
| Zoomed crops | design.md §6.5 | `preview.crops`, `scripts/render.py --crops` |
| Comparison with ChemDraw | — | `scripts/compare.py mech_basic_x*` (our arrows on ChemDraw's CurvePoints, our dots at its centers; templates/m3/) |

## 11. Implemented as (M4)

| Rule | Section | Implementation |
|---|---|---|
| ↔ arrows: reaction heads at both ends, shaft notch to notch | §4.1 | `arrows.resonance_arrow`; `tests/test_resonance_vs_chemdraw.py` (exact, 0.1 pt) |
| Square brackets: bar and lips, stroke exempt, lips × k | §4.1 | `brackets.square_bracket`, `brackets.bracket_pair`; `tests/test_resonance_vs_chemdraw.py` |
| Resonance placement: ↔ length and gaps, bracket gaps | §4.3 | `layout._place` (`put_bracket`), `layout._finish_row`; `scripts/compare.py resonance_x* --spec examples/m4/m4_resonance_x*.yaml` |
| Resonance bookkeeping | design.md §6.3 | `electrons.check_resonance` (`skeleton`, `check_step` per list) |
| Rows: units, greedy breaks, continuation rows start with an arrow, row gap, per-row numbers and vertical rule | §3.4 | `layout._place` (`Unit`), `layout._breaks`, `layout.build_scheme`, `layout._shift_row_records`; `scripts/compare.py scheme_wrap_x* --spec examples/m4/m4_wrap_x*.yaml` |
| Curved arrows avoid reaction and ↔ arrows; a bond target may take its far side (ring bonds) | §4.2 | `layout.arrow_parts`, `mechanism.Planner.build(far=…)` |
| Checks: ↔ arrows and bracket segments join the overlap checks; the step after a bracket may cross it; an arrow on a reagent over the reaction arrow may cross that arrow (user decision 2026-11-05) | design.md §6.4 | `geometry_checks.elements`, `check_mechanism` |

## 12. Implemented as (Lewis)

Unmeasured user decisions (2026-10-26), round defaults on the 0.5 pt grid, **approved by the user on 2026-10-28**: charges_lewis_x\* is the user's ChemDraw edit of our draft (`cdxml_gate.py --lewis-draft`), and it keeps every Lewis atom and dot where we put them (`tests/test_cdxml_species.py::test_charges_lewis_reference_keeps_our_lewis`). That approves the positions; it does not measure them (Open 17).

| Rule | Value | Implementation |
|---|---|---|
| Layout: square lattice, one bond length at every order, tree center as root, largest branch straight, H last; non-bonded lattice neighbors avoided where possible | — | `chem.lewis_layout` |
| `lewis.bond_length` (atom center to center) | 26.5 / 19.5 / 15.5: the smallest 0.5 pt value where 3 columns fit between two Cl or Br labels (first letter centered on both) with the inner dot edges ≥ 2 × `lone_pair.gap` from each element ink box (Br 26.31 / 19.18 / 15.35, Cl 25.50 / 18.53 / 14.86) | `presets/*.yaml` (`placement`, `nonlinear`); `tests/test_lewis.py` rederives it |
| Shared pairs: `order` columns `dot_spacing` apart along the bond, centered in the gap between the ink boxes; each column two dots `dot_spacing` apart across it | lone_pair values | `draw_mol.Drawer.bond_dots` |
| Bond lines (`bonds: lines`, shared only; user decision 2026-10-28) | the same lattice, labels and lone pairs; bonds drawn as in any labeled structure, no shared-pair dots | `chem.Molecule.lewis_lines`, `draw_mol.Drawer.draw` / `lewis_dot_sets` |
| Assigned electrons (`formal_charge`, `oxidation_state`): columns from the owner's ink box, the first `lone_pair.gap` out, then every `dot_spacing`; an odd last electron on the axis | lone_pair values | `draw_mol.Drawer.bond_dots`; Pauling values in `chem.ELECTRONEGATIVITY` |
| Lone pairs: every atom; sides facing a bond or a non-bonded lattice neighbor stay free | §4.2 rules | `draw_mol.Drawer.label_dots`, `lewis_facing` |
| Oxidation-state view: each atom's nonzero oxidation state drawn as a charge text (`+`, `2−`, `4+`), display only (formal charges stay in the chemistry); formal-charge view: the formal charges, as on any species | formal charge + Σ over bonds (+order to a more electronegative partner, −order to a less) | `chem.oxidation_states`, `chem.from_smiles` (user decision 2026-10-28, as the user drew them in charges_lewis_x\*) |
| Charges on a Lewis atom: the upper right corner first (every lattice bond is 45° from each corner, so only the 45° limit applies), then the clearance search | — | `draw_mol.Drawer.charge_center` (user review 2026-10-28: the user put every sign at the upper right, "aim for consistency") |
| Checks: shared dots ≥ 2 × gap from labels; assigned dots on the owner's half; dot sets ≥ `dot_spacing` apart | — | `geometry_checks.check_lewis` |

## 13. CDXML export (phase 1: M1 species, phase 2: M2 schemes §13.1, phase 3: M3 mechanisms §13.2, phase 4a: abbreviations, resonance, wrap §13.3, phase 4b: charges and Lewis §13.4; source: every reference's .cdxml)

`acs_chem_svg.cdxml_writer` writes what ChemDraw writes for its own structures, so ChemDraw can
edit a figure. These are file conventions, not drawn geometry. They are measured from the CDXML
files with `measure_refs.py --cdxml-symbols` and `tests/test_cdxml_species.py` (reference diff).

| Item | ChemDraw's files | Ours |
|---|---|---|
| Header: root attributes, color and font tables, page | the same at each scale in every reference (`test_header_matches_other_references`) | `presets/cdxml/<preset>.xml`, copied from m1_set_x\<k\> by `scripts/cdxml_headers.py` without CreationProgram, Name, BoundingBox and the page id |
| CreationProgram | "ChemDraw 26.1.0.6801" | "acs-chem-svg \<version\>" (user decision 2026-10-28) |
| Page | 540 × 720 per page; a drawing wider than that makes ChemDraw offer to expand, and the re-save has `BoundingBox="0 0 1080 720" WidthPages="2"` | whole pages, enough for the drawing |
| Labeled atom | `Element` (none for C), `NumHydrogens` (also 0), `NeedsClean="yes"`; one `<s face="96">` with the text in CDXML order (NH2 even when drawn H₂N) | the same |
| Label justification | `LabelJustification="Right"` + `Justification="Right"` when the bonds leave to the right and the text has two or more characters (H₂N, HO, –Cl, –Br); a one-character label (–F, –H, –I) stays Left | from `Label.h_side` and the text length |
| NH with H above / below (nh_x\*) | `LabelJustification="Left"`, no `Justification`, `LabelAlignment="Above"` / `"Below"`, `LineStarts="2 4"` / `"2 3"`; text `NH`; p.x at N's glyph origin; p.y the usual for Below, 8.37 / 6.74 / 5.30 higher for Above (−4.78 / −3.87 / −3.15 from the atom; `cdxml.h_above_rise`) | the same (`measure_refs.py --nh`, `tests/test_nh_vs_chemdraw.py`) |
| `LabelAlignment` | on 140 of 144 labels of two or more characters, never on one character (the other 4: a charged O, in mech_basic and labels_x1.0) | on labels of two or more characters |
| Label `<t p>` x | the element's glyph origin; right-justified, the end of the element's advance (–Cl: +1.11 at ×1.0 = half the l; H₂N: +3.61 = half the N) | the same, within 0.010 pt |
| Label `<t p>` y | atom y + half the ink top of the label's full-size glyphs: 0.3684 × size for O, C, S, Cl (an overshooting glyph anywhere: NO2, TMS, OMe, as §13.3); 0.3586 for N, H, F, I, Br, B, T; 0.3599 for Ph (the h). The drawn baseline is 0.3585 × size for all, so for O the stored p sits 0.095 pt below the drawn text | the same, from the font's glyph boxes, within 0.010 pt |
| Label BoundingBox | the label's ink box, offset by the same stored − drawn amount as p | the same, within 0.010 pt |
| Fragment BoundingBox | the drawn ink (half line widths, second lines of double bonds, charge signs) with the stored label boxes; dots: §13.2 | the same, within 0.010 pt (where a sign sets an edge, within our sign's placement, ≤ 0.39 pt) |
| Double bonds | no `DoublePosition` in any reference: ChemDraw places them itself | none |
| Charge | `<graphic GraphicType="Symbol" SymbolType="Plus|Minus">` + `<represent attribute="Charge" object=node>`; BoundingBox = (drawn sign center, center x − dx, center y), the center within 0.004 pt of the drawn bars | `cdxml.charge_box_dx` |
| `cdxml.charge_box_dx` | 7.5 / 4.8 / 2.7 = 0.075 × LabelSize² in 19 of 20 symbols (labels_x1.0, m1_set_x\*, mech_basic_x\*, resonance_x\*); one in resonance_x1.0 is −7.5 | 7.5 / 4.8 / 2.7 (`nonlinear`) |
| Ids | sequential; font ids are their own namespace (clearance_x1.0 has font 21 and bond 21) | sequential from 1 |

**DTD:** `reference_renders/cdxml/CDXML.dtd`, Revvity's file, fetched from the DOCTYPE URL
(2026-10-28; ChemDraw.app ships only hotkeys.dtd). It is local only and never committed (user
decision 2026-10-28), so the DTD tests skip without it. It dates from 2003. libxml2 rejects one of its enumerations (`AS` lists
`+ - ?`, which aren't name tokens), so `tests/cdxml_helpers.py` loosens it to CDATA. Every
reference then fails on only two attributes the DTD doesn't declare: `ShowPrimeLabels` on
`<CDXML>` and `IsNickname` on `<t>`. The test allows exactly those.

**Not written in phase 1:** abbreviations (written since phase 4a, §13.3) and Lewis structures (refused), lone-pair and radical
dots (the node keeps `Radical`; written since phase 3, §13.2), circled signs and `charge_labels` texts (a plain Plus/Minus with
the formal charge), atom maps (no `AtomID`). Each is reported on stderr.

**Gate (2026-10-28, `scripts/cdxml_gate.py --diff` on the user's ChemDraw 26 re-saves of all
seven gate files):** ChemDraw opens them without warnings (apart from the page size, now fixed)
and keeps our label p, justification, label and fragment boxes and charge boxes, all within
0.01 pt. It draws the re-saved references with the same ink as its originals: baselines and
anchors within 0.005 pt, bond ends at labels within 0.041. The charge centers are ours (up to
0.9 pt from its originals: placement defaults). It adds only its own bookkeeping (AS, AtomID, BS,
BondOrdering, BondCircularOrdering, Geometry) and no DoublePosition. It recomputed one label:
alanine's OH on a bond vertical to 1e-4, which we had set HO (a tie-tolerance bug in
`draw_mol.h_side`, fixed). Editing works: atoms drag, Clean Structure runs, and charges are
recognized as charges (a charge moved too far from its atom raises ChemDraw's warning). The
charge does not follow a dragged atom, as with ChemDraw's own. **Decision:** later phases write
labels explicitly, as here.

### 13.1 Schemes (phase 2: M2; source: scheme_basic_x\*, scheme_multi_x\*, scheme_wrap_x\*, equilibrium_x\*, arrows_large_x1.0)

Measured with `measure_refs.py --cdxml-scheme` (captions against the glyphs drawn in each
reference SVG, `_caption_glyphs` / `text_lines`; arrows against Head3D / Tail3D) and held by
`tests/test_cdxml_scheme.py` (`test_scheme_reference_diff`: the 9 scheme_basic / scheme_multi /
equilibrium references rebuilt at ChemDraw's coordinates and compared object by object).

| Item | ChemDraw's files | Ours |
|---|---|---|
| Page order | fragments, `<t>`, `<graphic>`, `<scheme>`, `<arrow>` (scheme_wrap_x\* also has the user's `<group>`s first); `<curve>`s after the graphics (§13.2) | the same, no groups |
| Caption `<t p>`, Left | first glyph origin on the first baseline, within 0.005 pt (single lines, plus signs, numbers, stage lines; formula fragments' labels too) | the same, from the richtext runs |
| Caption `<t p>`, Center | `CaptionJustification="Center" Justification="Center"`; p.x = the first line's advance center (every line's center is the same, within 0.005), p.y its baseline. ChemDraw centers only its stacked conditions | stacked or not, every `above:` / `below:` side is one centered `<t>`; stage lines are one Left `<t>` each (user decision 2026-10-28: every condition is a caption, including formulas ChemDraw stores as a one-node fragment) |
| Caption BoundingBox | the ink box within 0.005 pt for one line. Stacked: the right edge is the advance end of the lines that end in a full-size glyph ("reflux", "(4:1)"), else the ink; the bottom is 0.04–0.32 below the ink (larger at ×0.6), unexplained | the ink, the stacked right-edge rule, ink bottom (`KNOWN`, 0.35) |
| Formula-face superscript (`face="96"` "2. H₃O+") | box top 0.354 × size above the drawn "+" (3.54 / 2.83 / 2.12) | we write explicit runs (face 64): the ink (`KNOWN`) |
| `<s>` runs | one per style: Helvetica (21, no face), Helvetica Oblique (307, face 2: *n*), Helvetica Bold (310, face 1: numbers); formulas as one run with face 96, which ChemDraw subscripts itself; size the nominal caption size | one per richtext run, face 1 / 2 plus 32 (sub) / 64 (sup), nominal size. The header lists fonts 307 and 310, copied from scheme_multi_x1.0 |
| Lines | one `<t>` with `\n` between lines; `LineHeight="auto"`; `LineStarts` = each line's end offset in UTF-8 bytes with its newline ("…(4:1)\n90 °C, 12 h": 22 34, ° counting 2) | the same for a plain stack; a stack with a sub- or superscript is one `<t>` per line (gate below) |
| Line pitch | drawn 10.05 (no subscripts), 10.825 (subscript in the upper line), 10.975 (both) at ×1.0, for its formula face. For explicit script runs (ours) ChemDraw sets 11.775 whatever the scripts (gate) | ours (`text.line_pitch_*`) are fit to its formula face, so script stacks are split |
| `Warning` | "Chemical Interpretation is not possible…" / "ChemDraw can't interpret this label." on most captions | not written (ChemDraw's own note) |
| Arrow | `<graphic GraphicType="Line" SupersededBy=arrow ArrowType="FullHead|Equilibrium" HeadSize>` (BoundingBox = Head3D, Tail3D) + `<arrow>` with the same Z; steps cite the `<graphic>` | the same |
| Head attributes | `HeadSize / ArrowheadCenterSize / ArrowheadWidth` 1000 / 875 / 250 in every scheme and equilibrium reference at every scale, 1500 / 1313 / 375 in arrows_large_x1.0; ×LineWidth/100 = head_length and head_inset (`test_cdxml_head_attrs_match_preset_heads`) | `cdxml.arrow_heads` per `arrow.heads` variant |
| Equilibrium | `ArrowheadHead="HalfLeft" ArrowheadTail="HalfLeft"`, `ArrowShaftSpacing="300"` (= shaft spacing / LineWidth × 100), unbalanced `ArrowEquilibriumRatio="300"` (= long / short × 100) | derived from `equilibrium.shaft_spacing` and `ratios.equilibrium_long_per_short` |
| Head3D / Tail3D | forward: the drawn tip and tail. Balanced: our drawn ends moved in by `equilibrium.stagger` (box x reaches 0.90 beyond them). Favoring products: the ends. Favoring reactants: Head3D on the **left**, and the step's reactants are the species on the right | the same (the spec's `from` becomes ChemDraw's products there; the equilibrium is the same) |
| Axis points | Center3D the midpoint; MajorAxisEnd3D = center + (Head3D − Tail3D); MinorAxisEnd3D = center + that turned 90° (y down), within 0.01 pt in every arrow but arrows_large_x1.0's balanced one (minor axis 28.8 long: a copied arrow's) | the same |
| Arrow BoundingBox | x: the drawn ends. y: ± (ArrowheadWidth × LineWidth / 100, plus ArrowShaftSpacing × LineWidth / 200) with 0.5 more on the top of a forward arrow (both sides of ↔) and both sides of an equilibrium arrow whose Head3D is right, none when it is left; scheme_basic_x0.6 and scheme_wrap_x0.6 pad 0.18 instead | `cdxml.arrow_box_pad` 0.5 (`KNOWN` for scheme_basic_x0.6) |
| `<scheme>` | one per drawn sequence (scheme_multi: 3, equilibrium: 3); a wrapped sequence stays one (scheme_wrap: 5 steps, 1 scheme) | one per `new_row` |
| `<step>` | `ReactionStepReactants` / `Products` (fragment or group ids), `ReactionStepArrows` (the `<graphic>`), `ReactionStepObjectsAboveArrow` / `BelowArrow` (`<t>`, formula fragments, groups); empty lists left out | the same, with `<t>` for every condition; omitted products are not cited |

**Not written in phase 2:** abbreviations, so `m2_scheme_multi_x*`, `m2_sequence_x*` and
`m4_wrap_x*` were refused, and the scheme_multi_x\* round trip left out its OMe / OBn species;
Lewis structures; curved arrows, resonance steps and brackets. M3 is written since phase 3
(§13.2); abbreviations, resonance and those specs since phase 4a (§13.3).

**Gate (2026-10-28, `scripts/cdxml_gate.py --phase 2 --diff` on the user's ChemDraw 26 re-saves
of all seven files in `build/cdxml_gate/phase2/` (re-saves in `reference_renders/resaves/phase2/`)):**
- ChemDraw opens every file and keeps every `<scheme>` and `<step>`, each citing the same
  objects (it renumbers all ids, so `--diff` pairs objects by position). The exception is
  scheme_multi_x1.0's first row, whose round-trip export leaves out the four abbreviated
  structures (1, 3, 4, 5): with empty steps there, ChemDraw regroups that row into its own
  schemes.
- It keeps our arrows: every head and axis attribute, Head3D and Tail3D. It recomputes the
  `<arrow>` BoundingBox top 0.24 lower (pad 0.26, not the 0.5 stored in its references), so the
  box is bookkeeping it redoes.
- It keeps our caption p, runs, faces and LineStarts, and adds its `Warning`. It reads formula
  captions as chemistry: K₂CO₃, LiAlH₄, PCC, CH₂Cl₂ and H₂SO₄ become formula fragments, as in its
  own files, and H₃O⁺ becomes a hydronium fragment; the steps cite them in place of our `<t>`.
  It draws single lines exactly as its originals (glyph origins within 0.010 pt, sub- and
  superscripts included: our face 32 / 64 runs draw like its formula face). A superscript gets
  its 0.354 × size box top here too (3.54 on "3. H₃O+").
- **Line pitch:** a plain stack keeps 10.05 ("Acetone, reflux / 65%"). Stacks with explicit
  script runs are drawn 11.775 apart (PhB(OH)₂ / Pd(PPh₃)₄ 10.975, dioxane 10.825), 0.8–0.95 pt
  off. By the user's rule (more than 0.1 pt: one `<t>` per line), the writer splits stacks with a
  sub- or superscript into one `<t>` per line on our baselines; plain stacks stay one `<t>`.
- It re-justified a lone water from OH₂ to H₂O (Right, p at the end of the O advance). The
  renderer now draws a bondless O with H as H₂O (`draw_mol.LONE_H_LEFT`), and the export matches
  the re-save within 0.01 pt.
- Edit check (user): selecting an arrow does not show reactants and products, which is
  ChemDraw's standard behavior, not a fault of the files.

### 13.2 Mechanisms (phase 3: M3; source: mech_basic_x1.0 / x0.8 / x0.6)

Measured with `measure_refs.py --cdxml-symbols` (dot symbols against the circles drawn in each
reference SVG; each fragment's box against its outermost dots) and `measure_refs.py` (each
curve's anchors against its drawn stroke), and held by `tests/test_cdxml_mech.py`
(`test_mech_reference_diff`: the three references rebuilt at ChemDraw's coordinates, with its
CurvePoints and lone pairs on its atoms, and compared object by object).

| Item | ChemDraw's files | Ours |
|---|---|---|
| Curved arrow | `<curve CurveType ArrowheadHead ArrowheadType="Solid" CurvePoints>`, no BoundingBox; page-level after the `<graphic>`s (mech_basic's two fishhooks sit in a user `<group>` with tBu–Cl) | a free `<curve>`, no group; no step cites it (user decision 2026-10-28) |
| CurvePoints | (in-handle, anchor, out-handle) per anchor: 2e⁻ six points (2p0 − c1, p0, c1, c2, p3, 2p3 − c2), curve 40 to 0.01; the hand-drawn C–Br arrow (curve 91) has c2 = p3 | the same from `MechArrow` p0, c1, c2, p3: anchors and inner handles exact; end handles within 0.02 (we mirror its rounded points) |
| Last anchor | the SVG stroke end (the head's notch) in all 12 curves, within 0.002 pt; the tip is 4.05 pt beyond it (head_inset) (`test_curve_end_is_stroke_end`) | `MechArrow.p3`, as it is |
| 2e⁻ arrow | `CurveType="8" ArrowheadHead="Full"` | the same |
| Fishhook | `CurveType="32" ArrowheadHead="HalfRight"`, nine points: the six plus the second anchor and its handles again, a zero-length loop at the end (its SVG draws it as a second segment; with a 4.5 pt end handle it reaches 1.3 pt, under the head). All six reference fishhooks turn left on screen (y down) and are HalfRight | the six points only: with our longer end handles (up to 17.8 pt) the repeated anchor drew a loop up to 5.1 pt long and showed as a stray node in ChemDraw (user, gate 2026-10-28). `mechanism.fishhook_side` +1 → HalfRight (on all six, `test_fishhook_side_matches_reference`), −1 → HalfLeft (by symmetry; gate) |
| Lone pair | `<graphic GraphicType="Symbol" SymbolType="LonePair">` + `<represent attribute="Radical" object=node>`; BoundingBox = the two dot centers (to 0.001 pt), counterclockwise on screen around the atom (all 30), so their distance is `lone_pair.dot_spacing` (3.75 / 2.40 / 2.40) | the same from `Drawer.dot_sets()`; our dot centers are our placement (within 1.01 pt of its) |
| Radical dot | `SymbolType="Electron"`, the same `represent attribute="Radical"`; the node keeps `Radical="Doublet"`; BoundingBox = (dot center, center x − dx, center y) | `cdxml.electron_box_dx` |
| `cdxml.electron_box_dx` | 7.5 / 4.8 / 4.8 (all 6): not `charge_box_dx` at ×0.6 (2.7), nonlinear like the dot radius | 7.5 / 4.8 / 4.8 (`nonlinear`) |
| Fragment BoundingBox with dots | reaches 0.97 / 0.67 / 0.67 beyond the outermost lone-pair dot center (dot radius + 0.137) and 0.90 / 0.60 / 0.60 beyond a radical dot (+ 0.067; one per file, Cl·'s on the left) | `cdxml.lone_pair_box_reach`, `cdxml.electron_box_reach` (`nonlinear`), with the ink as before |
| Added reactant | (no reference) | a species a mechanism step adds, drawn above its arrow, is cited in `ReactionStepObjectsAboveArrow`, not as a reactant (user decision 2026-10-28) |
| `<step>` | cites the `<group>` holding tBu–Cl and its fishhooks as the reactant | cites the fragment |

Round trip, beyond the tolerances above: the tBu radical's dot is 2.8–5.0 pt from ChemDraw's
(ours takes the largest gap on the unlabeled C); the O⁻ and Br⁻ signs, which our dots push
out, are within 1.14 pt; mech_basic_x0.8 / x0.6 pad the arrow box 0.18, as scheme_basic_x0.6;
and the Cl· typed as text (no face 96) has its label box 0.1 higher.

**Not written in phase 3:** resonance steps and brackets (written since phase 4a, §13.3); Lewis structures'
shared and assigned electrons (refused); a second unpaired electron on an atom (warns). Atom
maps are not written either, but a mechanism doesn't warn about them: every species is mapped
for its arrows (user decision 2026-10-28).

**Gate (2026-10-28, `scripts/cdxml_gate.py --phase 3 --diff` on the user's ChemDraw 26 re-saves
of all six files in `build/cdxml_gate/phase3/` (re-saves in `reference_renders/resaves/phase3/`)):**
- ChemDraw keeps every `<curve>` (CurveType, ArrowheadHead, CurvePoints), every LonePair and
  Electron symbol with its box and `represent attribute="Radical"`, and every `<step>`, the
  added water above its arrow included.
- It draws our heads where we do: stroke ends within 0.005 pt, tips within 0.015, head shapes
  within 0.014 of our `MechArrow` (spec cases) or its original reference (round trips), barbs on
  the same side, HalfLeft included (m3_radical_propagation_x1.0). Every dot is drawn at the
  center its symbol declares (0.000).
- It recomputes the arrow box top 0.24 lower, as in phase 2, and some fragment boxes with dots
  and signs by up to 0.28: bookkeeping it redoes.
- **Fishhooks:** written with the references' nine points, our fishhooks showed an extra node
  near the end and looked distorted in ChemDraw (the repeated anchor's loop, above). The writer
  now writes six. Re-saved, ChemDraw keeps them unchanged and draws them as before, and on the
  mech_basic round trips exactly as its own nine-point originals (stroke ends 0.000, tips within
  0.011); the user found them clean to edit.
- Edit check (user): curves move as expected; dots stay behind when their atom is dragged and
  Clean Structure reshapes the structure, both as with ChemDraw's own drawings.

### 13.3 Abbreviations, resonance and wrapped schemes (phase 4a; source: abbrev_x\*, resonance_x\*, scheme_wrap_x\*, scheme_multi_x\*)

Measured with `measure_refs.py --cdxml-scheme` (now also resonance_x\*: ↔ arrows, and each
bracket `<graphic>` against the bracket drawn in the SVG) and the nickname nodes of abbrev_x\*,
and held by `tests/test_cdxml_species.py`: `test_reference_diff` on abbrev_x\* (all 21 keys in both
forms, three scales: nickname atoms, label p and box within 0.02 pt), `test_resonance_reference_diff`
and `test_scheme_reference_diff` on scheme_wrap_x\*.

| Item | ChemDraw's files | Ours |
|---|---|---|
| Abbreviation node | `<n NodeType="Fragment">` for formula-like groups (OMe, CO2Me, NO2, CF3, CN, NHAc, …), `NodeType="Nickname"` with `IsNickname="yes"` on the `<t>` for Ph, Bn, Et, Ac, TMS, TES, TBS, TIPS; `NodeType="Unspecified"` for PPh3 (metals_x\*, a typed label) and PMe3 (by analogy, its fragment PPh3's cut to the ipso carbons, unmeasured); no `Element`; `NeedsClean="yes"` | the same, NodeType and IsNickname from the template |
| Inner fragment | the group's atoms at real coordinates (Element, NumHydrogens, Charge: NO2's N+ and O−), plus a `NodeType="ExternalConnectionPoint"` 14.35 / 11.50 / 8.60 pt from the attachment atom toward the outer neighbor (0.05 pt short of BondLength) and bonded to it. The right-justified form is the left one turned 180°, to 0.01 pt | `presets/cdxml/<preset>.nicknames.xml`: each key's left-justified inner fragment from abbrev_x\<k\>, attachment atom at the origin (`scripts/cdxml_headers.py`), turned and moved rigidly onto the atom and its neighbor |
| Inner node attributes | abbrev_x\*: id, p, Element, NumHydrogens, Charge (a methyl C has NumHydrogens="3"). scheme_multi_x1.0's (expanded and contracted by the user): also Z, NeedsClean and ExternalConnectionNum, and no NumHydrogens on carbon | abbrev_x\*'s |
| Label `<t>` | the key text (OMe) in both forms, one face-96 `<s>`; Left: p.x the first glyph origin; Right (drawn MeO): `Justification="Right"`, p.x the end of the drawn text's advance (+3.89 = half the O); p.y as §13 (NO2, TMS: 3.68 for the O / S; Ph 3.60); BoundingBox the drawn text's ink in the stored frame | the same, from `Label.runs`, within 0.02 pt |
| RDKit | `Chem.MolsFromCDXML` expands both node types | round trips hold (`test_all_spec_round_trips`) |
| ↔ arrow | `<graphic ArrowType="Resonance" HeadSize>` + `<arrow ArrowheadHead="Full" ArrowheadTail="Full">`, 1000 / 875 / 250 at every scale; Head3D the right tip, Tail3D the left; axis points as §13.1 | `CdxArrow` kind `resonance` |
| ↔ BoundingBox | ± ArrowheadWidth × LineWidth / 100 plus 0.18 / 0.67 / 0.67 on both sides (resonance_x1.0 / x0.8 / x0.6; arrows_large_x1.0's large ↔: 0.5) | `cdxml.resonance_box_pad` (`approx`) |
| `<step>` | one per ↔, contributor k → k + 1, in the sequence's `<scheme>` | the same; a bracket continues its sequence's `<scheme>` (m4_nitration: step 3, two ↔ steps and step 5 in one) |
| Bracket | `<graphic GraphicType="Bracket" BracketType="Square" LipSize="60">`; BoundingBox the bar's x and its lips' stroke centers, within 0.005 pt: "[" (bottom, top), "]" (top, bottom). resonance_x1.0's "]" also has `<objecttag Name="bracketusage">`, and every file a `<bracketedgroup>` | free bracket graphics after the arrow graphics, no `<bracketedgroup>`, no objecttag (user decision 2026-10-28); `ResonanceLayout.bar_x` / `bar_y` |
| LipSize | 60 at every scale, while the drawn lips are 4.01 / 3.21 / 2.42 (0.073 × the bracket height in all three), so not lip_length / LineWidth | `cdxml.bracket_lip_size` (`exempt`) |
| Wrapped sequence | one `<scheme>` for both rows (scheme_wrap: 5 steps); its user `<group>`s hold a species with its number, or a formula caption with a caption cited elsewhere and the step's own arrow; steps cite the groups | one `<scheme>`, no groups (`test_wrapped_sequence_is_one_scheme`); `compare.cited_ids` reads a cited group as what it holds |

Round trip, beyond the tolerances: resonance_x0.8 / x0.6's C+ signs on unlabeled carbons are
1.52 pt below ChemDraw's (hand-placed, `charge.unlabeled_distance`), and resonance_x1.0's one
mirrored sign box (dx −7.5, §13) is ChemDraw's.

**Not written in phase 4a:** Lewis structures, circled signs, `charge_labels` texts, 2+ / 2−
and a second unpaired electron: written since phase 4b (§13.4).

**Gate (2026-10-28, `scripts/cdxml_gate.py --phase 4a --diff` on the user's ChemDraw 26
re-saves of all eight files in `build/cdxml_gate/phase4a/` (re-saves in `reference_renders/resaves/phase4a/`): the round trips of abbrev_x1.0,
resonance_x1.0 and scheme_wrap_x1.0 and the specs m1_abbrev, m4_resonance, m4_wrap, m4_nitration
and m2_sequence, x1.0):**
- ChemDraw keeps every node (336 on abbrev_x1.0), every nickname's NodeType, IsNickname and inner
  atoms and bonds (86 nicknames; inner atoms unmoved except one OBn in m2_sequence, 3.75 pt,
  the one the user expanded and contracted), every ↔ arrow and bracket, and every `<step>` and
  `<scheme>`; a wrapped sequence stays one scheme.
- Label p and boxes within 0.02 pt; arrow boxes 0.24 lower (FullHead, as in phase 2) and 0.08
  taller (↔): bookkeeping it redoes.
- It adds a `<bracketedgroup>` with two `<bracketattachment>`s around each bracket pair by
  itself, as in every reference: leaving it out costs nothing.
- scheme_wrap_x1.0's fourth step: the reference cites a hand-made group holding a caption above
  the arrow and one below as "below"; written as two captions, ChemDraw files the upper one
  above. Spec exports cite each caption on its own side and are unchanged.
- Edit check (user): nicknames expand and contract, contributors drag and brackets move as with
  ChemDraw's own drawings; no problems reported.

**General abbreviations (phase 4c; source: abbrev_2_x\*, ligands_x\*; 2026-11-03).**
`tests/test_cdxml_abbrev.py::test_abbrev2_reference_diff` holds abbrev_2_x1.0 and ligands_x1.0 (the
×0.8 / ×0.6 files' CDXML and SVG disagree, Open 32); `test_ring_nts_cdxml`,
`test_abbreviation_round_trip`, `test_text_label_cdxml`, `test_italic_label_runs`.

| Item | ChemDraw's files | Ours |
|---|---|---|
| NodeType | its own nicknames `Nickname` (Boc, Ts, Cbz, Fmoc, Me, *t*-Bu, *i*-Pr, *n*-Bu, Cy, Mes, N₃, as abbrev_x\*'s Ph … TIPS); formula labels and compositions `Fragment` (OTf, NHBoc, SnBu₃, ring NTs); labels typed on Pd `Unspecified` with the expansion inside (PMe₃, PCy₃, P(*t*-Bu)₃; ChemDraw read IPr as I + Pr); labels it can't interpret `Unspecified` with **no** inner fragment (L, Lₙ, dppf, XPhos, P(*o*-tol)₃, Bpin, *p*-Tol) | the template's; a generated one: registry `nickname: true` → `Nickname` + IsNickname, else `Fragment`; a text label `Unspecified`, no fragment. **Deliberate difference:** left-justified, P(*o*-tol)₃, Bpin, *p*-Tol keep their chemistry as a `Fragment` (RDKit reads them back); right-justified they are ChemDraw's `Unspecified` form (below) |
| Templates | — | `presets/cdxml/<preset>.nicknames.xml` from abbrev_x, abbrev_2_x and ligands_x (left-justified forms; OTs's right-justified one turned 180°) and metals_x: 55 keys. PMe₃'s is now ChemDraw's (ligands_x\*: methyls on the cross), not PPh₃'s cut down. Z, CrossingBonds, BondCircularOrdering, NeedsClean and AS / BS / AtomID are dropped |
| Generated fragment | — | any group, formula or composition without a template (`cdxml_writer.generated_nickname`): RDKit 2D coordinates of its SMILES scaled to BondLength, attachment atom at the origin, ExternalConnectionPoint BondLength − 0.05 out on −x; heteroatoms Element and NumHydrogens, charges, bond orders. Gate 4c passed (below) |
| Two-bond label (ring NTs) | one ExternalConnectionPoint per bond, BondLength − 0.05 along it (14.35 at ×1.0), some with `ExternalConnectionNum` 1 / 2; the inner `<fragment ConnectionOrder="ecp1 ecp2">`; the node `BondOrdering="bond1 bond2"` (outer bonds in the same order) | the same, ExternalConnectionNum always; the group turned (or mirrored) to fit both bonds by least squares, each ECP then put exactly on its bond |
| Label `<t>` | key text; `<s>` runs by style: plain and count digits face 96, an italic prefix 98 (*t* \| -Bu), a typed letter subscript 32 + 2 (Lₙ's n: 34); NTs below / above: `LabelAlignment` Below / Above, `LineStarts` "2 4" / "3 5", p.y 3.59 below the node / 4.98 above it (`cdxml.tail_above_rise` 8.57 / 6.90 / 5.43) | the same (`cdxml_writer.label_runs`) |
| p of N₃ and Ts | right-justified N₃: p.x the end of the ₃ centered on the atom (not of the drawn text, 0.70 further); Ts left-justified: p and box 0.55 right of the drawn glyphs (T's kerned advance centered) | `Label.internal_dx` |

**Gate (passed 2026-11-03):** `scripts/cdxml_gate.py --phase 4c` writes abbrev_2_x1.0, ligands_x1.0,
m1_abbrev2_x1.0 (a template, generated fragments for MEM and the spec's own group, OTBS, ring
NTs, Ln and dppf, *t*-Bu, N₃) and label_variants_x1.0; the user re-saved them in ChemDraw 26.1.
`--diff`: every node paired; every nickname's inner atoms and bonds kept (abbrev_2 80, ligands
19, m1_abbrev2 8; largest inner move 0.000); ring NTs keeps both ECPs, ConnectionOrder and
BondOrdering; text and Unspecified labels keep NodeType, no fragment and their `<s>` faces;
label p and boxes move ≤ 0.07 pt (rounding). ChemDraw adds AS / BS / AtomID /
BondCircularOrdering, `Warning="ChemDraw can't interpret this label."` on L, Lₙ, dppf, XPhos,
(*o*-tol)₃P, and `Warning="An atom in this label has an invalid valence."` on PMe₃, PCy₃,
P(*t*-Bu)₃ (P with four bonds) and IPr / PrI: the same warnings as in its own ligands_x1.0
(16 + 20). It also gives the left-justified IPr an inner fragment (I + Pr, 5 atoms, 1 ECP), as in
its own file: ChemDraw's chemistry for IPr is iodopropyl, which a text label can't prevent.

**Stored text of right-justified labels (user review of the phase 4c files, 2026-11-03).**
ChemDraw reverses a right-justified label it reads (OTf, SnBu₃: the key is stored). For the rest,
the review showed two failures:
- Bpin stored as pinB, a Fragment with its inner fragment, Right: ChemDraw showed **Bnip**
  (it reversed the unknown text glyph by glyph). Its own typed pinB survives only as an
  Unspecified node with its Warning and no inner fragment.
- IPr written Left-justified at the drawn start: ChemDraw showed **PrI** (it re-justified the
  label from its bond and reversed I + Pr).

Resolved by `label_variants_x1.0.cdxml` (`cdxml_gate.py --phase 4c`, four candidate forms per
row; user, 2026-11-03): only row B shows IPr and only row C shows pinB. So, right-justified:
- a text label ChemDraw reads as a formula (IPr) stores the text whose token reversal is ours
  (PrI; `abbrev.chemdraw_stored`), Unspecified, Right;
- a label ChemDraw can't read (Bpin, P(*o*-tol)₃, *p*-Tol, a spec's own group) is ChemDraw's own
  form: Unspecified, `Warning="ChemDraw can't interpret this label."`, no inner fragment, the
  displayed text (pinB, (*o*-tol)₃P, XyzO). It has no chemistry in ChemDraw (no Expand Label);
  left-justified, the same labels keep their generated or template fragment.
`abbrev.chemdraw_tokens` approximates ChemDraw's reading: the registry's groups less *o*-tol and
*p*-Tol (measured as unknown). (o-tol)₃P in row C is assumed to behave as pinB (same form).

Also from that review: `cdxml_to_molfile.py` reads a label typed as its `reversed:` text (pinB)
as its key, so abbrev_2_x\*'s left Bpin is a template; and a label atom never carries RDKit's
radicals (the SnBu₃ template's Sn had drawn a carbene lone pair).

### 13.4 Charges and Lewis structures (phase 4b; source: charges_lewis_x1.0 / x0.8 / x0.6)

charges_lewis_x\* is the user's ChemDraw 26 edit (2026-10-28) of `cdxml_gate.py
--lewis-draft`. Row 1 drew our charge structures without charges, and the user added them with
ChemDraw's tools. Rows 2–4 drew our Lewis structures in ChemDraw's objects, and the user kept
them. Measured with `measure_refs.py --charges-lewis`.

| Object | ChemDraw (charges_lewis_x\*) | Written |
|---|---|---|
| Circled ⊕ / ⊖ | `<graphic SymbolType="CirclePlus\|CircleMinus">`, `<represent attribute="Charge">`, BoundingBox (center, center x − 7.5 / 6.0 / 4.5, center y); drawn ring 3.59 / 2.89 / 2.24, bars 3.93 / 3.27 / 2.60 | the same; dx `cdxml.circled_box_dx` |
| Symbol size | ChemDraw draws a symbol from its box dx: Plus / Minus bars 3.93 at dx 7.5, 3.27 at 6.0, 2.60 at 4.5, 2.73 at 4.8, 1.80 at 2.7; Electron dot radius dx / 9 (0.833, 0.667, 0.533, 0.5). A symbol the user adds gets 0.75 × LabelSize; m1_set_x\*'s and mech_basic_x\*'s have their own dx | plain signs and dots keep `charge_box_dx` / `electron_box_dx` (their drawn sizes match §2.2 and §4.2's dots) |
| 2+ / 2− | no symbol: the user typed "Fe2+" into the label (superscript run, face 64; `Charge="2"` kept), but "O2−" lost the atom ("ChemDraw can't interpret this label": no Element, no Charge; RDKit reads C) | `Charge` ±2 and one free `<t>` "2+" / "2−" at digit size with `InterpretChemically="no"`, its ink centered where the drawing centers digits and bars; no symbol (user decisions 2026-10-28 and 2026-10-29, after the gates below) |
| "+3" (charge_labels) | the user typed "Fe+2" into the label, `Charge="2"` | `Charge` from the SMILES; one free `<t>` "+3", `InterpretChemically="no"` |
| Two unpaired electrons (CH₂) | two Electron symbols bound to the node; ChemDraw sets `Radical="Singlet"`; RDKit reads two radical electrons | the same |
| Singlet carbene / nitrene (`[CH2]`, `Cl[C]Cl` not in `triplet:`) | no reference | its pair a LonePair bound to the node, plus `Radical="Singlet"` so ChemDraw adds no H; unmeasured, ungated (Open 23) |
| Lewis, dots for bonds | bonds `Visible="no"`: ChemDraw hides them and keeps the chemistry. Shared pairs: free LonePair symbols, kept | the same |
| Lewis, bond lines | ordinary bonds, lone pairs bound LonePairs | the same |
| Assigned electrons | an Electron near a node becomes a radical: ChemDraw set `Radical` on HCN's C (Singlet), N and H (Doublet), and RDKit reads `[H][C]#[N]` | pairs: free LonePairs; an odd electron: a filled oval, `<graphic GraphicType="Oval" OvalType="Circle Filled">`, BoundingBox (center, center x + r, center y), Center3D / MajorAxisEnd3D / MinorAxisEnd3D, r = `lone_pair.dot_radius` (unmeasured form, user decision 2026-10-28) |
| Oxidation states | the user drew "+" with the charge tool (bound: every H got `Charge="1"`, RDKit reads `[H+]`) and "2−" as free text | free `<t>`s ("+", "2−"), `InterpretChemically="no"`; `Charge` unchanged |

ChemDraw flags many Lewis atoms "This atom has an unusual number of lone pairs" (CO₂, nitrate,
HCN, in both forms), but keeps them. Molecules RDKit reads from the reference (row 1, the user's
charges): C[N+](=O)[O-], C[CH+]C, C=C[O-], [Fe+2] (twice), C (the broken O²⁻) and [CH2] with two
radicals.

**Gate, round 1 (2026-10-28, `scripts/cdxml_gate.py --phase 4b --diff` on the user's
ChemDraw 26 re-saves of m1_charges, m1_lewis, m1_lewis_lines and m1_lewis_assigned, x1.0, as
first written: digit signs and Lewis electrons as free symbols):**
- **ChemDraw ignores `<represent>` and binds every symbol to the nearest atom on open** (ties,
  such as a shared pair centered on its bond, stay free with "Chemical symbol is not associated
  with an atom"), then recomputes each atom's `Charge` and `Radical` from its symbols:
  - nitromethane's circled ⊕, 9.66 pt from the =O and 10.04 from the N, became the O's
    (RDKit: `CN(=[O+])[O-]`);
  - the free signs of Fe²⁺, O²⁻ and "+3" were bound, and their charges became ±1 (`[Fe+]`,
    `[O-]`);
  - every oxidation-state "+" made its H a H⁺, the "2−" signs made C and O −1, and HCN's
    single assigned electrons made radicals again.
- Lewis structures in both forms kept their chemistry (RDKit SMILES unchanged); off-center
  shared pairs were bound to the nearer atom ("unusual number of lone pairs" warnings only).
  Labels and boxes within 0.01 pt; every symbol unmoved.
- Fixed: every sign nearest its own atom (§2.2), digits and oxidation states as free text, odd
  electrons as ovals.

**Gate, round 2 (2026-10-29, the user's re-saves of m1_charges and m1_lewis_assigned as
fixed):**
- Every bound sign stayed on its atom (nitromethane's ⊕ included; 9.48 pt from N, 12.68 from
  =O) and every symbol and label stayed put.
- The ovals work: ChemDraw keeps all four at our centers and draws them at radius 0.83–0.84
  (`lone_pair.dot_radius` 0.833). It rewrites them as BoundingBox (axis end, center) and
  `OvalType="Filled"`. HCN's formal-charge view now keeps its chemistry.
- **Still wrong (open):**
  - ChemDraw drops `Charge` from an atom with no sign symbol: Fe²⁺ came back Fe and O²⁻ came
    back O (RDKit `[Fe]`, `O`).
  - It read each free "+" text next to an H as a charge: the text disappeared and the H became
    H⁺. The "2−" texts (U+2212) it left alone.
- Round 3 (2026-10-29, re-saves of m1_charges and m1_lewis_assigned): `InterpretChemically="no"`
  on the free charge texts (ChemDraw writes it on two scheme_wrap_x1.0 captions) **does not
  help**: ChemDraw still converted all six oxidation-state "+" texts beside H into charges
  (the text gone, the H a H⁺). "2+", "2−" and "+3" texts survive unchanged; O²⁻'s O gets "An atom
  in this label has an invalid valence" once its charge is dropped. Still open: a lone "+" or
  "−" text next to an atom.
- Not fixed, by decision (user, 2026-10-29: "SVG is the primary output, and CDXML is for the
  appearance, not the chemistry"): an atom whose charge is drawn as text (2+ / 2−, a
  `charge_labels` text, an oxidation state on a formally charged atom) has no sign symbol, so
  ChemDraw drops its `Charge`. ChemDraw's own form, a superscript label run ("Fe" + "2+"), keeps
  the charge but sits lower than our placement, and a true minus in it ("O2−") loses the atom.

## 14. Chairs and crossings (source: chair_x1.0 / x0.8 / x0.6, ChemDraw 26.1.0.6801)

chair_x\* (the user, 2026-10-29), in CDXML order: two bare chairs from ChemDraw's Conformers
templates (`chair`, `flipped`); methyls on all six axial positions; on all six equatorial
positions; β-D-glucopyranose; an axial OAc; then seven crossing cases (single over single at 30°,
60° and 90°, single over double, double over single, a wedge over a single, a crossing next to a
CH₃ label; x0.8 / x0.6 have the last one without its label). Each crossing case is two
fragments; the front bond was brought to front. `measure_refs.py --chair`, `--crossing`.
Templates: `templates/chair/` (`cdxml_to_molfile.py`; crossing fragments merged, front bonds by
Z in `<reference>__<n>.front.json`).

### 14.1 Chairs (exact: ChemDraw's templates)

| Item | x1.0 | x0.8 | x0.6 | Notes |
|---|---|---|---|---|
| `chair.vertices` (pt from vertex 1, y down) | (0, 0) (−7.20, 12.47) (6.65, 8.51) (20.57, 12.47) (27.78, 0) (13.93, 3.96) | ×0.8 within 0.01 | ×0.6 within 0.01 | fragment 1; bonds 14.40–14.48 (x1.0): three pairs of parallel bonds, the two "middle" bonds 14.47 / 14.48 |
| Flipped chair | the exact mirror of `chair.vertices`: vertex k mirrors vertex (4 − k) mod 6 (0-based) | same | same | 2026-10-29 (user: the flip should be an exact mirror; mirroring fragment 1, which the axial and equatorial fragments 3–4 use). ChemDraw's own flipped template (fragment 2; β-D-glucopyranose, fragment 5, is drawn on it) is **not an exact mirror**: (0, 0) (7.20, 12.47) (21.13, 8.82) (34.97, 12.47) (27.78, 0) (13.85, 3.65), its middle vertices 0.31 pt (x0.8: 0.25, x0.6: 0.19) farther apart in y than the mirror's (8.51, 3.96), bonds 14.31–14.40 against the chair's 14.40–14.48. The same at all three sizes, so it is the template, not the drawing. Not used: our flipped chairs differ from ChemDraw's by up to 0.31 pt at those two vertices (and their substituents) |
| Axial bonds | vertical, 14.400 long | 11.520 | 8.640 | fragment 3; up at vertices 1, 3, 5 of `chair` (the right tip, 5, up; the left tip, 2, down) |
| Equatorial bonds, direction per vertex | −165°, 165°, −120°, 15°, −15°, 60° (± 0.05) | same | same | fragment 4 (deg ccw from +x, y up). At the middle vertices (3, 6) exactly parallel to the ring bond one removed; at the other four ChemDraw's 15° snap, **0.85–0.98° off** that parallel (15.9° bonds). β-D-glucopyranose (flipped template) agrees with the mirror at all five substituted vertices (0.24–0.37° off parallel there). Length = bond length (14.39–14.40) |
| Ring O of a pyranose | flipped vertex 5 | 5 | 5 | fragment 5 (`chair.heteroatom_vertex`): the back right |
| Bond ends at labels on chair bonds | §2.1 model | | | OH at 15°, −60°, 90°, ring O at 120° / 14.7°, OAc at 90°: within 0.08 pt (`tests/test_chair_vs_chemdraw.py`) |
| Depth | — | — | — | ChemDraw's Z (CrossingBonds pairs): an axial bond up from a lower ring atom is in front of the upper ring bond it crosses (axial OAc too); the lower ring bonds are in front of an axial bond down from an upper ring atom. So: lower on the page = nearer, and a substituent takes its ring atom's depth |

Not measured, chosen (`chair.py`): which flip (A-values, the bulkier groups equatorial; Eliel
& Wilen, Table 11.7: chemistry data, not style); where the ring goes on the templates
(crowding, then a ring heteroatom on vertex 5 of the flipped template, carbocycles on the
unflipped one, then fewest crossings, then the first substituted ring atom farthest right);
second-shell atoms on a 120° zigzag (iPr, Et, CH₂OH); RDKit for the rest. β-D-glucopyranose's
C3 / C4 HO labels overlap at the default bond length, in ChemDraw's drawing (fragment 5) as in
ours; `validate.py` reports it.

### 14.2 Crossings (exact)

| Item | Value | Notes |
|---|---|---|
| Which bond is cut | the lower Z | all 30 crossings (10 per preset): the higher-Z bond is drawn whole, the other has a gap |
| Gap edge | a line parallel to the front bond's outline on each side | the rear line's two corners at a cut end lie the same distance from the front axis (spread 0.000 for lines, 0.04 along a wedge's slanted side), so the end is oblique except at 90° |
| Gap edge distance | front half-extent + `margin_width` × `ratios.crossing_clearance_per_margin` (1.0; measured 0.991–1.005) | half-extent: lw/2 for a single bond (30°, 60°, 74°, 90°); the outer line's **center** (s/2, not s/2 + lw/2) for a centered double in front; the outline half-width at the crossing for a wedge. A double behind: both lines cut on the same edges |
| Next to a label | the rear line keeps the stub between its label clip and the gap | chair_x1.0 H₃C case: a 0.52 pt stub left of the gap (ChemDraw's; we draw it too) |

Implemented as: `draw_mol.Drawer.crossings` (axes that cross, no shared atom; front: nearer
`Atom.z`, then `Molecule.front` (spec `front:`), then the later bond, as ChemDraw puts what is
drawn later on top), `strip_edges`, `cut_pieces` (per line, so each line of a rear double bond
gets its own oblique ends), `main_segments` / `chains` (a cut main line ends one chain and
starts another; the joins stay at the atoms), `cut_polygons` (a piece in a chain is stroked to
its inner corner plus a cap polygon to the outer corner; a loose piece is one polygon). Not
measured: a rear wedge (clipped on the strip edges) and a rear hashed bond (a hash inside the
strip is dropped); an offset double or a triple in front (its outermost line center on each
side). `geometry_checks.check_labels` reports bond lines of unrelated bonds that touch without
a gap. CDXML: bonds get Z in depth order, rear first (`cdxml_writer.z_order`), no
CrossingBonds written; our export of chair_x1.0's templates has the same front bond as
ChemDraw's file at all 10 crossings. Gate (`cdxml_gate.py --phase 5`): passed (re-saves of
2026-10-29 in `reference_renders/resaves/phase5/`; closed by the user 2026-11-05). chair_x1.0 pairs
all 83 nodes and has the same front bond at all 10 crossings. The m1_chair_x1.0 re-save does not
line up with the current export (cause not traced): menthol is in a different chair pose and the other species sit
0.3–0.6 pt away, so its diff pairs 3 of 41 nodes and compares no crossings. A chair's CDXML has no wedges (ChemDraw's own chairs have none either), so
ChemDraw reads its ring stereocenters as undefined.

| Check | Worst | Test |
|---|---|---|
| M1 items on the chair templates (bond lengths, line widths and offsets, labels, ends at labels, wedge) | 0.076 / 0.069 / 0.075 pt | `tests/test_chair_vs_chemdraw.py::test_exact_geometry_within_tolerance` |
| Cut ends, distance from the front axis, both corners | 0.014 / 0.006 / 0.010 pt | `tests/test_chair_vs_chemdraw.py::test_crossing_cuts_within_tolerance` |

## 15. Catalytic cycles (source: arcs_x1.0 / x0.8 / x0.6, arcs_eq_x1.0 / x0.8 / x0.6, ChemDraw 26.1.0.6801)

arcs_x\* (the user, redrawn 2026-10-30): ChemDraw arc arrows (the arc tool) at 60°, 120° and
180° (59.6–60.1° as drawn), each tail and tip level, drawn counterclockwise and clockwise (not
flipped: a flipped copy keeps a stale Center3D and its MinorAxisEnd3D turns to −r); x1.0 has them
with large and standard heads, x0.8 and x0.6 with standard heads. ChemDraw has no cycle layout and no arc equilibrium arrow, so everything else in
a cycle is our placement (§15.2). Method: `measure_refs.py --arc` (`arc_summary`).

### 15.1 Arc arrows (exact)

| Item | x1.0 | x0.8 | x0.6 | Notes |
|---|---|---|---|---|
| Circle | the center on the chord's bisector, R = chord / (2 sin(\|AngularSize\| / 2)), on the page's left of tail → head for a positive AngularSize (`measure_refs.arc_circle`): Center3D and the `<graphic>`'s stored center within 0.07 (Head3D, Tail3D and AngularSize are rounded) | same | same | AngularSize = the sweep, + counterclockwise on the page, − clockwise |
| Tip | Head3D, on the circle | same | same | |
| Notch | on the circle (within 0.005), the chord to the tip = `head_inset` (7.878–7.879 large, 5.250–5.251 standard; counterclockwise) | 5.247–5.255 | 5.228–5.230 | the arc from notch to tip is 0.005–0.014 longer |
| Head axis | the chord notch → tip (within 0.01°): the straight head, not bent, symmetric about the chord | same | same | 2.1–5.9° off the tangent at the tip (x0.6: up to 7.2°); the back corners stand 1.43–2.43 off the circle, the outer one farther |
| Head shape | the straight head (§3.1) of the same variant: length 9.00 / 6.00, half-width 2.31–2.32 / 1.55–1.57, the same back edge and 0.05 hairline (counterclockwise) | 6.00, 1.57 | 5.98, 1.56 | `arrow.heads` |
| Shaft | width 0.60, a stroked path between Tail3D (within 0.007) and the notch, always stroked counterclockwise: from the tail on a counterclockwise arc, from the notch on a clockwise one | same | same | `measure_refs.arc_shaft` |
| Clockwise heads | the same shape (inset / length 0.875) on the chord (within 0.01°), the tip on Head3D, but scaled: length / measured 0.980, 0.983, 0.972 (standard) and 1.031, 0.991, 1.028 (large) at 60°, 120°, 180° | 0.979, 1.011, 0.973 | 1.028, 1.010, 1.021 | no pattern found (not radius, sweep or chord; x0.6's counterclockwise heads also moved from 1.000 to 0.996 between two saves of the same arcs): ours is the measured head both ways |
| Shaft cubics | ≤ 90° pieces from the tail; the 90° ones with handle 0.551915 r, the last (to the notch) 1–3.5 % shorter than an exact arc's, no closed form found. Up to **0.33 pt off the circle** (180°; 0.15 at 60°, 0.015 at 120°) | 0.24 | 0.15 | ChemDraw's approximation of its own circle: ours is the circle (`geometry.arc_segments`, equal pieces ≤ 90°, error < 0.03 % of r) |
| CDXML `<arrow>` | `AngularSize` (− clockwise; one decimal, none when whole), `Center3D` = the center, `MajorAxisEnd3D` = center + (r, 0), `MinorAxisEnd3D` = center + (0, r) (y down) in both directions, Head3D, Tail3D; head attributes as §13.1 | same | same | `cdxml_writer._Writer.arc` (§15.3) |
| CDXML `<graphic>` | `GraphicType="Arc"`, `SupersededBy` the arrow, `BoundingBox` = Head3D then the true center, `AngularSize` 90 for the 59.6–60.1° arcs (120 and 180 for the others; −90, −120, −180 clockwise) | same | same | the 60° value is unexplained; ours writes the whole-degree sweep, at least 90 (fits all three; Open 24) |
| `<arrow>` BoundingBox | 0.46–0.47 out from the ink on the tail side and top; the bottom 1.3–5.9 **inside** the ink for 60° and 120° arcs | | | not modelled (ChemDraw recomputes boxes, §13.1 gate) |

Implemented as `arrows.arc_arrow`: the shaft is the arc from the tail to the notch, and the head
is the straight `head_marker` with its origin at the notch, turned in the marker frame by the
angle from the arc's end tangent to the chord (`arc_turn`, one `<marker>` id per turn, so SVG's
`orient="auto"` puts it on the chord).

Arc equilibrium arrows (arcs_eq_x\*, the user, 2026-10-30: balanced and favoring products at
60°, 120° and 180° both ways, standard heads at every scale and large at x1.0; and headless
arcs). `measure_refs.py --arc-eq` (`arc_eq_arrows`); `arrows.arc_equilibrium_arrow`:

| Item | Measured | Ours |
|---|---|---|
| Shafts | two concentric arcs at R ± 0.900 (`ArrowShaftSpacing` 300), the one pointing at Head3D on the left of travel (outside for a clockwise arrow) | the same |
| Balanced | **no stagger**: each shaft's butt on its end of the arc (Tail3D / Head3D, within 0.01), its notch a head inset from the other end, its tip at that end ± 0.3 (+ clockwise, − counterclockwise) | each shaft end to end, tip at the end |
| Favoring products | the long shaft as a balanced one; the short one centered, L / 3 butt to tip (`ArrowEquilibriumRatio` 300), at 1.13–1.31 inside the circle (clockwise) or 0.61–0.83 outside (counterclockwise), not ±0.90 | short shaft L / ratio, centered, at ±spacing / 2 |
| Half heads | the straight half head on its own shaft's chord, its tip and notch on the shaft's **center line** (tip offset = the shaft's, ±0.90) | the same (straight half heads have them on the inner edge, §3.6) |
| CDXML | one `<arrow>`: `ArrowheadHead` and `ArrowheadTail` `HalfLeft`, `ArrowShaftSpacing` 300, `ArrowEquilibriumRatio` 300 when unbalanced, after `AngularSize`; Head3D / Tail3D the ends on the center circle; `<graphic ArrowType="HalfHead">` | the same; favoring reactants written reversed (Head3D at its start), as straight ones |
| Headless arc | `<arrow>` without `ArrowheadHead` (keeps the size attributes), `<graphic>` without `ArrowType` and its AngularSize the **other sign** from the arrow's; stored slightly elliptical (axes 0.10–0.23 apart) | a circle; the same attributes |
| `<graphic>` AngularSize | the arc tool's: the smallest of 90 / 120 / 180 / 270 that holds the sweep (113.7° → 120, 74.3° → 90) | the same |
| Glitches (left out) | x1.0's large balanced arrows are stored as ellipses (minor axis 60.16 vs major 80.84 at 60°) and drawn short or off the circle (the user: a ChemDraw bug) | — |

| Check | Worst | Test |
|---|---|---|
| Shaft offset, butt and notch along the arc | 0.45 pt (x0.6's 60.6° short shafts: 1.16) | `tests/test_arc_vs_chemdraw.py::test_arc_equilibrium_vs_chemdraw` |
| Head points to ChemDraw's head ink | 0.65 pt | the same |
| CDXML attributes, rewritten on ChemDraw's circles | 0.1 pt (headless: 0.25) | `tests/test_cdxml_cycle.py::test_reference_arcs_rewritten` |

| Check | Worst | Test |
|---|---|---|
| Shaft start and notch (counterclockwise; clockwise: the start, and ChemDraw's notch on the circle) | 0.006 / 0.006 / 0.022 pt | `tests/test_arc_vs_chemdraw.py::test_arc_shaft_ends_and_circle` |
| Head vertices and controls (tip, corners, back edge, notch), counterclockwise | 0.053 / 0.053 / 0.073 pt (the mean half-width, §3.1) | `tests/test_arc_vs_chemdraw.py::test_arc_head_on_the_chord` |
| Clockwise: tip; head on its chord; shape, scaled ≤ 3.5 % | within 0.1 pt; 0.05°; inset / length within 0.005 | the same test |

### 15.2 Cycle layout (ours, unmeasured; `cycle.py`, presets' `cycle:` keys)

Rules settled on the spike renders (`scripts/cycle_spike.py`, retired 2026-11-03 once `cycle.py` replaced it; user decisions 2026-10-29: equal
free arcs, labels inside the circle, tails as tangent arcs). Values at x1.0; x0.8 / x0.6 are x1.0
× k rounded to the 0.5 grid (`tail_join_in`, `tail_join_out` and `tail_turn` unchanged; the spike's 0.38 / 0.62 joins rounded to the 0.1 grid).

| Rule | Key (x1.0) |
|---|---|
| Species box = its ink plus its compound number, centered on the circle; the first at 12 o'clock, then round in `direction` | — |
| The arc arrows end `clear` from a species' ink (Euclidean, so a box corner is rounded); boxes keep `clear` from each other | `cycle.clear` 9.0 |
| Equal free arcs: each θk is set so every arrow is the same length (a wide species takes more of the circle) | — |
| Radius: the smallest one where the arrow length ≥ `scheme.arrow_min_length` + `tail_room` × (tails on the step's busier side), for the most demanding step | `cycle.tail_room` 24.0 |
| Step labels inside, on the inward line through the arc's middle, every ink corner ≥ `label_gap` inside the circle; the radius grows until labels keep `label_clear` from species and the center text, and `label_pair_clear` (one caption em; 4.0 read as one line, user decision 2026-10-30) from each other | `cycle.label_gap` 6.0, `cycle.label_clear` 4.0, `cycle.label_pair_clear` 10.0 |
| Tails: an arc of `tail_radius` tangent to the circle from outside, turning `tail_turn`; entering ones join at `tail_join_in` of the arc (earlier ones `tail_room` further back), leaving ones at `tail_join_out` (later ones further on); a leaving tail is an arc arrow (§15.1). Tails don't change the radius | `cycle.tail_radius` 20.0, `cycle.tail_turn` 80°, `cycle.tail_join_in` / `_out` 0.4 / 0.6 |
| Tail text or structure `tail_text_gap` beyond the free end, along the tail; consecutive texts share a tail (", ") | `cycle.tail_text_gap` 2.5 |
| A crowded tail (path or item within `label_clear` of a species, label, earlier tail or the arrows) tries radius ×1.3 / ×1.6, turns 65° / 100° / 50° / 120° and joins ±0.06 / −0.12; still crowded: a `validate.py` problem | — |

### 15.3 CDXML export (`render.py --cdxml`, `cdxml_writer.cycle_records`)

Every arrow and tail is a ChemDraw arc arrow as §15.1: an `<arrow>` with its legacy
`<graphic GraphicType="Arc">`, on the drawn circle. Rewritten on ChemDraw's own circles, the
arcs_x\* and arcs_eq_x\* arrows and graphics come back within 0.1 pt (headless 0.25),
AngularSize within 0.1° (`tests/test_cdxml_cycle.py`):

| Item | Written as |
|---|---|
| Forward arrow, leaving tail | an arc with `ArrowheadHead="Full"`; a tail on its own circle (§15.2) |
| Equilibrium arrow | ChemDraw's arc equilibrium arrow (§15.1), one `<arrow>` on the center circle |
| Entering tail | a headless arc (§15.1), drawn from the free end to the join |
| `<arrow>` BoundingBox | the drawn ink (ChemDraw's is not modelled, §15.1) |
| Captions | step labels, tail texts, the center text and compound numbers, as §13.1 (Left) |
| `<scheme>` | one, with a `<step>` per arrow: reactants the step's species and its entering species, products the next and its leaving species, `ReactionStepArrows` the arc graphic, the label and tail texts above the arrow; an equilibrium favoring reactants with reactants and products swapped (Head3D at its start) |

A dative bond is a plain single bond (§13), so RDKit reads the four-bonded P back as PH and
ChemDraw notes "An atom in this label has an invalid valence" on the P. `Order="dative"` clears
the note but ChemDraw draws a small arrow on the bond (gate 2026-10-30, `dative_test.cdxml`):
not used (user decision 2026-10-30).

Gate (2026-10-30, `cdxml_gate.py --phase 6`, both cycle examples): ChemDraw keeps every arc
(Head3D / Tail3D within 0.04, its recomputed boxes within 0.26), the arc equilibria, the headless
tails and the centered Pd labels, and draws them as ours (the user, by eye). It discards our
`<scheme>` and builds its own from the drawing: a step per arrow, leaving tails included,
entering (headless) tails ignored, labels above or below by position; PhBr, H2O and H+ captions
become formula fragments (as K2CO3, §13.1). It re-justified proline's NH with H above (both
ring bonds leave N downward); we now draw it so (nh_x\*, §2).

Re-save (2026-11-01, the gallery's Suzuki and enamine cycles): the PMe3 nicknames keep their
inner atoms and bonds, with the same invalid-valence note as PPh3 (dative as single, above). The NH
keeps H above (`LabelAlignment="Above"`, `LineStarts`, `<t p>` within 0.01). Arcs again within
0.05, the `<scheme>` rebuilt as before. New: ChemDraw also turns the bold Roman-numeral compound
numbers (I, II, III) into formula fragments, reading "I" as iodine with an invalid-valence note.
They are still drawn as before, so this is accepted, as with the other captions.

Held by `tests/test_cycle.py` (equal arcs within 0.05 pt, θ0 = 90°, no species within `clear` of
an arrow, the radius growing with count, size and labels, tail text not moving the radius).
