# Renderer spike: RDKit ACS1996 mode vs ChemDraw (M0, design.md §4.1)

To reproduce, run `uv run python scripts/spike_rdkit.py` (needs the local-only `reference_renders/`; `spike/` is gitignored). It writes `spike/metrics.json` and
these images:
- `spike/m1_x{1.0,0.8,0.6}_compare.png`: ChemDraw | RDKit | overlay.
- `spike/m1_x1.0_zoom.png`: per molecule at 10 px/pt.

## Setup

- **Versions:** RDKit 2026.09.1 `MolDraw2DSVG(-1, -1)`, with `SetACS1996Mode(opts, MeanBondLength(mol))`.
- **Scaled presets:** `fixedBondLength × k`, `fixedFontSize = 10k`.
- **Same input as ChemDraw:** each M1 molecule is rebuilt from the ChemDraw CDXML with
  identical atom coordinates, bond orders, charges, H counts, Kekulé bonds, and wedge/hash flags.
  `prepareMolsBeforeDrawing = False` stops RDKit from re-wedging or aromatizing.
- **Same measurements:** RDKit's SVG is mapped onto the ChemDraw atom positions by a
  least-squares fit (scale + translation; residual 0.000 pt). Then the bond-assignment and
  measurement code from `measure_refs.py` runs on both outputs.
- **Native values:** the RDKit numbers below are what RDKit actually drew (mapped value ÷ fit
  scale).
- **Precision:** RDKit writes coordinates to 0.1 px, so its numbers are good to about ±0.05 pt.

## Numeric diff (ChemDraw → RDKit)

| Parameter | ×1.0 CD | ×1.0 RDKit | ×0.8 CD | ×0.8 RDKit | ×0.6 CD | ×0.6 RDKit |
|---|---|---|---|---|---|---|
| Bond length | 14.40 | **14.08** | 11.52 | 11.52 | 8.64 | 8.64 |
| Line width | 0.60 | 0.60 | 0.60 | 0.61 | 0.60 | 0.60 |
| Ring double-bond spacing | 2.592 | 2.53 | 2.074 | 2.07 | 1.556 | 1.56 |
| Ring inner-line inset | 1.497 | 1.45 | 1.198 | 1.21 | 0.899 | 0.88 |
| Acyclic C=O spacing | 2.592 | 2.5–2.6 | 2.074 | 2.0–2.1 | 1.555 | 1.5–1.6 |
| Wedge wide end (at H label) | 3.00 | 2.81* | 2.40 | 2.41* | 1.80 | 1.78* |
| Hash count / pitch / last length | 3 / 2.56 / 3.00 | 3 / 2.72 / 2.52 | 2 / 3.98 / 2.40 | **3** / 2.20 / 2.12 | 2 / 2.86 / 1.80 | **3** / 1.59 / 1.49 |
| Bond end − ChemDraw bond end at labels (median, range) | — | **−0.72** (−1.09…−0.69) | — | −0.49 (−0.68…−0.36) | — | −0.15 (−0.35…−0.11) |
| Gap, bond end → element ink box (median) | 1.65 | **1.10** | 1.38 | 1.00 | 1.01 | 0.90 |
| O glyph ink height | 7.58 (Helvetica) | 7.30 (Roboto) | 6.06 | 5.80 | 4.55 | 4.40 |
| Right N=O bond in O=N⁺=O | 4.8 pt drawn | **1.8 pt stub** | — | **1.4 pt** | — | **0.7 pt** |

\* RDKit fills the wedge and strokes it with a 0.6 pt outline. The width shown includes the
stroke.

## Visual differences

- **Label clearance.** RDKit stops bonds a fixed padding from the label box. At ×1.0 the gap
  is about 0.55 pt smaller than ChemDraw's (1.10 vs 1.65), and the bonds visibly crowd O, OH and N.
  ChemDraw clips against the glyph hull plus MarginWidth (acs-style.md §2.1).
- **Charges.** RDKit puts the charge inline in the label run ("N⁺"). That squeezes the
  following bond to a stub (O=N⁺:O in the images). ChemDraw places ⊕ above-right as a separate
  0.5 pt-bar symbol that bonds don't clip against.
- **H placement.** Both write H₂N on the left and OH on the right. RDKit's subscript sits lower
  and smaller (see the H₂N overlay), and its glyphs are Roboto, not Helvetica.
- **Double bonds.** Spacing and centered/offset placement agree within RDKit's 0.1 px rounding.
  Ring insets agree (≤ 0.05 pt).
- **Wedges and hashes.** The wedge shape is close. Hash count and pitch differ at ×0.8 and ×0.6
  (3 lines vs ChemDraw's 2), and RDKit's last hash is shorter (2.52 vs 3.00 at ×1.0).
- **Bond length at ×1.0.** RDKit's flexicanvas rounds the canvas to whole pixels and rescales
  to fit, so ACS1996 mode drew 14.08 pt bonds, not 14.40. Separately, single molecules came out
  between 11.1 pt (propyne) and 14.0 pt (aspirin). The scale is therefore not under our control.
- **Joins.** RDKit strokes each bond separately with butt caps and miter joins between segments
  of the same bond. Ring corners look fine at this size.

## design.md §4.1 needs

| Need | RDKit MolDraw2DSVG | Custom renderer |
|---|---|---|
| Exact atom coordinates in the final SVG frame | `GetDrawCoords()` gives them, but the frame scale floats with the canvas (14.08 vs 14.40) | exact by construction |
| Label bounding boxes in the SVG frame | not exposed; must be recovered by parsing glyph paths (class `atom-N`) | from fontmetrics, exact |
| Editable `<text>` (design.md §1.5, §4.4) | **no**: labels are glyph outlines (0 `<text>` elements) | yes |
| Helvetica | ACS1996 mode uses built-in Roboto; a custom font needs a FreeType-loadable file | yes (Helvetica.ttc / Heros) |
| Clearance, charge placement, hash rule, wedge width | RDKit's own conventions; only coarse options | implement measured ChemDraw rules |
| Named groups, pt units, no transforms (§4.4) | post-processing needed | native |
| Mechanism overlays sharing one coordinate system | would need a second pass on parsed output | native scene graph |

## Recommendation

**Build the custom renderer (design.md §4.1 stands)**, and keep RDKit for chemistry: parsing,
sanitization, 2D coordinates, ring info, stereo perception and alignment.

RDKit gets the easy geometry right: line width, double-bond spacing and ring inset are within
0.05 pt at ×0.8 and ×0.6. But it misses the conventions this project exists for:
- **Label clearance:** 0.5–1.1 pt too short at ×1.0.
- **Inline charges** that destroy bonds.
- **Hash count** at the scaled presets.
- **Font:** Roboto instead of Helvetica.

It also can't satisfy three hard requirements:
- editable `<text>` output,
- a fixed 14.40 pt bond scale,
- exposed label boxes in the output frame.

Wrapping RDKit would mean post-processing its SVG to fix every one of these. That is more work
than drawing the bonds ourselves from the measured rules.
