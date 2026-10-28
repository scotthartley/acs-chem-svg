"""Regression gallery: approved renders in tests/gallery/ (CLAUDE.md, design.md §10).

Cases:
- <reference>.png: our drawing of each ChemDraw reference's templates, on ChemDraw's page
  (scheme references also with our arrows and captions at ChemDraw's positions).
- m1_smiles_x<k>.png: examples/m1_smiles_x<k>.yaml (RDKit layouts).
- m1_lewis_x<k>.png, m1_lewis_assigned_x1.0.png: examples/m1_lewis_*.yaml (Lewis structures).
- m1_abbrev_x1.0.png: examples/m1_abbrev_x1.0.yaml (abbreviations in both forms, from SMILES).
- m2_*.png: examples/m2_*.yaml schemes with our own layout (m2_sequence_x<k>: the M2 acceptance).
- mech_basic_x<k>.png: the mechanism references' templates with our curved arrows on ChemDraw's
  CurvePoints and our dots at its centers.
- m3_*.png: examples/m3_*.yaml mechanisms with our own layout (the M3 acceptance).
- resonance_x<k>.png: the resonance references' templates with our ↔ arrows and brackets on
  ChemDraw's.
- equilibrium_x<k>.png, arrows_large_x1.0.png: our equilibrium, forward and ↔ arrows on
  ChemDraw's, with the captions (no structures: the references' "B" is a stand-in atom).
- m4_*.png: examples/m4_*.yaml with our own layout (m4_nitration_x<k>: the M4 acceptance;
  m4_perkin_x1.0: omit / over).

- cdxml/m1_smiles_x1.0.cdxml, cdxml/m1_charges_x1.0.cdxml, cdxml/m2_scheme_basic_x1.0.cdxml,
  cdxml/m2_equilibrium_x1.0.cdxml, cdxml/m3_sn2_x1.0.cdxml, cdxml/m3_radical_propagation_x1.0.cdxml,
  cdxml/m3_carbonyl_addition_x1.0.cdxml, cdxml/m3_proton_transfer_x1.0.cdxml,
  cdxml/m1_abbrev_x1.0.cdxml, cdxml/m4_resonance_x1.0.cdxml, cdxml/m4_wrap_x1.0.cdxml,
  cdxml/m4_nitration_x1.0.cdxml: render.py --cdxml's
  output for those specs (cdxml_writer), compared as text after lxml C14N.

All PNGs are rendered at 4 px/pt with a transparent background. tests/test_gallery.py re-renders
every case and compares pixels (the CDXML cases, canonical text). Rewriting the gallery needs the user's approval of the new
images (CLAUDE.md), so --update is only for after that approval.

Usage: uv run python scripts/gallery.py [--review DIR] [--update]
"""

from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import measure_refs as mr  # noqa: E402
from compare import render_reference  # noqa: E402

from acs_chem_svg.cdxml_writer import spec_cdxml  # noqa: E402
from acs_chem_svg.layout import render_spec  # noqa: E402
from acs_chem_svg.preview import rasterize  # noqa: E402
from acs_chem_svg.svg_writer import to_svg  # noqa: E402

GALLERY = mr.REPO / "tests" / "gallery"
PX_PER_PT = 4.0
REFERENCES = ["m1_set_x1.0", "m1_set_x0.8", "m1_set_x0.6", "double_bonds_x1.0", "labels_x1.0",
              "scheme_basic_x1.0", "scheme_basic_x0.8", "scheme_basic_x0.6",
              "scheme_multi_x1.0", "scheme_multi_x0.8", "scheme_multi_x0.6",
              "mech_basic_x1.0", "mech_basic_x0.8", "mech_basic_x0.6",
              "resonance_x1.0", "resonance_x0.8", "resonance_x0.6",
              "abbrev_x1.0", "abbrev_x0.8", "abbrev_x0.6",
              "equilibrium_x1.0", "equilibrium_x0.8", "equilibrium_x0.6", "arrows_large_x1.0"]
SPECS = ["examples/m1_smiles_x1.0.yaml", "examples/m1_smiles_x0.8.yaml", "examples/m1_smiles_x0.6.yaml",
         "examples/m1_lewis_x1.0.yaml", "examples/m1_lewis_x0.8.yaml", "examples/m1_lewis_x0.6.yaml",
         "examples/m1_lewis_assigned_x1.0.yaml",
         "examples/m1_charges_x1.0.yaml", "examples/m1_charges_x0.8.yaml", "examples/m1_charges_x0.6.yaml",
         "examples/m1_abbrev_x1.0.yaml",
         "examples/m2_sequence_x1.0.yaml", "examples/m2_sequence_x0.8.yaml", "examples/m2_sequence_x0.6.yaml",
         "examples/m2_scheme_basic_x1.0.yaml", "examples/m2_scheme_basic_x0.8.yaml", "examples/m2_scheme_basic_x0.6.yaml",
         "examples/m2_scheme_multi_x1.0.yaml", "examples/m2_scheme_multi_x0.8.yaml", "examples/m2_scheme_multi_x0.6.yaml",
         "examples/m2_suzuki_x1.0.yaml", "examples/m2_stages_x1.0.yaml", "examples/m2_stages3_x1.0.yaml",
         "examples/m2_equilibrium_x1.0.yaml", "examples/m2_equilibrium_x0.8.yaml", "examples/m2_equilibrium_x0.6.yaml",
         "examples/m3_sn2_x1.0.yaml", "examples/m3_sn2_x0.8.yaml", "examples/m3_sn2_x0.6.yaml",
         "examples/m3_e2_x1.0.yaml", "examples/m3_carbonyl_addition_x1.0.yaml",
         "examples/m3_radical_initiation_x1.0.yaml", "examples/m3_radical_propagation_x1.0.yaml",
         "examples/m3_radical_propagation_2_x1.0.yaml", "examples/m3_proton_transfer_x1.0.yaml",
         "examples/m4_nitration_x1.0.yaml", "examples/m4_nitration_x0.8.yaml", "examples/m4_nitration_x0.6.yaml",
         "examples/m4_resonance_x1.0.yaml", "examples/m4_resonance_x0.8.yaml", "examples/m4_resonance_x0.6.yaml",
         "examples/m4_wrap_x1.0.yaml", "examples/m4_wrap_x0.8.yaml", "examples/m4_wrap_x0.6.yaml",
         "examples/m4_perkin_x1.0.yaml",
         "examples/benchmarks/nitration_accepted.yaml"]


CDXML_SPECS = ["examples/m1_smiles_x1.0.yaml", "examples/m1_charges_x1.0.yaml", "examples/m2_scheme_basic_x1.0.yaml",
               "examples/m2_equilibrium_x1.0.yaml", "examples/m3_sn2_x1.0.yaml", "examples/m3_radical_propagation_x1.0.yaml",
               "examples/m3_carbonyl_addition_x1.0.yaml", "examples/m3_proton_transfer_x1.0.yaml",
               "examples/m1_abbrev_x1.0.yaml", "examples/m4_resonance_x1.0.yaml", "examples/m4_wrap_x1.0.yaml",
               "examples/m4_nitration_x1.0.yaml"]


def spec_svg(path: Path) -> str:
    """The same layout function as scripts/render.py."""
    return to_svg(render_spec(path))[0]


def cases() -> dict[str, str]:
    """name -> SVG text for every gallery case. The ChemDraw-reference cases need
    reference_renders/ (local only, not in the public repo); without it they are left out."""
    out = {ref: render_reference(ref).svg for ref in REFERENCES if (mr.REF_DIR / f"{ref}.svg").exists()}
    for s in SPECS:
        out[Path(s).stem] = spec_svg(mr.REPO / s)
    return out


def cdxml_cases() -> dict[str, str]:
    """name -> CDXML text (render.py --cdxml) for every CDXML case."""
    return {Path(s).stem: spec_cdxml(mr.REPO / s, warn=lambda m: None) for s in CDXML_SPECS}


def canonical(text: str) -> bytes:
    """C14N of a CDXML text, insignificant whitespace dropped (the DOCTYPE isn't part of it)."""
    from lxml import etree
    root = etree.fromstring(text.encode(), etree.XMLParser(load_dtd=False, no_network=True, remove_blank_text=True))
    return etree.tostring(root, method="c14n")


def png(svg: str) -> bytes:
    return rasterize(svg, PX_PER_PT)


def alpha(png_bytes: bytes):
    import numpy as np
    from PIL import Image
    return np.asarray(Image.open(io.BytesIO(png_bytes)).convert("RGBA"))[:, :, 3].astype(float) / 255


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--update", action="store_true", help="overwrite tests/gallery (needs user approval)")
    ap.add_argument("--review", type=Path, help="write new or changed renders here, for the user to review")
    args = ap.parse_args(argv)
    GALLERY.mkdir(parents=True, exist_ok=True)
    for name, svg in cases().items():
        path = GALLERY / f"{name}.png"
        data = png(svg)
        if args.update:
            path.write_bytes(data)
            print(f"wrote {path.relative_to(mr.REPO)}")
        else:
            same = path.exists() and path.read_bytes() == data
            print(f"{name:22} {'identical' if same else 'DIFFERS or missing'}")
            if not same and args.review:
                args.review.mkdir(parents=True, exist_ok=True)
                (args.review / f"{name}.png").write_bytes(data)
    (GALLERY / "cdxml").mkdir(exist_ok=True)
    for name, text in cdxml_cases().items():
        path = GALLERY / "cdxml" / f"{name}.cdxml"
        if args.update:
            path.write_text(text)
            print(f"wrote {path.relative_to(mr.REPO)}")
        else:
            same = path.exists() and canonical(path.read_text()) == canonical(text)
            print(f"{'cdxml/' + name:22} {'identical' if same else 'DIFFERS or missing'}")
            if not same and args.review:
                (args.review / "cdxml").mkdir(parents=True, exist_ok=True)
                (args.review / "cdxml" / f"{name}.cdxml").write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
