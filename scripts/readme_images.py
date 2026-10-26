"""Render the README's example images into docs/images/ (PNG on white, so they read in GitHub's
dark mode too).

Usage: uv run python scripts/readme_images.py
"""

from __future__ import annotations

import io
import sys
from pathlib import Path

from PIL import Image

from acs_chem_svg.layout import render_spec
from acs_chem_svg.preview import rasterize
from acs_chem_svg.svg_writer import to_svg

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "images"
SPECS = {
    "synthesis": "examples/m2_sequence_x1.0.yaml",
    "nitration": "examples/benchmarks/nitration_accepted.yaml",
}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, spec in SPECS.items():
        svg, _ = to_svg(render_spec(ROOT / spec))
        im = Image.open(io.BytesIO(rasterize(svg, 3.0))).convert("RGBA")
        flat = Image.new("RGBA", im.size, "white")
        flat.alpha_composite(im)
        path = OUT / f"{name}.png"
        flat.convert("RGB").save(path, optimize=True)
        print(path.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
