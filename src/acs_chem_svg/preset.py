"""Typed access to presets/*.yaml (design.md §3, references/acs-style.md).

Values are used as written; nothing is derived from `scale` here. The `exempt`,
`nonlinear`, `assumed_exempt`, `approx` and `placement` lists document how each value
relates to the x1.0 preset (tests/test_m0.py checks them); the renderer does not read them.
`placement` keys are round defaults (0.5 pt grid) for hand-placed positions.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict

PRESET_DIR = Path(__file__).resolve().parents[2] / "presets"


class _Section(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class Lengths(_Section):
    bond_length: float
    line_width: float
    bond_spacing: float
    ring_inner_inset: float
    bold_width: float
    wedge_wide_end: float
    hash_line_width: float
    hash_spacing: float
    margin_width: float
    label_outline_offset_x: float
    label_outline_offset_y: float


class Ratios(_Section):
    bond_spacing_fraction: float
    ring_inner_inset_tan: float
    wedge_end_per_bold: float
    label_clearance_per_margin: float
    label_band_per_size: float
    label_baseline_per_size: float
    subscript_scale: float
    subscript_drop_per_size: float
    superscript_scale: float
    superscript_raise_per_size: float
    arrow_back_control_per_half_width: float
    equilibrium_long_per_short: float


class Angles(_Section):
    chain_angle: float


class Text(_Section):
    label_size: float
    caption_size: float
    number_weight: str
    line_pitch_extra: float
    line_pitch_sub_upper: float
    line_pitch_sub_lower: float


class Charge(_Section):
    bar_length: float
    bar_thickness: float
    unlabeled_distance: float
    labeled_offset_x: float
    labeled_offset_y: float
    above_label: float
    circle_radius: float  # circled signs (unmeasured)
    digit_scale: float  # 2+, +3: digit size / label_size (unmeasured)
    digit_gap: float  # digit ink to the sign bars (unmeasured)


class Head(_Section):
    head_length: float
    head_half_width: float
    head_inset: float
    shaft_end_from_tip: float


class Heads(_Section):
    standard: Head  # ChemDraw's default (HeadSize 1000)
    large: Head  # HeadSize 1500 (arrows_large_x1.0)


class Arrow(_Section):
    shaft_width: float
    head_outline: float
    default_head: Literal["standard", "large"]
    heads: Heads


class Equilibrium(_Section):
    shaft_spacing: float
    stagger: float
    back_control_shift: float


class CurvedArrow(_Section):
    stroke_width: float
    head_length: float
    head_half_width: float
    head_inset: float
    fishhook_barb: float
    fishhook_offset: float
    tail_gap: float
    bond_tail_gap: float
    head_gap: float
    curvature: float
    loop_curvature: float


class LonePair(_Section):
    dot_radius: float
    dot_spacing: float
    gap: float


class Radical(_Section):
    distance: float


class Lewis(_Section):
    bond_length: float


class Bracket(_Section):
    stroke_width: float
    lip_length: float
    gap_x: float
    gap_above: float
    gap_below: float


class Resonance(_Section):
    arrow_length: float
    arrow_gap: float


class Scheme(_Section):
    arrow_min_length: float
    conditions_pad: float
    multistep_indent: float
    conditions_above_baseline: float
    conditions_below_baseline: float
    center_above_axis: float
    arrow_gap: float
    plus_gap: float
    number_gap: float
    row_gap: float
    wrap_width: float


class Preset(_Section):
    name: str
    scale: float
    font_family: str
    exempt: list[str]
    nonlinear: list[str]
    assumed_exempt: list[str]
    approx: list[str]
    placement: list[str]
    lengths: Lengths
    ratios: Ratios
    angles: Angles
    text: Text
    charge: Charge
    arrow: Arrow
    equilibrium: Equilibrium
    curved_arrow: CurvedArrow
    lone_pair: LonePair
    radical: Radical
    lewis: Lewis
    bracket: Bracket
    resonance: Resonance
    scheme: Scheme

    def head(self, variant: str | None = None) -> Head:
        """Reaction, equilibrium and ↔ arrowhead: `variant` ("standard" | "large"), or the
        preset's arrow.default_head."""
        v = variant or self.arrow.default_head
        if v not in Heads.model_fields:
            raise ValueError(f"unknown arrow head {v!r}; use one of: {', '.join(Heads.model_fields)}")
        return getattr(self.arrow.heads, v)


def available() -> list[str]:
    return sorted(p.stem for p in PRESET_DIR.glob("*.yaml"))


@lru_cache(maxsize=None)
def load_preset(name: str = "acs1996") -> Preset:
    """Load presets/<name>.yaml, or a path to a preset file."""
    path = Path(name) if str(name).endswith((".yaml", ".yml")) else PRESET_DIR / f"{name}.yaml"
    if not path.exists():
        raise FileNotFoundError(f"Unknown preset {name!r}; available: {', '.join(available())}")
    return Preset.model_validate(yaml.safe_load(path.read_text()))
