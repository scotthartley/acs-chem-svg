"""A new style's preset values from one ChemDraw file, the style probe (acs-style.md §16).

Every preset value has one of three sources:
- settings: the CDXML root of the file (bond length, line width, label size, ...);
- probe: a specimen drawn in the file, measured with the same methods as the ACS presets;
- derived: carried over from acs1996 (dimensionless ratios, placement defaults scaled by the bond
  length ratio, chair vertices), and not re-verified.

`sources` maps each role to a file. new_preset.py build gives every role the one probe file;
`references()` gives the ACS references instead, so the methods can be checked against the committed
presets before any probe exists (tests/test_preset_derivation.py).
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass
from pathlib import Path

from lxml import etree

from acs_chem_svg.draw_mol import element_x
from acs_chem_svg.fontmetrics import load_face
from acs_chem_svg.geometry import length, sub

from .cdxml import charges_lewis, cdxml_symbols
from .core import REF_DIR, label_glyphs, load
from .labels import _tail_line, nh_summary
from .measurements import measure
from .mech import mech_summary
from .scheme import equilibrium_summary

BASE_PRESET = "acs1996"
BASE_BOND_LENGTH = 14.4  # acs1996's lengths.bond_length: scale = bond length / this
SUPPORTED_FONT = "Helvetica"  # fontmetrics, svg_writer and preview draw only this (acs-style.md §16, open)
RESTYLE_TOL = 0.1  # pt: specimens vs the file's settings; M1's tolerance
RULE_TOL = 0.05  # pt: a measured value vs its rule, in the report
ROLES = ("settings", "stereo", "charges", "circled", "mech", "arrows", "large", "equilibrium", "brackets", "nh", "tail")
# dimensionless placement defaults: not scaled
UNSCALED_PLACEMENT = ("curved_arrow.curvature", "curved_arrow.loop_curvature", "cycle.tail_join_in", "cycle.tail_join_out",
                      "cycle.tail_turn")
# Choices, not measurements: value when new_preset.py isn't given one (None: acs1996's). Large heads are a
# user decision for acs1996 alone, so the default is ChemDraw's own head.
OPTIONS = {"text.number_weight": None, "arrow.default_head": "standard", "scheme.wrap_width": None}
# Keys the settings and the specimens supply (read_settings, measure_specimens); the rest are derived or options.
SETTINGS_KEYS = ("font_family", "lengths.bond_length", "lengths.line_width", "lengths.bold_width", "lengths.margin_width",
                 "lengths.hash_spacing", "lengths.bond_spacing", "lengths.ring_inner_inset", "ratios.bond_spacing_fraction",
                 "angles.chain_angle", "text.label_size", "text.caption_size")
PROBE_KEYS = ("lengths.hash_line_width", "lengths.wedge_wide_end", "charge.bar_length", "charge.bar_thickness", "cdxml.charge_box_dx",
              "curved_arrow.stroke_width", "curved_arrow.head_length",
              "curved_arrow.head_half_width", "curved_arrow.head_inset", "curved_arrow.fishhook_barb", "curved_arrow.fishhook_offset",
              "lone_pair.dot_radius", "lone_pair.dot_spacing", "cdxml.electron_box_dx",
              *(f"arrow.heads.{v}.{k}" for v in ("standard", "large") for k in ("head_length", "head_half_width", "head_inset", "shaft_end_from_tip")),
              "arrow.shaft_width", "arrow.head_outline", "cdxml.arrow_heads.standard", "cdxml.arrow_heads.large",
              "equilibrium.shaft_spacing", "equilibrium.stagger", "equilibrium.back_control_shift",
              "bracket.stroke_width", "cdxml.bracket_lip_size",
              "lengths.label_h_line_pitch", "cdxml.h_above_rise", "lengths.label_tail_above_pitch", "cdxml.tail_above_rise")
# carried over unchanged, unmeasured: dimensionless, or constant in pt at every ACS scale
COPIED_PREFIXES = ("ratios.",)
COPIED = ("lengths.label_outline_offset_x", "lengths.label_outline_offset_y", "charge.digit_scale", "text.line_pitch_extra",
          "cdxml.arrow_box_pad", "cdxml.resonance_box_pad", "chair.axial_up", "chair.equatorial_deg", "chair.heteroatom_vertex")
# acs1996 x k, 3 decimals. ChemDraw draws a bracket's lip at 0.0726 x its bar's height, whatever the style (resonance_x*,
# gate re-saves): a probe's lip would report the probe's bracket, so the default keeps acs1996's proportion to the bond.
SCALED = ("text.line_pitch_sub_upper", "text.line_pitch_sub_lower", "bracket.lip_length")


class MissingSpecimen(ValueError):
    """The file lacks a specimen a value is measured from, or its specimens don't follow its settings."""


class UnsupportedFont(ValueError):
    """The file's label font isn't the one the renderer draws."""


@dataclass(frozen=True)
class Source:
    name: str
    ref_dir: Path

    @property
    def cdxml(self) -> Path:
        return self.ref_dir / f"{self.name}.cdxml"

    @property
    def label(self) -> str:
        return f"{self.name}.cdxml"


@dataclass(frozen=True)
class Value:
    value: object
    source: str  # settings | probe | derived | option
    note: str


def single_file(path: Path) -> dict[str, Source]:
    """Every role from one re-saved probe: <path>.cdxml and <path>.svg side by side."""
    path = Path(path)
    for ext in (".cdxml", ".svg"):
        if not path.with_suffix(ext).exists():
            raise MissingSpecimen(f"{path.with_suffix(ext)} not found: save the probe from ChemDraw as CDXML and as SVG")
    src = Source(path.stem, path.parent)
    return dict.fromkeys(ROLES, src)


def references(scale: str, ref_dir: Path = REF_DIR) -> dict[str, Source]:
    """Each role from the ACS reference that holds its specimens at `scale` ("1.0", "0.8", "0.6").
    Large heads exist only at x1.0 (arrows_large_x1.0)."""
    fam = {"settings": "m1_set", "stereo": "m1_set", "charges": "m1_set", "circled": "charges_lewis", "mech": "mech_basic",
           "arrows": "scheme_basic", "equilibrium": "equilibrium", "brackets": "resonance", "nh": "nh", "tail": "abbrev_2"}
    out = {role: Source(f"{f}_x{scale}", ref_dir) for role, f in fam.items()}
    out["large"] = Source("arrows_large_x1.0", ref_dir)
    return out


def _parse(src: Source):
    return etree.parse(str(src.cdxml), etree.XMLParser(load_dtd=False, no_network=True)).getroot()


def _median(values, what: str, src: Source, nd: int = 3) -> float:
    vals = [v for v in values if v is not None]
    if not vals:
        raise MissingSpecimen(f"{src.label}: no {what} found")
    return round(statistics.median(vals), nd)


# --------------------------------------------------------------------------- settings


def read_settings(sources: dict[str, Source]) -> dict[str, Value]:
    """The values the CDXML root states, and those that follow from them exactly."""
    src = sources["settings"]
    root = _parse(src)

    def num(key: str) -> float:
        v = root.get(key)
        if v is None:
            raise MissingSpecimen(f"{src.label}: root attribute {key} missing")
        return float(v)

    def val(v, attr: str) -> Value:
        return Value(v, "settings", f"{src.label} {attr}")
    fonts = {f.get("id"): f.get("name") for f in root.find("fonttable")}
    font = fonts.get(root.get("LabelFont"))
    if font != SUPPORTED_FONT:
        raise UnsupportedFont(f"{src.label}: label font {font!r}; only {SUPPORTED_FONT} is supported so far "
                              "(fontmetrics, svg_writer and preview hardcode it; acs-style.md §16)")
    bl, bs = num("BondLength"), num("BondSpacing")
    spacing = bl * bs / 100
    return {
        "font_family": val(font, "LabelFont, fonttable"),
        "lengths.bond_length": val(bl, "BondLength"),
        "lengths.line_width": val(num("LineWidth"), "LineWidth"),
        "lengths.bold_width": val(num("BoldWidth"), "BoldWidth"),
        "lengths.margin_width": val(num("MarginWidth"), "MarginWidth"),
        "lengths.hash_spacing": val(num("HashSpacing"), "HashSpacing"),
        "lengths.bond_spacing": val(round(spacing, 3), "BondLength x BondSpacing / 100"),
        "lengths.ring_inner_inset": val(round(spacing * math.tan(math.radians(30)), 3), "BondLength x BondSpacing / 100 x tan 30"),
        "ratios.bond_spacing_fraction": val(bs / 100, "BondSpacing / 100"),
        "angles.chain_angle": val(num("ChainAngle"), "ChainAngle"),
        "text.label_size": val(num("LabelSize"), "LabelSize"),
        "text.caption_size": val(num("CaptionSize"), "CaptionSize"),
    }


def check_restyled(sources: dict[str, Source], settings: dict[str, Value]) -> None:
    """The specimens must follow the file's settings: a probe whose settings changed but whose
    drawing did not would measure the wrong style."""
    src = sources["stereo"]
    m = measure(load(src.name, src.ref_dir))
    for key, got, what in (("lengths.bond_length", m["bond_length"], "bond length"),
                           ("lengths.line_width", m["line_width"], "line width")):
        want = settings[key].value
        if not got:
            raise MissingSpecimen(f"{src.label}: no bonds to measure the {what} on")
        if abs(got["median"] - want) > RESTYLE_TOL:
            raise MissingSpecimen(f"{src.label}: specimens drawn with {what} {got['median']}, settings say {want}: "
                                  "the settings were not applied to the drawing")
    doc = load(src.name, src.ref_dir)
    sizes = [g.size for a in doc.atoms.values() if a.label for g in label_glyphs(doc, a)[0]]
    want = settings["text.label_size"].value
    if sizes and abs(statistics.median(sizes) - want) > RESTYLE_TOL:
        raise MissingSpecimen(f"{src.label}: labels drawn at {statistics.median(sizes)} pt, settings say {want}: "
                              "the settings were not applied to the drawing")


# --------------------------------------------------------------------------- specimens


def _stereo(src: Source) -> dict[str, tuple]:
    st = measure(load(src.name, src.ref_dir))["stereo"]
    hashes = [s["line_thickness"]["median"] for s in st if s.get("line_thickness")]
    wedges = [s["wide_end_width"] for s in st if "Hash" not in s["display"]]
    return {"lengths.hash_line_width": (_median(hashes, "hashed bond", src), "hash polygon thickness"),
            "lengths.wedge_wide_end": (_median(wedges, "wedge bond", src), "wedge polygon width at the drawn end")}


def _charges(src: Source) -> dict[str, tuple]:
    """Plain signs only: a circled sign's bars are its own (charges_lewis)."""
    doc = load(src.name, src.ref_dir)
    m = measure(doc)
    rings = [s.center for s in doc.shapes if s.kind == "stroke" and s.r and s.center]
    plain = [c for c in m["charges"] if not any(length(sub(c["center"], r)) < 0.3 for r in rings)]
    dx = [s["dx"] for s in cdxml_symbols(src.name, src.ref_dir)["symbols"]]
    return {"charge.bar_length": (_median([c["bar_length"] for c in plain], "plain charge sign", src), "filled bars (measure charges)"),
            "charge.bar_thickness": (_median([c["bar_thickness"] for c in plain], "plain charge sign", src), "filled bars (measure charges)"),
            "cdxml.charge_box_dx": (_median(dx, "Plus / Minus symbol", src), "symbol BoundingBox dx (cdxml_symbols)")}


def circled_ring(src: Source) -> tuple[float, float]:
    """(box dx, ring radius) of the file's circled signs. They set no preset value: applying a style
    rewrites every sign's box to 0.075 x LabelSize^2, where a sign drawn fresh has 0.75 x LabelSize
    (acs-style.md §16). The ring still tests ring_radius_of_box on the new style."""
    sym = [s for s in charges_lewis(src.name, src.ref_dir)["symbols"] if s["symbol"] in ("CirclePlus", "CircleMinus")]
    return (_median([s["dx"] for s in sym], "circled sign", src), _median([s.get("ring_radius") for s in sym], "circled sign", src))


def _mech(src: Source) -> dict[str, tuple]:
    m = mech_summary(src.name, src.ref_dir)
    electrons = cdxml_symbols(src.name, src.ref_dir).get("electrons") or []

    def med(st, what):
        if not st:
            raise MissingSpecimen(f"{src.label}: no {what} found")
        return st["median"]
    full, half = m["full_head"], m["fishhook_head"]
    return {"curved_arrow.stroke_width": (med(m["curved_arrow_stroke"], "curved arrow"), "stroked cubic (mech_summary)"),
            "curved_arrow.head_length": (med(full["length"], "two-electron arrow"), "4-vertex head (mech_summary)"),
            "curved_arrow.head_half_width": (med(full["half_width"], "two-electron arrow"), "4-vertex head (mech_summary)"),
            "curved_arrow.head_inset": (med(full["inset"], "two-electron arrow"), "4-vertex head (mech_summary)"),
            "curved_arrow.fishhook_barb": (med(half["barb"], "fishhook arrow"), "3-vertex half head (mech_summary)"),
            "curved_arrow.fishhook_offset": (med(half["stroke_end_off_axis_toward_barb"], "fishhook arrow"), "stroke end off the head's axis (mech_summary)"),
            "lone_pair.dot_radius": (med(m["dot_radius"], "lone pair"), "circle paths (mech_summary)"),
            "lone_pair.dot_spacing": (med(m["dot_spacing"], "lone pair"), "center to center (mech_summary)"),
            "cdxml.electron_box_dx": (_median([e["dx"] for e in electrons], "radical (Electron symbol)", src), "symbol BoundingBox dx (cdxml_symbols)")}


def _forward(src: Source) -> dict[int, list[dict]]:
    """Forward arrows by HeadSize."""
    out: dict[int, list[dict]] = {}
    for a in measure(load(src.name, src.ref_dir))["arrows"]:
        if not a.get("equilibrium") and len(a["heads"]) == 1 and a.get("shaft"):
            out.setdefault(int(a["cdxml"]["HeadSize"]), []).append(a)
    return out


def _head(arrows: list[dict], prefix: str, src: Source) -> dict[str, tuple]:
    heads = [next(iter(a["heads"].values())) for a in arrows]
    return {f"{prefix}.head_length": (_median([h["length"] for h in heads], "arrowhead", src), "fill path"),
            f"{prefix}.head_half_width": (_median([(h["half_width_left"] + h["half_width_right"]) / 2 for h in heads], "arrowhead", src),
                                          "mean of ChemDraw's two sides"),
            f"{prefix}.head_inset": (_median([h["inset_from_tip"] for h in heads], "arrowhead", src), "tip to notch"),
            f"{prefix}.shaft_end_from_tip": (_median([a["shaft"]["u_end_from_tip"] for a in arrows], "arrow shaft", src), "shaft polygon")}


def _arrows(std_src: Source, large_src: Source) -> dict[str, tuple]:
    """The smallest HeadSize among `std_src`'s forward arrows is the standard head, the largest
    among `large_src`'s the large one."""
    by_size = _forward(std_src)
    large_by_size = by_size if large_src == std_src else _forward(large_src)
    if not by_size or not large_by_size:
        raise MissingSpecimen(f"{std_src.label}: no forward arrow with a shaft found")
    std, large = min(by_size), max(large_by_size)
    if std == large:
        raise MissingSpecimen(f"{std_src.label}: every forward arrow has HeadSize {std}: the probe needs a standard and "
                              "a large head (did the settings overwrite the large one?)")
    arrows = by_size[std]

    def attrs(group: list[dict]) -> list[int]:
        c = group[0]["cdxml"]
        return [int(c["HeadSize"]), int(c["ArrowheadCenterSize"]), int(c["ArrowheadWidth"])]
    out = _head(arrows, "arrow.heads.standard", std_src) | _head(large_by_size[large], "arrow.heads.large", large_src)
    out["arrow.shaft_width"] = (_median([a["shaft"]["width"] for a in arrows], "arrow shaft", std_src), "shaft polygon")
    outline = [a["head_outline_stroke"]["median"] for a in arrows if a.get("head_outline_stroke")]
    out["arrow.head_outline"] = (_median(outline, "arrowhead outline", std_src), "hairline over the fill")
    out["cdxml.arrow_heads.standard"] = (attrs(arrows), "HeadSize, ArrowheadCenterSize, ArrowheadWidth")
    out["cdxml.arrow_heads.large"] = (attrs(large_by_size[large]), "HeadSize, ArrowheadCenterSize, ArrowheadWidth")
    return out


def _equilibrium(src: Source) -> dict[str, tuple]:
    recs = [a for a in equilibrium_summary(load(src.name, src.ref_dir))["arrows"] if "spacing (center to center)" in a]
    balanced = [a for a in recs if a["kind"] == "balanced"]
    return {"equilibrium.shaft_spacing": (_median([a["spacing (center to center)"] for a in recs], "equilibrium arrow", src), "shaft center to center (equilibrium_summary)"),
            "equilibrium.stagger": (_median([a["stagger_tail"] for a in balanced], "balanced equilibrium arrow", src), "butt behind the other tip (equilibrium_summary)"),
            "equilibrium.back_control_shift": (_median([a["back_control_shift (forward - back)"] for a in recs], "equilibrium arrow", src), "half-head control shift (equilibrium_summary)")}


def _brackets(src: Source) -> dict[str, tuple]:
    br = measure(load(src.name, src.ref_dir))["brackets"]
    lips = [int(g.get("LipSize")) for g in _parse(src).iter("graphic") if g.get("GraphicType") == "Bracket" and g.get("LipSize")]
    return {"bracket.stroke_width": (_median([b["stroke_width"] for b in br], "bracket", src), "stroked segments"),
            "cdxml.bracket_lip_size": (int(_median(lips, "bracket LipSize", src)), "LipSize")}


def _first_line_y(nodes) -> float:
    """The usual first-line <t> p.y relative to its node: the commonest value (0.01 pt) over labels
    not set above. A median would drift with the many labels stored another way (abbrev_2_x1.0)."""
    ys = [round(float(n.find("t").get("p").split()[1]) - float(n.get("p").split()[1]), 2) for n in nodes
          if n.find("t") is not None and n.find("t").get("LabelAlignment") != "Above"]
    return statistics.mode(ys)


def _nh(src: Source) -> dict[str, tuple]:
    nh = nh_summary(src.name, src.ref_dir)
    if not nh["line_pitch"] or "Above" not in nh["t_dy"]:
        raise MissingSpecimen(f"{src.label}: no NH with its H above N found")
    usual = _first_line_y(_parse(src).iter("n"))
    return {"lengths.label_h_line_pitch": (nh["line_pitch"]["median"], "H above / below N, baseline to baseline (nh_summary)"),
            "cdxml.h_above_rise": (round(usual - nh["t_dy"]["Above"]["median"], 2), "<t> p.y against the usual first line (nh_summary)")}


def _tail(src: Source) -> dict[str, tuple]:
    doc = load(src.name, src.ref_dir)
    pitches = [_tail_line(doc, label_glyphs(doc, a)[0][0], "Above")["pitch"] for a in doc.atoms.values()
               if a.nickname and a.alignment == "Above"]
    root = _parse(src)
    # a composed label's tail above its element: its <t> has LineStarts, unlike an NH's two-letter case
    above = [float(n.find("t").get("p").split()[1]) - float(n.get("p").split()[1]) for n in root.iter("n")
             if n.get("NodeType") in ("Fragment", "Nickname") and n.find("t") is not None
             and n.find("t").get("LabelAlignment") == "Above"]
    return {"lengths.label_tail_above_pitch": (_median(pitches, "ring NTs with Ts above N", src), "baseline to baseline (_tail_line)"),
            "cdxml.tail_above_rise": (round(_first_line_y(root.iter("n")) - _median(above, "ring NTs with Ts above N", src), 2),
                                      "<t> p.y against the usual first line")}


def measure_specimens(sources: dict[str, Source]) -> dict[str, Value]:
    """Every probe value. Raises MissingSpecimen if a specimen isn't found: nothing falls back."""
    groups = [(sources["stereo"], _stereo(sources["stereo"])), (sources["charges"], _charges(sources["charges"])),
              (sources["mech"], _mech(sources["mech"])),
              (sources["arrows"], _arrows(sources["arrows"], sources["large"])),
              (sources["equilibrium"], _equilibrium(sources["equilibrium"])), (sources["brackets"], _brackets(sources["brackets"])),
              (sources["nh"], _nh(sources["nh"])), (sources["tail"], _tail(sources["tail"]))]
    return {k: Value(v, "probe", f"{src.label}: {note}") for src, g in groups for k, (v, note) in g.items()}


# --------------------------------------------------------------------------- derivation


def lewis_bond_length(label_size: float, dot_spacing: float, dot_radius: float, gap: float) -> float:
    """lewis.bond_length (the base preset's comment, tests/test_lewis.py): the smallest 0.5 pt value where
    three shared-pair columns fit between two Cl or Br labels with 2 x lone_pair.gap to each box."""
    face = load_face("Regular")
    need = 0.0
    for el in ("Cl", "Br"):
        b = face.text_ink_bounds(el, label_size)
        ex = element_x(face, el, label_size, "right")
        need = max(need, (ex + b[2]) - (ex + b[0]) + 2 * dot_spacing + 2 * dot_radius + 4 * gap)
    return math.ceil(2 * need) / 2


def ring_radius_of_box(dx: float) -> float:
    """ChemDraw's circled-sign ring from its box dx: 7.5 -> 3.59, 6 -> 2.89, 4.8 -> 2.39, 4.5 -> 2.24
    (charges_lewis_x*, style_probe_x0.8), all within 0.03 pt of this line."""
    return round(0.25 + 4 / 9 * dx, 3)


def bar_length_of_box(dx: float) -> float:
    """ChemDraw's sign bars from the box dx: exact on 7.5, 6, 4.8, 4.5 and 2.7 (charges_lewis_x*)."""
    return round(0.6 + 4 / 9 * dx, 3)


def scale_of(bond_length: float) -> float:
    return round(bond_length / BASE_BOND_LENGTH, 3)


def rounded_placement(v0: float, k: float) -> float:
    return math.floor(2 * k * v0 + 0.5) / 2


def derive(base: dict[str, object], placement: list[str], known: dict[str, Value]) -> dict[str, Value]:
    """The values `known` (settings, specimens, options) doesn't hold, from the flat base preset and
    the bond length ratio. Raises if a value has no rule."""
    k = known["lengths.bond_length"].value / BASE_BOND_LENGTH
    out: dict[str, Value] = {}
    dr, sp = known["lone_pair.dot_radius"].value, known["lone_pair.dot_spacing"].value
    circled = round(0.75 * known["text.label_size"].value, 3)

    def d(v, why: str) -> Value:
        return Value(v, "derived", f"derived from {BASE_PRESET} {why}, unmeasured")
    rules = {"cdxml.lone_pair_box_reach": d(round(dr + 0.137, 2), "lone_pair.dot_radius + 0.137 (three ACS scales)"),
             "cdxml.electron_box_reach": d(round(dr + 0.067, 2), "lone_pair.dot_radius + 0.067 (three ACS scales)"),
             "cdxml.circled_box_dx": d(circled, "rule: a circled sign drawn fresh, 0.75 x LabelSize (three ACS scales)"),
             "charge.circle_radius": d(ring_radius_of_box(circled), "rule: 0.25 + 4/9 x cdxml.circled_box_dx (charges_lewis_x*, the probe)")}
    for key, v0 in base.items():
        if key in known or key in ("name", "scale"):
            continue
        if rules.get(key):
            out[key] = rules[key]
        elif key == "lewis.bond_length":
            continue  # after lone_pair.gap
        elif key == "chair.vertices":
            out[key] = d([[round(k * x, 2), round(k * y, 2)] for x, y in v0], f"x {k:.3f} (a template, exact in the bond length)")
        elif key in UNSCALED_PLACEMENT:
            out[key] = d(v0, "(dimensionless placement default)")
        elif key in placement:
            out[key] = d(rounded_placement(v0, k), f"x {k:.3f}, on the 0.5 pt grid (placement default)")
        elif key.startswith(COPIED_PREFIXES) or key in COPIED:
            out[key] = d(v0, "(not scaled at any ACS scale)")
        elif key in SCALED:
            out[key] = d(round(k * v0, 3), f"x {k:.3f}")
        else:
            raise ValueError(f"no source for preset value {key!r}: measure it, or give it a rule")
    gap = out["lone_pair.gap"].value
    out["lewis.bond_length"] = d(lewis_bond_length(known["text.label_size"].value, sp, dr, gap),
                                 "rule: three shared-pair columns between Cl or Br labels (tests/test_lewis.py)")
    return out


def rule_checks(values: dict[str, Value], ring: tuple[float, float] | None = None) -> list[tuple[str, float, float, str]]:
    """Measured values against rules fitted to the ACS scales (acs-style.md §16): (key, measured, rule,
    rule text). `ring`: circled_ring's result. They set nothing; they test the rules on a new style."""
    v = {k: x.value for k, x in values.items()}
    dr, sp, size = v["lone_pair.dot_radius"], v["lone_pair.dot_spacing"], v["text.label_size"]
    checks = [] if ring is None else [("circled ring (the probe's own box)", ring[1], ring_radius_of_box(ring[0]), "0.25 + 4/9 x box dx")]
    return checks + [
            ("lone_pair.dot_spacing", sp, round(4.5 * dr, 3), "4.5 x dot_radius"),
            ("cdxml.electron_box_dx", v["cdxml.electron_box_dx"], round(2 * sp, 3), "2 x dot_spacing"),
            ("cdxml.charge_box_dx", v["cdxml.charge_box_dx"], round(0.075 * size ** 2, 3), "0.075 x LabelSize^2"),
            ("lengths.wedge_wide_end", v["lengths.wedge_wide_end"], round(1.5 * v["lengths.bold_width"], 3), "1.5 x BoldWidth"),
            ("charge.bar_length", v["charge.bar_length"], bar_length_of_box(v["cdxml.charge_box_dx"]), "0.6 + 4/9 x box dx"),
            ("arrow.heads.standard.head_length", v["arrow.heads.standard.head_length"],
             round(v["cdxml.arrow_heads.standard"][0] * v["lengths.line_width"] / 100, 3), "HeadSize x LineWidth / 100")]
