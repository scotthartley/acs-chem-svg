"""A new style preset from one ChemDraw file (style-spec.md §16).

    uv run python scripts/new_preset.py probe [--out build/style_probe]
    uv run python scripts/new_preset.py build --name NAME --probe PATH/style_probe.cdxml
        [--wrap-width PT] [--default-head standard|large] [--number-weight bold|regular] [--dry-run]

`probe` writes style_probe.cdxml (presets/probe/style_probe.yaml, at acs1996). In ChemDraw, open it and
apply the target style's settings to the whole drawing. Then delete every charge sign, lone pair and
radical dot, and add a few fresh ones with ChemDraw's tools: applying a style rewrites their sizes, which
are measured. Save it as CDXML and as SVG, same name. Tidy nothing else: lost bold and italics in the
captions don't matter.
`build` measures that file and writes presets/NAME.yaml, presets/cdxml/NAME.xml and an empty
presets/cdxml/NAME.nicknames.xml (abbreviations then use generated fragments), and prints each value's
source and its departure from scaled ACS. It stops if a specimen is missing.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml
from lxml import etree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import measure_refs as mr  # noqa: E402
from cdxml_headers import header  # noqa: E402
from chemdraw_refs.probe import (  # noqa: E402
    BASE_PRESET, MissingSpecimen, OPTIONS, RULE_TOL, UNSCALED_PLACEMENT, UnsupportedFont, Value, check_restyled, circled_ring,
    derive, measure_specimens, read_settings, restyled_signs, rule_checks, scale_of, single_file,
)

from chem_svg_cdxml.cdxml_writer import HEADER_DIR  # noqa: E402
from chem_svg_cdxml.figure import spec_cdxml  # noqa: E402
from chem_svg_cdxml.preset import PRESET_DIR, Preset, load_preset  # noqa: E402

PROBE_SPEC = PRESET_DIR / "probe" / "style_probe.yaml"
LISTS = ("exempt", "nonlinear", "assumed_exempt", "approx", "placement", "derived")
DEPARTURE_PT = 0.5  # pt: the report flags values this far from scaled ACS (the placement tolerance)
PARSER = etree.XMLParser(load_dtd=False, no_network=True)


def flat(d: dict, prefix: str = "") -> dict:
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            out.update(flat(v, f"{prefix}{k}."))
        else:
            out[f"{prefix}{k}"] = v
    return out


def base_values() -> tuple[dict, dict]:
    """The base preset as loaded: its nested dict, and flat without the lists."""
    raw = yaml.safe_load((PRESET_DIR / f"{BASE_PRESET}.yaml").read_text())
    return raw, flat({k: v for k, v in raw.items() if k not in LISTS})


# --------------------------------------------------------------------------- probe


def large_head_on_second_forward_arrow(cdxml: str) -> str:
    """Turn the second forward arrow's heads large: a spec has one head variant, and the probe needs
    a specimen of each. The <arrow> and its legacy <graphic> both state the size."""
    root = etree.fromstring(cdxml.encode(), PARSER)
    tree = root.getroottree()
    forward = [a for a in root.iter("arrow") if a.get("ArrowheadHead") == "Full" and not a.get("ArrowheadTail")]
    if len(forward) < 2:
        raise ValueError("style_probe.yaml: expected two forward arrows")
    a = forward[1]
    size, center, width = load_preset(BASE_PRESET).cdxml.arrow_heads["large"]
    a.set("HeadSize", str(size))
    a.set("ArrowheadCenterSize", str(center))
    a.set("ArrowheadWidth", str(width))
    for g in root.iter("graphic"):
        if g.get("SupersededBy") == a.get("id"):
            g.set("HeadSize", str(size))
    return etree.tostring(tree, xml_declaration=True, encoding="UTF-8", doctype=tree.docinfo.doctype, pretty_print=True).decode()


def write_probe(out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "style_probe.cdxml"
    path.write_text(large_head_on_second_forward_arrow(spec_cdxml(PROBE_SPEC)))
    return path


# --------------------------------------------------------------------------- build


def collect(sources: dict, options: dict) -> dict[str, Value]:
    """Every preset value with its source (style-spec.md §16)."""
    raw, base = base_values()
    settings = read_settings(sources)
    check_restyled(sources, settings)
    known = settings | measure_specimens(sources)
    for key, default in OPTIONS.items():
        v = options.get(key)
        if v is not None:
            known[key] = Value(v, "option", "new_preset.py option")
        elif default is not None:
            known[key] = Value(default, "option", "default: ChemDraw's own")
        else:
            known[key] = Value(base[key], "derived", f"derived from {BASE_PRESET}, unmeasured")
    return known | derive(base, raw["placement"], known)


def _inline(v) -> str:
    return yaml.safe_dump(v, default_flow_style=True, width=10 ** 6).strip().removesuffix("...").strip()


def preset_yaml(name: str, values: dict[str, Value], source_label: str) -> str:
    """The preset's text, in the base preset's key order, each value with its source as a comment."""
    raw, _ = base_values()
    head = [f"# {name}: generated by scripts/new_preset.py from {source_label}; rebuild it rather than editing it.",
            "# Settings and probe values are measured; `derived` values are carried over from acs1996 (style-spec.md §16).",
            f"name: {name}",
            f"scale: {scale_of(values['lengths.bond_length'].value)}",
            f"font_family: {values['font_family'].value}  # {values['font_family'].note}"]
    for key in LISTS:
        items = {"placement": raw["placement"], "derived": sorted(k for k, v in values.items() if v.source == "derived")}.get(key, [])
        head += [f"{key}: {_inline(items)}" if not items else f"{key}:"] + [f"  - {i}" for i in items]

    def section(d: dict, prefix: str, indent: str) -> list[str]:
        out = []
        for k, v in d.items():
            key = f"{prefix}{k}"
            if isinstance(v, dict):
                out.append(f"{indent}{k}:")
                out += section(v, f"{key}.", indent + "  ")
            else:
                val = values[key]
                out.append(f"{indent}{k}: {_inline(val.value)}  # {val.source}: {val.note}")
        return out
    body = []
    for sec, content in raw.items():
        if isinstance(content, dict):
            body += ["", f"{sec}:"] + section(content, f"{sec}.", "  ")
    return "\n".join(head + body) + "\n"


def report(values: dict[str, Value], ring: tuple[float, float] | None = None) -> list[str]:
    """Each value, its source and its departure from scaled ACS; rules that fail on this style."""
    raw, base = base_values()
    k = values["lengths.bond_length"].value / base["lengths.bond_length"]
    constant = set(raw["exempt"]) | set(raw["assumed_exempt"]) | set(UNSCALED_PLACEMENT)
    rows = [f"bond length ratio k = {k:.3f}; ACS = what acs1996 predicts: constant if the ACS presets keep it constant, else x k",
            f"{'value':42s} {'new':>9s} {'ACS':>9s} {'diff':>7s}  source"]
    flagged = []
    for key, val in values.items():
        v0 = base.get(key)
        if isinstance(val.value, (int, float)) and isinstance(v0, (int, float)):
            expect = v0 if key in constant or key.startswith(("ratios.", "angles.")) else k * v0
            diff = val.value - expect
            far = abs(diff) > DEPARTURE_PT
            if far:
                flagged.append(key)
            rows.append(f"{key:42s} {val.value:9.3f} {expect:9.3f} {diff:+7.3f}  {val.source}{' <-' if far else ''}")
        else:
            rows.append(f"{key:42s} {_inline(val.value)[:40]}  {val.source}")
    counts = {s: sum(v.source == s for v in values.values()) for s in ("settings", "probe", "derived", "option")}
    rows.append("sources: " + ", ".join(f"{n} {s}" for s, n in counts.items()))
    rows.append(f"{len(flagged)} value(s) more than {DEPARTURE_PT} pt from ACS: {', '.join(flagged) or 'none'}")
    bad = [c for c in rule_checks(values, ring) if abs(c[1] - c[2]) > RULE_TOL]
    rows.append("rules fitted to the ACS scales: " + ("all hold" if not bad else "do NOT hold on this style:"))
    rows += [f"  {key}: measured {got}, {text} = {rule}" for key, got, rule, text in bad]
    if restyled_signs(values):
        rows.append("warning: the charge signs have the size applying a style gives, not a fresh sign's. If they weren't redrawn in the "
                    "style, charge bars and dots here are off: replace them and rebuild.")
    return rows


def write_preset(name: str, values: dict[str, Value], probe: Path | None, source_label: str,
                 preset_dir: Path = PRESET_DIR, header_dir: Path = HEADER_DIR, header_bytes: bytes | None = None) -> Path:
    path = preset_dir / f"{name}.yaml"
    text = preset_yaml(name, values, source_label)
    Preset.model_validate(yaml.safe_load(text))  # the schema, before anything is written
    if header_bytes is None:  # caption fonts the probe lost in the re-save: the committed base template's
        header_bytes = header(probe, font_source=HEADER_DIR / f"{BASE_PRESET}.xml")
    path.write_text(text)
    header_dir.mkdir(parents=True, exist_ok=True)
    (header_dir / f"{name}.xml").write_bytes(header_bytes)
    root = etree.Element("nicknames", {"source": "none: generated fragments (cdxml_writer.generated_nickname)"})
    (header_dir / f"{name}.nicknames.xml").write_bytes(etree.tostring(root, xml_declaration=True, encoding="UTF-8", pretty_print=True))
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("probe", help="write style_probe.cdxml")
    p.add_argument("--out", type=Path, default=mr.REPO / "build" / "style_probe")
    b = sub.add_parser("build", help="measure a re-saved probe and write the preset")
    b.add_argument("--name", required=True)
    b.add_argument("--probe", type=Path, required=True)
    b.add_argument("--wrap-width", type=float)
    b.add_argument("--default-head", choices=["standard", "large"])
    b.add_argument("--number-weight", choices=["bold", "regular"])
    b.add_argument("--dry-run", action="store_true", help="print the report, write nothing")
    args = ap.parse_args(argv)
    if args.cmd == "probe":
        path = write_probe(args.out)
        print(f"wrote {path}\nIn ChemDraw: open it and apply the target style's settings to the whole drawing. Delete every charge\n"
              "sign, lone pair and radical dot, then add fresh ones with ChemDraw's tools: three lone pairs (horizontal and vertical), a\n"
              "radical dot away from the others, and one each of +, -, circled + and circled -. Applying a style rewrites their sizes.\n"
              "Tidy nothing else. Save as CDXML and as SVG, same name. Then: new_preset.py build --name NAME --probe <saved>.cdxml")
        return 0
    if args.name in mr.ACS_PRESETS:
        print(f"error: {args.name} is a measured ACS preset; pick another name", file=sys.stderr)
        return 1
    options = {"scheme.wrap_width": args.wrap_width, "arrow.default_head": args.default_head, "text.number_weight": args.number_weight}
    try:
        sources = single_file(args.probe)
        values, ring = collect(sources, options), circled_ring(sources["circled"])
    except (MissingSpecimen, UnsupportedFont) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    print("\n".join(report(values, ring)))
    if not args.dry_run:
        path = write_preset(args.name, values, args.probe, args.probe.name)
        print(f"wrote {path.relative_to(mr.REPO)}, presets/cdxml/{args.name}.xml and .nicknames.xml")
    return 0


if __name__ == "__main__":
    sys.exit(main())
