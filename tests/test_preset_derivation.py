"""New style presets (scripts/new_preset.py, acs-style.md §16): hold-out against the committed ACS
presets. The derivation rules run on x0.8 and x0.6 with no reference; the measuring methods run on
their references (skipped without reference_renders/), each role from the file that holds its
specimen; real probes of those settings, when saved, run the same comparison on one file."""

import shutil
import sys
from pathlib import Path

import pytest
import yaml
from lxml import etree

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import new_preset as np  # noqa: E402
from chemdraw_refs import probe  # noqa: E402

from acs_chem_svg import cdxml_writer  # noqa: E402
from acs_chem_svg import preset as preset_mod  # noqa: E402
from acs_chem_svg.figure import render_spec, spec_cdxml  # noqa: E402
from acs_chem_svg.preset import Preset  # noqa: E402
from acs_chem_svg.spec import load_spec  # noqa: E402

TARGETS = {"0.8": "acs1996_x0.8", "0.6": "acs1996_x0.6"}
MEASURED_TOL = 0.1  # M1's tolerance
PLACEMENT_TOL = 0.5
HAS_REFS = probe.REF_DIR.exists()
needs_refs = pytest.mark.skipif(not HAS_REFS, reason="no reference_renders/")
PROBES = REPO / "reference_renders" / "probes"
# Today's misses of the derivation rules:
KNOWN_DERIVED_MISSES = {"cdxml.resonance_box_pad"}  # 0.18 / 0.67 / 0.67 (not drawn; acs1996 yaml, approx)
# A measured value the committed presets don't follow: user decision 2026-10-24 scales ChemDraw's charge
# bars with k in the ACS family; a generated preset keeps what ChemDraw draws (user decision 2026-11-06),
# 2.733 / 1.798 against 3.148 / 2.361.
KNOWN_MEASURED_MISSES = {"charge.bar_length"}


def _committed(name):
    raw = yaml.safe_load((np.PRESET_DIR / f"{name}.yaml").read_text())
    return raw, np.flat({k: v for k, v in raw.items() if k not in np.LISTS})


def _diffs(values, committed):
    """key -> |generated - committed|; lists compare element by element, the largest difference."""
    out = {}
    for key, val in values.items():
        got, want = val.value, committed[key]
        if key == "chair.vertices":
            out[key] = max(abs(a - b) for p, q in zip(got, want) for a, b in zip(p, q))
        elif isinstance(want, bool) or not isinstance(want, (int, float)):
            out[key] = 0.0 if got == want else float("inf")
        else:
            out[key] = abs(got - want)
    return out


def _known(name):
    """The settings and probe values as the committed preset has them: the derivation rules alone."""
    raw, flat = _committed(name)
    known = {k: probe.Value(flat[k], "probe", "hold-out: committed value") for k in (*probe.SETTINGS_KEYS, *probe.PROBE_KEYS)}
    for key in probe.OPTIONS:
        known[key] = probe.Value(flat[key], "option", "hold-out: committed value")
    return known


@pytest.mark.parametrize("name", TARGETS.values())
def test_derivation_rules_hold_out(name):
    """Derived values within 0.5 pt (placement) or 0.1 pt of the committed preset, the known misses excepted."""
    base_raw, base = np.base_values()
    known = _known(name)
    derived = probe.derive(base, base_raw["placement"], known)
    assert known.keys() | derived.keys() == _committed(name)[1].keys() - {"name", "scale"}
    d = _diffs(derived, _committed(name)[1])
    off = {k for k, v in d.items() if v > (PLACEMENT_TOL if k in base_raw["placement"] else MEASURED_TOL)}
    assert off == KNOWN_DERIVED_MISSES


@pytest.mark.parametrize("name", TARGETS.values())
def test_rules_hold_on_acs(name):
    """The rules the report checks a new style against, on the values they were fitted to."""
    _, flat = _committed(name)
    values = {k: probe.Value(v, "probe", "") for k, v in flat.items()}
    bad = [c for c in probe.rule_checks(values) if abs(c[1] - c[2]) > probe.RULE_TOL]
    assert [c[0] for c in bad] == ["charge.bar_length"]  # the k-scaled bars: KNOWN_MEASURED_MISSES


@pytest.mark.parametrize("scale", ["1.0", "0.8", "0.6"])
@needs_refs
def test_sign_sizes_follow_their_box(scale):
    """ChemDraw's sign bars and rings depend on the symbol's box dx alone (acs-style.md §16); the
    plain signs and the circled ones have different boxes at x0.8 and x0.6."""
    syms = np.mr.charges_lewis(f"charges_lewis_x{scale}")["symbols"]
    assert len(syms) > 4
    for s in syms:
        if s.get("bar_length"):
            assert s["bar_length"] == pytest.approx(probe.bar_length_of_box(s["dx"]), abs=0.003), s
        if s.get("ring_radius"):
            assert s["ring_radius"] == pytest.approx(probe.ring_radius_of_box(s["dx"]), abs=0.03), s


def _check_measured(values, name):
    base_raw, _ = np.base_values()
    assert values.keys() == _committed(name)[1].keys() - {"name", "scale"}
    d = _diffs(values, _committed(name)[1])
    off = {k for k, v in d.items() if v > (PLACEMENT_TOL if k in base_raw["placement"] else MEASURED_TOL)}
    assert off == KNOWN_DERIVED_MISSES | KNOWN_MEASURED_MISSES


@needs_refs
@pytest.mark.parametrize("scale,name", TARGETS.items())
def test_methods_on_the_references(scale, name):
    """collect() on the references, role by role: each measured value within 0.1 pt of the committed preset."""
    values = np.collect(probe.references(scale), {"text.number_weight": "bold"})
    assert {k for k, v in values.items() if v.source == "settings"} == set(probe.SETTINGS_KEYS)
    assert {k for k, v in values.items() if v.source == "probe"} == set(probe.PROBE_KEYS)
    _check_measured(values, name)


@pytest.mark.parametrize("scale,name", TARGETS.items())
def test_real_probe_hold_out(scale, name):
    """A probe saved from ChemDraw with x<scale>'s settings: the same comparison, from one file."""
    path = PROBES / f"style_probe_x{scale}.cdxml"
    if not path.exists():
        pytest.skip(f"no {path.name}: one ChemDraw round on the probe")
    sources = probe.single_file(path)
    values = np.collect(sources, {})
    _check_measured(values, name)
    bad = [c for c in probe.rule_checks(values, probe.circled_ring(sources["circled"])) if abs(c[1] - c[2]) > probe.RULE_TOL]
    assert bad == []


# --------------------------------------------------------------------------- the probe


def test_probe_holds_every_specimen(tmp_path):
    """The probe file: what the build measures must be there before anyone opens it in ChemDraw."""
    path = np.write_probe(tmp_path)
    root = etree.parse(str(path), np.PARSER).getroot()
    arrows = list(root.iter("arrow"))
    forward = [a.get("HeadSize") for a in arrows if a.get("ArrowheadHead") == "Full" and not a.get("ArrowheadTail")]
    assert forward[:2] == ["1000", "1500"]  # standard, then large
    assert any(a.get("ArrowShaftSpacing") for a in arrows)  # equilibrium
    assert any(a.get("ArrowheadTail") == "Full" for a in arrows)  # ↔
    symbols = {g.get("SymbolType") for g in root.iter("graphic")}
    assert {"Plus", "Minus", "CirclePlus", "CircleMinus", "LonePair", "Electron"} <= symbols
    assert sum(g.get("GraphicType") == "Bracket" for g in root.iter("graphic")) == 2
    assert {c.get("ArrowheadHead") for c in root.iter("curve")} >= {"Full", "HalfLeft", "HalfRight"}
    assert {"WedgeBegin", "WedgedHashBegin"} <= {b.get("Display") for b in root.iter("b")}
    labels = {("".join(s.text for s in n.find("t").iter("s")), n.find("t").get("LabelAlignment")) for n in root.iter("n")
              if n.find("t") is not None}
    assert {("NH", "Above"), ("NH", "Below"), ("NTs", "Above")} <= labels
    fonts = {f.get("id"): f.get("name") for f in root.find("fonttable")}
    assert {fonts[s.get("font")] for s in root.iter("s")} >= {"Helvetica", "Helvetica Bold", "Helvetica Oblique"}


def _variant(tmp_path, attrs):
    """m1_set_x0.8 with its root attributes changed: a settings file whose drawing didn't follow."""
    root = etree.parse(str(probe.REF_DIR / "m1_set_x0.8.cdxml"), np.PARSER)
    for k, v in attrs.items():
        root.getroot().set(k, v)
    root.write(str(tmp_path / "m.cdxml"))
    shutil.copy(probe.REF_DIR / "m1_set_x0.8.svg", tmp_path / "m.svg")
    return {r: probe.Source("m", tmp_path) for r in probe.ROLES}


@needs_refs
def test_refuses_specimens_that_ignore_the_settings(tmp_path):
    src = _variant(tmp_path, {"BondLength": "20"})
    with pytest.raises(probe.MissingSpecimen, match="settings were not applied"):
        probe.check_restyled(src, probe.read_settings(src))


@needs_refs
def test_refuses_another_font(tmp_path):
    root = etree.parse(str(probe.REF_DIR / "m1_set_x0.8.cdxml"), np.PARSER)
    table = root.getroot().find("fonttable")
    etree.SubElement(table, "font", {"id": "999", "charset": "x-mac-roman", "name": "Arial"})
    root.getroot().set("LabelFont", "999")
    root.write(str(tmp_path / "m.cdxml"))
    with pytest.raises(probe.UnsupportedFont, match="Arial"):
        probe.read_settings({"settings": probe.Source("m", tmp_path)})


@needs_refs
def test_refuses_missing_specimens():
    only = probe.references("0.8")
    with pytest.raises(probe.MissingSpecimen, match="no circled|no Plus|circled sign"):
        probe.circled_ring(only["settings"])  # m1_set has no circled sign
    with pytest.raises(probe.MissingSpecimen, match="standard and"):
        probe._arrows(only["arrows"], only["arrows"])  # one head size only
    with pytest.raises(probe.MissingSpecimen, match="NH"):
        probe._nh(only["settings"])  # m1_set has no NH
    with pytest.raises(probe.MissingSpecimen):
        probe.single_file(REPO / "no_such_probe.cdxml")


@needs_refs
def test_header_from_any_file_without_the_references():
    """A preset's header template from a file by path, its caption fonts from the committed base
    template: the same bytes as the committed template of that scale."""
    got = np.header(probe.REF_DIR / "m1_set_x0.8.cdxml", font_source=np.HEADER_DIR / "acs1996.xml")
    assert got == (np.HEADER_DIR / "acs1996_x0.8.xml").read_bytes()


# --------------------------------------------------------------------------- a generated preset in use


SPECS = ["examples/m1/m1_smiles_x1.0.yaml", "examples/m1/m1_abbrev2_x1.0.yaml", "examples/m1/m1_charges_x1.0.yaml",
         "examples/m1/m1_chair_x1.0.yaml", "examples/m1/m1_lewis_x1.0.yaml", "examples/m2/m2_equilibrium_x1.0.yaml",
         "examples/m3/m3_radical_propagation_x1.0.yaml", "examples/m4/m4_resonance_x1.0.yaml", "examples/m4/m4_wrap_x1.0.yaml",
         "examples/m6/m6_suzuki_cycle_x1.0.yaml"]


@pytest.fixture
def generated(tmp_path, monkeypatch):
    """A preset written by new_preset from x0.8's values, in a temporary presets/: no ChemDraw file needed."""
    base_raw, base = np.base_values()
    known = _known("acs1996_x0.8")
    values = known | probe.derive(base, base_raw["placement"], known)
    pdir = tmp_path / "presets"
    (pdir / "cdxml").mkdir(parents=True)
    shutil.copy(np.PRESET_DIR / "acs1996.yaml", pdir)
    name = "holdout_x08"
    np.write_preset(name, values, None, "hold-out", pdir, pdir / "cdxml",
                    header_bytes=(np.HEADER_DIR / "acs1996_x0.8.xml").read_bytes())
    monkeypatch.setattr(preset_mod, "PRESET_DIR", pdir)
    monkeypatch.setattr(cdxml_writer, "HEADER_DIR", pdir / "cdxml")
    cdxml_writer.load_nicknames.cache_clear()
    yield name, pdir
    cdxml_writer.load_nicknames.cache_clear()
    preset_mod.load_preset.cache_clear()


@pytest.mark.parametrize("spec", ["examples/m1/m1_abbrev2_x1.0.yaml", "examples/m3/m3_radical_propagation_x1.0.yaml",
                                  "examples/m4/m4_resonance_x1.0.yaml", "examples/m1/m1_charges_x1.0.yaml"])
def test_preset_built_from_the_real_probe(tmp_path, monkeypatch, spec):
    """write_preset on a real re-saved probe, its header included: it draws and exports, and the CDXML
    states the probe's settings with the caption fonts present."""
    path = PROBES / "style_probe_x0.8.cdxml"
    if not path.exists():
        pytest.skip(f"no {path.name}: one ChemDraw round on the probe")
    pdir = tmp_path / "presets"
    (pdir / "cdxml").mkdir(parents=True)
    shutil.copy(np.PRESET_DIR / "acs1996.yaml", pdir)
    np.write_preset("from_probe", np.collect(probe.single_file(path), {}), path, path.name, pdir, pdir / "cdxml")
    monkeypatch.setattr(preset_mod, "PRESET_DIR", pdir)
    monkeypatch.setattr(cdxml_writer, "HEADER_DIR", pdir / "cdxml")
    for f in (preset_mod.load_preset, cdxml_writer.load_nicknames):
        f.cache_clear()
    raw = yaml.safe_load((REPO / spec).read_text())
    raw["preset"] = "from_probe"
    out = tmp_path / Path(spec).name
    out.write_text(yaml.safe_dump(raw, allow_unicode=True))
    render_spec(out)
    root = etree.fromstring(spec_cdxml(out, warn=lambda m: None).encode(), np.PARSER)
    assert root.get("BondLength") == "11.52"
    assert {f.get("name") for f in root.find("fonttable")} >= {"Helvetica", "Helvetica Bold", "Helvetica Oblique"}
    for f in (preset_mod.load_preset, cdxml_writer.load_nicknames):
        f.cache_clear()


def test_generated_preset_lists_its_derived_values(generated):
    name, pdir = generated
    p = Preset.model_validate(yaml.safe_load((pdir / f"{name}.yaml").read_text()))
    assert "scheme.arrow_gap" in p.derived and "lengths.bond_length" not in p.derived
    text = (pdir / f"{name}.yaml").read_text()
    assert "derived from acs1996" in text and "hold-out: committed value" in text


@pytest.mark.parametrize("spec", SPECS)
def test_generated_preset_renders_and_exports(generated, tmp_path, spec):
    """Every spec kind draws at a generated preset and exports CDXML, abbreviations by generated fragments."""
    name, _ = generated
    raw = yaml.safe_load((REPO / spec).read_text())
    raw["preset"] = name
    path = tmp_path / Path(spec).name
    path.write_text(yaml.safe_dump(raw, allow_unicode=True))
    load_spec(path)
    render_spec(path)
    root = etree.fromstring(spec_cdxml(path, warn=lambda m: None).encode(), np.PARSER)
    assert float(root.get("BondLength")) == pytest.approx(11.52)
