"""scripts/render.py writes CDXML by default, and a spec it can't export still gets SVG + PNG."""

from pathlib import Path

import render
from render import main

from acs_chem_svg.cdxml_writer import CdxmlUnsupported

ROOT = Path(__file__).resolve().parents[1]
SCHEME = ROOT / "examples" / "m3" / "m3_carbonyl_addition_x1.0.yaml"
SPECIES = ROOT / "examples" / "m1" / "m1_charges_x1.0.yaml"


def test_cdxml_is_written_by_default(tmp_path):
    for spec in (SCHEME, SPECIES):
        paths = render.render(spec, tmp_path)
        assert [p.suffix for p in paths] == [".svg", ".png", ".cdxml"]
        assert paths[2].read_text().lstrip().startswith("<?xml")


def test_no_cdxml_writes_none_and_the_same_svg(tmp_path):
    for spec in (SCHEME, SPECIES):
        on = render.render(spec, tmp_path / "on")
        off = render.render(spec, tmp_path / "off", with_cdxml=False)
        assert [p.suffix for p in off] == [".svg", ".png"]
        assert not list((tmp_path / "off").glob("*.cdxml"))
        assert on[0].read_text() == off[0].read_text()


def test_flags(tmp_path, capsys):
    assert main([str(SPECIES), "--out", str(tmp_path / "a"), "--no-cdxml"]) == 0
    assert not list((tmp_path / "a").glob("*.cdxml"))
    assert main([str(SPECIES), "--out", str(tmp_path / "b"), "--cdxml"]) == 0  # old commands still work
    assert list((tmp_path / "b").glob("*.cdxml"))


def test_unsupported_cdxml_warns_and_still_renders(tmp_path, monkeypatch, capsys):
    def refuse(*args, **kwargs):
        raise CdxmlUnsupported("test: not exportable")
    monkeypatch.setattr(render, "to_cdxml", refuse)
    assert main([str(SCHEME), "--out", str(tmp_path)]) == 0
    assert sorted(p.suffix for p in tmp_path.iterdir()) == [".png", ".svg"]
    assert capsys.readouterr().err.count("warning: no CDXML: test: not exportable") == 1
