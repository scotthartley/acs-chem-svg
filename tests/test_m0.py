"""M0 smoke tests: font faces load, presets parse and follow the documented scale rules."""

from pathlib import Path

import pytest
import yaml

from acs_chem_svg import fontmetrics

PRESETS = Path(__file__).resolve().parents[1] / "presets"
# dimensionless placement defaults (control offset / chord): the same in every preset
RATIO_PLACEMENT = ("curved_arrow.curvature", "curved_arrow.loop_curvature")


def _flat(d, prefix=""):
    out = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flat(v, key + "."))
        else:
            out[key] = v
    return out


def _load(name):
    return yaml.safe_load((PRESETS / f"{name}.yaml").read_text())


@pytest.mark.parametrize("sub", ["Regular", "Bold"])
def test_fallback_faces_load(sub):
    face = fontmetrics.load_fallback(sub)
    assert face.subfamily == sub
    assert face.advance("H", 10) == pytest.approx(7.22, abs=0.01)


@pytest.mark.skipif(not fontmetrics.SYSTEM_HELVETICA.exists(), reason="no system Helvetica")
@pytest.mark.parametrize("sub", ["Regular", "Bold"])
def test_system_helvetica_faces_load(sub):
    face = fontmetrics.load_system_helvetica(sub)
    assert (face.family, face.subfamily) == ("Helvetica", sub)
    if sub == "Regular":  # labels (Regular) hang from cap height / 2 (acs-style.md §2)
        assert face.cap_height(10) / 2 == pytest.approx(3.586, abs=0.01)


@pytest.mark.skipif(not fontmetrics.SYSTEM_HELVETICA.exists(), reason="no system Helvetica")
def test_kerning_matches_resvg_behaviour():
    h = fontmetrics.load_system_helvetica("Regular")
    assert h.advance("11", 10) == pytest.approx(2 * h.advance("1", 10, kern=False) - 0.737, abs=0.002)


def test_presets_parse_and_scale():
    base = _load("acs1996")
    lists = ("exempt", "nonlinear", "assumed_exempt", "approx", "placement")
    flat_base = _flat({k: v for k, v in base.items() if k not in lists})
    for key in base["placement"]:
        grid = 10 if key in RATIO_PLACEMENT else 2  # round ratios on a 0.1 grid; pt on a 0.5 pt grid
        assert round(grid * flat_base[key], 9).is_integer(), key
    for name, k in (("acs1996_x0.8", 0.8), ("acs1996_x0.6", 0.6)):
        p = _load(name)
        assert p["scale"] == k
        assert all(p[x] == base[x] for x in lists)
        flat = _flat({kk: v for kk, v in p.items() if kk not in lists})
        assert flat.keys() == flat_base.keys()
        for key, v0 in flat_base.items():
            v = flat[key]
            if key in ("name", "scale", "font_family") or key in p["approx"] or not isinstance(v0, (int, float)):
                continue
            if (key.startswith(("ratios.", "angles.")) or key in p["exempt"] or key in p["assumed_exempt"]
                    or key in RATIO_PLACEMENT):
                assert v == v0, key
            elif key in p["placement"]:  # round defaults, each preset rounded from its own measurement
                assert (2 * v).is_integer(), key
                if key not in p["nonlinear"]:
                    assert abs(v - k * v0) <= 0.5, key
            elif key in p["nonlinear"]:
                assert v < v0, key
            else:
                assert v == pytest.approx(k * v0, abs=0.012), key
