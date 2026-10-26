"""Rich-text parsing and line layout (richtext.py, acs-style.md §3)."""

import pytest

from acs_chem_svg import richtext
from acs_chem_svg.preset import load_preset

P = load_preset("acs1996")


def runs(s):
    return [(r.text, r.style, r.script) for r in richtext.parse(s)]


@pytest.mark.parametrize("text, expected", [
    ("K2CO3", [("K", "regular", "base"), ("2", "regular", "sub"), ("CO", "regular", "base"), ("3", "regular", "sub")]),
    ("Pd(PPh3)4", [("Pd(PPh", "regular", "base"), ("3", "regular", "sub"), (")", "regular", "base"),
                   ("4", "regular", "sub")]),
    ("CH2Cl2", [("CH", "regular", "base"), ("2", "regular", "sub"), ("Cl", "regular", "base"), ("2", "regular", "sub")]),
    ("C12H25", [("C", "regular", "base"), ("12", "regular", "sub"), ("H", "regular", "base"), ("25", "regular", "sub")]),
    ("H3O+", [("H", "regular", "base"), ("3", "regular", "sub"), ("O", "regular", "base"), ("+", "regular", "sup")]),
    ("Cl−", [("Cl", "regular", "base"), ("−", "regular", "sup")]),
    ("*n*-BuLi", [("n", "italic", "base"), ("-BuLi", "regular", "base")]),
])
def test_auto_sub_sup_italic(text, expected):
    assert runs(text) == expected


@pytest.mark.parametrize("text", ["0 °C", "−78 °C", "(4:1)", "1,4-dioxane", "12 h", "90 °C, 12 h", "acetone, reflux",
                                  "THF, 0 °C", "n-BuLi", "65%"])
def test_untouched(text):
    assert runs(text) == [(text, "regular", "base")]


def test_escapes():
    assert runs("SN\\2") == [("SN2", "regular", "base")]
    assert runs("\\*") == [("*", "regular", "base")]
    assert runs("NH\\4+") == [("NH4+", "regular", "base")]


def test_unclosed_italic():
    with pytest.raises(ValueError, match="unclosed"):
        richtext.parse("*n-BuLi")


def test_bold_base_style():
    assert runs("11") == [("11", "regular", "base")]
    assert [(r.text, r.style) for r in richtext.parse("4a", style="bold")] == [("4a", "bold")]


def test_layout_offsets():
    """Sub- and superscripts: 0.75 x size, dropped 0.23 / raised 0.3586 x size; runs are
    consecutive (each starts at the previous origin + advance)."""
    ln = richtext.layout_line("H3O+", 10.0, 50.0, 10.0, P)
    h, three, o, plus = ln.runs
    assert (three.size, three.y) == (pytest.approx(7.5), pytest.approx(52.3))
    assert (plus.size, plus.y) == (pytest.approx(7.5), pytest.approx(50 - 3.586))
    assert three.x == pytest.approx(h.x + h.advance())
    assert plus.x == pytest.approx(o.x + o.advance())
    assert ln.advance == pytest.approx(sum(r.advance() for r in ln.runs))


def test_italic_run_face():
    ln = richtext.layout_line("*n*-BuLi", 0.0, 0.0, 10.0, P)
    assert ln.runs[0].style == "italic" and ln.runs[0].face.subfamily in ("Oblique", "Italic")
    assert ln.runs[1].x == pytest.approx(5.562, abs=0.002)  # scheme_multi_x1.0: n advance, no kern across


def test_line_pitch_rule():
    plain = richtext.layout_line("90 °C, 12 h", 0, 0, 10.0, P)
    sub = richtext.layout_line("1,4-dioxane/H2O (4:1)", 0, 0, 10.0, P)
    sub2 = richtext.layout_line("Pd(PPh3)4", 0, 0, 10.0, P)
    assert richtext.line_pitch(plain, plain, 10.0, P) == pytest.approx(10.05)  # scheme_basic_x1.0
    assert richtext.line_pitch(sub, plain, 10.0, P) == pytest.approx(10.825)  # scheme_multi_x1.0
    assert richtext.line_pitch(sub, sub2, 10.0, P) == pytest.approx(10.975)
