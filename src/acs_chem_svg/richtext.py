"""Rich text for conditions, compound numbers and abbreviation labels (acs-style.md §3).

`parse` turns a plain string into styled runs:
- digits right after a letter, `)` or a subscripted digit are subscripts (K2CO3, Pd(PPh3)4,
  CH2Cl2, C12H25);
- a `+` or `−` (U+2212) right after a formula token, and before a space, punctuation or the
  end, is a superscript (H3O+, Cl−). ASCII `-` is always a hyphen (n-BuLi, 1,4-dioxane);
- `*…*` is italic (*n*-BuLi);
- `\\` escapes the next character, which is then set plain (SN\\2, a literal \\*).
Everything else is untouched: °C, −78, (4:1), 1,4-, 12 h.

`layout_line` sets runs as consecutive TextRuns, as ChemDraw does: each style or script change
starts a new run at the previous run's origin + advance, with no kerning across the boundary
(scheme_multi_x*: K|2|CO|3, 1. |n|-BuLi). Pair kerning applies within a run (bold "11",
acs-style.md §5). Sub- and superscript sizes and offsets are preset ratios.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

from .fontmetrics import FontFace, load_face
from .geometry import union_box
from .preset import Preset
from .scene import TextRun

STYLES = ("regular", "italic", "bold")
SCRIPTS = ("base", "sub", "sup")
_SUP_SIGNS = "+−"
_SUP_FOLLOW = " ,;:)]/"


@dataclass(frozen=True)
class Run:
    text: str
    style: str = "regular"  # regular | italic | bold
    script: str = "base"  # base | sub | sup


def parse(s: str, style: str = "regular") -> list[Run]:
    """Split `s` into runs (see the module docstring). `style` is the base style ("bold" for
    compound numbers); `*…*` switches to italic."""
    chars: list[tuple[str, str, str]] = []  # (char, style, script)
    italic = False
    prev_kind = ""  # "letter", "paren", "sub", "" (anything else) for the previous character
    i = 0
    while i < len(s):
        ch = s[i]
        st = "italic" if italic else style
        if ch == "\\" and i + 1 < len(s):
            chars.append((s[i + 1], st, "base"))
            prev_kind = ""
            i += 2
            continue
        if ch == "*":
            italic = not italic
            i += 1
            continue
        nxt = s[i + 1] if i + 1 < len(s) else ""
        if ch.isdigit() and prev_kind in ("letter", "paren", "sub"):
            chars.append((ch, st, "sub"))
            prev_kind = "sub"
        elif ch in _SUP_SIGNS and prev_kind in ("letter", "paren", "sub") and (nxt == "" or nxt in _SUP_FOLLOW):
            chars.append((ch, st, "sup"))
            prev_kind = ""
        else:
            chars.append((ch, st, "base"))
            prev_kind = "letter" if ch.isalpha() and ch.isascii() else "paren" if ch in ")]" else ""
        i += 1
    if italic:
        raise ValueError(f"unclosed '*' (italic) in {s!r}")
    runs: list[Run] = []
    for ch, st, sc in chars:
        if runs and runs[-1].style == st and runs[-1].script == sc:
            runs[-1] = Run(runs[-1].text + ch, st, sc)
        else:
            runs.append(Run(ch, st, sc))
    return runs


def has_sub(runs: list[Run]) -> bool:
    return any(r.script == "sub" for r in runs)


@lru_cache(maxsize=None)
def face_for(style: str) -> FontFace:
    return load_face({"regular": "Regular", "italic": "Oblique", "bold": "Bold"}[style])


@dataclass
class Line:
    """One laid-out line: glyph origin x0 of the first run, advance, and its TextRuns."""
    runs: list[TextRun]
    x0: float
    baseline: float
    advance: float
    has_sub: bool
    text: str = ""
    parts: list[Run] = field(default_factory=list)

    def ink(self):
        return union_box(r.bbox() for r in self.runs)

    def shift(self, dx: float, dy: float = 0.0) -> None:
        for r in self.runs:
            r.x += dx
            r.y += dy
        self.x0 += dx
        self.baseline += dy


def run_metrics(run: Run, size: float, preset: Preset) -> tuple[float, float]:
    """(font size, baseline offset, y down) of a run set in text of `size`."""
    r = preset.ratios
    if run.script == "sub":
        return r.subscript_scale * size, r.subscript_drop_per_size * size
    if run.script == "sup":
        return r.superscript_scale * size, -r.superscript_raise_per_size * size
    return size, 0.0


def layout_line(text: str | list[Run], x: float, baseline: float, size: float, preset: Preset,
                source: str = "", style: str = "regular") -> Line:
    """Set one line with its first glyph origin at (x, baseline)."""
    parts = parse(text, style) if isinstance(text, str) else text
    runs, cx = [], x
    for p in parts:
        sz, dy = run_metrics(p, size, preset)
        face = face_for(p.style)
        runs.append(TextRun(cx, baseline + dy, sz, p.text, face, weight="bold" if p.style == "bold" else "normal",
                            style="italic" if p.style == "italic" else "normal", source=source))
        cx += face.advance(p.text, sz)
    plain = "".join(p.text for p in parts)
    return Line(runs, x, baseline, cx - x, has_sub(parts), plain, parts)


def line_pitch(upper: Line, lower: Line, size: float, preset: Preset) -> float:
    """Baseline to baseline (acs-style.md §3): size + line_pitch_extra, plus fitted increments
    when the upper or the lower line has a subscript (approx, Open 9)."""
    t = preset.text
    return (size + t.line_pitch_extra + (t.line_pitch_sub_upper if upper.has_sub else 0.0)
            + (t.line_pitch_sub_lower if lower.has_sub else 0.0))


def layout_block(lines: list[str | list[Run]], x: float, first_baseline: float, size: float, preset: Preset,
                 align: str = "left", source: str = "") -> list[Line]:
    """Stack lines downward from `first_baseline`. align: "left" (origins at x), "center"
    (advance centers at x)."""
    out: list[Line] = []
    base = first_baseline
    for text in lines:
        ln = layout_line(text, 0.0, 0.0, size, preset, source)
        if out:
            base += line_pitch(out[-1], ln, size, preset)
        ln.shift(x - (ln.advance / 2 if align == "center" else 0.0), base)
        out.append(ln)
    return out
