"""Glyph metrics for label and caption layout (M0 scope).

Faces are loaded with fontTools. System Helvetica is preferred; the bundled
TeX Gyre Heros (fonts/, GUST license) is the fallback. Every FontFace records
the file it came from so renders can report which font was measured.

All public measurements are in pt for a given font size in pt. Runs are
kerned (GPOS pairs, else the legacy 'kern' table), matching resvg.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

from fontTools.pens.boundsPen import BoundsPen
from fontTools.ttLib import TTCollection, TTFont

# System Helvetica's tables make fontTools log harmless parse notes ("1 extra bytes in
# post.stringData array", "'created' timestamp seems very low") that use sessions mistook for
# the fallback-font warning (M5 r0 trials). Real problems still raise.
logging.getLogger("fontTools").setLevel(logging.ERROR)

REPO_ROOT = Path(__file__).resolve().parents[2]
SYSTEM_HELVETICA = Path("/System/Library/Fonts/Helvetica.ttc")
FALLBACK_FILES = {
    "Regular": REPO_ROOT / "fonts" / "texgyreheros-regular.otf",
    "Bold": REPO_ROOT / "fonts" / "texgyreheros-bold.otf",
    "Oblique": REPO_ROOT / "fonts" / "texgyreheros-italic.otf",  # italic runs in conditions (n-BuLi)
}


@dataclass(eq=False)  # hashed by identity for the lru_caches below
class FontFace:
    path: Path
    index: int | None  # face index within a .ttc, else None
    family: str
    subfamily: str
    font: TTFont = field(repr=False)

    @property
    def is_fallback(self) -> bool:
        return self.path != SYSTEM_HELVETICA

    @property
    def upm(self) -> int:
        return self.font["head"].unitsPerEm

    def _glyph_name(self, ch: str) -> str:
        cmap = self.font.getBestCmap()
        try:
            return cmap[ord(ch)]
        except KeyError:
            raise KeyError(f"{ch!r} (U+{ord(ch):04X}) not in {self.path.name}") from None

    def has_char(self, ch: str) -> bool:
        return ord(ch) in self.font.getBestCmap()

    def advance(self, text: str, size: float, kern: bool = True) -> float:
        """Advance width of a run in pt, with pair kerning if `kern`."""
        hmtx = self.font["hmtx"]
        names = [self._glyph_name(ch) for ch in text]
        units = sum(hmtx[g][0] for g in names)
        if kern:
            units += sum(self._gpos_pair(a, b) for a, b in zip(names, names[1:]))
        return units * size / self.upm

    def kerning(self, left: str, right: str, size: float) -> float:
        """GPOS pair adjustment in pt between two characters (0 if none)."""
        return self._gpos_pair(self._glyph_name(left), self._glyph_name(right)) * size / self.upm

    def _gpos_pair(self, g1: str, g2: str) -> int:
        """Pair adjustment in font units: GPOS PairPos if present, else the legacy
        'kern' table. resvg applies both (scripts/check_rasterizer.py)."""
        if "GPOS" not in self.font:
            if "kern" not in self.font:
                return 0
            return sum(t.kernTable.get((g1, g2), 0) for t in self.font["kern"].kernTables)
        for lookup in self.font["GPOS"].table.LookupList.Lookup:
            subtables = lookup.SubTable
            if lookup.LookupType == 9:
                subtables = [s.ExtSubTable for s in subtables]
            for st in subtables:
                if getattr(st, "LookupType", lookup.LookupType) != 2 or g1 not in st.Coverage.glyphs:
                    continue
                if st.Format == 1:
                    for pv in st.PairSet[st.Coverage.glyphs.index(g1)].PairValueRecord:
                        if pv.SecondGlyph == g2:
                            return getattr(pv.Value1, "XAdvance", 0) or 0
                else:
                    c1 = st.ClassDef1.classDefs.get(g1, 0)
                    c2 = st.ClassDef2.classDefs.get(g2, 0)
                    v = st.Class1Record[c1].Class2Record[c2].Value1
                    xa = (getattr(v, "XAdvance", 0) or 0) if v else 0
                    if xa:
                        return xa
        return 0

    def ink_bounds(self, ch: str, size: float) -> tuple[float, float, float, float] | None:
        """Ink box of one glyph as (xmin, ymin, xmax, ymax) in pt, y up, origin at
        the glyph origin on the baseline. None for blank glyphs."""
        return _ink_bounds_units(self, self._glyph_name(ch), size)

    def text_ink_bounds(self, text: str, size: float) -> tuple[float, float, float, float] | None:
        """Ink box of a kerned run set from x=0 on the baseline, y up, in pt."""
        x = 0.0
        box = None
        for i, ch in enumerate(text):
            b = self.ink_bounds(ch, size)
            if b is not None:
                b = (b[0] + x, b[1], b[2] + x, b[3])
                box = b if box is None else (
                    min(box[0], b[0]), min(box[1], b[1]), max(box[2], b[2]), max(box[3], b[3]))
            x += self.advance(ch, size)
            if i + 1 < len(text):
                x += self.kerning(ch, text[i + 1], size)
        return box

    def outline(self, ch: str, size: float, steps: int = 8) -> list[list[tuple[float, float]]]:
        """Glyph contours flattened to closed polylines, in pt, y up, origin at
        the glyph origin on the baseline."""
        return _outline_units(self, self._glyph_name(ch), size, steps)

    def cap_height(self, size: float) -> float:
        os2 = self.font["OS/2"]
        h = getattr(os2, "sCapHeight", 0) or 0
        if not h:
            h = self.ink_bounds("H", size)[3] * self.upm / size
        return h * size / self.upm

    def describe(self) -> str:
        idx = f"#{self.index}" if self.index is not None else ""
        return f"{self.family} {self.subfamily} ({self.path}{idx})"


@lru_cache(maxsize=None)
def _ink_bounds_units(face: FontFace, glyph: str, size: float):
    gs = face.font.getGlyphSet()
    pen = BoundsPen(gs)
    gs[glyph].draw(pen)
    if pen.bounds is None:
        return None
    s = size / face.upm
    return tuple(v * s for v in pen.bounds)


@lru_cache(maxsize=None)
def _outline_units(face: FontFace, glyph: str, size: float, steps: int):
    from fontTools.pens.recordingPen import DecomposingRecordingPen

    gs = face.font.getGlyphSet()
    pen = DecomposingRecordingPen(gs)
    gs[glyph].draw(pen)
    s = size / face.upm
    contours, cur, last = [], [], (0.0, 0.0)

    def bez(pts, t):  # de Casteljau
        while len(pts) > 1:
            pts = [((1 - t) * a[0] + t * b[0], (1 - t) * a[1] + t * b[1]) for a, b in zip(pts, pts[1:])]
        return pts[0]

    for op, args in pen.value:
        if op == "moveTo":
            cur, last = [args[0]], args[0]
        elif op == "lineTo":
            cur.append(args[0])
            last = args[0]
        elif op in ("qCurveTo", "curveTo"):
            pts = list(args)
            if op == "qCurveTo" and len(pts) > 2:  # implied on-curve points
                segs, p0 = [], last
                for i in range(len(pts) - 2):
                    mid = ((pts[i][0] + pts[i + 1][0]) / 2, (pts[i][1] + pts[i + 1][1]) / 2)
                    segs.append([p0, pts[i], mid])
                    p0 = mid
                segs.append([p0, pts[-2], pts[-1]])
            else:
                segs = [[last, *pts]]
            for seg in segs:
                cur.extend(bez(seg, k / steps) for k in range(1, steps + 1))
            last = pts[-1]
        elif op in ("closePath", "endPath"):
            if cur:
                contours.append(cur)
            cur = []
    return [[(x * s, y * s) for x, y in c] for c in contours]


def _names(font: TTFont) -> tuple[str, str]:
    name = font["name"]
    return str(name.getDebugName(1)), str(name.getDebugName(2))


def list_system_helvetica_faces() -> list[tuple[int, str, str]]:
    """(index, family, subfamily) for every face in the system Helvetica.ttc."""
    if not SYSTEM_HELVETICA.exists():
        return []
    coll = TTCollection(str(SYSTEM_HELVETICA))
    return [(i, *_names(f)) for i, f in enumerate(coll.fonts)]


@lru_cache(maxsize=None)
def load_system_helvetica(subfamily: str) -> FontFace:
    coll = TTCollection(str(SYSTEM_HELVETICA))
    for i, f in enumerate(coll.fonts):
        fam, sub = _names(f)
        if fam == "Helvetica" and sub == subfamily:
            return FontFace(SYSTEM_HELVETICA, i, fam, sub, f)
    raise LookupError(f"No Helvetica {subfamily} face in {SYSTEM_HELVETICA}")


@lru_cache(maxsize=None)
def load_fallback(subfamily: str) -> FontFace:
    path = FALLBACK_FILES[subfamily]
    f = TTFont(str(path))
    fam, sub = _names(f)
    return FontFace(path, None, fam, sub, f)


def load_face(subfamily: str = "Regular", allow_fallback: bool = True) -> FontFace:
    """Helvetica `subfamily` ("Regular", "Bold" or "Oblique"), else the bundled fallback."""
    try:
        return load_system_helvetica(subfamily)
    except (OSError, LookupError):
        if not allow_fallback:
            raise
    face = load_fallback(subfamily)
    warnings.warn(f"System Helvetica unavailable; measuring with {face.describe()}", stacklevel=2)
    return face
