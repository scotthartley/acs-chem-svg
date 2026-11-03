"""Abbreviation labels as rules (acs-style.md §3.5; abbrev_x*, abbrev_2_x*, ligands_x*, metals_x*).

The registry (data/abbreviations.yaml, plus a spec's `abbreviations:` block) holds three kinds of
entry: **groups** (Me, Boc, TBS: named substituents that also merge into a heteroatom's label),
**formulas** (OMe, CO2Me, PPh3: labels on their attachment atom) and **text** labels (L, Lₙ,
dppf: no chemistry). Everything else is derived:

- `tokenize`: the key split into tokens, longest match first against the group names, then
  element symbols (case-sensitive), then italic-prefixed names (s-Bu). Parentheses make a group;
  digits after a token or `)` are its count.
- `reverse`: the right-justified text (bond leaving the label to the right): the tokens in
  reverse order, also inside parentheses (OTf -> TfO, B(OH)2 -> (HO)2B, P(t-Bu)3 -> (t-Bu)3P).
  A single token stays as it is (Boc, t-Bu, N3, dppf). An entry's `reversed:` overrides it (Bpin -> pinB).
- `anchor`: the glyphs that stand for the attachment atom, at the centered end: the element
  symbol when the text starts (reversed: ends) with that element's token, else the first (last)
  character (Ph: P; TMS reversed: S; N3 right-justified: the 3).

ChemDraw reads a label it can't interpret in its own way; two of ours differ by user decision
(2026-11-03): IPr stays IPr (ChemDraw: PrI) and P(o-tol)3 reverses to (o-tol)3P (ChemDraw:
unchanged).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field, replace
from functools import lru_cache
from pathlib import Path

import yaml

DATA = Path(__file__).resolve().parent / "data" / "abbreviations.yaml"
KINDS = ("group", "formula", "text")
_SECTIONS = {"groups": "group", "formulas": "formula", "text": "text"}
_ELEMENTS = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
    "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm "
    "Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr"
).split()
ELEMENTS = frozenset(_ELEMENTS)
_PREFIXED = re.compile(r"([a-z]{1,4})-([A-Z][a-z]*)")  # s-Bu, sec-Bu: an italic prefix on a name


@dataclass(frozen=True)
class Entry:
    key: str  # plain ASCII: the spec's name for it and the CDXML label text
    kind: str  # group | formula | text
    smiles: str | None = None  # one `*` = the bond to the rest; None for a text label
    text: str = ""  # richtext markup (default: the key)
    reversed: str | None = None  # right-justified markup, where the rule is wrong
    nickname: bool = False  # ChemDraw writes it as NodeType="Nickname"

    @property
    def markup(self) -> str:
        return self.text or self.key


@dataclass(frozen=True)
class Token:
    """One unit of a label: an element (O), a named group (Me, t-Bu), a parenthesized group
    ((OH), (t-Bu)) or an unparsed literal (text labels). `count` is its trailing digits."""
    key: str  # plain text as in the label key, without the count
    markup: str  # richtext markup, without the count
    count: int = 1
    element: str | None = None  # set when the token is an element symbol
    inner: tuple["Token", ...] = ()  # a parenthesized group's tokens

    def render(self) -> str:
        return self.markup + (str(self.count) if self.count > 1 else "")


@dataclass(frozen=True)
class Registry:
    entries: dict[str, Entry] = field(default_factory=dict)

    def get(self, key: str) -> Entry | None:
        return self.entries.get(key)

    def __contains__(self, key: str) -> bool:
        return key in self.entries

    def keys(self, kind: str | None = None) -> list[str]:
        return [k for k, e in self.entries.items() if kind is None or e.kind == kind]

    def with_entries(self, extra: dict[str, dict] | None) -> "Registry":
        """This registry plus a spec's `abbreviations:` block ({key: {kind, smiles, text,
        reversed, nickname}}; kind defaults to formula, or text without a SMILES). Same keys replace."""
        if not extra:
            return self
        out = dict(self.entries)
        for key, d in extra.items():
            out[key] = entry_from(key, d)
        return Registry(out)


def entry_from(key: str, d: dict | None, kind: str | None = None) -> Entry:
    d = dict(d or {})
    kind = d.pop("kind", None) or kind or ("formula" if d.get("smiles") else "text")
    if kind not in KINDS:
        raise ValueError(f"abbreviation {key!r}: kind must be one of {', '.join(KINDS)}, not {kind!r}")
    unknown = set(d) - {"smiles", "text", "reversed", "nickname"}
    if unknown:
        raise ValueError(f"abbreviation {key!r}: unknown field(s) {', '.join(sorted(unknown))}")
    smiles = d.get("smiles")
    if kind == "text" and smiles:
        raise ValueError(f"abbreviation {key!r}: a text label has no SMILES")
    if kind != "text":
        if not smiles or smiles.count("*") != 1:
            raise ValueError(f"abbreviation {key!r}: the SMILES needs exactly one '*' (the attachment bond)")
        from rdkit import Chem
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            raise ValueError(f"abbreviation {key!r}: can't parse SMILES {smiles!r}")
        star = next(a for a in m.GetAtoms() if a.GetAtomicNum() == 0)
        if star.GetDegree() != 1:
            raise ValueError(f"abbreviation {key!r}: '*' must have exactly one bond")
    return Entry(key, kind, smiles, d.get("text") or "", d.get("reversed"), bool(d.get("nickname")))


@lru_cache(maxsize=1)
def default() -> Registry:
    raw = yaml.safe_load(DATA.read_text())
    entries = {}
    for section, kind in _SECTIONS.items():
        for key, d in (raw.get(section) or {}).items():
            if key in entries:
                raise ValueError(f"{DATA.name}: {key!r} is defined twice")
            entries[key] = entry_from(key, d, kind)
    return Registry(entries)


# --------------------------------------------------------------------------- grammar


def _markup_of(name: str, reg: Registry) -> str:
    e = reg.get(name)
    if e is not None:
        return e.markup
    m = _PREFIXED.fullmatch(name)
    return f"*{m.group(1)}*-{m.group(2)}" if m else name


def tokenize(label: str, reg: Registry | None = None) -> list[Token]:
    """`label` (a plain key) as tokens; ValueError naming the part that is neither a group nor an
    element. A text entry, or any key that is a group, is one token."""
    reg = reg or default()
    e = reg.get(label)
    if e is not None and e.kind in ("text", "group"):
        return [Token(label, e.markup)]
    names = sorted((k for k in reg.keys() if reg.get(k).kind == "group"), key=len, reverse=True)
    toks, _ = _tokens(label, 0, names, reg, label)
    return toks


def _tokens(s: str, i: int, names: list[str], reg: Registry, label: str, close: bool = False):
    out: list[Token] = []
    while i < len(s):
        ch = s[i]
        if ch == ")":
            if not close:
                raise ValueError(f"label {label!r}: unmatched ')'")
            return out, i + 1
        if ch == "(":
            inner, i = _tokens(s, i + 1, names, reg, label, close=True)
            key = "(" + "".join(t.key + (str(t.count) if t.count > 1 else "") for t in inner) + ")"
            tok = Token(key, "(" + "".join(t.render() for t in inner) + ")", inner=tuple(inner))
        else:
            tok, i = _one(s, i, names, reg, label)
        m = re.match(r"\d+", s[i:])
        if m:
            tok = replace(tok, count=int(m.group()))
            i += m.end()
        out.append(tok)
    if close:
        raise ValueError(f"label {label!r}: unclosed '('")
    return out, i


def _one(s: str, i: int, names: list[str], reg: Registry, label: str) -> tuple[Token, int]:
    for n in names:
        if s.startswith(n, i):
            return Token(n, _markup_of(n, reg)), i + len(n)
    if s[i:i + 2] in ELEMENTS:
        return Token(s[i:i + 2], s[i:i + 2], element=s[i:i + 2]), i + 2
    if s[i] in ELEMENTS:
        return Token(s[i], s[i], element=s[i]), i + 1
    m = _PREFIXED.match(s, i)
    if m:
        return Token(m.group(), _markup_of(m.group(), reg)), m.end()
    rest = re.match(r"[A-Za-z-]+", s[i:])
    part = rest.group() if rest else s[i]
    raise ValueError(f"label {label!r}: can't read {part!r} (neither a defined group nor an element). "
                     f"Define it as a group in the spec's `abbreviations:` block, or give the label a "
                     f"`reversed:` text")


def _reversed_tokens(toks: list[Token] | tuple[Token, ...]) -> list[Token]:
    out = []
    for t in reversed(toks):
        if t.inner:
            inner = _reversed_tokens(t.inner)
            t = replace(t, inner=tuple(inner), markup="(" + "".join(x.render() for x in inner) + ")",
                        key="(" + "".join(x.key + (str(x.count) if x.count > 1 else "") for x in inner) + ")")
        out.append(t)
    return out


def markup(label: str, side: str = "right", reg: Registry | None = None) -> str:
    """The label's richtext markup as drawn: side "right" (H side, i.e. the text runs right of the
    atom, left-justified) is the key's own text; side "left" (right-justified) is its reverse."""
    reg = reg or default()
    e = reg.get(label)
    if side != "left":
        return e.markup if e else "".join(t.render() for t in tokenize(label, reg))
    if e is not None and e.reversed:
        return e.reversed
    return "".join(t.render() for t in _reversed_tokens(tokenize(label, reg)))


def reverse(label: str, reg: Registry | None = None) -> str:
    """The right-justified text, plain (markup stripped): OTf -> TfO, B(OH)2 -> (HO)2B."""
    return plain(markup(label, "left", reg))


# groups ChemDraw doesn't know (abbrev_2_x*, ligands_x*: P(o-tol)3 and p-Tol are Unspecified,
# "can't interpret this label"; pin likewise, so Bpin is not a token sequence)
CHEMDRAW_UNKNOWN = frozenset({"o-tol", "p-Tol"})


def chemdraw_tokens(label: str) -> list[Token] | None:
    """`label` as ChemDraw would read it, approximated by the built-in registry's groups less
    CHEMDRAW_UNKNOWN (a spec's own groups are unknown to it too); None when it can't read it.
    ChemDraw reverses a right-justified label it reads, and keeps one it can't as stored."""
    reg = default()
    names = sorted((k for k in reg.keys() if reg.get(k).kind == "group" and k not in CHEMDRAW_UNKNOWN),
                   key=len, reverse=True)
    try:
        toks, _ = _tokens(label, 0, names, reg, label)
    except ValueError:
        return None
    if any(t.key in CHEMDRAW_UNKNOWN or any(u.key in CHEMDRAW_UNKNOWN for u in t.inner) for t in toks):
        return None  # read through _PREFIXED (o-tol): ChemDraw doesn't know it
    return toks


def chemdraw_stored(text: str) -> str:
    """The right-justified text to store so that ChemDraw's token reversal shows `text` (phase 4c
    label variants, user 2026-11-03: IPr stored as PrI shows IPr). `text` must be ChemDraw-readable."""
    toks = chemdraw_tokens(text)
    if toks is None:
        raise ValueError(f"ChemDraw can't read {text!r}")
    return "".join(t.render() for t in _reversed_tokens(toks))


def reads_as_formula(label: str) -> bool:
    """True when a text label also reads as two or more element / group tokens (IPr: I + Pr).
    ChemDraw interprets such a label on open and reverses it when right-justified (ligands_x*:
    PrI), against the user decision that text labels never reverse (2026-11-03)."""
    toks = chemdraw_tokens(label)
    return toks is not None and len(toks) > 1


def plain(markup_text: str) -> str:
    from . import richtext
    return "".join(r.text for r in richtext.parse(markup_text))


def anchor(label: str, element: str, side: str = "right", reg: Registry | None = None) -> str:
    """The glyphs that stand for the attachment atom (§3.5): `element` when the drawn text starts
    with its token (side "left", right-justified: ends with it, count 1), else the first (last)
    character of the drawn text."""
    reg = reg or default()
    text = plain(markup(label, side, reg))
    try:
        toks = tokenize(label, reg)
    except ValueError:
        toks = []
    if side == "left":
        toks = _reversed_tokens(toks)
        end = toks[-1] if toks else None
        ok = end is not None and end.element == element and end.count == 1 and text.endswith(element)
        return element if ok else text[-1:]
    first = toks[0] if toks else None
    ok = first is not None and first.element == element and text.startswith(element)
    return element if ok else text[:1]


def attachment_element(entry: Entry) -> str | None:
    """Element of the atom the entry attaches through (OMe: O, CO2Me: C); None for text labels."""
    if not entry.smiles:
        return None
    from rdkit import Chem
    m = Chem.MolFromSmiles(entry.smiles)
    star = next(a for a in m.GetAtoms() if a.GetAtomicNum() == 0)
    return star.GetNeighbors()[0].GetSymbol()


def compose(element: str, num_h: int, groups: list[str], reg: Registry | None = None) -> str:
    """The key of a heteroatom merged with its collapsed groups (§3.5 composition): the element,
    its H count, then the groups in order, identical ones counted (NMe2, SnBu3); a counted group
    whose text has a hyphen in parentheses (P(t-Bu)3). NHBoc, NTs, OTBS."""
    out = element + ("H" + (str(num_h) if num_h > 1 else "") if num_h else "")
    seen: list[tuple[str, int]] = []
    for g in groups:
        if seen and seen[-1][0] == g:
            seen[-1] = (g, seen[-1][1] + 1)
        elif any(s == g for s, _ in seen):
            seen = [(s, n + 1) if s == g else (s, n) for s, n in seen]
        else:
            seen.append((g, 1))
    for g, n in seen:
        out += (f"({g})" if "-" in g and n > 1 else g) + (str(n) if n > 1 else "")
    return out


# --------------------------------------------------------------------------- per-atom label info


@dataclass(frozen=True)
class LabelInfo:
    """What drawing and CDXML export need of one abbreviation label, fixed when the molecule is
    condensed (so a spec's own entries travel with the atom): key, kind (group | formula | text |
    composed), the markup of both forms, ChemDraw's NodeType class, and the collapsed atoms as
    SMILES with one `*` per outer bond ([*:1], [*:2] in bond order for composed labels)."""
    key: str
    kind: str
    text: str  # left-justified markup (bond leaving to the left)
    reversed: str  # right-justified markup
    element: str  # attachment element ("*" for a text label)
    nickname: bool = False
    smiles: str | None = None
    anchor_text: str = ""  # anchor glyphs of the left-justified form (anchor())
    anchor_reversed: str = ""  # ... and of the right-justified form

    def markup(self, side: str) -> str:
        return self.reversed if side == "left" else self.text

    def anchor(self, side: str) -> str:
        return self.anchor_reversed if side == "left" else self.anchor_text

    def to_prop(self) -> str:
        import json
        return json.dumps(self.__dict__)

    @staticmethod
    def from_prop(s: str) -> "LabelInfo":
        import json
        return LabelInfo(**json.loads(s))


def info(key: str, reg: Registry | None = None, element: str | None = None) -> LabelInfo:
    """LabelInfo of a registry key, or of a key the grammar reads (a composed label)."""
    reg = reg or default()
    e = reg.get(key)
    el = (attachment_element(e) if e else None) or element or "*"
    return LabelInfo(key, e.kind if e else "composed", markup(key, "right", reg), markup(key, "left", reg), el,
                     e.nickname if e else False, e.smiles if e else None,
                     anchor(key, el, "right", reg), anchor(key, el, "left", reg))

