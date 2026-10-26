"""Scheme spec: YAML validated with pydantic (references/scheme-schema.md).

    preset: acs1996                      # or acs1996_x0.8, acs1996_x0.6
    species:
      ester:  {smiles: "COC(=O)c1ccc(O)cc1", abbreviate: [OMe]}
      bnbr:   {smiles: "BrCc1ccccc1"}
      ether:  {smiles: "COC(=O)c1ccc(OCc2ccccc2)cc1", abbreviate: [OMe, OBn]}
      alanine:
        molfile: templates/m1/m1_set_x1.0__3.mol
        orientation: {rotate: 30, flip_h: false, flip_v: false}
    steps:                               # M2: one linear row
      - {from: [ester, bnbr], to: ether, above: K2CO3, below: "acetone, reflux"}
      - {to: alcohol, stages: ["*n*-BuLi, THF, −78 °C", "H3O+"]}   # from: the previous `to`
    numbers: {ester: 1, bnbr: 2, ether: "3a"}
    align: auto                          # or false; per species: align: false | align_to: <name>

A spec without `steps` is an M1 spec: its species are drawn side by side, unchanged from M1.

M3 mechanisms: a step lists curved arrows, drawn on its reactants and checked by electrons.py.

    steps:
      - from: [hydroxide, mebr]
        to: [methanol, bromide]
        arrows:
          - {type: two_electron, from: "lp(1)", to: "between(1,2)"}
          - {from: "bond(2-3)", to: "atom(3)"}         # type defaults to two_electron
    species:
      ether: {smiles: "...", show_lone_pairs: [4]}     # also draw these atoms' lone pairs
    overrides:
      arrows: [{step: 1, index: 0, bow: left, curvature: 0.6, tail_side: top}]

Endpoints use atom-map numbers: from lp(n), bond(a-b) / bond(a=b) / bond(a#b) (the written
order is checked), electron(n) (a radical); to atom(n) or between(a,b) (a new or stronger bond).

M4: resonance steps, several sequences, and wrapping into rows.

    layout: {width: 504}                 # pt; default: the preset's scheme.wrap_width
    steps:
      - {from: [benzene, nitronium], to: sigma_a, arrows: [...]}
      - resonance: [sigma_a, sigma_b, sigma_c]     # the first is the previous step's product
        arrows: [[...], [...]]                     # one list per contributor except the last
      - {from: [sigma_c, bisulfate], to: [nitrobenzene, acid]}
      - {new_row: true, from: [...], to: ...}     # a new sequence, on a new row

A step's `from` starts with one or more of the previous step's products, in order; the rest
of those products are by-products left behind, and further reactants are added ("+ X").
A resonance step may follow a step with several products: the first product opens the bracket
and the rest are drawn after it ("+ HX").

Generic atoms (X, R, Nu): write a real element in the SMILES as a stand-in and relabel it by
atom-map number; the electron bookkeeping keeps the stand-in (Cl for a radical X·).

    species:
      xrad: {smiles: "[Cl:9]"}
      hx:   {smiles: "[H:7][Cl:9]"}
    labels: {9: X}                       # drawn X, HX
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Annotated, Literal, Union

import yaml
from pydantic import BaseModel, ConfigDict, Discriminator, Field, Tag, field_validator, model_validator
from rdkit import Chem

from . import chem
from .preset import available, load_preset

REPO_ROOT = Path(__file__).resolve().parents[2]
MAX_STAGES = 2  # multistep arrows with more stages are unmeasured (acs-style.md §6 Open)
LABEL_TEXT = re.compile(r"^[A-Z][a-z]?'*$")  # a generic atom label: X, R, Nu, R'


@dataclass(frozen=True)
class Endpoint:
    """A parsed arrow endpoint: kind lp | bond | electron | atom | between; atoms are map
    numbers; order is the bond order written in bond(a-b) (1), bond(a=b) (2), bond(a#b) (3)."""
    kind: str
    atoms: tuple[int, ...]
    order: int | None = None
    text: str = ""


_ENDPOINT = re.compile(r"^\s*(lp|atom|electron)\(\s*(\d+)\s*\)\s*$|^\s*bond\(\s*(\d+)\s*([-=#])\s*(\d+)\s*\)\s*$"
                       r"|^\s*between\(\s*(\d+)\s*,\s*(\d+)\s*\)\s*$")
SOURCES = ("lp", "bond", "electron")
TARGETS = ("atom", "between")


def parse_endpoint(text: str) -> Endpoint:
    m = _ENDPOINT.match(text)
    if not m:
        raise ValueError(f"arrow endpoint {text!r}: expected lp(n), bond(a-b), bond(a=b), bond(a#b), electron(n), "
                         "atom(n) or between(a,b) with atom-map numbers")
    if m.group(1):
        return Endpoint(m.group(1), (int(m.group(2)),), text=text.strip())
    if m.group(3):
        a, b = int(m.group(3)), int(m.group(5))
        if a == b:
            raise ValueError(f"arrow endpoint {text!r}: a bond needs two different atoms")
        return Endpoint("bond", (a, b), {"-": 1, "=": 2, "#": 3}[m.group(4)], text.strip())
    a, b = int(m.group(6)), int(m.group(7))
    if a == b:
        raise ValueError(f"arrow endpoint {text!r}: between needs two different atoms")
    return Endpoint("between", (a, b), text=text.strip())


class CurvedArrow(BaseModel):
    model_config = ConfigDict(extra="forbid", populate_by_name=True)
    type: Literal["two_electron", "fishhook"] = "two_electron"
    from_: str = Field(alias="from")
    to: str

    @field_validator("from_")
    @classmethod
    def _source(cls, v):
        e = parse_endpoint(v)
        if e.kind not in SOURCES:
            raise ValueError(f"arrow from {v!r}: an arrow starts at lp(n), bond(a-b) or electron(n)")
        return v

    @field_validator("to")
    @classmethod
    def _target(cls, v):
        e = parse_endpoint(v)
        if e.kind not in TARGETS:
            raise ValueError(f"arrow to {v!r}: an arrow ends at atom(n) or between(a,b)")
        return v

    @property
    def source(self) -> Endpoint:
        return parse_endpoint(self.from_)

    @property
    def target(self) -> Endpoint:
        return parse_endpoint(self.to)

    @property
    def electrons(self) -> int:
        return 2 if self.type == "two_electron" else 1


class ArrowOverride(BaseModel):
    """Semantic adjustments to one curved arrow (design.md §5): no raw coordinates."""
    model_config = ConfigDict(extra="forbid")
    step: int  # 1-based
    index: int  # 0-based, in the step's `arrows` list
    bow: Literal["left", "right"] | None = None  # side of the tail -> head direction the curve bulges to
    curvature: float | None = Field(None, gt=0, le=3)  # control offset / chord
    tail_side: Literal["top", "right", "bottom", "left"] | None = None  # which lone pair (lp tails)


class Overrides(BaseModel):
    model_config = ConfigDict(extra="forbid")
    arrows: list[ArrowOverride] = Field(default_factory=list)


class Orientation(BaseModel):
    model_config = ConfigDict(extra="forbid")
    rotate: float | Literal["auto"] | None = None  # degrees, counter-clockwise on screen
    flip_h: bool = False
    flip_v: bool = False


class Species(BaseModel):
    model_config = ConfigDict(extra="forbid")
    smiles: str | None = None
    molfile: str | None = None
    orientation: Orientation | None = None
    kekule: list[str] | None = None  # "m=n" atom-map pairs that are double bonds
    # explicit labels on terminal carbons (CH3): atom ids are atom-map numbers if set, else 1-based atom indices.
    # Molfile atom aliases (CH3, H3C, ...) are honored as well.
    show_carbons: Literal["none", "terminal"] | list[int | str] = "none"
    abbreviate: list[str] | None = None  # collapse these groups to labels (chem.ABBREVIATIONS): OMe, OBn, ...
    align: bool = True  # take the orientation of the species it is drawn from (steps); false = own layout
    align_to: str | None = None  # align to this species instead
    show_lone_pairs: list[int] = Field(default_factory=list)  # atom-map numbers whose lone pairs are drawn (M3)

    @model_validator(mode="after")
    def _one_source(self):
        if (self.smiles is None) == (self.molfile is None):
            raise ValueError("give exactly one of smiles or molfile")
        return self

    @field_validator("kekule")
    @classmethod
    def _pairs(cls, v):
        for p in v or []:
            parts = p.split("=")
            if len(parts) != 2 or not all(x.strip().isdigit() for x in parts):
                raise ValueError(f"kekule entry {p!r} must look like '1=2' (atom-map numbers)")
        return v

    @field_validator("abbreviate")
    @classmethod
    def _known_abbrevs(cls, v):
        unknown = [a for a in v or [] if a not in chem.ABBREVIATIONS]
        if unknown:
            raise ValueError(f"unknown abbreviation {', '.join(unknown)}; known: {', '.join(chem.ABBREVIATIONS)}")
        return v


Names = str | list[str]


def _names(v: Names | None) -> list[str]:
    return [] if v is None else [v] if isinstance(v, str) else list(v)


class Step(BaseModel):
    model_config = ConfigDict(extra="forbid", populate_by_name=True)
    from_: Names | None = Field(None, alias="from")  # default: the previous step's `to`
    new_row: bool = False  # start a new sequence on a new row (`from` required)
    to: Names
    above: Names | None = None  # conditions above the arrow, one string per line
    below: Names | None = None
    stages: list[str] | None = None  # "1. ..." above and "2. ..." below, left-justified
    arrows: list[CurvedArrow] = Field(default_factory=list)  # M3 curved arrows, drawn on the reactants

    @model_validator(mode="before")
    @classmethod
    def _species_keys(cls, data):
        if isinstance(data, dict):
            misplaced = [k for k in ("align", "align_to", "orientation") if k in data]
            if misplaced:
                raise ValueError(f"{'/'.join(misplaced)} {'belongs' if len(misplaced) == 1 else 'belong'} on the species (species.<name>.{misplaced[0]}), "
                                 "not on a step")
        return data

    @model_validator(mode="after")
    def _stages(self):
        if self.stages is not None:
            if self.above is not None or self.below is not None:
                raise ValueError("give either stages or above/below, not both")
            if not 1 <= len(self.stages) <= MAX_STAGES:
                raise ValueError(f"stages: {len(self.stages)} given; 1 to {MAX_STAGES} are supported "
                                 "(more stages are unmeasured, acs-style.md §6 Open)")
        return self

    @property
    def sources(self) -> list[str]:
        return _names(self.from_)

    @property
    def products(self) -> list[str]:
        return _names(self.to)

    @property
    def arrow_lists(self) -> list[list[CurvedArrow]]:
        return [self.arrows]

    @property
    def all_arrows(self) -> list[CurvedArrow]:
        return self.arrows

    def lines(self) -> tuple[list[str], list[str]]:
        """(above, below) text lines; stages are numbered."""
        if self.stages:
            numbered = [f"{i}. {t}" for i, t in enumerate(self.stages, 1)]
            return numbered[:1], numbered[1:]
        return _names(self.above), _names(self.below)


class ResonanceStep(BaseModel):
    """Resonance contributors in square brackets, joined by double-headed arrows. Curved arrows:
    list k is drawn on contributor k and gives contributor k + 1."""
    model_config = ConfigDict(extra="forbid")
    resonance: list[str] = Field(min_length=2)
    arrows: list[list[CurvedArrow]] = Field(default_factory=list)

    @model_validator(mode="after")
    def _lists(self):
        if self.arrows and len(self.arrows) != len(self.resonance) - 1:
            raise ValueError(f"resonance arrows: {len(self.arrows)} list(s) given for {len(self.resonance)} contributors; "
                             "give one list per contributor except the last (list k turns contributor k into k + 1)")
        return self

    @property
    def sources(self) -> list[str]:
        return self.resonance[:1]

    @property
    def products(self) -> list[str]:
        return self.resonance[-1:]

    @property
    def arrow_lists(self) -> list[list[CurvedArrow]]:
        return self.arrows

    @property
    def all_arrows(self) -> list[CurvedArrow]:
        return [a for lst in self.arrows for a in lst]


AnyStep = Annotated[Union[Annotated[Step, Tag("reaction")], Annotated[ResonanceStep, Tag("resonance")]],
                    Discriminator(lambda v: "resonance" if (isinstance(v, dict) and "resonance" in v)
                                  or isinstance(v, ResonanceStep) else "reaction")]


class Layout(BaseModel):
    model_config = ConfigDict(extra="forbid")
    width: float | None = Field(None, gt=0)  # pt; rows wrap to stay within it (default: preset scheme.wrap_width)


@dataclass(frozen=True)
class ArrowSet:
    """One set of curved arrows applied together: a reaction step's arrows on its reactants, or
    one resonance list on contributor k (giving contributor k + 1). `offset` is the index of
    arrows[0] in the step's flattened list (overrides.arrows index, group names)."""
    step: int
    sources: tuple[str, ...]
    products: tuple[str, ...]
    arrows: tuple[CurvedArrow, ...]
    offset: int = 0
    contributor: int | None = None  # resonance: k (0-based)


def kept_products(src: list[str], prev: list[str]) -> int:
    """How many of the previous step's products `src` starts with (the rest are by-products)."""
    k = 0
    while k < min(len(src), len(prev)) and src[k] == prev[k]:
        k += 1
    return k


class Spec(BaseModel):
    model_config = ConfigDict(extra="forbid")
    preset: str = "acs1996"
    species: dict[str, Species] = Field(min_length=1)
    steps: list[AnyStep] | None = None
    layout: Layout = Field(default_factory=Layout)
    numbers: dict[str, int | str] | None = None
    align: Literal["auto"] | bool = "auto"
    labels: dict[int, str] = Field(default_factory=dict)  # atom-map number -> drawn symbol (X); chemistry keeps the SMILES element
    overrides: Overrides = Field(default_factory=Overrides)

    @field_validator("preset")
    @classmethod
    def _known(cls, v):
        if v not in available():
            raise ValueError(f"unknown preset {v!r}; available: {', '.join(available())}")
        return v

    @field_validator("labels")
    @classmethod
    def _label_text(cls, v):
        for m, t in v.items():
            if not LABEL_TEXT.match(t):
                raise ValueError(f"labels: {m}: {t!r} must be a plain symbol: a capital letter, an optional "
                                 "lower-case one and primes (X, R, Nu, R')")
        return v

    @model_validator(mode="after")
    def _scheme(self):
        names = set(self.species)
        for key, n in [("numbers", k) for k in self.numbers or {}] + [
                ("align_to", sp.align_to) for sp in self.species.values() if sp.align_to]:
            if n not in names:
                raise ValueError(f"{key}: unknown species {n!r}")
        if not self.steps:
            if self.numbers:
                raise ValueError("numbers need steps (a scheme)")
            if self.overrides.arrows:
                raise ValueError("overrides.arrows need steps with arrows")
            return self
        prev: list[str] = []  # the previous step's products, as written
        prev_ids: list[str] = []  # ... and as occurrence ids
        seen: dict[str, int] = {}

        def occurrence(n: str, aligned: bool) -> str:
            """A species' id at this place in the scheme: its name the first time, then a clone
            `name~k` (declared in `species:` to set its own align/orientation, else copied)."""
            seen[n] = seen.get(n, 0) + 1
            if seen[n] == 1:
                return n
            cid = f"{n}~{seen[n]}"
            if cid not in self.species:
                clone = self.species[n].model_copy(deep=True)
                if aligned and clone.align:
                    clone.orientation = None  # it takes the orientation of its own predecessor
                self.species[cid] = clone
            return cid

        for i, st in enumerate(self.steps, 1):
            if isinstance(st, ResonanceStep):
                for n in st.resonance:
                    if n not in names:
                        raise ValueError(f"step {i}: unknown species {n!r}")
                if prev and prev[0] != st.resonance[0]:
                    raise ValueError(f"step {i}: the first contributor {st.resonance[0]!r} must be the previous step's "
                                     f"first product (it has {prev}); the bracket opens with it and the other "
                                     "products follow the bracket")
                written = st.products
                st.resonance = (prev_ids[:1] if prev else [occurrence(st.resonance[0], False)]) + \
                    [occurrence(n, True) for n in st.resonance[1:]]
                prev, prev_ids = written, st.products
                continue
            if st.new_row and not st.sources:
                raise ValueError(f"step {i}: a step with new_row needs `from` (it starts a new sequence)")
            src = st.sources or prev
            if not src:
                raise ValueError(f"step {i}: the first step needs `from`")
            for n in src + st.products:
                if n not in names:
                    raise ValueError(f"step {i}: unknown species {n!r}")
            kept = 0 if st.new_row or not prev else kept_products(src, prev)
            if prev and not st.new_row and not kept:
                raise ValueError(f"step {i}: `from` {src} must start with the previous step's first product "
                                 f"{prev[0]!r} (then any more of {prev} in order; the rest are by-products): the steps "
                                 "form one linear chain; start another with new_row: true")
            if not st.products:
                raise ValueError(f"step {i}: `to` is empty")
            written = st.products
            st.from_ = prev_ids[:kept] + [occurrence(n, False) for n in src[kept:]]
            st.to = [occurrence(n, True) for n in written]
            prev, prev_ids = written, st.products
        order = self.order()
        unused = sorted(names - set(order))
        if unused:
            raise ValueError(f"species not used in any step: {', '.join(unused)}")
        for o in self.overrides.arrows:
            if not 1 <= o.step <= len(self.steps) or not 0 <= o.index < len(self.steps[o.step - 1].all_arrows):
                raise ValueError(f"overrides.arrows: step {o.step} has no arrow with index {o.index} "
                                 "(steps are 1-based, arrow indices 0-based; a resonance step counts its lists in order)")
        return self

    def arrow_override(self, step: int, index: int) -> ArrowOverride:
        for o in self.overrides.arrows:
            if (o.step, o.index) == (step, index):
                return o
        return ArrowOverride(step=step, index=index)

    def order(self) -> list[str]:
        """Species in drawing order, left to right and row by row."""
        if not self.steps:
            return list(self.species)
        out: list[str] = []
        prev: list[str] = []
        for st in self.steps:
            if isinstance(st, ResonanceStep):
                out += st.resonance[1:] if prev else st.resonance
            elif st.new_row or not prev:
                out += st.sources + st.products
            else:
                out += st.sources[kept_products(st.sources, prev):] + st.products
            prev = st.products
        return out

    def kept(self, i: int) -> int:
        """Step i (1-based): how many of the previous step's products its reactants keep (0 for
        the first step and new rows)."""
        st = self.steps[i - 1]
        if i == 1 or isinstance(st, ResonanceStep) or st.new_row:
            return 0
        return kept_products(st.sources, self.steps[i - 2].products)

    def arrow_sets(self) -> list[ArrowSet]:
        """Every set of curved arrows, in step order (ArrowSet)."""
        out = []
        for i, st in enumerate(self.steps or [], 1):
            if isinstance(st, ResonanceStep):
                off = 0
                for k, lst in enumerate(st.arrows):
                    if lst:
                        out.append(ArrowSet(i, (st.resonance[k],), (st.resonance[k + 1],), tuple(lst), off, k))
                    off += len(lst)
            elif st.arrows:
                out.append(ArrowSet(i, tuple(st.sources), tuple(st.products), tuple(st.arrows)))
        return out


def load_spec(path: Path | str) -> Spec:
    return Spec.model_validate(yaml.safe_load(Path(path).read_text()))


def resolve_molfile(ref: str, base_dir: Path) -> Path:
    for root in (base_dir, REPO_ROOT):
        p = (root / ref).expanduser()
        if p.exists():
            return p
    raise FileNotFoundError(f"molfile {ref!r} not found relative to {base_dir} or {REPO_ROOT}")


def _parsed(sp: Species, base_dir: Path, name: str = ""):
    try:
        if sp.smiles is not None:
            return chem.parse_smiles(sp.smiles, sp.kekule), False
        return chem.parse_molfile(resolve_molfile(sp.molfile, base_dir).read_text(), sp.kekule), True
    except ValueError as e:
        raise ValueError(f"species {name}: {e}" if name else str(e)) from None


def parsed_species(spec: Spec, base_dir: Path) -> dict[str, Chem.Mol]:
    """Each species as parsed (Kekulé, explicit H atoms kept, abbreviations not collapsed):
    the full chemistry, for electron bookkeeping."""
    return {n: _parsed(sp, base_dir, n)[0] for n, sp in spec.species.items()}


def check_labels(spec: Spec, base_dir: Path) -> None:
    """Every `labels` map number must be on an atom of some species."""
    maps = {a.GetAtomMapNum() for m in parsed_species(spec, base_dir).values() for a in m.GetAtoms()}
    unknown = sorted(set(spec.labels) - maps)
    if unknown:
        raise ValueError(f"labels: no atom has map number {', '.join(map(str, unknown))}")


def relabel(spec: Spec, mol: chem.Molecule) -> chem.Molecule:
    """The drawing's atoms with `labels` symbols in place of their stand-in elements."""
    if not spec.labels or not any(a.map in spec.labels for a in mol.atoms):
        return mol
    return replace(mol, atoms=[replace(a, element=spec.labels[a.map]) if a.map in spec.labels else a
                               for a in mol.atoms])


def molecules(spec: Spec, base_dir: Path) -> list[chem.Molecule]:
    """M1: every species laid out on its own."""
    bl = load_preset(spec.preset).lengths.bond_length
    out = []
    for name, sp in spec.species.items():
        o = sp.orientation or Orientation()
        if sp.smiles is not None:
            rot = "auto" if o.rotate is None else o.rotate
            mol = chem.from_smiles(sp.smiles, name, bl, sp.kekule, rot, o.flip_h, o.flip_v, sp.abbreviate)
        else:
            rot = 0.0 if o.rotate is None else o.rotate
            mol = chem.from_molfile(resolve_molfile(sp.molfile, base_dir), name, bl, sp.kekule,
                                    rot, o.flip_h, o.flip_v, sp.abbreviate)
        out.append(relabel(spec, chem.show_carbons(mol, sp.show_carbons)))
    check_labels(spec, base_dir)
    return out


def reference_for(spec: Spec, name: str) -> str | None:
    """The species `name` takes its orientation from: align_to, else the step source it shares
    the most atoms with (chem.best_reference), else none (its own layout)."""
    sp = spec.species[name]
    if sp.align_to:
        return sp.align_to
    if spec.align is False or not sp.align or not spec.steps:
        return None
    for st in spec.steps:
        if isinstance(st, ResonanceStep):
            if name in st.resonance[1:]:
                return st.resonance[st.resonance.index(name) - 1]
        elif name in st.products:
            return st.sources[0] if len(st.sources) == 1 else "|".join(st.sources)
    return None


def scheme_prepared(spec: Spec, base_dir: Path) -> dict[str, chem.Prepared]:
    """M2: species of a scheme, with each product aligned to the species it is drawn from
    (chem.prepare). Species are built in drawing order so references exist first."""
    bl = load_preset(spec.preset).lengths.bond_length
    built: dict[str, chem.Prepared] = {}
    for name in spec.order():
        sp = spec.species[name]
        mol, from_file = _parsed(sp, base_dir, name)
        ref = reference_for(spec, name) if not from_file or sp.align_to else None  # molfiles keep their layout
        cands = [built[r] for r in (ref.split("|") if ref else []) if r in built]
        if ref and not cands:
            raise ValueError(f"{name}: align_to {ref!r} must be drawn before it")
        if cands and sp.orientation is not None:
            raise ValueError(f"{name}: orientation is set but the species is aligned to {ref}; "
                             "set align: false to orient it on its own")
        o = sp.orientation or Orientation()
        built[name] = chem.prepare(mol, name, bl, sp.abbreviate, from_file, o.rotate, o.flip_h, o.flip_v, cands)
    return built


def scheme_molecules(spec: Spec, base_dir: Path) -> dict[str, chem.Molecule]:
    check_labels(spec, base_dir)
    return {n: relabel(spec, chem.show_carbons(p.molecule(), spec.species[n].show_carbons))
            for n, p in scheme_prepared(spec, base_dir).items()}
