"""Electron bookkeeping for curved arrows (design.md §6.3).

The model is a step's reactants taken together. Each atom has N nonbonding electrons, and
each bond holds E electrons (2 × order). Atoms are keyed by (species, atom index); arrows
address them by atom-map number. Implicit H atoms stay put: only H written as atoms ([H:5])
can move.

All of a step's arrows apply at once. A two_electron arrow moves 2 electrons and a fishhook
moves 1.
- Sources: lp(n) takes from n's lone pairs, bond(a-b) from the bond (the written order must
  match), and electron(n) from n's unpaired electron.
- Destinations: atom(n) adds nonbonding electrons; between(a,b) adds bond electrons (a new
  bond if there was none).

Afterwards, formal charge = V − N − E/2, summed over the atom's bonds (implicit H count 2
each). Unpaired electrons are N mod 2 on atoms an arrow touched (electrons pair up: no carbenes
or triplets); other atoms keep theirs.

Errors:
- a source lacks the electrons (no lone pair, no such bond, wrong order written, no radical);
- a fishhook draws from a lone pair;
- a bond ends with an odd electron count;
- a second-row atom ends above an octet, or H above 2;
- the result does not match the declared products.

The comparison removes H atoms on both sides (mapped ones included), drops stereo, and
compares canonical SMILES with atom maps, charges and radicals. Passing proves internal
consistency, not mechanistic correctness.

Resonance steps (check_resonance): every contributor has the first one's atom maps, the same
connectivity (bonds ignoring order, H counts) and the same total charge; and arrow list k,
applied to contributor k, gives contributor k + 1.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem

from .chem import _ORDER, lone_pairs, valence_electrons
from .spec import CurvedArrow, ResonanceStep, Spec, parsed_species

SECOND_ROW = {"Li", "Be", "B", "C", "N", "O", "F", "Ne"}
_BOND_TYPE = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
_ORDER_WORD = {1: "single", 2: "double", 3: "triple"}

Key = tuple[str, int]  # (species, atom index)


@dataclass
class AtomState:
    element: str
    map: int
    n: int  # nonbonding electrons
    h: int  # implicit / counted H (bonded, 2 electrons each)
    radicals: int
    lone_pairs: int


@dataclass
class Model:
    atoms: dict[Key, AtomState] = field(default_factory=dict)
    bonds: dict[frozenset, int] = field(default_factory=dict)  # electrons per bond
    by_map: dict[int, Key] = field(default_factory=dict)

    def bond_electrons(self, k: Key) -> int:
        return sum(e for pair, e in self.bonds.items() if k in pair) + 2 * self.atoms[k].h

    def charge(self, k: Key) -> int:
        a = self.atoms[k]
        return valence_electrons(Chem.Atom(a.element)) - a.n - self.bond_electrons(k) // 2

    def label(self, m: int) -> str:
        return f"{m} ({self.atoms[self.by_map[m]].element})"


@dataclass
class StepCheck:
    step: int
    sources: tuple[str, ...] = ()  # the species the arrows are drawn on
    products: tuple[str, ...] = ()
    problems: list[str] = field(default_factory=list)
    lp_tails: set[int] = field(default_factory=set)  # map numbers of atoms with an lp(n) tail
    targets: set[int] = field(default_factory=set)  # map numbers in atom(n) / between(a,b) heads
    gained_lone_pairs: set[int] = field(default_factory=set)  # map numbers whose lone pairs increase
    result: str = ""  # canonical SMILES of what the arrows give
    expected: str = ""  # ... and of the declared products
    suggested: list[tuple[str, str]] = field(default_factory=list)  # (product name or "?", SMILES) per fragment


def duplicate_maps(mols: dict[str, Chem.Mol]) -> list[str]:
    """Map numbers used more than once across these species, as 'n (species, species)'."""
    seen: dict[int, list[str]] = {}
    for name, mol in mols.items():
        for a in mol.GetAtoms():
            if a.GetAtomMapNum():
                seen.setdefault(a.GetAtomMapNum(), []).append(name)
    return [f"{m} ({', '.join(names)})" for m, names in sorted(seen.items()) if len(names) > 1]


def build_model(mols: dict[str, Chem.Mol]) -> Model:
    out = Model()
    for name, mol in mols.items():
        for a in mol.GetAtoms():
            k = (name, a.GetIdx())
            rad = a.GetNumRadicalElectrons()
            lp = lone_pairs(a)
            out.atoms[k] = AtomState(a.GetSymbol(), a.GetAtomMapNum(), 2 * lp + rad,
                                     a.GetTotalNumHs(includeNeighbors=False), rad, lp)
            if a.GetAtomMapNum():
                out.by_map[a.GetAtomMapNum()] = k
        for b in mol.GetBonds():
            pair = frozenset({(name, b.GetBeginAtomIdx()), (name, b.GetEndAtomIdx())})
            out.bonds[pair] = 2 * _ORDER[b.GetBondType()]
    return out


def apply_arrows(model: Model, arrows: list[CurvedArrow], where: str) -> tuple[Model, list[str], set[Key]]:
    """The model after all arrows, problems found, and the atoms the arrows touched."""
    problems: list[str] = []
    demand_lp, demand_rad, demand_bond = Counter(), Counter(), Counter()
    dn, de = Counter(), Counter()
    touched: set[Key] = set()

    def key(m: int, tag: str) -> Key | None:
        if m not in model.by_map:
            problems.append(f"{tag}: no atom with map number {m} in the step's reactants")
            return None
        return model.by_map[m]

    for k, arr in enumerate(arrows):
        tag = f"{where} arrows[{k}] ({arr.from_} -> {arr.to})"
        src, dst, q = arr.source, arr.target, arr.electrons
        if src.kind == "lp":
            if arr.type == "fishhook":
                problems.append(f"{tag}: a fishhook moves one electron, so it cannot start at a lone pair; "
                                "start it at a bond or electron(n)")
                continue
            a = key(src.atoms[0], tag)
            if a is None:
                continue
            demand_lp[a] += q
            dn[a] -= q
            touched.add(a)
        elif src.kind == "electron":
            a = key(src.atoms[0], tag)
            if a is None:
                continue
            demand_rad[a] += q
            dn[a] -= q
            touched.add(a)
        else:  # bond
            a, b = (key(m, tag) for m in src.atoms)
            if a is None or b is None:
                continue
            pair = frozenset({a, b})
            if pair not in model.bonds:
                problems.append(f"{tag}: atoms {model.label(src.atoms[0])} and {model.label(src.atoms[1])} "
                                "are not bonded")
                continue
            order = model.bonds[pair] // 2
            if order != src.order:
                problems.append(f"{tag}: atoms {src.atoms[0]} and {src.atoms[1]} are joined by a "
                                f"{_ORDER_WORD[order]} bond, not {_ORDER_WORD[src.order]} "
                                f"(write bond({src.atoms[0]}{'-=#'[order - 1]}{src.atoms[1]}))")
                continue
            demand_bond[pair] += q
            de[pair] -= q
            touched |= {a, b}
        if dst.kind == "atom":
            a = key(dst.atoms[0], tag)
            if a is not None:
                dn[a] += q
                touched.add(a)
        else:
            a, b = (key(m, tag) for m in dst.atoms)
            if a is not None and b is not None:
                de[frozenset({a, b})] += q
                touched |= {a, b}

    for a, need in demand_lp.items():
        have = 2 * model.atoms[a].lone_pairs
        if need > have:
            m = model.atoms[a].map
            problems.append(f"{where}: lp({m}): atom {model.label(m)} has "
                            f"{model.atoms[a].lone_pairs or 'no'} lone pair{'s' if model.atoms[a].lone_pairs != 1 else ''}"
                            f"{f', and {need // 2} arrows start there' if need > 2 else ''}")
    for a, need in demand_rad.items():
        if need > model.atoms[a].radicals:
            m = model.atoms[a].map
            problems.append(f"{where}: electron({m}): atom {model.label(m)} has "
                            f"{model.atoms[a].radicals or 'no'} unpaired electron{'s' if model.atoms[a].radicals != 1 else ''}")
    for pair, need in demand_bond.items():
        if need > model.bonds[pair]:
            ms = sorted(model.atoms[k].map for k in pair)
            problems.append(f"{where}: bond {ms[0]}-{ms[1]} holds {model.bonds[pair]} electrons, "
                            f"but the arrows take {need}")
    if problems:
        return model, problems, touched

    new = Model({k: AtomState(**vars(v)) for k, v in model.atoms.items()}, dict(model.bonds), dict(model.by_map))
    for a, d in dn.items():
        new.atoms[a].n += d
    for pair, d in de.items():
        new.bonds[pair] = new.bonds.get(pair, 0) + d
    for pair in [p for p, e in new.bonds.items() if e == 0]:
        del new.bonds[pair]
    odd = [p for p, e in new.bonds.items() if e % 2]
    for pair, e in new.bonds.items():
        ms = "-".join(str(new.atoms[k].map or "?") for k in sorted(pair))
        if e % 2:
            problems.append(f"{where}: bond {ms} ends with {e} electrons; a bond needs an even number "
                            f"(fishhooks come in pairs){_pairing_hint(model, pair, odd, de)}")
        elif e > 6:
            problems.append(f"{where}: bond {ms} ends as a {e // 2}-fold bond")
    for k in touched:
        a = new.atoms[k]
        a.radicals = a.n % 2
        a.lone_pairs = (a.n - a.radicals) // 2
        total = a.n + new.bond_electrons(k)
        name = f"atom {a.map or '?'} ({a.element})"
        if a.element == "H" and total > 2:
            problems.append(f"{where}: {name} ends with {total} electrons (H holds at most 2)")
        elif a.element in SECOND_ROW and total > 8:
            problems.append(f"{where}: {name} ends with {total} electrons (octet exceeded)")
    return new, problems, touched


def _pairing_hint(model: Model, pair: frozenset, odd: list[frozenset], de: Counter) -> str:
    """The fishhook that would even out an odd bond: from a half-broken bond to the half-formed
    bond next to it, else from the bond to one of its atoms, or to the new bond."""
    def maps(p):
        return sorted(model.atoms[k].map for k in p)

    def bond(p):
        a, b = maps(p)
        return f"bond({a}{'-=#'[model.bonds[p] // 2 - 1]}{b})"

    if not all(maps(pair)):
        return ""
    broken = de[pair] < 0
    partner = next((q for q in odd if q != pair and q & pair and all(maps(q)) and (de[q] < 0) != broken), None)
    if partner is not None:
        src, dst = (pair, partner) if broken else (partner, pair)
        return f"; add {{type: fishhook, from: \"{bond(src)}\", to: \"between({','.join(map(str, maps(dst)))})\"}}"
    a, b = maps(pair)
    if broken:
        return f"; add a fishhook from {bond(pair)} to atom({a}) or atom({b})"
    return f"; add a fishhook (from a bond or an electron) to between({a},{b})"


def to_mol(model: Model) -> Chem.Mol:
    rw = Chem.RWMol()
    index = {}
    for k, a in model.atoms.items():
        at = Chem.Atom(a.element)
        at.SetFormalCharge(model.charge(k))
        at.SetNumExplicitHs(a.h)
        at.SetNoImplicit(True)
        at.SetNumRadicalElectrons(a.radicals)
        at.SetAtomMapNum(a.map)
        index[k] = rw.AddAtom(at)
    for pair, e in model.bonds.items():
        a, b = sorted(pair)
        rw.AddBond(index[a], index[b], _BOND_TYPE[e // 2])
    return rw.GetMol()


def combine(mols: list[Chem.Mol]) -> Chem.Mol:
    out = Chem.Mol()
    for m in mols:
        out = Chem.CombineMols(out, m)
    return out


def canonical(mol: Chem.Mol) -> str:
    """Canonical SMILES with maps, charges and radicals; H atoms removed (mapped ones too)
    and stereo dropped, so explicit and implicit H compare equal."""
    m = Chem.RWMol(mol)
    m.UpdatePropertyCache(strict=False)
    for a in m.GetAtoms():
        a.SetNumExplicitHs(a.GetTotalNumHs())  # keep each atom's H count as given (implicit H of `C` too):
        a.SetNoImplicit(True)  # ... radicals must not gain H
    Chem.SanitizeMol(m)
    Chem.RemoveStereochemistry(m)
    ps = Chem.RemoveHsParameters()
    ps.removeMapped = True
    m = Chem.RemoveHs(m, ps)
    return Chem.MolToSmiles(m)


_TEMP_MAP = 99999  # marks added H atoms while writing SMILES, so they stay atoms


def product_smiles(model: Model) -> list[str]:
    """What the arrows give, one SMILES per molecule, to paste into a spec: map numbers,
    charges, radicals and mapped H atoms ([H:17]) kept, unmapped H implicit except on atoms that
    hold a mapped H, the Kekulé form the arrows made (no aromaticity perception)."""
    mol = Chem.RWMol(to_mol(model))
    Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_SETAROMATICITY)
    ps = Chem.RemoveHsParameters()
    ps.removeMapped = False
    mol = Chem.RemoveHs(mol, ps, sanitize=False)
    # an atom holding a mapped H atom gets its other H as atoms too: [H][O:1][H:2] draws H-O-H, not HO-H
    hosts = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() != "H" and a.GetTotalNumHs()
             and any(n.GetSymbol() == "H" and n.GetAtomMapNum() for n in a.GetNeighbors())]
    if hosts:
        n0 = mol.GetNumAtoms()
        mol = Chem.AddHs(mol, onlyOnAtoms=hosts)
        for a in list(mol.GetAtoms())[n0:]:
            a.SetAtomMapNum(_TEMP_MAP)
    return [Chem.MolToSmiles(f, kekuleSmiles=True).replace(f"[H:{_TEMP_MAP}]", "[H]")
            for f in Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)]


def _maps(mol: Chem.Mol) -> set[int]:
    return {a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomMapNum()}


def name_fragments(smiles: list[str], products: dict[str, Chem.Mol]) -> list[tuple[str, str]]:
    """Pair each fragment with the declared product sharing the most atom maps (each product
    used once; "?" when none does), in the declared products' order, unmatched ones last."""
    frags = [(s, _maps(Chem.MolFromSmiles(s, sanitize=False))) for s in smiles]
    free = dict(products)
    named: dict[str, str] = {}
    rest = []
    for s, maps in sorted(frags, key=lambda f: -len(f[1])):
        best = max(free, key=lambda n: len(maps & _maps(free[n])), default=None)
        if best is not None and maps & _maps(free[best]):
            named[best] = s
            del free[best]
        else:
            rest.append(("?", s))
    return [(n, named[n]) for n in products if n in named] + rest


def check_step(n: int, reactants: dict[str, Chem.Mol], products: dict[str, Chem.Mol],
               arrows: list[CurvedArrow], where: str | None = None, what: str = "the products") -> StepCheck:
    out = StepCheck(n, tuple(reactants), tuple(products))
    where = where or f"step {n}"
    for side, mols in (("reactants", reactants), ("products", products)):
        dup = duplicate_maps(mols)
        if dup:
            out.problems.append(f"{where}: map number {'; '.join(dup)} is used more than once in the step's {side}")
    if out.problems:
        return out
    model = build_model(reactants)
    out.lp_tails = {a.source.atoms[0] for a in arrows if a.source.kind == "lp"}
    out.targets = {m for a in arrows for m in a.target.atoms}
    new, problems, _ = apply_arrows(model, arrows, where)
    out.problems += problems
    if problems:
        return out
    try:
        out.result = canonical(to_mol(new))
        out.suggested = name_fragments(product_smiles(new), products)
    except Exception as e:  # RDKit sanitization: an impossible structure
        out.problems.append(f"{where}: the arrows give a structure RDKit rejects: {e}")
        return out
    out.expected = canonical(combine(list(products.values())))
    if out.result != out.expected:
        msg = f"{where}: the arrows give {out.result}, but {what} {'is' if what.startswith('the next') else 'are'} {out.expected}"
        have, want = _maps(to_mol(new)), set().union(*(_maps(m) for m in products.values()))
        if have - want:
            msg += f"; map number(s) {', '.join(map(str, sorted(have - want)))} of the reactants are missing from {what}"
        if want - have:
            msg += f"; map number(s) {', '.join(map(str, sorted(want - have)))} of {what} are not in the reactants"
        for k, arr in enumerate(arrows):  # r1 sn2: lp(1) -> atom(2) written for a new 1-2 bond
            a, b = arr.source.atoms[0], arr.target.atoms[0]
            if arr.source.kind in ("lp", "electron") and arr.target.kind == "atom" and a != b:
                msg += (f"; arrows[{k}] ({arr.from_} -> {arr.to}) puts the electrons on atom {b} as a lone pair, "
                        f"with no bond: to make the {a}-{b} bond write to: \"between({a},{b})\"")
        out.problems.append(msg)
    before = {a.map: a.lone_pairs for a in model.atoms.values() if a.map}
    out.gained_lone_pairs = {a.map for a in new.atoms.values() if a.map and a.lone_pairs > before[a.map]}
    return out


def skeleton(mol: Chem.Mol) -> str:
    """Canonical SMILES of the connectivity alone: all bonds single, no charges or radicals,
    each atom's H count kept (H atoms made implicit), atom maps kept."""
    m = Chem.RWMol(mol)
    ps = Chem.RemoveHsParameters()
    ps.removeMapped = True
    m = Chem.RWMol(Chem.RemoveHs(m, ps, sanitize=False))
    m.UpdatePropertyCache(strict=False)
    for a in m.GetAtoms():
        a.SetNumExplicitHs(a.GetTotalNumHs())
        a.SetNoImplicit(True)
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)
        a.SetNumRadicalElectrons(0)
    for b in m.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    return Chem.MolToSmiles(m, canonical=True)


def check_resonance(n: int, contributors: dict[str, Chem.Mol], arrow_lists: list[list[CurvedArrow]]) -> list[StepCheck]:
    """Step n's contributors (in order) and its arrow lists: one StepCheck per non-empty list
    (contributor k -> k + 1); structural problems go on the first check, or on a check of
    their own when there are no arrows."""
    names = list(contributors)
    first, m0 = names[0], contributors[names[0]]
    head = StepCheck(n, (first,), tuple(names[1:]))
    maps0 = sorted(a.GetAtomMapNum() for a in m0.GetAtoms() if a.GetAtomMapNum())
    q0 = Chem.GetFormalCharge(m0)
    for k, nm in enumerate(names[1:], 2):
        m = contributors[nm]
        tag = f"step {n}: contributor {nm} ({k} of {len(names)})"
        maps = sorted(a.GetAtomMapNum() for a in m.GetAtoms() if a.GetAtomMapNum())
        if maps != maps0:
            head.problems.append(f"{tag} has atom maps {maps}, but {first} has {maps0}; contributors keep every atom")
        elif skeleton(m) != skeleton(m0):
            head.problems.append(f"{tag} has a different connectivity from {first} ({skeleton(m)} vs {skeleton(m0)}): "
                                 "resonance moves electrons, not atoms")
        elif Chem.GetFormalCharge(m) != q0:
            head.problems.append(f"{tag} has total charge {Chem.GetFormalCharge(m):+d}, but {first} has {q0:+d}")
    out = []
    for k, arrows in enumerate(arrow_lists):
        if not arrows:
            continue
        a, b = names[k], names[k + 1]
        out.append(check_step(n, {a: contributors[a]}, {b: contributors[b]}, arrows,
                              where=f"step {n} ({a} -> {b})", what=f"the next contributor {b}"))
    if not out or head.problems:
        if out:
            out[0].problems = head.problems + out[0].problems
        else:
            out.append(head)
    return out


def check_spec(spec: Spec, base_dir: Path) -> list[StepCheck]:
    """Bookkeeping for every step with arrows, and every resonance step."""
    if not spec.steps or not (spec.arrow_sets() or any(isinstance(st, ResonanceStep) for st in spec.steps)):
        return []
    mols = parsed_species(spec, base_dir)
    out = []
    for i, st in enumerate(spec.steps, 1):
        if isinstance(st, ResonanceStep):
            out += check_resonance(i, {s: mols[s] for s in st.resonance}, st.arrows)
        elif st.arrows:
            out.append(check_step(i, {s: mols[s] for s in st.sources}, {s: mols[s] for s in st.products}, st.arrows))
    return out


def lone_pair_atoms(spec: Spec, checks: list[StepCheck]) -> dict[str, set[int]]:
    """species -> map numbers whose lone pairs are drawn: lp-tail and arrow-target atoms on the
    species an arrow set is drawn on, atoms that gain a lone pair on its drawn products, and
    each species' show_lone_pairs. (An omitted product keeps its name, which may be drawn
    elsewhere: its gained pairs belong to this step only.)"""
    out = {n: set(sp.show_lone_pairs) for n, sp in spec.species.items()}
    for c in checks:
        for n in c.sources:
            out[n] |= c.lp_tails | c.targets
        for n in c.products:
            if n not in getattr(spec.steps[c.step - 1], "omit", ()):
                out[n] |= c.gained_lone_pairs
    return out
