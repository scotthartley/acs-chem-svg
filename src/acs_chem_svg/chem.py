"""Chemistry via RDKit: parse, sanitize, 2D coordinates, rings, stereo (design.md §4.1).

The output is a plain `Molecule` in pt, y down, so drawing never touches RDKit.

- A SMILES carbon with one H written as an atom ([H:6][CH2:5]...) gets all its H's as atoms
  (unmapped), so a drawn C–H never hides the others. Molfiles keep their own H atoms.
- Kekulé bonds are kept as written. Aromatic SMILES are kekulized by RDKit, or
  per `kekule` (pairs of atom-map numbers that are double bonds).
- SMILES get CoordGen 2D coordinates. A molfile keeps its layout. Either way the median
  bond is scaled to the preset bond length and y is flipped to the SVG frame.
- Wedges come from the molfile when present, otherwise from RDKit's WedgeMolBonds.
- Carbons are unlabeled unless asked for (ChemDraw's ShowTerminalCarbonLabels is off): a
  molfile atom alias on a carbon (CH3, H3C, CH2, ..., as ChemDraw imports it), or
  `show_carbons` ("terminal" or atom ids). Only terminal (degree <= 1) carbons can be
  labeled: where ChemDraw puts the H of a CH2/CH with bonds on both sides is unmeasured
  (acs-style.md §6 Open 8).
- Abbreviations (M2): named groups (OMe, OBn, ...) collapse to one labeled atom at the
  attachment atom's position (RDKit rdAbbreviations). A molfile alias that names one does the same.
- Alignment (M2): a scheme product takes the coordinates of the atoms it shares with the
  species it is drawn from (atom maps, else the MCS), and only the rest is laid out.
- Lewis structures (`lewis`): every H becomes an atom and every atom is labeled; lewis_layout
  places the (acyclic) tree on a square lattice at the preset's lewis.bond_length.
"""

from __future__ import annotations

import itertools
import math
import re
import statistics
from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem, rdBase, Geometry
from rdkit.Chem import rdAbbreviations, rdCoordGen, rdDepictor, rdFMCS

from .geometry import Point, add, rotate

_SANITIZE = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
_ORDER = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3}
_MOLFILE_STEREO = {1: "wedge", 6: "hash"}
_CARBON_ALIAS = re.compile(r"C(H\d?)?|H\d?C")  # CH3, H3C, CH2, CH, C
STANDARD_BOND = 1.5  # A; molfile bond length that maps exactly onto the preset bond length
_RDKIT_DIR = {Chem.BondDir.BEGINWEDGE: "wedge", Chem.BondDir.BEGINDASH: "hash"}
# label -> group SMILES (* = attachment). Drawn left-justified with the first letter on the
# atom (OMe, OBn: scheme_multi_x*); the label text goes through richtext (CO2Me -> CO₂Me).
ABBREVIATIONS = {
    "OMe": "*OC", "OEt": "*OCC", "OBn": "*OCc1ccccc1", "OPh": "*Oc1ccccc1", "OAc": "*OC(C)=O",
    "SMe": "*SC", "NHAc": "*NC(C)=O", "CO2Me": "*C(=O)OC", "CO2Et": "*C(=O)OCC", "CO2H": "*C(=O)O",
    "CF3": "*C(F)(F)F", "CN": "*C#N", "NO2": "*[N+](=O)[O-]", "Ph": "*c1ccccc1", "Bn": "*Cc1ccccc1",
    "Et": "*CC", "Ac": "*C(C)=O",
    "TMS": "*[Si](C)(C)C", "TES": "*[Si](CC)(CC)CC", "TBS": "*[Si](C)(C)C(C)(C)C",
    "TIPS": "*[Si](C(C)C)(C(C)C)C(C)C",
}
# The right-justified form (bond leaving the label to the right), as ChemDraw displays each
# nickname (abbrev_x1.0 SVG; the CDXML keeps the key text). Formula-like ones are reversed;
# Ph, Bn, Et, Ac and the silyl groups keep their text.
ABBREVIATIONS_REVERSED = {
    "OMe": "MeO", "OEt": "EtO", "OBn": "BnO", "OPh": "PhO", "OAc": "AcO",
    "SMe": "MeS", "NHAc": "AcHN", "CO2Me": "MeO2C", "CO2Et": "EtO2C", "CO2H": "HO2C",
    "CF3": "F3C", "CN": "NC", "NO2": "O2N", "Ph": "Ph", "Bn": "Bn",
    "Et": "Et", "Ac": "Ac",
    "TMS": "TMS", "TES": "TES", "TBS": "TBS", "TIPS": "TIPS",
}
MCS_TIMEOUT = 2  # s
# Pauling electronegativities (CRC Handbook of Chemistry and Physics, "Electronegativity of the
# Elements"): chemistry constants for Lewis oxidation-state views, not style numbers.
ELECTRONEGATIVITY = {
    "H": 2.20, "Li": 0.98, "Be": 1.57, "B": 2.04, "C": 2.55, "N": 3.04, "O": 3.44, "F": 3.98,
    "Na": 0.93, "Mg": 1.31, "Si": 1.90, "P": 2.19, "S": 2.58, "Cl": 3.16, "K": 0.82, "Ca": 1.00,
    "Br": 2.96, "Rb": 0.82, "Sr": 0.95, "I": 2.66, "Cs": 0.79, "Ba": 0.89,
}
LEWIS_MODES = ("shared", "formal_charge", "oxidation_state")
LEWIS_SEARCH_LIMIT = 20000  # lewis_layout: placements tried before settling for the best found
ALIGN_TOL = 1e-4  # A; matched atoms must hold their reference position this closely


@dataclass
class Atom:
    idx: int
    id: str  # atom-map number if set, else idx + 1; used in scene sources
    element: str
    x: float
    y: float
    charge: int = 0
    num_h: int = 0  # implicit + explicit-count H (H atoms are separate Atoms)
    map: int | None = None
    show_label: bool = False  # explicit carbon label (CH3); heteroatoms are always labeled
    abbrev: str | None = None  # abbreviation label text (OMe); `element` is its attachment atom
    radicals: int = 0  # unpaired electrons (RDKit GetNumRadicalElectrons)
    lone_pairs: int = 0  # (valence electrons - charge - bond orders - H - radicals) / 2
    charge_label: str | None = None  # drawn charge text in place of `charge` (spec charge_labels: +3); display only

    @property
    def p(self) -> Point:
        return (self.x, self.y)


@dataclass
class Bond:
    idx: int
    a: int  # for stereo bonds, the narrow end (stereocenter)
    b: int
    order: int = 1
    stereo: str | None = None  # "wedge" | "hash"
    ring: bool = False

    @property
    def id(self) -> str:
        return str(self.idx + 1)

    def other(self, i: int) -> int:
        return self.b if i == self.a else self.a


@dataclass
class Molecule:
    name: str
    atoms: list[Atom]
    bonds: list[Bond]
    rings: list[list[int]] = field(default_factory=list)  # SSSR, atom indices in ring order
    lewis: str | None = None  # Lewis structure: "shared" | "formal_charge" | "oxidation_state" (draw_mol)
    circled_charges: bool = False  # ±1 signs drawn in a ring (spec charges.circled / circled_charges)

    def neighbors(self, i: int) -> list[tuple[int, Bond]]:
        return [(b.other(i), b) for b in self.bonds if i in (b.a, b.b)]

    def degree(self, i: int) -> int:
        return sum(1 for b in self.bonds if i in (b.a, b.b))

    def by_id(self, atom_id: str) -> Atom:
        return next(a for a in self.atoms if a.id == atom_id)


# --------------------------------------------------------------------------- input


def _valence_hint(mol: Chem.Mol, e: Exception) -> str:
    """RDKit's valence message, with the atom's map number, and a hint when the atom both
    writes H in its bracket and has an [H:n] atom attached (r1 nitration: [H:17][CH:11])."""
    msg = str(e).strip()
    m = re.search(r"atom # (\d+)", msg)
    if not m:
        return msg
    a = mol.GetAtomWithIdx(int(m.group(1)))
    where = f"atom {a.GetSymbol()} with map number {a.GetAtomMapNum()}" if a.GetAtomMapNum() else \
        f"the unmapped {a.GetSymbol()} (atom {a.GetIdx() + 1} in the SMILES)"
    msg = f"{msg}: {where} has too many bonds"
    h_atoms = [n for n in a.GetNeighbors() if n.GetSymbol() == "H"]
    if h_atoms and a.GetNumExplicitHs():
        hm = h_atoms[0].GetAtomMapNum()
        msg += (f"; it has an H atom attached ([H{f':{hm}' if hm else ''}]) and also writes H in its bracket: "
                f"drop the bracket H ([CH:11] -> [C:11])")
    return msg


def _sanitize(mol: Chem.Mol, kekule: list[str] | None) -> Chem.Mol:
    aromatic = {b.GetIdx() for b in mol.GetBonds() if b.GetIsAromatic() or b.GetBondType() == Chem.BondType.AROMATIC}
    block = rdBase.BlockLogs()  # the error is raised with a better message; don't also log it
    try:
        Chem.SanitizeMol(mol, _SANITIZE)
    except Chem.AtomValenceException as e:
        raise ValueError(_valence_hint(mol, e)) from None
    finally:
        del block
    Chem.Kekulize(mol, clearAromaticFlags=True)
    if kekule:
        by_map = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}
        for i in aromatic:
            mol.GetBondWithIdx(i).SetBondType(Chem.BondType.SINGLE)
        for pair in kekule:
            m1, m2 = (int(t) for t in pair.split("="))
            try:
                b = mol.GetBondBetweenAtoms(by_map[m1], by_map[m2])
            except KeyError as e:
                raise ValueError(f"kekule {pair!r}: no atom with map number {e.args[0]}") from None
            if b is None:
                raise ValueError(f"kekule {pair!r}: atoms {m1} and {m2} are not bonded")
            b.SetBondType(Chem.BondType.DOUBLE)
        Chem.SanitizeMol(mol, _SANITIZE)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)  # E/Z for CoordGen, R/S for wedges
    return mol


def parse_smiles(smiles: str, kekule: list[str] | None = None) -> Chem.Mol:
    params = Chem.SmilesParserParams()
    params.removeHs = False
    params.sanitize = False
    mol = Chem.MolFromSmiles(smiles, params)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES {smiles!r}")
    mol = _sanitize(mol, kekule)
    partial = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "C" and a.GetTotalNumHs(includeNeighbors=False)
               and any(n.GetAtomicNum() == 1 for n in a.GetNeighbors())]
    return Chem.AddHs(mol, onlyOnAtoms=partial) if partial else mol


def parse_molfile(text: str, kekule: list[str] | None = None) -> Chem.Mol:
    mol = Chem.MolFromMolBlock(text, sanitize=False, removeHs=False)
    if mol is None:
        raise ValueError("RDKit could not parse the molfile")
    return _sanitize(mol, kekule)


# --------------------------------------------------------------------------- coordinates


def _coords(mol: Chem.Mol) -> list[Point]:
    """Conformer coordinates with y flipped to the SVG frame (y down)."""
    conf = mol.GetConformer()
    return [(conf.GetAtomPosition(i).x, -conf.GetAtomPosition(i).y) for i in range(mol.GetNumAtoms())]


def _set_coords(mol: Chem.Mol, pts: list[Point]) -> None:
    conf = mol.GetConformer()
    for i, (x, y) in enumerate(pts):
        conf.SetAtomPosition(i, (x, -y, 0.0))


def _bond_angles(mol, pts):
    return [math.degrees(math.atan2(-(pts[b.GetEndAtomIdx()][1] - pts[b.GetBeginAtomIdx()][1]),
                                    pts[b.GetEndAtomIdx()][0] - pts[b.GetBeginAtomIdx()][0]))
            for b in mol.GetBonds()]


def _off(angle: float, step: float, phase: float = 0.0) -> float:
    """Distance (deg) from `angle` to the nearest phase + n*step."""
    return abs(((angle - phase + step / 2) % step) - step / 2)


def auto_rotation(mol: Chem.Mol, pts: list[Point]) -> float:
    """Rotation (deg) in ChemDraw's habit: zigzag bonds at ±30°, rings with vertical sides
    (the most bonds at 30° + n·60°), then the most bonds on the 30° grid, then a layout
    wider than tall, then the smallest turn. All-collinear molecules lie horizontal."""
    angles = _bond_angles(mol, pts)
    if not angles:
        return 0.0
    if all(_off(a - angles[0], 180) < 0.5 for a in angles):
        return -angles[0]
    best = None
    for a in angles:
        for target in range(-180, 180, 30):
            r = target - a
            zig = sum(1 for b in angles if _off(b + r, 60, 30) < 0.5)
            grid = sum(1 for b in angles if _off(b + r, 30) < 0.5)
            rp = [rotate(p, r) for p in pts]
            w = max(p[0] for p in rp) - min(p[0] for p in rp)
            h = max(p[1] for p in rp) - min(p[1] for p in rp)
            key = (-zig, -grid, w < h - 1e-6, _off(r, 360))
            if best is None or key < best[0]:
                best = (key, r)
    return best[1]


def layout_2d(mol: Chem.Mol) -> None:
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)


def orient(mol: Chem.Mol, rotate_deg: float | str | None = None, flip_h: bool = False,
           flip_v: bool = False) -> None:
    """Rotate (counter-clockwise on screen) and mirror the conformer in place, before
    wedges are assigned. `rotate_deg="auto"` uses auto_rotation()."""
    pts = _coords(mol)
    if rotate_deg == "auto":
        rotate_deg = auto_rotation(mol, pts)
    if rotate_deg:
        pts = [rotate(p, float(rotate_deg)) for p in pts]
    if flip_h:
        pts = [(-x, y) for x, y in pts]
    if flip_v:
        pts = [(x, -y) for x, y in pts]
    _set_coords(mol, pts)


# --------------------------------------------------------------------------- output


def _stereo(mol: Chem.Mol) -> dict[int, tuple[str, int]]:
    """bond idx -> (kind, narrow-end atom idx)."""
    out = {}
    for b in mol.GetBonds():
        props = b.GetPropsAsDict(includePrivate=True, includeComputed=False)
        kind = _MOLFILE_STEREO.get(props.get("_MolFileBondStereo"))
        if kind:
            out[b.GetIdx()] = (kind, b.GetBeginAtomIdx())
    if out:
        return out
    if not any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()):
        return out
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    wm = Chem.Mol(mol)
    Chem.WedgeMolBonds(wm, wm.GetConformer())
    for b in wm.GetBonds():
        kind = _RDKIT_DIR.get(b.GetBondDir())
        if kind:
            out[b.GetIdx()] = (kind, b.GetBeginAtomIdx())
    return out


def layout_scale(mol: Chem.Mol, bond_length: float) -> float:
    """pt per layout unit: the median bond becomes `bond_length`. A 1.5 A layout (RDKit,
    templates/) maps exactly, so a template keeps ChemDraw's page positions; the median of
    slightly uneven bonds would drift far atoms by ~0.03 pt."""
    pts = _coords(mol)
    lengths = [math.dist(pts[b.GetBeginAtomIdx()], pts[b.GetEndAtomIdx()]) for b in mol.GetBonds()]
    med = statistics.median(lengths) if lengths else STANDARD_BOND  # a lone atom (Br⁻ template): 1.5 A units
    return bond_length / (STANDARD_BOND if abs(med / STANDARD_BOND - 1) < 0.002 else med)


def to_molecule(mol: Chem.Mol, name: str, bond_length: float, scale: float | None = None) -> Molecule:
    """`scale` (pt per layout unit) defaults to layout_scale(); aligned species pass their
    reference's, so shared atoms land on the same pt coordinates."""
    k = layout_scale(mol, bond_length) if scale is None else scale
    pts = [(x * k, y * k) for x, y in _coords(mol)]
    stereo = _stereo(mol)
    atoms = []
    for a in mol.GetAtoms():
        i = a.GetIdx()
        mp = a.GetAtomMapNum() or None
        show, abbrev = _alias_label(a)
        element = abbrev_element(abbrev) if abbrev and a.GetAtomicNum() == 0 else a.GetSymbol()
        atoms.append(Atom(i, str(mp) if mp else str(i + 1), element, pts[i][0], pts[i][1],
                          a.GetFormalCharge(), 0 if abbrev else a.GetTotalNumHs(includeNeighbors=False), mp,
                          show, abbrev, a.GetNumRadicalElectrons(), 0 if abbrev else lone_pairs(a)))
    bonds = []
    for b in mol.GetBonds():
        order = _ORDER.get(b.GetBondType())
        if order is None:
            raise ValueError(f"Unsupported bond type {b.GetBondType()} (bond {b.GetIdx() + 1})")
        a, e = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        kind = None
        if b.GetIdx() in stereo:
            kind, narrow = stereo[b.GetIdx()]
            a, e = narrow, (e if narrow == a else a)
        bonds.append(Bond(b.GetIdx(), a, e, order, kind, b.IsInRing()))
    rings = [list(r) for r in Chem.GetSymmSSSR(mol)]
    return Molecule(name, atoms, bonds, rings)


def valence_electrons(a: Chem.Atom) -> int:
    return Chem.GetPeriodicTable().GetNOuterElecs(a.GetAtomicNum())


def bond_order_sum(a: Chem.Atom) -> int:
    """Bond orders to every neighbor, plus one per implicit or counted H."""
    return sum(_ORDER[b.GetBondType()] for b in a.GetBonds()) + a.GetTotalNumHs(includeNeighbors=False)


def lone_pairs(a: Chem.Atom) -> int:
    """Nonbonding electron pairs: (V - formal charge - bond orders - radicals) / 2, for a
    Kekulé mol (H has none; dummy atoms none)."""
    if a.GetAtomicNum() in (0, 1) and a.GetFormalCharge() >= 0:
        return 0
    n = valence_electrons(a) - a.GetFormalCharge() - bond_order_sum(a) - a.GetNumRadicalElectrons()
    return max(0, n // 2)


def _alias_label(a: Chem.Atom) -> tuple[bool, str | None]:
    """(explicit carbon label, abbreviation). An abbreviation comes from rdAbbreviations
    (atomLabel) or a molfile alias naming one (OMe on the O, as cdxml_to_molfile writes ChemDraw
    nicknames). A carbon alias that is its own label (CH3, H3C, ...) shows the carbon. Any other
    alias is a condensed label we don't draw."""
    if a.HasProp("atomLabel") and a.GetAtomicNum() == 0:
        return False, a.GetProp("atomLabel")
    if not a.HasProp("molFileAlias"):
        return False, None
    alias = a.GetProp("molFileAlias")
    if alias in ABBREVIATIONS:
        if abbrev_element(alias) != a.GetSymbol():
            raise ValueError(f"atom {a.GetIdx() + 1}: alias {alias!r} attaches through "
                             f"{abbrev_element(alias)}, not {a.GetSymbol()}")
        return False, alias
    if a.GetSymbol() == "C" and _CARBON_ALIAS.fullmatch(alias):
        n_h = int(re.sub(r"\D", "", alias) or ("1" if "H" in alias else "0"))
        if n_h != a.GetTotalNumHs(includeNeighbors=False):
            raise ValueError(f"atom {a.GetIdx() + 1}: alias {alias!r} but the carbon has "
                             f"{a.GetTotalNumHs(includeNeighbors=False)} H")
        if a.GetDegree() > 1:
            raise NotImplementedError(f"atom {a.GetIdx() + 1}: alias {alias!r} on a non-terminal carbon "
                                      "(unmeasured, acs-style.md §6 Open 8)")
        return True, None
    raise ValueError(f"atom {a.GetIdx() + 1}: alias {alias!r} is a nickname or condensed label we don't "
                     f"draw; known abbreviations: {', '.join(ABBREVIATIONS)}")


def show_carbons(m: Molecule, which: str | list[str | int] | None) -> Molecule:
    """Label carbons: "terminal" (degree <= 1) or a list of atom ids (the atom-map number if
    set, else the 1-based atom index) of terminal carbons. None / "none" leaves them as they are."""
    if which in (None, "none"):
        return m
    ids = None if isinstance(which, str) else {str(x) for x in which}
    if which != "terminal" and ids is None:
        raise ValueError(f"show_carbons must be 'none', 'terminal' or a list of atom ids, not {which!r}")
    if ids is not None:
        unknown = ids - {a.id for a in m.atoms if a.element == "C"}
        if unknown:
            raise ValueError(f"show_carbons: no carbon with id {', '.join(sorted(unknown))} in {m.name}")
        inner = sorted(a.id for a in m.atoms if a.id in ids and m.degree(a.idx) > 1)
        if inner:
            raise NotImplementedError(f"show_carbons: carbon {', '.join(inner)} in {m.name} is not terminal; "
                                      "CH2/CH labels are unmeasured (acs-style.md §6 Open 8)")
    for a in m.atoms:
        if a.element == "C" and m.degree(a.idx) <= 1 and (which == "terminal" or a.id in ids):
            a.show_label = True
    return m


def abbrev_element(label: str) -> str:
    """Element of the atom an abbreviation attaches through (OMe: O, CO2Me: C, Ph: C)."""
    star = Chem.MolFromSmiles(ABBREVIATIONS[label]).GetAtomWithIdx(0)
    return star.GetNeighbors()[0].GetSymbol()


def condense(mol: Chem.Mol, labels: list[str] | None, name: str = "") -> Chem.Mol:
    """Collapse each named group (ABBREVIATIONS) to one dummy atom carrying atomLabel, at the
    attachment atom (its index, coordinates and atom-map number are kept). Groups are found by RDKit
    rdAbbreviations on an aromatic copy (its patterns are aromatic); the Kekulé bonds drawn
    are those of `mol`. Every label must be found."""
    if not labels:
        return mol
    defs = rdAbbreviations.ParseAbbreviations("\n".join(f"{lab} {ABBREVIATIONS[lab]}" for lab in labels))
    ar = Chem.Mol(mol)
    Chem.SetAromaticity(ar)
    found = Chem.GetMolSubstanceGroups(rdAbbreviations.LabelMolAbbreviations(ar, defs, maxCoverage=1.0))
    missing = [lab for lab in labels if lab not in {sg.GetProp("LABEL") for sg in found}]
    if missing:
        raise ValueError(f"{name}: abbreviate {', '.join(missing)}: no such terminal group in the structure")
    rw = Chem.RWMol(mol)
    drop = []
    for sg in found:
        group = set(sg.GetAtoms())
        attach = next(i for i in group
                      if any(n.GetIdx() not in group for n in mol.GetAtomWithIdx(i).GetNeighbors()))
        dummy = Chem.Atom(0)
        dummy.SetNoImplicit(True)
        dummy.SetProp("atomLabel", sg.GetProp("LABEL"))
        dummy.SetAtomMapNum(mol.GetAtomWithIdx(attach).GetAtomMapNum())  # aligned species keep the label in place
        rw.ReplaceAtom(attach, dummy)
        drop += [i for i in group if i != attach]
    for i in sorted(drop, reverse=True):
        rw.RemoveAtom(i)
    out = rw.GetMol()
    Chem.SanitizeMol(out, _SANITIZE)
    return out


def _orient_molfile(mol: Chem.Mol, rotate_deg, flip_h: bool, flip_v: bool) -> None:
    if rotate_deg or flip_h or flip_v:
        # mirror images of wedged drawings invert stereo, so wedges must be re-derived
        if (flip_h ^ flip_v) and any(b.HasProp("_MolFileBondStereo") for b in mol.GetBonds()):
            Chem.ReapplyMolBlockWedging(mol)
            Chem.AssignChiralTypesFromBondDirs(mol)
            for b in mol.GetBonds():
                b.SetBondDir(Chem.BondDir.NONE)
                b.ClearProp("_MolFileBondStereo")
        orient(mol, rotate_deg, flip_h, flip_v)


def from_smiles(smiles: str, name: str, bond_length: float, kekule: list[str] | None = None,
                rotate_deg: float | str | None = "auto", flip_h: bool = False, flip_v: bool = False,
                abbreviate: list[str] | None = None, lewis: str | None = None) -> Molecule:
    """`lewis` (a LEWIS_MODES value): a Lewis structure with every H an atom and every atom
    labeled, laid out by lewis_layout with `bond_length` as the lattice step (lewis.bond_length)."""
    if lewis is not None:
        mol = Chem.AddHs(parse_smiles(smiles))
        Chem.RemoveStereochemistry(mol)
        lewis_layout(mol, bond_length, name)
        orient(mol, 0.0 if rotate_deg in (None, "auto") else rotate_deg, flip_h, flip_v)  # a multiple of 90
        _set_coords(mol, [(round(x, 6) + 0.0, round(y, 6) + 0.0) for x, y in _coords(mol)])
        m = to_molecule(mol, name, bond_length, scale=1.0)
        for a in m.atoms:
            a.show_label = True
        m.lewis = lewis
        return m
    mol = condense(parse_smiles(smiles, kekule), abbreviate, name)
    layout_2d(mol)
    orient(mol, rotate_deg, flip_h, flip_v)
    return to_molecule(mol, name, bond_length)


def from_molfile(path_or_text: str | Path, name: str, bond_length: float, kekule: list[str] | None = None,
                 rotate_deg: float | str | None = None, flip_h: bool = False, flip_v: bool = False,
                 abbreviate: list[str] | None = None) -> Molecule:
    text = Path(path_or_text).read_text() if isinstance(path_or_text, Path) else path_or_text
    mol = condense(parse_molfile(text, kekule), abbreviate, name)
    _orient_molfile(mol, rotate_deg, flip_h, flip_v)
    return to_molecule(mol, name, bond_length)


# --------------------------------------------------------------------------- Lewis layout

_RIGHT, _LEFT, _UP, _DOWN = (1, 0), (-1, 0), (0, -1), (0, 1)  # lattice steps, y down


def _atom_name(a: Chem.Atom) -> str:
    return f"{a.GetSymbol()} (map {a.GetAtomMapNum()})" if a.GetAtomMapNum() else \
        f"{a.GetSymbol()} (atom {a.GetIdx() + 1} in the SMILES)"


def lewis_layout(mol: Chem.Mol, bond_length: float, name: str = "") -> None:
    """Lewis structure layout (user decision 2026-10-26): every bond on a square lattice at
    `bond_length`, whatever its order. Acyclic molecules with at most 4 neighbors per atom.

    - Root: the heavy-atom tree center (least eccentricity; ties: more neighbors, then lower
      index). Its heavy neighbors go right, left, up, down; H atoms fill what is left.
    - Each other atom: its largest heavy branch continues straight, the others go up then down
      (right then left on a vertical bond); H atoms fill the remaining directions, straight first.
    - Two atoms never share a lattice point. Non-bonded atoms one step apart would look bonded:
      the search (preferred orders first, every order of an atom's heavy branches) keeps the
      placement with the fewest such pairs; draw_mol keeps their facing sides free of lone pairs.
      A methyl next to a three-way center (acetic acid) cannot avoid one.
    Writes the conformer (pt, SVG frame flipped as _set_coords expects)."""
    where = f"lewis {name}: " if name else "lewis: "
    if mol.GetRingInfo().NumRings():
        a = next(a for a in mol.GetAtoms() if a.IsInRing())
        raise ValueError(f"{where}rings are not supported (atom {_atom_name(a)} is in a ring); "
                         "Lewis structures are drawn for acyclic molecules only")
    for a in mol.GetAtoms():
        if a.GetDegree() > 4:
            raise ValueError(f"{where}atom {_atom_name(a)} has {a.GetDegree()} neighbors; a Lewis structure "
                             "places at most 4 (one per side)")
    n = mol.GetNumAtoms()
    nbrs = [[x.GetIdx() for x in a.GetNeighbors()] for a in mol.GetAtoms()]
    heavy = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1} or set(range(n))

    def ecc(i: int) -> int:
        dist, todo = {i: 0}, [i]
        while todo:
            j = todo.pop(0)
            for k in nbrs[j]:
                if k in heavy and k not in dist:
                    dist[k] = dist[j] + 1
                    todo.append(k)
        return max(dist.values())
    root = min(sorted(heavy), key=lambda i: (ecc(i), -len(nbrs[i])))

    parent, order, todo = {root: None}, [], [root]
    while todo:  # breadth-first; atoms not reached are in other fragments
        i = todo.pop(0)
        order.append(i)
        for k in nbrs[i]:
            if k not in parent:
                parent[k] = i
                todo.append(k)
    if len(order) != n:
        raise ValueError(f"{where}draw each molecule of a salt or mixture as its own species")

    def size(i: int) -> tuple[int, int]:  # (heavy atoms, all atoms) in the branch from i
        out, stack = [0, 0], [i]
        while stack:
            j = stack.pop()
            out[0] += j in heavy
            out[1] += 1
            stack += [k for k in nbrs[j] if parent.get(k) == j]
        return out[0], out[1]
    kids = {i: sorted((k for k in nbrs[i] if parent[k] == i), key=lambda k: (k not in heavy, [-x for x in size(k)], k))
            for i in order}

    def options(i: int, into: tuple[int, int] | None):
        """Each way to give atom i's children directions: heavy children over every ordering
        (preferred first), H children in the remaining directions in order."""
        dirs = [_RIGHT, _LEFT, _UP, _DOWN] if into is None else \
            [into] + ([_UP, _DOWN] if into[1] == 0 else [_RIGHT, _LEFT])
        hv = [k for k in kids[i] if k in heavy]
        hs = [k for k in kids[i] if k not in heavy]
        for perm in itertools.permutations(dirs, len(hv)):
            rest = [d for d in dirs if d not in perm]
            yield list(zip(hv, perm)) + list(zip(hs, rest))

    pos: dict[int, tuple[int, int]] = {root: (0, 0)}
    occupied = {(0, 0): root}
    best: list = [None, None]  # (cost, positions)
    tried = [0]

    def search(k: int, cost: int) -> None:
        if best[0] is not None and (cost >= best[0] or best[0] == 0):
            return
        if k == len(order):
            best[0], best[1] = cost, dict(pos)
            return
        tried[0] += 1
        if tried[0] > LEWIS_SEARCH_LIMIT and best[0] is not None:
            return
        i = order[k]
        pi = pos[i]
        into = None if parent[i] is None else tuple(a - b for a, b in zip(pi, pos[parent[i]]))
        for opt in options(i, into):
            new = {c: (pi[0] + d[0], pi[1] + d[1]) for c, d in opt}
            if any(p in occupied for p in new.values()):
                continue
            added = 0
            for c, p in new.items():
                for d in (_RIGHT, _LEFT, _UP, _DOWN):
                    j = occupied.get((p[0] + d[0], p[1] + d[1]))
                    if j is not None and j != i:
                        added += 1
            for c, p in new.items():
                pos[c], occupied[p] = p, c
            search(k + 1, cost + added)
            for c, p in new.items():
                del pos[c], occupied[p]
    search(0, 0)
    if best[1] is None:
        raise ValueError(f"{where}no square-lattice layout keeps every atom on its own point (too many branches "
                         "around one center, as in neopentane); draw it with bond lines instead")
    if mol.GetNumConformers() == 0:
        mol.AddConformer(Chem.Conformer(n), assignId=True)
    _set_coords(mol, [(best[1][i][0] * bond_length, best[1][i][1] * bond_length) for i in range(n)])


# --------------------------------------------------------------------------- scheme alignment (M2)


@dataclass
class Prepared:
    """An RDKit mol with its final 2D layout, and the pt-per-unit scale to draw it at."""
    mol: Chem.Mol
    name: str
    bond_length: float
    scale: float
    aligned_to: str | None = None
    matched: list[tuple[int, int]] = field(default_factory=list)  # (reference idx, idx)

    def molecule(self) -> Molecule:
        return to_molecule(self.mol, self.name, self.bond_length, self.scale)


def _mcs_copy(mol: Chem.Mol) -> Chem.Mol:
    """Copy for MCS: aromaticity perceived (so Kekulé forms don't matter) and each atom's
    isotope set to its element, or to a per-label code for abbreviations (OMe != OBn)."""
    labels = sorted(ABBREVIATIONS)
    c = Chem.Mol(mol)
    for a in c.GetAtoms():
        if a.GetAtomicNum():
            a.SetIsotope(a.GetAtomicNum())
        else:
            lab = a.GetProp("atomLabel") if a.HasProp("atomLabel") else ""
            a.SetIsotope(200 + labels.index(lab) if lab in labels else 199)
    Chem.FastFindRings(c)
    Chem.SetAromaticity(c)
    return c


def correspondence(ref: Chem.Mol, mol: Chem.Mol) -> list[tuple[int, int]]:
    """(reference atom idx, atom idx) pairs: shared atom-map numbers if at least two, else the
    maximum common substructure (same elements and abbreviations, same bond orders with
    aromatic rings matched as aromatic, complete rings only). Fewer than 3 atoms: none."""
    maps_r = {a.GetAtomMapNum(): a.GetIdx() for a in ref.GetAtoms() if a.GetAtomMapNum()}
    maps_m = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}
    common = sorted(maps_r.keys() & maps_m.keys())
    if len(common) >= 2:
        return [(maps_r[k], maps_m[k]) for k in common]
    r, m = _mcs_copy(ref), _mcs_copy(mol)
    res = rdFMCS.FindMCS([r, m], atomCompare=rdFMCS.AtomCompare.CompareIsotopes,
                         bondCompare=rdFMCS.BondCompare.CompareOrder, ringMatchesRingOnly=True,
                         completeRingsOnly=True, timeout=MCS_TIMEOUT)
    if res.numAtoms < 3:
        return []
    patt = Chem.MolFromSmarts(res.smartsString)
    return list(zip(r.GetSubstructMatch(patt), m.GetSubstructMatch(patt)))


def _scale_conformer(mol: Chem.Mol, k: float) -> None:
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x * k, p.y * k, 0.0))


def _max_deviation(mol: Chem.Mol, ref: Chem.Mol, pairs) -> float:
    c, cr = mol.GetConformer(), ref.GetConformer()
    return max(math.dist((c.GetAtomPosition(m).x, c.GetAtomPosition(m).y),
                         (cr.GetAtomPosition(r).x, cr.GetAtomPosition(r).y)) for r, m in pairs)


def align_coords(mol: Chem.Mol, ref: Chem.Mol, pairs: list[tuple[int, int]]) -> None:
    """Lay out `mol` with the matched atoms fixed at the reference's coordinates. CoordGen
    (coordMap) first; it holds them only to ~0.005 A, so if that misses ALIGN_TOL use RDKit's
    GenerateDepictionMatching2DStructure, which holds them exactly. RDKit lays new atoms out at
    1.5 A bonds, so the reference is scaled to 1.5 A for it and the result scaled back."""
    conf = ref.GetConformer()
    params = rdCoordGen.CoordGenParams()
    params.SetCoordMap({m: Geometry.Point2D(conf.GetAtomPosition(r).x, conf.GetAtomPosition(r).y)
                        for r, m in pairs})
    rdCoordGen.AddCoords(mol, params)
    if _max_deviation(mol, ref, pairs) <= ALIGN_TOL:
        return
    pts = _coords(ref)
    med = statistics.median(math.dist(pts[b.GetBeginAtomIdx()], pts[b.GetEndAtomIdx()]) for b in ref.GetBonds())
    tmp = Chem.Mol(ref)
    _scale_conformer(tmp, STANDARD_BOND / med)
    rdDepictor.GenerateDepictionMatching2DStructure(mol, tmp, list(pairs), -1, True)
    _scale_conformer(mol, med / STANDARD_BOND)
    dev = _max_deviation(mol, ref, pairs)
    if dev > ALIGN_TOL:
        raise RuntimeError(f"alignment left matched atoms {dev:.4f} A from their reference positions")


def stable_core(ref: Chem.Mol, mol: Chem.Mol, pairs: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """The matched atoms worth pinning: those whose neighbor count is unchanged, and of them only
    the largest connected group. An atom that changes its neighbors (the sp3 C of a sigma
    complex becoming sp2 in nitrobenzene) and what hangs off it are laid out afresh, so their
    bonds get proper angles. All pairs if that leaves fewer than 3 atoms."""
    keep = {m: r for r, m in pairs if ref.GetAtomWithIdx(r).GetDegree() == mol.GetAtomWithIdx(m).GetDegree()}
    if len(keep) == len(pairs):
        return pairs
    groups, seen = [], set()
    for start in keep:
        if start in seen:
            continue
        group, todo = set(), [start]
        while todo:
            a = todo.pop()
            if a in group:
                continue
            group.add(a)
            todo += [n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in keep]
        seen |= group
        groups.append(group)
    best = max(groups, key=len, default=set())
    if len(best) < 3:
        return pairs
    return [(r, m) for r, m in pairs if m in best]


LINEAR_TOL = 10.0  # degrees an sp atom may deviate from 180° after alignment (a choice, not a style value)


def bent_linear_atoms(mol: Chem.Mol) -> list[int]:
    """Atoms that must be drawn straight (two neighbors and a triple bond or two double bonds:
    O=N+=O, C≡N, allenes) but whose layout bends them by more than LINEAR_TOL."""
    pts = _coords(mol)
    out = []
    for a in mol.GetAtoms():
        bonds = a.GetBonds()
        orders = sorted(b.GetBondType() for b in bonds)
        if len(bonds) != 2 or not (Chem.BondType.TRIPLE in orders
                                   or orders == [Chem.BondType.DOUBLE, Chem.BondType.DOUBLE]):
            continue
        p, (q, r) = pts[a.GetIdx()], (pts[n.GetIdx()] for n in a.GetNeighbors())
        v, w = (q[0] - p[0], q[1] - p[1]), (r[0] - p[0], r[1] - p[1])
        cos = (v[0] * w[0] + v[1] * w[1]) / max(math.hypot(*v) * math.hypot(*w), 1e-12)
        if 180.0 - math.degrees(math.acos(max(-1.0, min(1.0, cos)))) > LINEAR_TOL:
            out.append(a.GetIdx())
    return out


def prepare(mol: Chem.Mol, name: str, bond_length: float, abbreviate: list[str] | None = None,
            from_file: bool = False, rotate_deg: float | str | None = None, flip_h: bool = False,
            flip_v: bool = False, refs: list[Prepared] | None = None) -> Prepared:
    """Condense abbreviations, then lay out: aligned to the best of `refs` (most matched atoms)
    when one shares at least 3 atoms, else on its own (SMILES: CoordGen + orientation, default
    auto; molfile: its own layout, oriented as asked). An alignment that bends an sp atom
    (nitronium pinned to bent HNO3 atoms) is dropped for the species' own layout."""
    mol = condense(mol, abbreviate, name)
    best = None
    for ref in refs or []:
        pairs = correspondence(ref.mol, mol)
        if pairs and (best is None or len(pairs) > len(best[1])):
            best = (ref, pairs)
    if best:
        ref, pairs = best
        pairs = stable_core(ref.mol, mol, pairs)  # `matched` holds the pinned atoms
        align_coords(mol, ref.mol, pairs)
        if not bent_linear_atoms(mol):
            return Prepared(mol, name, bond_length, ref.scale, ref.name, pairs)
    if from_file:
        _orient_molfile(mol, rotate_deg or 0.0, flip_h, flip_v)
    else:
        layout_2d(mol)
        orient(mol, "auto" if rotate_deg is None else rotate_deg, flip_h, flip_v)
    return Prepared(mol, name, bond_length, layout_scale(mol, bond_length))


def translate(m: Molecule, d: Point) -> None:
    for a in m.atoms:
        a.x, a.y = add(a.p, d)


def centroid(m: Molecule) -> Point:
    return (statistics.fmean(a.x for a in m.atoms), statistics.fmean(a.y for a in m.atoms))
