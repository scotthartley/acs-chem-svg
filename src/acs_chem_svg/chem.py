"""Chemistry via RDKit: parse, sanitize, 2D coordinates, rings, stereo (design.md §4.1).

The output is a plain `Molecule` in pt, y down, so drawing never touches RDKit.

- A SMILES carbon with one H written as an atom ([H:6][CH2:5]...) gets all its H's as atoms
  (unmapped), so a drawn C–H never hides the others. Molfiles keep their own H atoms.
- Kekulé bonds are kept as written. Aromatic SMILES are kekulized by RDKit, or
  per `kekule` (pairs of atom-map numbers that are double bonds).
- SMILES get CoordGen 2D coordinates. A molfile keeps its layout. A four-coordinate d-block
  metal with four separate ligands (Br[Pd](Ph)(PMe3)PMe3, cisplatin) is drawn as a cross, its
  ligands laid out one by one; @SP1-3 in the SMILES choose cis or trans (_square_planar). A
  two-coordinate one with two separate ligands (Pd(PMe3)2, L-Au-Cl) is drawn straight. Either way the median
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
import zlib
import re
import statistics
from dataclasses import dataclass, field, replace
from pathlib import Path

from rdkit import Chem, rdBase, Geometry
from rdkit.Chem import rdAbbreviations, rdCoordGen, rdDepictor, rdFMCS

from . import abbrev
from .geometry import Point, add, rotate, segments_cross

_SANITIZE = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
# a dative bond (P->[Pd]) is drawn as a single bond; bond_order_sum gives its pair to the donor
_ORDER = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2, Chem.BondType.TRIPLE: 3, Chem.BondType.DATIVE: 1}
_MOLFILE_STEREO = {1: "wedge", 6: "hash"}
_CARBON_ALIAS = re.compile(r"C(H\d?)?|H\d?C")  # CH3, H3C, CH2, CH, C
STANDARD_BOND = 1.5  # A; molfile bond length that maps exactly onto the preset bond length
_RDKIT_DIR = {Chem.BondDir.BEGINWEDGE: "wedge", Chem.BondDir.BEGINDASH: "hash"}
# Abbreviations come from the registry (abbrev.py, data/abbreviations.yaml): label -> group
# SMILES (* = attachment) for groups and formulas, and the right-justified text by the reversal
# rule (OMe -> MeO, CO2Me -> MeO2C, TMS unchanged; acs-style.md §3.5). Views of abbrev.default().
ABBREVIATIONS = {k: e.smiles for k, e in abbrev.default().entries.items() if e.smiles}
ABBREVIATIONS_REVERSED = {k: abbrev.reverse(k) for k in ABBREVIATIONS}
ABBREV_PROP = "_abbrevInfo"  # a collapsed dummy atom's abbrev.LabelInfo (JSON)
# composition (§3.5): a heteroatom merges with its collapsed groups when it keeps this many other bonds
COMPOSE_MAX_SKELETON = 2
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
_D_BLOCK = set(range(21, 31)) | set(range(39, 49)) | {57} | set(range(72, 81))  # Sc-Zn, Y-Cd, La, Hf-Hg
_SP_ORDER = {1: (0, 1, 2, 3), 2: (0, 2, 1, 3), 3: (0, 1, 3, 2)}  # @SP1-3: bond order round the square
_ARM_STEP = {4: 90.0, 2: 180.0}  # _square_planar: degrees between a d-block metal's arms, by its degree (user decision 2026-10-30 for 2)
CROSS_CROWDED_PER_BOND = 0.6  # _square_planar: atoms of two ligands nearer than this x bond length crowd (a choice, as chair.CROWDED_PER_BOND)
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
    abbrev: str | None = None  # abbreviation label key (OMe, NTs); `element` is its attachment atom
    abbrev_info: abbrev.LabelInfo | None = None  # its texts, anchors and CDXML data (None: from the registry by key)
    radicals: int = 0  # unpaired electrons (RDKit GetNumRadicalElectrons)
    lone_pairs: int = 0  # (valence electrons - charge - bond orders - H - radicals) / 2
    singlet: bool = False  # a carbene / nitrene: RDKit's two radicals are this lone pair (nonbonding()); always drawn
    charge_label: str | None = None  # drawn charge text in place of `charge` (spec charge_labels: +3); display only
    z: float = 0.0  # depth, larger = nearer the viewer (chair.py); the nearer of two crossing bonds is drawn whole

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
    lewis_lines: bool = False  # Lewis structure with bond lines instead of shared-pair dots
    circled_charges: bool = False  # ±1 signs drawn in a ring (spec charges.circled / circled_charges)
    front: list[int] = field(default_factory=list)  # bond indices drawn whole where they cross a bond at equal depth (spec front)

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


def _arm_metal(mol: Chem.Mol, i: int) -> bool:
    """A d-block atom _square_planar lays out: four (a cross) or two (a line) separate ligands."""
    a = mol.GetAtomWithIdx(i)
    return a.GetAtomicNum() in _D_BLOCK and a.GetDegree() in _ARM_STEP and _cross_ligands(mol, i) is not None


def _cross_angles(mol: Chem.Mol, pts: list[Point]) -> list[float]:
    """Directions (deg) of the bonds of metals drawn on arms (_square_planar): a cross, bonds
    90° apart, or a line, 180° apart."""
    out = []
    for a in mol.GetAtoms():
        if not _arm_metal(mol, a.GetIdx()):
            continue
        c = pts[a.GetIdx()]
        dirs = [math.degrees(math.atan2(-(pts[n.GetIdx()][1] - c[1]), pts[n.GetIdx()][0] - c[0])) for n in a.GetNeighbors()]
        if all(_off(d - dirs[0], _ARM_STEP[a.GetDegree()]) < 1.0 for d in dirs):
            out += dirs
    return out


def _off(angle: float, step: float, phase: float = 0.0) -> float:
    """Distance (deg) from `angle` to the nearest phase + n*step."""
    return abs(((angle - phase + step / 2) % step) - step / 2)


def auto_rotation(mol: Chem.Mol, pts: list[Point], near: float | None = None) -> float:
    """Rotation (deg) in ChemDraw's habit: zigzag bonds at ±30°, rings with vertical sides
    (the most bonds at 30° + n·60°), then the most bonds on the 30° grid, then a layout
    wider than tall, then the smallest turn. A metal drawn as a cross or a line keeps its arms
    on the axes (only quarter turns are tried). All-collinear molecules lie horizontal.
    `near`: only rotations within 30° of it (orient_to's fitted angle; 45° for a cross or a
    line, so a quarter turn is always in reach), the closest last."""
    angles = _bond_angles(mol, pts)
    if not angles:
        return 0.0 if near is None else near
    cross = _cross_angles(mol, pts)  # a metal cross sits on the axes: only turns that keep it there
    turn = (lambda r: _off(r, 360)) if near is None else (lambda r: _off(r - near, 360))
    if all(_off(a - angles[0], 180) < 0.5 for a in angles):
        return min((-angles[0], 180 - angles[0]), key=turn)
    best = None
    for a in cross or angles:
        for target in range(-180, 180, 90 if cross else 30):
            r = target - a
            if near is not None and turn(r) > (45 if cross else 30) + 1e-6:
                continue
            zig = sum(1 for b in angles if _off(b + r, 60, 30) < 0.5)
            grid = sum(1 for b in angles if _off(b + r, 30) < 0.5)
            rp = [rotate(p, r) for p in pts]
            w = max(p[0] for p in rp) - min(p[0] for p in rp)
            h = max(p[1] for p in rp) - min(p[1] for p in rp)
            key = (-zig, -grid, w < h - 1e-6, turn(r))
            if best is None or key < best[0]:
                best = (key, r)
    return best[1]


def _dative_as_single(mol: Chem.Mol) -> Chem.Mol:
    """mol itself, or a copy with its dative bonds made single, to lay out: CoordGen lays a
    dative bond out as a single bond anyway, but logs "unrecognized bond type" for each."""
    if not any(b.GetBondType() == Chem.BondType.DATIVE for b in mol.GetBonds()):
        return mol
    m = Chem.RWMol(mol)
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DATIVE:
            b.SetBondType(Chem.BondType.SINGLE)
    return m


def _take_conformer(mol: Chem.Mol, laid: Chem.Mol) -> None:
    if laid is not mol:
        mol.RemoveAllConformers()
        mol.AddConformer(Chem.Conformer(laid.GetConformer()), assignId=True)


def _cross_ligands(mol: Chem.Mol, i: int) -> dict[int, set[int]] | None:
    """Donor atom -> the atoms of its ligand, for a metal atom i whose four neighbors are on four
    separate ligands; None for a chelate (two donors on one ligand) or a ring through i."""
    out: dict[int, set[int]] = {}
    for nb in mol.GetAtomWithIdx(i).GetNeighbors():
        if any(nb.GetIdx() in comp for comp in out.values()):
            return None
        comp, todo = set(), [nb.GetIdx()]
        while todo:
            a = todo.pop()
            if a not in comp:
                comp.add(a)
                todo += [n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() != i]
        out[nb.GetIdx()] = comp
    return out


def _cross_order(mol: Chem.Mol, i: int, pts: list[Point]) -> list[int]:
    """The metal's four neighbors in order round the cross: the SMILES's square-planar stereo
    (@SP1 U, @SP2 4, @SP3 Z shapes, OpenSMILES), else counter-clockwise as CoordGen laid them."""
    a = mol.GetAtomWithIdx(i)
    bonded = [b.GetOtherAtomIdx(i) for b in a.GetBonds()]
    if a.GetChiralTag() == Chem.ChiralType.CHI_SQUAREPLANAR and a.HasProp("_chiralPermutation"):
        perm = _SP_ORDER.get(a.GetUnsignedProp("_chiralPermutation"))
        if perm:
            return [bonded[k] for k in perm]
    c = pts[i]
    return sorted(bonded, key=lambda j: math.atan2(-(pts[j][1] - c[1]), pts[j][0] - c[0]))


def _cross_crowding(placed: list[dict[int, Point]], step: float) -> float:
    """How far atoms of different ligands come inside CROSS_CROWDED_PER_BOND x `step` of each
    other, summed (screen coordinates, the layout's own units)."""
    lim = CROSS_CROWDED_PER_BOND * step
    return sum(max(0.0, lim - math.dist(p, q))
               for a, b in itertools.combinations(placed, 2) for p in a.values() for q in b.values())


def _square_planar(mol: Chem.Mol, laid: Chem.Mol) -> None:
    """Redraw each four-coordinate d-block metal with four separate ligands as a cross: each
    ligand laid out alone with CoordGen (as a substituent on the lone metal atom) and turned so
    its bond to the metal lies on its arm, the arms 90° apart, the first nearest CoordGen's own
    direction. Each ligand is then kept as laid out or mirrored about its arm, whichever crowds
    the others least (_cross_crowding; ties keep it). CoordGen draws such a metal unevenly
    (Br–Pd 0.88, Pd–P 1.13 bond lengths) and ignores a coordMap round it. `laid` is mol's layout
    copy (_dative_as_single), its conformer set. A two-coordinate metal with two separate
    ligands is laid out the same way on two arms 180° apart: a line (Pd(PMe3)2; CoordGen bends
    it when the ligands are abbreviations)."""
    conf = laid.GetConformer()
    for a in laid.GetAtoms():
        i = a.GetIdx()
        if not _arm_metal(laid, i):
            continue
        ligands = _cross_ligands(laid, i)
        step_deg = _ARM_STEP[a.GetDegree()]
        pts = _coords(laid)
        c = pts[i]
        order = _cross_order(mol, i, pts)
        d0 = pts[order[0]]
        start = round(math.atan2(-(d0[1] - c[1]), d0[0] - c[0]) / (math.pi / 2)) * math.pi / 2  # on an axis
        options, step = [], 0.0
        for k, donor in enumerate(order):
            keep = sorted(ligands[donor] | {i})
            sub = Chem.RWMol(laid)
            for j in sorted(set(range(laid.GetNumAtoms())) - set(keep), reverse=True):
                sub.RemoveAtom(j)
            rdDepictor.Compute2DCoords(sub)
            sp = _coords(sub)
            mi, di = keep.index(i), keep.index(donor)
            step = math.dist(sp[mi], sp[di])
            arm = math.degrees(start) + k * step_deg
            turn = arm - math.degrees(math.atan2(-(sp[di][1] - sp[mi][1]), sp[di][0] - sp[mi][0]))
            rel = {old: rotate((sp[new][0] - sp[mi][0], sp[new][1] - sp[mi][1]), turn)
                   for new, old in enumerate(keep) if old != i}
            # mirrored about the arm: rotate onto +x, flip y, rotate back
            flip = {j: rotate((lambda q: (q[0], -q[1]))(rotate(v, -arm)), arm) for j, v in rel.items()}
            options.append([{j: add(c, v) for j, v in o.items()} for o in (rel, flip)])
        best = min(itertools.product((0, 1), repeat=len(options)),
                   key=lambda pick: (round(_cross_crowding([o[f] for o, f in zip(options, pick)], step), 6), sum(pick)))
        for o, f in zip(options, best):
            for j, (x, y) in o[f].items():
                conf.SetAtomPosition(j, (x, -y, 0.0))


def layout_2d(mol: Chem.Mol) -> None:
    rdDepictor.SetPreferCoordGen(True)
    m = _dative_as_single(mol)
    if m is mol and any(a.GetAtomicNum() in _D_BLOCK for a in mol.GetAtoms()):
        m = Chem.Mol(mol)  # _square_planar reads mol's stereo after laying out its copy
    rdDepictor.Compute2DCoords(m)
    _square_planar(mol, m)
    _take_conformer(mol, m)


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


def to_molecule(mol: Chem.Mol, name: str, bond_length: float, scale: float | None = None,
                depth: dict[int, float] | None = None, flat: set[int] | None = None) -> Molecule:
    """`scale` (pt per layout unit) defaults to layout_scale(); aligned species pass their
    reference's, so shared atoms land on the same pt coordinates. `depth`: atom idx -> Atom.z
    (chairs). `flat`: atoms whose stereo the drawing already shows (a chair's ring): no wedge
    starts there."""
    k = layout_scale(mol, bond_length) if scale is None else scale
    pts = [(x * k, y * k) for x, y in _coords(mol)]
    stereo = {bi: st for bi, st in _stereo(mol).items() if st[1] not in (flat or ())}
    atoms = []
    for a in mol.GetAtoms():
        i = a.GetIdx()
        mp = a.GetAtomMapNum() or None
        show, key = _alias_label(a)
        info = label_info(a) if key else None
        element = info.element if info and a.GetAtomicNum() == 0 else a.GetSymbol()
        if info is not None and info.kind == "text":  # no chemistry: drawn on its stand-in, which keeps its element
            info = replace(info, element=element)
        lp, rad = nonbonding(a)
        if key:  # a label stands for its whole group: RDKit's radicals on a template's Sn (SnBu3,
            rad = 0  # one drawn bond) are not electrons to draw, nor a carbene's pair
        atoms.append(Atom(i, str(mp) if mp else str(i + 1), element, pts[i][0], pts[i][1],
                          a.GetFormalCharge(), 0 if key else a.GetTotalNumHs(includeNeighbors=False), mp,
                          show, key, info, rad, 0 if key else lp, singlet=not key and rad != a.GetNumRadicalElectrons(),
                          z=(depth or {}).get(i, 0.0)))
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
    """Electrons the atom puts into its bonds: the bond orders to every neighbor, plus one per
    implicit or counted H. A dative bond's pair is all the donor's (2) and none of the acceptor's,
    so P in P->[Pd] keeps no lone pair and Pd keeps its d electrons."""
    def share(b: Chem.Bond) -> int:
        if b.GetBondType() == Chem.BondType.DATIVE:
            return 2 if b.GetBeginAtomIdx() == a.GetIdx() else 0
        return _ORDER[b.GetBondType()]
    return sum(share(b) for b in a.GetBonds()) + a.GetTotalNumHs(includeNeighbors=False)


def lone_pairs(a: Chem.Atom) -> int:
    """Nonbonding electron pairs: (V - formal charge - bond orders - radicals) / 2, for a
    Kekulé mol (H has none; dummy atoms none)."""
    if a.GetAtomicNum() in (0, 1) and a.GetFormalCharge() >= 0:
        return 0
    n = valence_electrons(a) - a.GetFormalCharge() - bond_order_sum(a) - a.GetNumRadicalElectrons()
    return max(0, n // 2)


TRIPLET_PROP = "_triplet"  # set on an atom listed in its species' triplet: (spec._parsed)


def nonbonding(a: Chem.Atom) -> tuple[int, int]:
    """(lone pairs, unpaired electrons). RDKit reads [CH2], Cl[C]Cl and [N] as two radicals;
    unless the atom is a declared triplet (TRIPLET_PROP), they are a singlet carbene's or
    nitrene's lone pair."""
    rad = a.GetNumRadicalElectrons()
    if rad == 2 and not a.HasProp(TRIPLET_PROP):
        return lone_pairs(a) + 1, 0
    return lone_pairs(a), rad


def _alias_label(a: Chem.Atom, reg: abbrev.Registry | None = None) -> tuple[bool, str | None]:
    """(explicit carbon label, abbreviation key). An abbreviation comes from condense (atomLabel)
    or a molfile alias naming a registry group or formula (OMe on the O, as cdxml_to_molfile
    writes ChemDraw nicknames). A carbon alias that is its own label (CH3, H3C, ...) shows the
    carbon. Any other alias is a condensed label we don't draw."""
    if a.HasProp("atomLabel") and a.GetAtomicNum() == 0:
        return False, a.GetProp("atomLabel")
    if not a.HasProp("molFileAlias"):
        return False, None
    alias = a.GetProp("molFileAlias")
    reg = reg or abbrev.default()
    e = reg.get(alias)
    if e is None and _CARBON_ALIAS.fullmatch(alias) is None:
        # (a template's stand-in atom has implicit H for the atoms the label holds: not checked)
        info = composed_info(alias, a.GetSymbol(), a.GetDegree(), None, reg)
        if info is not None:  # a composed label read from a template (ring NTs, abbrev_2_x*)
            if not a.HasProp(ABBREV_PROP):
                a.SetProp(ABBREV_PROP, info.to_prop())
            return False, alias
    if e is not None and e.kind == "text":  # a text label on its stand-in atom (L, dppf: ligands_x*)
        if not a.HasProp(ABBREV_PROP):
            a.SetProp(ABBREV_PROP, abbrev.info(alias, reg).to_prop())
        return False, alias
    if e is not None and e.smiles:
        if abbrev_element(alias, reg) != a.GetSymbol():
            raise ValueError(f"atom {a.GetIdx() + 1}: alias {alias!r} attaches through "
                             f"{abbrev_element(alias, reg)}, not {a.GetSymbol()}")
        if not a.HasProp(ABBREV_PROP):
            a.SetProp(ABBREV_PROP, abbrev.info(alias, reg).to_prop())
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
                     f"draw; define it under the spec's `abbreviations:` (acs-style.md §3.5)")


def composed_info(key: str, element: str, degree: int, n_h: int | None = None,
                  reg: abbrev.Registry | None = None) -> abbrev.LabelInfo | None:
    """LabelInfo of a composed label read back from its key (a molfile alias on its heteroatom:
    NTs on a ring N): the element's token, an H count (checked against `n_h` when given), then
    registry groups with counts; None when the key doesn't read that way."""
    reg = reg or abbrev.default()
    try:
        toks = abbrev.tokenize(key, reg)
    except ValueError:
        return None
    if not toks or toks[0].element != element or toks[0].count != 1:
        return None
    rest = toks[1:]
    h = 0
    if rest and rest[0].element == "H":
        h, rest = rest[0].count, rest[1:]
    groups = []
    for t in rest:
        name = t.inner[0].key if t.inner and len(t.inner) == 1 and t.inner[0].count == 1 else t.key
        e = reg.get(name)
        if e is None or e.kind != "group":
            return None
        groups += [e.smiles] * t.count
    if not groups or (n_h is not None and h != n_h):
        return None
    info = abbrev.info(key, reg, element=element)
    return replace(info, kind="composed", smiles=_composed_smiles(element, h, groups, degree))


def label_info(a: Chem.Atom) -> abbrev.LabelInfo:
    """The LabelInfo condense stored on a collapsed atom, else the registry's for its key."""
    if a.HasProp(ABBREV_PROP):
        return abbrev.LabelInfo.from_prop(a.GetProp(ABBREV_PROP))
    key = a.GetProp("atomLabel") if a.HasProp("atomLabel") else a.GetProp("molFileAlias")
    return abbrev.info(key, element=a.GetSymbol() if a.GetAtomicNum() else None)


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


def abbrev_element(label: str, reg: abbrev.Registry | None = None) -> str:
    """Element of the atom an abbreviation attaches through (OMe: O, CO2Me: C, Ph: C)."""
    e = (reg or abbrev.default()).get(label)
    if e is None or not e.smiles:
        raise ValueError(f"{label!r} is not a group or formula abbreviation")
    return abbrev.attachment_element(e)


def condense(mol: Chem.Mol, labels: list[str] | None, name: str = "",
             reg: abbrev.Registry | None = None) -> Chem.Mol:
    """Collapse each named group or formula (`labels`, registry keys) to one dummy atom carrying
    atomLabel and its abbrev.LabelInfo, at the attachment atom (its index, coordinates and
    atom-map number are kept). Groups are found by RDKit rdAbbreviations on an aromatic copy (its
    patterns are aromatic); the Kekulé bonds drawn are those of `mol`. Every label must be found.
    rdAbbreviations misses a group attached by a dative bond (PMe3 on Pd), so the copy has its
    dative bonds made single; the drawn bond stays dative. Then compose() merges heteroatoms with
    their collapsed groups (OTBS, NHBoc, ring NTs)."""
    reg = reg or abbrev.default()
    for a in mol.GetAtoms():  # molfile aliases naming a spec's own entries carry its texts
        if a.HasProp("molFileAlias") and reg.get(a.GetProp("molFileAlias")) is not None:
            a.SetProp(ABBREV_PROP, abbrev.info(a.GetProp("molFileAlias"), reg).to_prop())
    if not labels:
        return mol
    text = [lab for lab in labels if reg.get(lab) is not None and reg.get(lab).kind == "text"]
    if text:
        raise ValueError(f"{name}: abbreviate {', '.join(text)}: a text label has no structure to collapse; "
                         "draw it on a stand-in atom with `labels:` (acs-style.md §3.5)")
    defs = rdAbbreviations.ParseAbbreviations("\n".join(f"{lab} {reg.get(lab).smiles}" for lab in labels))
    ar = Chem.Mol(_dative_as_single(mol))
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
        dummy.SetProp(ABBREV_PROP, abbrev.info(sg.GetProp("LABEL"), reg).to_prop())
        dummy.SetAtomMapNum(mol.GetAtomWithIdx(attach).GetAtomMapNum())  # aligned species keep the label in place
        rw.ReplaceAtom(attach, dummy)
        drop += [i for i in group if i != attach]
    for i in sorted(drop, reverse=True):
        rw.RemoveAtom(i)
    out = rw.GetMol()
    Chem.SanitizeMol(out, _SANITIZE)
    return compose(out, reg)


def _composable(a: Chem.Atom) -> bool:
    """A heteroatom that can take its groups into its label: not C or H, not a d-block metal
    (Pd keeps its Ph and PPh3 bonds: metals_x*), uncharged, no unpaired electrons."""
    return (a.GetAtomicNum() not in (0, 1, 6) and a.GetAtomicNum() not in _D_BLOCK
            and not a.GetFormalCharge() and not a.GetNumRadicalElectrons())


def compose(mol: Chem.Mol, reg: abbrev.Registry | None = None) -> Chem.Mol:
    """§3.5 composition: a heteroatom (_composable) bonded to collapsed groups (registry kind
    group) and to one or two other atoms (the skeleton) becomes one labeled atom: element, H count,
    then the groups (abbrev.compose: NHBoc, OTBS, SnBu3, P(t-Bu)3, ring NTs). Its LabelInfo holds
    the merged atoms as SMILES with [*:k] for its k-th skeleton bond, by neighbor index."""
    reg = reg or abbrev.default()

    def group_key(n: Chem.Atom) -> str | None:
        if n.GetAtomicNum() or not n.HasProp("atomLabel") or not n.HasProp(ABBREV_PROP):
            return None
        info = abbrev.LabelInfo.from_prop(n.GetProp(ABBREV_PROP))
        return info.key if info.kind == "group" else None

    rw = Chem.RWMol(mol)
    drop = []
    for a in mol.GetAtoms():
        if not _composable(a):
            continue
        nbs = sorted(a.GetNeighbors(), key=lambda n: n.GetIdx())
        groups = [n for n in nbs if group_key(n)]
        skeleton = [n for n in nbs if not group_key(n)]
        if not groups or not 1 <= len(skeleton) <= COMPOSE_MAX_SKELETON:
            continue
        if any(mol.GetBondBetweenAtoms(a.GetIdx(), g.GetIdx()).GetBondType() != Chem.BondType.SINGLE for g in groups):
            continue
        n_h = a.GetTotalNumHs(includeNeighbors=False)
        keys = [group_key(g) for g in groups]
        key = abbrev.compose(a.GetSymbol(), n_h, keys, reg)
        smiles = _composed_smiles(a.GetSymbol(), n_h, [reg.get(k).smiles for k in keys], len(skeleton))
        info = replace(abbrev.info(key, reg, element=a.GetSymbol()), kind="composed", smiles=smiles)
        dummy = Chem.Atom(0)
        dummy.SetNoImplicit(True)
        dummy.SetProp("atomLabel", key)
        dummy.SetProp(ABBREV_PROP, info.to_prop())
        dummy.SetAtomMapNum(a.GetAtomMapNum())
        rw.ReplaceAtom(a.GetIdx(), dummy)
        drop += [g.GetIdx() for g in groups]
    if not drop:
        return mol
    for i in sorted(drop, reverse=True):
        rw.RemoveAtom(i)
    out = rw.GetMol()
    Chem.SanitizeMol(out, _SANITIZE)
    return out


def _composed_smiles(element: str, n_h: int, groups: list[str], n_outer: int) -> str:
    """The merged atoms of a composed label: the heteroatom with its H, each group's SMILES
    bonded by its `*`, and [*:1]..[*:n] for the skeleton bonds."""
    rw = Chem.RWMol()
    x = rw.AddAtom(Chem.Atom(element))
    rw.GetAtomWithIdx(x).SetNumExplicitHs(n_h)
    rw.GetAtomWithIdx(x).SetNoImplicit(True)
    for k in range(1, n_outer + 1):
        d = rw.AddAtom(Chem.Atom(0))
        rw.GetAtomWithIdx(d).SetAtomMapNum(k)
        rw.AddBond(x, d, Chem.BondType.SINGLE)
    for g in groups:
        gm = Chem.MolFromSmiles(g)
        off = rw.GetNumAtoms()
        star = next(a.GetIdx() for a in gm.GetAtoms() if a.GetAtomicNum() == 0)
        nb = gm.GetAtomWithIdx(star).GetNeighbors()[0].GetIdx()
        combo = Chem.RWMol(Chem.CombineMols(rw, gm))
        combo.AddBond(x, off + nb, Chem.BondType.SINGLE)
        combo.RemoveAtom(off + star)
        rw = combo
    m = rw.GetMol()
    Chem.SanitizeMol(m)
    return Chem.MolToSmiles(m)


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
                abbreviate: list[str] | None = None, lewis: str | None = None, lewis_lines: bool = False,
                chair: ChairSpec | None = None, registry: abbrev.Registry | None = None) -> Molecule:
    """`lewis` (a LEWIS_MODES value): a Lewis structure with every H an atom and every atom
    labeled, laid out by lewis_layout with `bond_length` as the lattice step (lewis.bond_length);
    `lewis_lines`: its bonds drawn as lines (shared only). `chair`: drawn as a chair (chair.py)."""
    if chair is not None:
        p = prepare(parse_smiles(smiles, kekule), name, bond_length, abbreviate, chair=chair, registry=registry)
        return p.molecule()
    if lewis is not None:
        mol = Chem.AddHs(parse_smiles(smiles))
        Chem.RemoveStereochemistry(mol)
        lewis_layout(mol, bond_length, name)
        orient(mol, 0.0 if rotate_deg in (None, "auto") else rotate_deg, flip_h, flip_v)  # a multiple of 90
        _set_coords(mol, [(round(x, 6) + 0.0, round(y, 6) + 0.0) for x, y in _coords(mol)])
        m = to_molecule(mol, name, bond_length, scale=1.0)
        for a in m.atoms:
            a.show_label = True
        m.lewis, m.lewis_lines = lewis, lewis_lines
        if lewis == "oxidation_state":  # each nonzero oxidation state drawn as a charge text (display only)
            for a, ox in zip(m.atoms, oxidation_states(m)):
                if ox:
                    a.charge_label = f"{abs(ox) if abs(ox) > 1 else ''}{'+' if ox > 0 else '-'}"
        return m
    mol = condense(parse_smiles(smiles, kekule), abbreviate, name, registry)
    layout_2d(mol)
    orient(mol, rotate_deg, flip_h, flip_v)
    return to_molecule(mol, name, bond_length)


@dataclass
class ChairSpec:
    """A chair species (chair.layout): the preset's templates and the spec's options (atom maps)."""
    templates: object  # preset.Chair
    ring: list[int] | None = None
    flip: bool | None = None
    axial: list[int] = field(default_factory=list)
    equatorial: list[int] = field(default_factory=list)
    start: int | None = None


def from_molfile(path_or_text: str | Path, name: str, bond_length: float, kekule: list[str] | None = None,
                 rotate_deg: float | str | None = None, flip_h: bool = False, flip_v: bool = False,
                 abbreviate: list[str] | None = None, registry: abbrev.Registry | None = None) -> Molecule:
    text = Path(path_or_text).read_text() if isinstance(path_or_text, Path) else path_or_text
    mol = condense(parse_molfile(text, kekule), abbreviate, name, registry)
    _orient_molfile(mol, rotate_deg, flip_h, flip_v)
    return to_molecule(mol, name, bond_length)


# --------------------------------------------------------------------------- Lewis layout

_RIGHT, _LEFT, _UP, _DOWN = (1, 0), (-1, 0), (0, -1), (0, 1)  # lattice steps, y down


def _atom_name(a: Chem.Atom) -> str:
    return f"{a.GetSymbol()} (map {a.GetAtomMapNum()})" if a.GetAtomMapNum() else \
        f"{a.GetSymbol()} (atom {a.GetIdx() + 1} in the SMILES)"


def oxidation_states(m: Molecule) -> list[int]:
    """Each atom's oxidation state: its formal charge plus, per bond, +order when the other atom
    is more electronegative (ELECTRONEGATIVITY), -order when less, 0 when equal; the bookkeeping
    of Lewis oxidation_state (draw_mol.Drawer.bond_dots: the more electronegative atom gets the
    bond's electrons). User decision 2026-10-28: the oxidation_state view draws them."""
    def en(a):
        if a.element not in ELECTRONEGATIVITY:
            raise ValueError(f"lewis oxidation_state: no electronegativity for {a.element} (atom {a.id}); "
                             f"known: {', '.join(ELECTRONEGATIVITY)}")
        return ELECTRONEGATIVITY[a.element]
    out = [a.charge for a in m.atoms]
    for b in m.bonds:
        ea, eb = en(m.atoms[b.a]), en(m.atoms[b.b])
        if ea != eb:
            hi, lo = (b.a, b.b) if ea > eb else (b.b, b.a)
            out[hi] -= b.order
            out[lo] += b.order
    return out


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
    depth: dict[int, float] | None = None  # chairs: atom idx -> Atom.z
    chair_ring: list[int] | None = None  # chairs: the ring, whose stereo the geometry shows

    def molecule(self) -> Molecule:
        return to_molecule(self.mol, self.name, self.bond_length, self.scale, self.depth,
                           set(self.chair_ring or ()))


def _mcs_copy(mol: Chem.Mol) -> Chem.Mol:
    """Copy for MCS: aromaticity perceived (so Kekulé forms don't matter) and each atom's
    isotope set to its element, or to a per-label code for abbreviations (OMe != OBn)."""
    c = Chem.Mol(mol)
    for a in c.GetAtoms():
        if a.GetAtomicNum():
            a.SetIsotope(a.GetAtomicNum())
        else:  # a stable per-label code below the map codes (1000 + map)
            lab = a.GetProp("atomLabel") if a.HasProp("atomLabel") else ""
            a.SetIsotope(200 + zlib.crc32(lab.encode()) % 700 if lab else 199)
    Chem.FastFindRings(c)
    Chem.SetAromaticity(c)
    return c


def correspondence(ref: Chem.Mol, mol: Chem.Mol) -> list[tuple[int, int]]:
    """(reference atom idx, atom idx) pairs. With at least two shared atom-map numbers: a maximum
    common substructure seeded by them (each shared map its own atom type, other atoms by element
    or abbreviation, any bond order, rings only to rings, so cyclohexene's C=C still matches
    norcarane's C–C and its unmapped CH2s count), or the map pairs alone if that MCS leaves a map
    out. Else the maximum common substructure (same elements and abbreviations, same bond orders
    with aromatic rings matched as aromatic, complete rings only). Fewer than 3 atoms: none."""
    maps_r = {a.GetAtomMapNum(): a.GetIdx() for a in ref.GetAtoms() if a.GetAtomMapNum()}
    maps_m = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}
    common = sorted(maps_r.keys() & maps_m.keys())
    r, m = _mcs_copy(ref), _mcs_copy(mol)
    if len(common) >= 2:
        by_map = [(maps_r[k], maps_m[k]) for k in common]
        for k in common:
            r.GetAtomWithIdx(maps_r[k]).SetIsotope(1000 + k)
            m.GetAtomWithIdx(maps_m[k]).SetIsotope(1000 + k)
        res = rdFMCS.FindMCS([r, m], atomCompare=rdFMCS.AtomCompare.CompareIsotopes,
                             bondCompare=rdFMCS.BondCompare.CompareAny, ringMatchesRingOnly=True,
                             timeout=MCS_TIMEOUT)
        if res.numAtoms <= len(common):
            return by_map
        patt = Chem.MolFromSmarts(res.smartsString)
        pairs = list(zip(r.GetSubstructMatch(patt), m.GetSubstructMatch(patt)))
        return pairs if set(by_map) <= set(pairs) else by_map
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


class AlignmentError(RuntimeError):
    """align_coords could not hold the matched atoms at the reference's coordinates."""


def align_coords(mol: Chem.Mol, ref: Chem.Mol, pairs: list[tuple[int, int]]) -> None:
    """Lay out `mol` with the matched atoms fixed at the reference's coordinates. CoordGen
    (coordMap) first; it holds them only to ~0.005 A, so if that misses ALIGN_TOL use RDKit's
    GenerateDepictionMatching2DStructure, which holds them exactly. RDKit lays new atoms out at
    1.5 A bonds, so the reference is scaled to 1.5 A for it and the result scaled back."""
    conf = ref.GetConformer()
    params = rdCoordGen.CoordGenParams()
    params.SetCoordMap({m: Geometry.Point2D(conf.GetAtomPosition(r).x, conf.GetAtomPosition(r).y)
                        for r, m in pairs})
    m = _dative_as_single(mol)
    rdCoordGen.AddCoords(m, params)
    _take_conformer(mol, m)
    if _max_deviation(mol, ref, pairs) <= ALIGN_TOL:
        return
    pts = _coords(ref)
    med = statistics.median(math.dist(pts[b.GetBeginAtomIdx()], pts[b.GetEndAtomIdx()]) for b in ref.GetBonds())
    tmp = Chem.Mol(_dative_as_single(ref))
    _scale_conformer(tmp, STANDARD_BOND / med)
    rdDepictor.GenerateDepictionMatching2DStructure(m, tmp, list(pairs), -1, True)
    _take_conformer(mol, m)
    _scale_conformer(mol, med / STANDARD_BOND)
    dev = _max_deviation(mol, ref, pairs)
    if dev > ALIGN_TOL:
        raise AlignmentError(f"alignment left matched atoms {dev:.4f} A from their reference positions")


def stable_core(ref: Chem.Mol, mol: Chem.Mol, pairs: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """The matched atoms worth pinning: those whose neighbor count is unchanged, and of them only
    the largest connected group. An atom that changes its neighbors (the sp3 C of a sigma
    complex becoming sp2 in nitrobenzene) and what hangs off it are laid out afresh, so their
    bonds get proper angles. None if that leaves fewer than 3 atoms (prepare then orients
    the species' own layout instead of pinning it)."""
    keep = {m: r for r, m in pairs if ref.GetAtomWithIdx(r).GetDegree() == mol.GetAtomWithIdx(m).GetDegree()}
    if len(keep) == len(pairs):
        return pairs if len(pairs) >= 3 else []
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
        return []
    return [(r, m) for r, m in pairs if m in best]


LINEAR_TOL = 10.0  # degrees an sp atom may deviate from 180° after alignment (a choice, not a style value)
BOND_TOL = 0.15  # fraction a pinned layout's bond may differ from its median (a choice, not a style value)
ANGLE_TOL = 15.0  # degrees a pinned acyclic 2- or 3-bond atom's narrowest gap may fall below 120° (a choice)
OVERLAID_PER_BOND = 0.25  # non-bonded atoms nearer than this x the median bond overlap (a choice; geometry_checks too)


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


def off_arm_metals(mol: Chem.Mol) -> list[int]:
    """Metals _square_planar draws on arms (a cross or a line) whose layout puts a bond more
    than LINEAR_TOL off them (a Pd(II) pinned to Pd(0)'s P-Pd-P, its new ligands placed by
    CoordGen), or whose order round the cross contradicts the SMILES's @SP (cis for trans)."""
    pts = _coords(mol)
    out = []
    for a in mol.GetAtoms():
        i = a.GetIdx()
        if not _arm_metal(mol, i):
            continue
        c = pts[i]
        dirs = [math.degrees(math.atan2(-(pts[n.GetIdx()][1] - c[1]), pts[n.GetIdx()][0] - c[0])) for n in a.GetNeighbors()]
        if any(_off(d - dirs[0], _ARM_STEP[a.GetDegree()]) > LINEAR_TOL for d in dirs):
            out.append(i)
            continue
        drawn = sorted((n.GetIdx() for n in a.GetNeighbors()),
                       key=lambda j: math.atan2(-(pts[j][1] - c[1]), pts[j][0] - c[0]))
        want = _cross_order(mol, i, pts)
        turns = [want[k:] + want[:k] for k in range(len(want))]
        if drawn not in turns + [t[::-1] for t in turns]:  # either way round (a mirror image)
            out.append(i)
    return out


def uneven_bonds(mol: Chem.Mol) -> list[int]:
    """Bonds longer or shorter than the median by more than BOND_TOL (a ring closed across
    pinned atoms, as in perkin's tetrahedral intermediate)."""
    pts = _coords(mol)
    lengths = [math.dist(pts[b.GetBeginAtomIdx()], pts[b.GetEndAtomIdx()]) for b in mol.GetBonds()]
    if not lengths:
        return []
    med = statistics.median(lengths)
    return [i for i, d in enumerate(lengths) if abs(d / med - 1) > BOND_TOL]


def crowded_atoms(mol: Chem.Mol) -> list[int]:
    """Acyclic atoms with 2 or 3 neighbors whose narrowest angle between bonds is under
    120° - ANGLE_TOL (a trigonal anion pinned to a tetrahedral drawing's 60° gap)."""
    pts = _coords(mol)
    ring = mol.GetRingInfo()
    out = []
    for a in mol.GetAtoms():
        if a.GetDegree() not in (2, 3) or ring.NumAtomRings(a.GetIdx()):
            continue
        p = pts[a.GetIdx()]
        dirs = sorted(math.degrees(math.atan2(-(pts[n.GetIdx()][1] - p[1]), pts[n.GetIdx()][0] - p[0])) % 360
                      for n in a.GetNeighbors())
        gaps = [b - a_ for a_, b in zip(dirs, dirs[1:])] + [360 - dirs[-1] + dirs[0]]
        if min(gaps) < 120 - ANGLE_TOL:
            out.append(a.GetIdx())
    return out


def overlaid_atoms(mol: Chem.Mol) -> list[int]:
    """Atoms of a layout that lays one group over another: non-bonded atoms nearer than
    OVERLAID_PER_BOND x the median bond, or the ends of two bonds with no shared atom that
    cross (a tolyl pinned across a PMe3 on the core shared with the reactant). The layout-time
    twin of geometry_checks.check_overlaid."""
    pts = _coords(mol)
    ends = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    if not ends:
        return []
    lim = OVERLAID_PER_BOND * statistics.median(math.dist(pts[i], pts[j]) for i, j in ends)
    bonded = {frozenset(e) for e in ends}
    out = set()
    for i, j in itertools.combinations(range(mol.GetNumAtoms()), 2):
        if frozenset((i, j)) not in bonded and math.dist(pts[i], pts[j]) < lim:
            out |= {i, j}
    for e1, e2 in itertools.combinations(ends, 2):
        if not set(e1) & set(e2) and segments_cross(*(pts[i] for i in e1 + e2)):
            out |= set(e1 + e2)
    return sorted(out)


def _fit_angle(pts: list[Point], ref_pts: list[Point], pairs) -> tuple[float, float]:
    """(rotation in deg, residual): the rigid rotation that best carries the matched atoms of
    `pts` onto `ref_pts` (2D Kabsch on centered points; screen frame, as geometry.rotate)."""
    p = [pts[m] for _, m in pairs]
    q = [ref_pts[r] for r, _ in pairs]
    cp = (statistics.fmean(x for x, _ in p), statistics.fmean(y for _, y in p))
    cq = (statistics.fmean(x for x, _ in q), statistics.fmean(y for _, y in q))
    p = [(x - cp[0], -(y - cp[1])) for x, y in p]  # y up, so the angle is counter-clockwise
    q = [(x - cq[0], -(y - cq[1])) for x, y in q]
    dot = sum(a[0] * b[0] + a[1] * b[1] for a, b in zip(p, q))
    cross = sum(a[0] * b[1] - a[1] * b[0] for a, b in zip(p, q))
    t = math.degrees(math.atan2(cross, dot))
    k = math.cos(math.radians(t))
    s_ = math.sin(math.radians(t))
    res = sum((a[0] * k - a[1] * s_ - b[0]) ** 2 + (a[0] * s_ + a[1] * k - b[1]) ** 2 for a, b in zip(p, q))
    return t, res


def _collinear(pts: list[Point]) -> bool:
    if len(pts) < 3:
        return True
    a = pts[0]
    b = max(pts, key=lambda p: math.dist(p, a))
    d = math.dist(a, b)
    return d < 1e-9 or all(abs((b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])) / d < 0.05 * d
                           for p in pts)


def orient_to(mol: Chem.Mol, ref: Chem.Mol, pairs: list[tuple[int, int]], relayout: bool = True) -> None:
    """Orient mode: the species' own layout (CoordGen, unless `relayout` is off: a molfile keeps
    its coordinates), turned (and mirrored, when at least 3 non-collinear atoms are matched and
    `relayout` is on) to best fit the matched atoms onto the reference's, then snapped to the
    best auto_rotation within 30° of that fit, so bonds sit on ChemDraw's grid and still point
    the way the reactant's did. Nothing is pinned."""
    if relayout:
        layout_2d(mol)
    pts, ref_pts = _coords(mol), _coords(ref)
    fits = [(pts, *_fit_angle(pts, ref_pts, pairs))]
    if relayout and not _collinear([ref_pts[r] for r, _ in pairs]):
        mirrored = [(-x, y) for x, y in pts]
        fits.append((mirrored, *_fit_angle(mirrored, ref_pts, pairs)))
    pts, t, _ = min(fits, key=lambda f: f[2])
    r = auto_rotation(mol, pts, near=t)
    _set_coords(mol, [rotate(p, r) for p in pts])


def prepare(mol: Chem.Mol, name: str, bond_length: float, abbreviate: list[str] | None = None,
            from_file: bool = False, rotate_deg: float | str | None = None, flip_h: bool = False,
            flip_v: bool = False, refs: list[Prepared] | None = None, chair: ChairSpec | None = None,
            registry: abbrev.Registry | None = None) -> Prepared:
    """Condense abbreviations, then lay out against the best of `refs`: the one whose stable
    core (stable_core) is largest, then with the most matched atoms. Pin mode: with a core of at
    least 3 atoms, those atoms keep the reference's coordinates (align_coords); with no core, all
    matched atoms (≥ 3) do, which suits a trigonal center turning tetrahedral (acetone to the
    cyanohydrin alkoxide). Either is rejected if it bends an sp atom (nitronium on bent HNO3
    atoms), stretches a bond (a ring closed across pinned atoms), crowds an acyclic atom's
    bonds (a trigonal anion on chloroform's tetrahedron), lays one group over another
    (overlaid_atoms: Suzuki's tolyl across a PMe3) or takes a metal off its cross or line
    (off_arm_metals), or can't be held at all (AlignmentError). Otherwise orient mode (orient_to): the species' own layout, turned the way the
    matched atoms point in the reference; nothing pinned. No reference: its own layout (SMILES:
    CoordGen + orientation, default auto; molfile: its own layout, oriented as asked). A `chair`
    species is laid out by chair.layout and aligns to nothing; nothing aligns to a chair either
    (its pseudo-3D coordinates would bend a flat drawing)."""
    mol = condense(mol, abbreviate, name, registry)
    if chair is not None:
        from . import chair as chair_mod
        m, z, pl = chair_mod.layout(mol, chair.templates, bond_length, chair.ring, chair.flip,
                                    chair.axial, chair.equatorial, chair.start)
        return Prepared(m, name, bond_length, bond_length / STANDARD_BOND, depth=z, chair_ring=pl.ring)
    best = None
    for ref in refs or []:
        if ref.chair_ring is not None:
            continue
        pairs = correspondence(ref.mol, mol)
        if not pairs:
            continue
        core = stable_core(ref.mol, mol, pairs)
        key = (len(core) >= 3, len(core), len(pairs))
        if best is None or key > best[0]:
            best = (key, ref, pairs, core)
    if best:
        _, ref, pairs, core = best
        for pins in [core] + ([pairs] if not core and len(pairs) >= 3 else []):
            if not pins:
                continue
            pinned = Chem.Mol(mol)
            try:
                align_coords(pinned, ref.mol, pins)
            except AlignmentError:  # RDKit can't hold them either (a ligand's P pinned as a cross)
                continue
            if not (bent_linear_atoms(pinned) or uneven_bonds(pinned) or crowded_atoms(pinned)
                    or overlaid_atoms(pinned) or off_arm_metals(pinned)):
                return Prepared(pinned, name, bond_length, ref.scale, ref.name, pins)  # `matched`: the pinned atoms
        orient_to(mol, ref.mol, pairs, relayout=not from_file)
        return Prepared(mol, name, bond_length, layout_scale(mol, bond_length), ref.name)
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
