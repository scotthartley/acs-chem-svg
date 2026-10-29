"""Chair conformers of single six-membered rings (acs-style.md §14.1).

A chair species is drawn on ChemDraw's chair template (preset `chair.vertices`) or its exact
mirror image (the flipped chair), with each ring substituent axial (vertical)
or equatorial (along the template's measured direction). The rest of the molecule is laid out by
RDKit with the ring and the first substituent atoms pinned.

- Faces: each substituent's side of the ring comes from its atom's chiral tag (faces): the side
  from which the ring, in `ring` order, runs counterclockwise is +1. A substituent on an atom
  without a tag has no face; it goes equatorial (of a geminal pair without tags, the bulkier).
- Conformer: of the two ring flips, the one with the smaller sum of A-values over its axial
  substituents (A_VALUES), unless the spec asks for the other (`flip`) or names axial or
  equatorial atoms.
- Placement: the ring goes onto a template in one of 24 ways (2 templates x 6 rotations x 2
  directions); half of them give each flip. Of those, the fewest crowded atom pairs, then a ring
  heteroatom on the sugar vertex (preset chair.heteroatom_vertex of the flipped template,
  β-D-glucopyranose's O; carbocycles: the unflipped template), then the fewest crossings, then
  the first substituted ring atom farthest right.
- Depth (Atom.z): each ring atom's screen y (lower on the page = nearer: chair_x*'s Z order puts
  an axial bond on a lower ring atom in front of the upper ring bonds it crosses, and the lower
  ring bonds in front of an axial bond from an upper atom); every other atom takes the depth of
  the ring atom it hangs off.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Geometry import Point2D

from .geometry import Point
from .preset import Chair

# Conformational free energies (kcal/mol) for an axial substituent on cyclohexane: chemistry
# data, not style numbers. Representative values from E. L. Eliel, S. H. Wilen, Stereochemistry of
# Organic Compounds (Wiley, 1994), Table 11.7; used only to rank the two flips.
A_VALUES = {
    "H": 0.0, "D": 0.006, "F": 0.25, "Cl": 0.53, "Br": 0.48, "I": 0.47,
    "OH": 0.87, "OMe": 0.60, "OAc": 0.71, "OBn": 0.6, "OPh": 0.65, "OEt": 0.6,
    "SH": 1.21, "SMe": 1.04, "NH2": 1.6, "NHMe": 1.3, "NMe2": 2.1, "NHAc": 1.6, "NO2": 1.1,
    "CN": 0.2, "C#CH": 0.41, "CH3": 1.74, "Et": 1.79, "iPr": 2.21, "tBu": 4.9, "vinyl": 1.68,
    "Ph": 2.8, "Bn": 1.68, "CHO": 0.73, "Ac": 1.2, "CO2H": 1.35, "CO2Me": 1.27, "CO2Et": 1.2,
    "CF3": 2.4, "TMS": 2.5, "TES": 2.5, "TBS": 2.5, "TIPS": 2.5, "CH2OH": 1.76,
}
CROWDED_PER_BOND = 0.6  # placement score: non-bonded atoms nearer than this x bond length are crowded (a choice)
CROWDED_LABELED_PER_BOND = 1.0  # ... both labeled (heteroatoms, H): their labels need more room (a choice)


@dataclass
class Placement:
    """The ring on a template: ring[k] sits on vertex (start + direction * k) % 6."""
    ring: list[int]
    flipped: bool
    start: int
    direction: int
    axial: set[int] = field(default_factory=set)  # substituent atom indices drawn axial
    pts: dict[int, Point] = field(default_factory=dict)  # pinned atoms: pt, y down, ring[0] at the origin


def template(c: Chair, flipped: bool) -> tuple[list[Point], list[bool], list[Point]]:
    """(vertices pt y down, axial up per vertex, equatorial unit vector per vertex, y down). The
    flipped template is the exact mirror: its vertex k mirrors vertex (4 - k) % 6 of the other,
    with its vertex 1 at the origin (preset.Chair)."""
    eq = [(math.cos(math.radians(d)), -math.sin(math.radians(d))) for d in c.equatorial_deg]
    if not flipped:
        return [tuple(v) for v in c.vertices], list(c.axial_up), eq
    m = [(4 - k) % 6 for k in range(6)]
    x0 = c.vertices[m[0]][0]
    return ([(x0 - c.vertices[j][0], c.vertices[j][1]) for j in m], [c.axial_up[j] for j in m],
            [(-eq[j][0], eq[j][1]) for j in m])


def orientation(pts: list[Point]) -> int:
    """+1 if the screen polygon (y down) runs counterclockwise seen from above the chair: its
    vertices read with y up (front = low on the page) turn counterclockwise."""
    area = sum(p[0] * (-q[1]) - q[0] * (-p[1]) for p, q in zip(pts, pts[1:] + pts[:1]))
    return 1 if area > 0 else -1


# --------------------------------------------------------------------------- chemistry


def chair_rings(mol: Chem.Mol) -> list[list[int]]:
    """Six-membered rings in ring order."""
    out = []
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) == 6:
            out.append(_ordered(mol, list(ring)))
    return out


def _ordered(mol: Chem.Mol, ring: list[int]) -> list[int]:
    order, members = [ring[0]], set(ring)
    while len(order) < 6:
        prev = order[-2] if len(order) > 1 else None
        nxt = [n.GetIdx() for n in mol.GetAtomWithIdx(order[-1]).GetNeighbors()
               if n.GetIdx() in members and n.GetIdx() != prev and n.GetIdx() not in order]
        order.append(min(nxt))
    return order


def ring_problem(mol: Chem.Mol, ring: list[int]) -> str | None:
    """Why `ring` can't be drawn as a chair, or None."""
    ri = mol.GetRingInfo()
    for i in ring:
        a = mol.GetAtomWithIdx(i)
        if a.GetIsAromatic() or a.GetHybridization() not in (Chem.HybridizationType.SP3,):
            if a.GetSymbol() in ("O", "N", "S") and not a.GetIsAromatic() and all(
                    b.GetBondType() == Chem.BondType.SINGLE for b in a.GetBonds()):
                continue  # ring O / N / S: no tetrahedral carbon, but single bonds only
            return (f"ring atom {_name(a)} is not sp3 (a chair needs six sp3 ring atoms: no double "
                    "bonds or aromatic atoms in or on the ring, so cyclohexene and benzene are drawn flat)")
        if ri.NumAtomRings(i) > 1:
            return (f"ring atom {_name(a)} is in more than one ring: fused, bridged and spiro chairs "
                    "(decalins, steroids) are not supported yet")
    return None


def _name(a: Chem.Atom) -> str:
    return f"{a.GetSymbol()} (map {a.GetAtomMapNum()})" if a.GetAtomMapNum() else f"{a.GetSymbol()} (atom {a.GetIdx() + 1})"


def substituents(mol: Chem.Mol, ring: list[int]) -> dict[int, list[int]]:
    """Ring atom -> its neighbors outside the ring (atoms of `mol`: implicit H not included)."""
    members = set(ring)
    return {i: [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() not in members]
            for i in ring}


def _ccw(v0, v1, v2, v3) -> bool:
    """Seen from neighbor 0 toward the center, do neighbors 1-3 run counterclockwise? (vectors
    from the center, 3D)."""
    a = [v2[k] - v1[k] for k in range(3)]
    b = [v3[k] - v1[k] for k in range(3)]
    n = (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
    return sum(n[k] * v0[k] for k in range(3)) > 0


def faces(mol: Chem.Mol, ring: list[int]) -> dict[int, int | None]:
    """Substituent atom -> +1 / -1: its side of the ring (+1: the side from which `ring` runs
    counterclockwise), or None for an atom without a chiral tag. From the chiral tags: the ring is
    put flat (a regular hexagon, counterclockwise in ring order) and the substituent tried on each
    side, keeping the side whose handedness matches the tag."""
    mh = Chem.AddHs(mol)  # every neighbor explicit, so the tag's neighbor order is the bond order
    pos = {r: (math.cos(math.pi / 3 * k), math.sin(math.pi / 3 * k), 0.0) for k, r in enumerate(ring)}
    out: dict[int, int | None] = {}
    for k, i in enumerate(ring):
        a = mh.GetAtomWithIdx(i)
        subs = [n.GetIdx() for n in a.GetNeighbors() if n.GetIdx() not in pos]
        tag = a.GetChiralTag()
        heavy = [s for s in subs if s < mol.GetNumAtoms()]
        if tag not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW) or len(subs) != 2:
            out.update({s: None for s in heavy})
            continue
        c = pos[i]
        out_dir = (c[0] * 0.5, c[1] * 0.5)

        def place(s, up):
            return (c[0] + out_dir[0], c[1] + out_dir[1], 1.0 if up else -1.0)
        trial = {subs[0]: place(subs[0], True), subs[1]: place(subs[1], False)}
        nbrs = [n.GetIdx() for n in a.GetNeighbors()]  # bond order
        vec = [tuple(q - p for q, p in zip(trial.get(j) or pos[j], c)) for j in nbrs]
        ok = _ccw(*vec) == (tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        first_up = 1 if ok else -1
        for s, f in ((subs[0], first_up), (subs[1], -first_up)):
            if s < mol.GetNumAtoms():
                out[s] = f
    return out


def a_value(mol: Chem.Mol, ring_atom: int, s: int) -> float:
    """A-value (kcal/mol) of substituent atom s on ring_atom: by abbreviation label, else by
    group (A_VALUES); unknown groups by their first atom and its neighbors."""
    a = mol.GetAtomWithIdx(s)
    if a.HasProp("atomLabel"):
        return A_VALUES.get(a.GetProp("atomLabel"), 1.5)
    sym = a.GetSymbol()
    if sym in ("H", "F", "Cl", "Br", "I"):
        return A_VALUES[sym]
    heavy = [n for n in a.GetNeighbors() if n.GetIdx() != ring_atom and n.GetAtomicNum() > 1]
    if sym == "O":
        if not heavy:
            return A_VALUES["OH"]
        n = heavy[0]
        if any(b.GetBondType() == Chem.BondType.DOUBLE for b in n.GetBonds()):
            return A_VALUES["OAc"]
        return A_VALUES["OMe"]
    if sym == "N":
        if a.GetFormalCharge() > 0:
            return A_VALUES["NO2"]
        return {0: A_VALUES["NH2"], 1: A_VALUES["NHMe"]}.get(len(heavy), A_VALUES["NMe2"])
    if sym == "S":
        return A_VALUES["SH"] if not heavy else A_VALUES["SMe"]
    if sym == "C":
        if a.GetIsAromatic():
            return A_VALUES["Ph"]
        orders = [b.GetBondType() for b in a.GetBonds()]
        if Chem.BondType.TRIPLE in orders:
            return A_VALUES["CN"] if any(n.GetSymbol() == "N" for n in heavy) else A_VALUES["C#CH"]
        if Chem.BondType.DOUBLE in orders:
            if any(n.GetSymbol() == "O" for n in heavy):
                return A_VALUES["CO2Me"] if len(heavy) > 1 else A_VALUES["CHO"]
            return A_VALUES["vinyl"]
        return {0: A_VALUES["CH3"], 1: A_VALUES["Et"], 2: A_VALUES["iPr"]}.get(len(heavy), A_VALUES["tBu"])
    return 1.5


# --------------------------------------------------------------------------- placement


@dataclass
class ChairChoice:
    ring: list[int]
    flip: bool | None = None  # None: the lower-energy flip; True: the other one
    axial: set[int] = field(default_factory=set)  # substituent atoms that must be axial
    equatorial: set[int] = field(default_factory=set)
    start: int | None = None  # ring atom on the right tip of the template


def _axial_set(mol, ring, subs, fc, axial_up, verts, sigma) -> set[int]:
    out = set()
    for k, i in enumerate(ring):
        v = verts[k]
        known = [s for s in subs[i] if fc.get(s) is not None]
        free = [s for s in subs[i] if fc.get(s) is None]
        for s in known:
            if (fc[s] == sigma) == axial_up[v]:
                out.add(s)
        if len(subs[i]) == 2 and len(free) == 2:  # geminal without stereo: the bulkier equatorial
            out.add(min(free, key=lambda s: (a_value(mol, i, s), -s)))
        elif len(subs[i]) == 2 and len(free) == 1:  # the other goes where the known one isn't
            if known[0] not in out:
                out.add(free[0])
    return out


def candidates(mol: Chem.Mol, c: Chair, choice: ChairChoice, bond_length: float,
               fc: dict[int, int | None]) -> list[Placement]:
    """Every placement of the ring (24), with its axial set and pinned positions (pt, y down)."""
    ring = choice.ring
    subs = substituents(mol, ring)
    out = []
    for flipped in (False, True):
        tv, up, eq = template(c, flipped)
        sigma_t = orientation(tv)
        for start, direction in itertools.product(range(6), (1, -1)):
            verts = [(start + direction * k) % 6 for k in range(6)]
            ax = _axial_set(mol, ring, subs, fc, up, verts, sigma_t * direction)
            pl = Placement(ring, flipped, start, direction, ax)
            o = tv[verts[0]]
            for k, i in enumerate(ring):
                p = (tv[verts[k]][0] - o[0], tv[verts[k]][1] - o[1])
                pl.pts[i] = p
                for s in subs[i]:
                    d = (0.0, -1.0 if up[verts[k]] else 1.0) if s in ax else eq[verts[k]]
                    pl.pts[s] = (p[0] + d[0] * bond_length, p[1] + d[1] * bond_length)
            out.append(pl)
    return out


def energy(mol: Chem.Mol, ring: list[int], pl: Placement) -> float:
    subs = substituents(mol, ring)
    return sum(a_value(mol, i, s) for i in ring for s in subs[i] if s in pl.axial)


def _rot(d: Point, deg: float) -> Point:
    c, s_ = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    return (d[0] * c - d[1] * s_, d[0] * s_ + d[1] * c)


def second_shell(mol: Chem.Mol, pl: Placement, bond_length: float) -> list[dict[int, Point]]:
    """Ways to place the atoms one bond beyond each substituent (chain_angle 120°: a zigzag on
    from the ring bond): one or two further heavy atoms go at ±60° off the substituent bond's line
    (one: either side, so two choices); three or more are left to RDKit (tBu)."""
    ring = set(pl.ring)
    options: list[list[dict[int, Point]]] = []
    for i in pl.ring:
        for n in mol.GetAtomWithIdx(i).GetNeighbors():
            s = n.GetIdx()
            if s in ring or n.HasProp("atomLabel"):
                continue
            further = [m.GetIdx() for m in n.GetNeighbors() if m.GetIdx() != i and m.GetAtomicNum() > 1]
            if not further or len(further) > 2 or any(m in pl.pts for m in further):
                continue
            d = ((pl.pts[s][0] - pl.pts[i][0]) / bond_length, (pl.pts[s][1] - pl.pts[i][1]) / bond_length)

            def at(deg):
                r = _rot(d, deg)
                return (pl.pts[s][0] + r[0] * bond_length, pl.pts[s][1] + r[1] * bond_length)
            if len(further) == 1:
                options.append([{further[0]: at(60)}, {further[0]: at(-60)}])
            else:
                options.append([{further[0]: at(60), further[1]: at(-60)}, {further[0]: at(-60), further[1]: at(60)}])
    out = []
    for combo in itertools.product(*options) if options else [()]:
        merged: dict[int, Point] = {}
        for o in combo:
            merged.update(o)
        out.append(merged)
    return out[:16]  # (2 choices per substituent; more than 4 such substituents: the first 16)


def _pinned_layout(mol: Chem.Mol, pts: dict[int, Point], bond_length: float) -> Chem.Mol:
    """A copy of mol laid out by RDKit with `pts` (pt, y down) pinned (1.5 A per bond length, y up)."""
    m = Chem.Mol(mol)
    k = 1.5 / bond_length
    cm = {i: Point2D(p[0] * k, -p[1] * k) for i, p in pts.items()}
    prefer = rdDepictor.GetPreferCoordGen()
    rdDepictor.SetPreferCoordGen(False)  # RDKit's depictor holds coordMap atoms exactly
    try:
        rdDepictor.Compute2DCoords(m, coordMap=cm)
    finally:
        rdDepictor.SetPreferCoordGen(prefer)
    return m


def _score(mol: Chem.Mol, m: Chem.Mol, ring: list[int]) -> tuple[float, int]:
    """(crowding, crossings) of a laid-out copy, in bond lengths (1.5 A). Crowding sums how far
    non-bonded atom pairs come inside CROWDED_PER_BOND (CROWDED_LABELED_PER_BOND when one is
    labeled), and labeled atoms inside CROWDED_PER_BOND of a bond not on them."""
    from .geometry import seg_point_dist, segments_cross
    conf = m.GetConformer()
    pts = [(conf.GetAtomPosition(i).x / 1.5, conf.GetAtomPosition(i).y / 1.5) for i in range(m.GetNumAtoms())]
    labeled = [a.GetSymbol() != "C" or a.HasProp("atomLabel") for a in m.GetAtoms()]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in m.GetBonds()]
    near = {frozenset(b) for b in bonds}
    crowd = 0.0
    members = set(ring)
    for i, j in itertools.combinations(range(len(pts)), 2):
        if frozenset((i, j)) in near or (i in members and j in members):
            continue
        lim = CROWDED_LABELED_PER_BOND if labeled[i] or labeled[j] else CROWDED_PER_BOND
        crowd += max(0.0, lim - math.dist(pts[i], pts[j]))
    for i in range(len(pts)):
        if labeled[i]:
            for a, b in bonds:
                if i not in (a, b):
                    crowd += max(0.0, CROWDED_PER_BOND - seg_point_dist(pts[i], pts[a], pts[b]))
    cross = sum(1 for (a, b), (c, d) in itertools.combinations(bonds, 2)
                if not {a, b} & {c, d} and segments_cross(pts[a], pts[b], pts[c], pts[d]))
    return round(crowd, 6), cross


def choose(mol: Chem.Mol, c: Chair, choice: ChairChoice, bond_length: float) -> tuple[Placement, Chem.Mol]:
    """The placement to draw (module docstring) and mol laid out on it."""
    ring = choice.ring
    fc = faces(mol, ring)
    cands = candidates(mol, c, choice, bond_length, fc)
    # substituents without a face take the one that puts them where the spec asks (axial /
    # equatorial), else equatorial in the preferred flip, so that `flip` moves them axial like the rest
    subs = substituents(mol, ring)
    known_ax = {s for s in choice.axial if fc.get(s) is not None}
    known_eq = {s for s in choice.equatorial if fc.get(s) is not None}
    ok = [p for p in cands if known_ax <= p.axial and not (known_eq & p.axial)] or cands
    p0 = min(ok, key=lambda p: energy(mol, ring, p))
    tv, up, _ = template(c, p0.flipped)
    sigma = orientation(tv) * p0.direction
    for k, i in enumerate(ring):
        v = (p0.start + p0.direction * k) % 6
        free = [s for s in subs[i] if fc.get(s) is None]
        want = {s: s in p0.axial for s in free}
        for s in free:
            if s in choice.axial or s in choice.equatorial:
                want[s] = s in choice.axial
                for t in free:  # a geminal partner goes the other way
                    if t != s and t not in choice.axial and t not in choice.equatorial:
                        want[t] = not want[s]
        for s in free:
            fc[s] = sigma if want[s] == up[v] else -sigma
    cands = candidates(mol, c, choice, bond_length, fc)
    if choice.axial or choice.equatorial:
        cands = [p for p in cands if choice.axial <= p.axial and not (choice.equatorial & p.axial)]
        if not cands:
            raise ValueError("chair: no ring flip puts the requested atoms axial / equatorial "
                             "(their faces come from the SMILES stereo)")
    else:
        energies = sorted({round(energy(mol, ring, p), 6) for p in cands})
        want = energies[-1] if choice.flip else energies[0]
        cands = [p for p in cands if round(energy(mol, ring, p), 6) == want]
    if choice.start is not None:
        tips = {False: 4, True: 3}  # right tip vertex (0-based) of each template
        k0 = ring.index(choice.start)
        cands = [p for p in cands if (p.start + p.direction * k0) % 6 == tips[p.flipped]] or cands
    hetero = [k for k, i in enumerate(ring) if mol.GetAtomWithIdx(i).GetSymbol() != "C"]
    substituted = [k for k, i in enumerate(ring) if substituents(mol, ring)[i]]

    def key(p: Placement, m: Chem.Mol):
        crowd, cross = _score(mol, m, ring)
        verts = [(p.start + p.direction * k) % 6 for k in range(6)]
        if hetero:
            pref = 0 if p.flipped and verts[hetero[0]] == c.heteroatom_vertex - 1 else 1
        else:
            pref = 0 if not p.flipped else 1
        right = -p.pts[ring[substituted[0]]][0] if substituted else 0.0
        return (crowd, pref, cross, round(right, 6))
    scored = []
    for n, p in enumerate(cands):
        for k, extra in enumerate(second_shell(mol, p, bond_length)):
            m = _pinned_layout(mol, {**p.pts, **extra}, bond_length)
            scored.append((key(p, m), n, k, p, m))
    _, _, _, p, m = min(scored, key=lambda t: t[:3])
    return p, m


def depths(mol: Chem.Mol, ring: list[int], pl: Placement) -> dict[int, float]:
    """Atom index -> depth (module docstring): ring atoms their screen y, the rest that of the ring
    atom they hang off."""
    z = {i: pl.pts[i][1] for i in ring}
    todo = list(ring)
    while todo:
        i = todo.pop()
        for n in mol.GetAtomWithIdx(i).GetNeighbors():
            j = n.GetIdx()
            if j not in z:
                z[j] = z[i]
                todo.append(j)
    return z


def layout(mol: Chem.Mol, c: Chair, bond_length: float, ring: list[int] | None = None, flip: bool | None = None,
           axial: list[int] = (), equatorial: list[int] = (), start: int | None = None,
           ) -> tuple[Chem.Mol, dict[int, float], Placement]:
    """mol laid out as a chair (1.5 A per bond length, y up), atom depths, and the placement.
    `ring`, `axial`, `equatorial`, `start`: atom-map numbers (ring: its six atoms, needed when
    there is more than one six-membered ring)."""
    by_map = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}

    def idx(maps, what):
        missing = [m for m in maps if m not in by_map]
        if missing:
            raise ValueError(f"chair: {what}: no atom has map number {', '.join(map(str, missing))}")
        return [by_map[m] for m in maps]
    rings = chair_rings(mol)
    if ring:
        want = set(idx(ring, "ring"))
        rings = [r for r in rings if set(r) == want]
        if not rings:
            raise ValueError(f"chair: ring {ring} is not a six-membered ring of the species")
    if not rings:
        raise ValueError("chair: the species has no six-membered ring")
    if len(rings) > 1:
        raise ValueError("chair: the species has more than one six-membered ring; name the one to draw "
                         "with chair: {ring: [six atom-map numbers]}")
    r = rings[0]
    problem = ring_problem(mol, r)
    if problem:
        raise ValueError(f"chair: {problem}")
    subs = {s for ss in substituents(mol, r).values() for s in ss}
    ax, eq = set(idx(axial, "axial")), set(idx(equatorial, "equatorial"))
    if (ax | eq) - subs:
        raise ValueError("chair: axial / equatorial atoms must be bonded to the ring (substituent atoms)")
    if ax & eq:
        raise ValueError("chair: an atom can't be both axial and equatorial")
    st = None
    if start is not None:
        (st,) = idx([start], "start")
        if st not in r:
            raise ValueError("chair: start must be a ring atom")
    pl, m = choose(mol, c, ChairChoice(r, flip, ax, eq, st), bond_length)
    return m, depths(mol, r, pl), pl
