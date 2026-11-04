"""Convert each fragment of a ChemDraw reference CDXML into a molfile template.

The output keeps ChemDraw's atom coordinates, wedges, charges and H counts, with each
CDXML atom id stored as the atom-map number (ids above 999, which V2000 can't hold, in
<reference>__<n>.ids.json, by atom index). Rendering these templates compares our
drawing conventions with ChemDraw's on an identical layout (scripts/compare.py,
tests/test_m1_vs_chemdraw.py). Coordinates are pt x 1.5 / BondLength (1.5 A bonds,
y up), so chem.from_molfile scales them back to the exact positions on ChemDraw's SVG page.

Writes templates/m1/<reference>__<n>.mol (templates/m2/ for scheme references, templates/m3/
for mechanism references, radicals kept; templates/m4/ for resonance references; templates/abbrev/
for the abbreviation references; templates/chair/ for chair_x*), with n =
fragment order in the CDXML (1-based). Fragments whose bonds cross (CDXML CrossingBonds: the
crossing cases of chair_x*, one bond per fragment) are merged into the first one's template, and
templates/chair/<reference>__<n>.front.json lists the bonds ChemDraw draws in front at each
crossing (the higher Z) as CDXML atom-id pairs (acs-style.md §14.2). Labeled carbons (CH3 as text) become molfile atom
aliases, which chem.from_molfile reads as explicit carbon labels. Nicknames that are known
abbreviations (OMe, OBn) become their attachment atom with the nickname as its alias; fragments
with other nicknames (CH3Cl as text) are skipped. Formula text (K2CO3) has no atoms: it is a
caption (measure_refs.Doc.text_frags).

Usage: uv run python scripts/cdxml_to_molfile.py [reference names]   (default: all below)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Geometry import Point3D

sys.path.insert(0, str(Path(__file__).resolve().parent))
import measure_refs as mr  # noqa: E402

from acs_chem_svg import abbrev, chem  # noqa: E402
from acs_chem_svg.chem import ABBREVIATIONS  # noqa: E402

OUT = mr.REPO / "templates" / "m1"
OUT_M2 = mr.REPO / "templates" / "m2"
OUT_M3 = mr.REPO / "templates" / "m3"
OUT_M4 = mr.REPO / "templates" / "m4"
OUT_ABBREV = mr.REPO / "templates" / "abbrev"
OUT_CHAIR = mr.REPO / "templates" / "chair"
ANUM = {v: k for k, v in mr.ELEMENTS.items()}
DIRS = {"WedgeBegin": Chem.BondDir.BEGINWEDGE, "WedgedHashBegin": Chem.BondDir.BEGINDASH}
M1_REFERENCES = ["m1_set_x1.0", "m1_set_x0.8", "m1_set_x0.6", "double_bonds_x1.0", "labels_x1.0",
                 "double_bonds_2", "stereo_x1.0", "stereo_x0.8", "stereo_x0.6",
                 "clearance_x1.0", "clearance_x0.8", "clearance_x0.6"]
M2_REFERENCES = ["scheme_basic_x1.0", "scheme_basic_x0.8", "scheme_basic_x0.6",
                 "scheme_multi_x1.0", "scheme_multi_x0.8", "scheme_multi_x0.6"]
M3_REFERENCES = ["mech_basic_x1.0", "mech_basic_x0.8", "mech_basic_x0.6"]
M4_REFERENCES = ["resonance_x1.0", "resonance_x0.8", "resonance_x0.6",
                 "scheme_wrap_x1.0", "scheme_wrap_x0.8", "scheme_wrap_x0.6"]
ABBREV_REFERENCES = ["abbrev_x1.0", "abbrev_x0.8", "abbrev_x0.6", "abbrev_2_x1.0", "abbrev_2_x0.8", "abbrev_2_x0.6",
                     "ligands_x1.0", "ligands_x0.8", "ligands_x0.6"]
CHAIR_REFERENCES = mr.CHAIR_REFS



def build_mol(doc: mr.Doc, pt_per_unit: float = 9.6, atom_ids=None, map_ids: bool = False,
              nicknames: bool = False):
    """RDKit mol with the CDXML atoms/bonds/coords; returns (mol, cdxml atom ids by index).

    Coordinates are pt / pt_per_unit with y up (9.6 turns 14.4 pt bonds into 1.5 A). Pass
    `atom_ids` to keep one fragment, and `map_ids` to store each CDXML id as the atom map.
    `nicknames` keeps nickname nodes (OMe) as their attachment atom; the caller sets the alias."""
    rw = Chem.RWMol()
    ids = []
    for aid, a in doc.atoms.items():
        if (a.nickname and not nicknames) or (atom_ids is not None and aid not in atom_ids):
            continue
        at = Chem.Atom(ANUM.get(a.element, 0))  # 0: a text label's stand-in (dppf: no inner fragment)
        at.SetFormalCharge(a.charge)
        if a.label and a.num_h is not None and not a.nickname:
            at.SetNumExplicitHs(a.num_h)
            at.SetNoImplicit(True)
        if a.radicals:  # tBu·, Cl· (mech_basic): no implicit H on the radical center
            at.SetNumRadicalElectrons(a.radicals)
            at.SetNoImplicit(True)
            if not a.label:
                at.SetNumExplicitHs(4 - a.radicals - sum(1 for b in doc.bonds if aid in (b.a, b.b)))
        if map_ids:
            at.SetAtomMapNum(int(aid))
        rw.AddAtom(at)
        ids.append(aid)
    idx = {aid: i for i, aid in enumerate(ids)}
    order = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    for b in doc.bonds:
        if b.a not in idx or b.b not in idx:
            continue
        a, e = (b.a, b.b) if not b.display.endswith("End") else (b.b, b.a)
        rw.AddBond(idx[a], idx[e], order[b.order])
        display = b.display.replace("End", "Begin")  # *End: the bond was reversed above
        if display in DIRS:
            rw.GetBondBetweenAtoms(idx[a], idx[e]).SetBondDir(DIRS[display])
    mol = rw.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for aid, i in idx.items():
        x, y = doc.atoms[aid].p
        conf.SetAtomPosition(i, Point3D(x / pt_per_unit, -y / pt_per_unit, 0.0))
    mol.AddConformer(conf, assignId=True)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                     ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    Chem.AssignChiralTypesFromBondDirs(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol, ids

def fragments(doc: mr.Doc) -> list[tuple[int, str, set[str]]]:
    """(n, fragment id, atom ids) for fragments drawn with real atoms and bonds. Fragments
    with a nickname other than a known abbreviation are skipped: we don't draw those labels.
    Fragments with crossing bonds (CDXML CrossingBonds) are merged into the first of them."""
    out = []
    _, pairs = mr._bond_z(doc.name) if (mr.REF_DIR / f"{doc.name}.cdxml").exists() else ({}, set())
    frag_of = {b.id: b.frag for b in doc.bonds}
    order = list(doc.fragments)
    merged = {}  # fragment id -> the fragment id it is drawn with
    for p in pairs:
        f1, f2 = sorted((frag_of[i] for i in p), key=order.index)
        if f1 != f2:
            merged[f2] = merged.get(f1, f1)
    for n, fid in enumerate(doc.fragments, 1):
        if fid in merged:
            continue
        group = {fid} | {f for f, to in merged.items() if to == fid}
        atoms = {aid for aid, a in doc.atoms.items() if a.frag in group}
        if atoms and not any(doc.atoms[a].nickname and not _drawable(doc, a) for a in atoms):
            out.append((n, fid, atoms))
    return out


def label_key(label: str | None) -> str | None:
    """A nickname's registry key: a label typed in its right-justified form, where that is an
    entry's `reversed:` override (pinB), is its key (Bpin); any other label is itself."""
    typed = {abbrev.plain(e.reversed): k for k, e in abbrev.default().entries.items() if e.reversed}
    return typed.get(label, label)


def _drawable(doc: mr.Doc, aid: str) -> bool:
    """A nickname we draw from a template alias: a registry group or formula, or a label composed
    of its element and registry groups (ring NTs)."""
    a = doc.atoms[aid]
    label = label_key(a.label)
    e = abbrev.default().get(label)
    if label in ABBREVIATIONS or (e is not None and e.kind == "text"):  # text: L, dppf on a stand-in
        return True
    degree = sum(1 for b in doc.bonds if aid in (b.a, b.b))
    return chem.composed_info(label, a.element, degree) is not None


def front_bonds(doc: mr.Doc, atoms: set[str]) -> list[list[str]]:
    """The bonds among `atoms` drawn in front where two cross: the higher CDXML Z of each
    CrossingBonds pair, as [atom id, atom id]."""
    z, pairs = mr._bond_z(doc.name)
    bonds = {b.id: b for b in doc.bonds}
    out = []
    for p in sorted(pairs, key=sorted):
        front = bonds[max(p, key=lambda i: z[i])]
        if front.a in atoms and front.b in atoms:
            out.append([front.a, front.b])
    return out


def ids_path(name: str, n: int) -> Path:
    """CDXML atom ids by atom index, for a template whose ids don't fit V2000 atom maps (> 999)."""
    return template_path(name, n).with_suffix(".ids.json")


def front_path(name: str, n: int) -> Path:
    return template_path(name, n).with_suffix(".front.json")


def with_wedges(block: str, doc: mr.Doc, ids: list[str]) -> str:
    """Write ChemDraw's wedge and hash flags into the V2000 bond block. RDKit's writer derives
    wedges from chirality and drops wedges on non-stereocenters, but ChemDraw wedges are
    drawing marks (stereo_x*: wedges into ring and chain carbons)."""
    flags = {}
    for b in doc.bonds:
        if "Wedge" in b.display and b.a in ids and b.b in ids:
            narrow, wide = (b.a, b.b) if b.display.endswith("Begin") else (b.b, b.a)
            flags[(ids.index(narrow) + 1, ids.index(wide) + 1)] = 6 if "Hash" in b.display else 1
    lines = block.splitlines()
    n_atoms, n_bonds = int(lines[3][:3]), int(lines[3][3:6])
    for k in range(4 + n_atoms, 4 + n_atoms + n_bonds):
        a, e, order = int(lines[k][:3]), int(lines[k][3:6]), lines[k][6:9]
        if (a, e) in flags or (e, a) in flags:
            flag = flags.get((a, e)) or flags[(e, a)]
            if (a, e) not in flags:
                a, e = e, a  # the narrow end must come first
            lines[k] = f"{a:3d}{e:3d}{order}{flag:3d}" + lines[k][12:]
        elif lines[k][9:12].strip() not in ("", "0"):
            lines[k] = lines[k][:9] + "  0" + lines[k][12:]
    return "\n".join(lines) + "\n"


def template_path(name: str, n: int) -> Path:
    out = (OUT_M2 if name.startswith("scheme_") else OUT_M3 if name.startswith("mech_")
           else OUT_M4 if name.startswith("resonance_") else OUT_ABBREV if name.startswith(("abbrev_", "ligands_"))
           else OUT_CHAIR if name.startswith("chair_") else OUT)
    return out / f"{name}__{n}.mol"


def convert(name: str) -> list[Path]:
    doc = mr.load(name)
    pt_per_unit = float(doc.settings["BondLength"]) / 1.5
    for a in doc.atoms.values():  # a label ChemDraw kept as text (p-Tol, Bpin: no inner fragment)
        if a.nickname:
            a.label = label_key(a.label)
        if a.nickname and a.element == "X" and a.label in ABBREVIATIONS:
            a.element = chem.abbrev_element(a.label)
        elif a.nickname and a.element == "X":  # a text label (dppf): a dummy stand-in atom
            a.element = "*"
    written = []
    for n, fid, atoms in fragments(doc):
        big = max(int(a) for a in atoms) > 999  # V2000 atom maps hold 3 digits: ids go to a sidecar
        mol, ids = build_mol(doc, pt_per_unit, atoms, map_ids=not big, nicknames=True)
        mol.GetConformer().Set3D(False)
        mol.SetProp("_Name", f"{name} fragment {fid}")
        for a, aid in zip(mol.GetAtoms(), ids):
            if (doc.atoms[aid].element == "C" or doc.atoms[aid].nickname) and doc.atoms[aid].label:
                a.SetProp("molFileAlias", doc.atoms[aid].label)  # explicit carbon label or abbreviation
        path = template_path(name, n)
        path.write_text(with_wedges(Chem.MolToMolBlock(mol), doc, ids))
        written.append(path)
        if big:
            ids_path(name, n).write_text(json.dumps(ids) + "\n")
            written.append(ids_path(name, n))
        elif ids_path(name, n).exists():
            ids_path(name, n).unlink()
        front = front_bonds(doc, atoms) if name.startswith("chair_") else []
        if front:
            front_path(name, n).write_text(json.dumps(front) + "\n")
            written.append(front_path(name, n))
    return written


def main(argv=None) -> int:
    names = (argv if argv is not None else sys.argv[1:]) or (M1_REFERENCES + M2_REFERENCES + M3_REFERENCES
                                                          + M4_REFERENCES + ABBREV_REFERENCES + CHAIR_REFERENCES)
    for d in (OUT, OUT_M2, OUT_M3, OUT_M4, OUT_ABBREV, OUT_CHAIR):
        d.mkdir(parents=True, exist_ok=True)
    for name in names:
        for p in convert(name):
            print(p.relative_to(mr.REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
