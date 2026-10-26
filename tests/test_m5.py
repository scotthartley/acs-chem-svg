"""M5: natural-language workflow aids (species reused by name, the --maps overlay, ...)."""

import re
from pathlib import Path

from render import render

from rdkit import Chem

from acs_chem_svg.spec import Spec

ROOT = Path(__file__).resolve().parents[1]

REUSE = {
    "species": {"a": {"smiles": "[CH3:1][Br:2]"}, "b": {"smiles": "[OH-:3]", "orientation": {"rotate": 90}},
                "c": {"smiles": "[CH3:1][OH:3]"}, "d": {"smiles": "[Br-:2]"}},
    "steps": [{"from": ["a", "b"], "to": ["c", "d"]}, {"new_row": True, "from": ["a", "b"], "to": "c"}],
}


def test_a_repeated_species_is_cloned_per_occurrence():
    spec = Spec.model_validate(REUSE)
    assert spec.order() == ["a", "b", "c", "d", "a~2", "b~2", "c~2"]
    assert spec.steps[1].sources == ["a~2", "b~2"] and spec.steps[1].products == ["c~2"]
    assert spec.species["b~2"].orientation.rotate == 90  # a reactant keeps its orientation


def test_a_clone_can_be_declared_to_set_its_own_options():
    spec = Spec.model_validate({**REUSE, "species": {**REUSE["species"], "c~2": {"smiles": "[CH3:1][OH:3]",
                                                                                 "align": False}}})
    assert spec.species["c~2"].align is False


def test_reuse_renders_the_same_pixels_as_numbered_copies(tmp_path):
    reuse = render(ROOT / "tests" / "specs" / "nitration_reuse_x1.0.yaml", tmp_path / "reuse")
    copies = render(ROOT / "examples" / "m4_nitration_x1.0.yaml", tmp_path / "copies")
    assert reuse[1].read_bytes() == copies[1].read_bytes()
    ids = [re.findall(r'<g id="([^"]*)"', p.read_text())[1:] for p in (reuse[0], copies[0])]  # [0]: the file stem
    assert ids[0] == ids[1]  # clones are named like the copies (species-nitronium_2)


def test_maps_overlay_leaves_the_svg_unchanged(tmp_path):
    from render import maps_svg

    from acs_chem_svg.layout import build_scheme
    from acs_chem_svg.spec import load_spec

    path = ROOT / "examples" / "m3_carbonyl_addition_x1.0.yaml"
    plain = render(path, tmp_path / "plain")
    with_maps = render(path, tmp_path / "maps", with_maps=True)
    assert plain[0].read_text() == with_maps[0].read_text()
    assert "atom-maps" not in with_maps[0].read_text()
    assert with_maps[-1].name == "m3_carbonyl_addition_x1.0_maps.png" and with_maps[-1].stat().st_size
    lay = build_scheme(load_spec(path), path.parent, path.stem)
    drawn = sorted(str(a.map) for p in lay.species.values() for a in p.mol.atoms if a.map)
    overlay = maps_svg(lay).split('<g id="atom-maps">')[1]
    assert sorted(re.findall(r">(\d+)</text>", overlay)) == drawn


def test_validate_prints_nothing_on_stderr():
    # r0 trials read fontTools' parse notes as the fallback-font warning
    import subprocess
    import sys
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "validate.py"), str(ROOT / "examples" / "m3_sn2_x1.0.yaml")],
                       capture_output=True, text=True)
    assert r.returncode == 0 and r.stderr == ""


def test_the_reported_bow_reproduces_the_default():
    # render.py prints each curved arrow's bow so overrides aren't blind (r0 carbonyl_addition)
    import yaml as _yaml

    from acs_chem_svg.layout import build_scheme
    for name in ("m3_sn2_x1.0", "m3_carbonyl_addition_x1.0", "m3_e2_x1.0"):
        path = ROOT / "examples" / f"{name}.yaml"
        lay = build_scheme(Spec.model_validate(_yaml.safe_load(path.read_text())), path.parent, name)
        raw = _yaml.safe_load(path.read_text())
        raw["overrides"] = {"arrows": [{"step": a.step, "index": a.index, "bow": a.bow} for a in lay.curved]}
        again = build_scheme(Spec.model_validate(raw), path.parent, name)
        assert [(a.p0, a.p3, a.bow) for a in again.curved] == [(a.p0, a.p3, a.bow) for a in lay.curved], name


def test_an_alignment_that_bends_nitronium_is_dropped():
    # r0 nitration (x0.8): nitronium aligned to bent HNO3 atoms was drawn bent
    from acs_chem_svg import chem
    ref = chem.prepare(chem.parse_smiles("[O:1]=[N+:2]([O-:3])[OH+:4][H:10]"), "nitric_h", 14.4)
    p = chem.prepare(chem.parse_smiles("[O:1]=[N+:2]=[O:3]"), "nitronium", 14.4, refs=[ref])
    assert p.aligned_to is None
    assert chem.bent_linear_atoms(p.mol) == []
    pinned = Chem.Mol(chem.parse_smiles("[O:1]=[N+:2]=[O:3]"))
    chem.align_coords(pinned, ref.mol, chem.correspondence(ref.mol, pinned))
    assert chem.bent_linear_atoms(pinned) == [1]


def test_silyl_abbreviations():
    # inbox 2026-10-25 (use session): "unknown abbreviation TMS"
    from acs_chem_svg.spec import molecules
    spec = Spec.model_validate({"species": {
        "tmsa": {"smiles": "C#C[Si](C)(C)C", "abbreviate": ["TMS"]},
        "tbs": {"smiles": "CC(C)(C)[Si](C)(C)OCCO", "abbreviate": ["TBS"]}}})
    tmsa, tbs = molecules(spec, ROOT)
    assert [a.abbrev for a in tmsa.atoms if a.abbrev] == ["TMS"] and len(tmsa.atoms) == 3
    assert [a.abbrev for a in tbs.atoms if a.abbrev] == ["TBS"]
