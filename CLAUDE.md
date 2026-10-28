# acs-chem-svg: development sessions

`docs/design.md` is the design brief (README.md is the introduction for new users; code and
references cite the brief as `design.md §n`). `references/acs-style.md` is the measured style spec. Where
they disagree, `acs-style.md` (the ChemDraw references) wins.

## Location

- The repo is at `~/Coding/productivity/acs-chem-svg`.
- Install per project, in each project folder that should have it (not globally, and never
  inside this repo: a link back to the repo makes tree walks such as pytest collection loop):
  `mkdir -p .claude/skills && ln -sfn ~/Coding/productivity/acs-chem-svg .claude/skills/acs-chem-svg`.
  `-n` matters: a plain `ln -s` rerun on an existing link nests a new link inside the repo.
  Development sessions here don't need the skill installed.
- Use sessions run `uv run --project ~/Coding/productivity/acs-chem-svg ...`.

- `reference_renders/`, `spike/` and `issues/` are local only (gitignored, not in the public
  repo). `reference_renders/cdxml/CDXML.dtd` is Revvity's DTD: never commit it. Tests that need `reference_renders/` skip without it, and `scripts/gallery.py` leaves out
  the ChemDraw-reference cases; keep the folder here so those run before every commit.

## Rules

- **Numbers come from measurement only.** Every style number in `presets/` or
  `references/acs-style.md` must trace to a file in `reference_renders/chemdraw/` plus a method in
  `scripts/measure_refs.py`, or be marked `README, unmeasured` (the design brief's values) / `unmeasured`. Never fill a value
  from memory. If the references can't answer it, add it to the Open list and ask for a new
  reference.
  Hand-placed placement keys (each preset's `placement:` list: scheme gaps, conditions offsets,
  number gap, charge positions) are round user-chosen defaults on a 0.5 pt grid; the measured
  value stays in the comment as provenance.
- **Exact vs approximate.** ChemDraw-generated geometry is exact: M1 holds it to 0.1 pt.
  Hand-placed positions (arrows, conditions, plus signs, numbers, lone pairs) are starting
  defaults only.
- **Triage `issues/inbox.md`** at the start of a session. Each item becomes one of:
  - a fix plus a regression test,
  - a dated `references/lessons.md` entry,
  - discarded (say why).
- **Before every commit,** run `uv run pytest`, `uv run python scripts/check_rasterizer.py`,
  and (from M1) the regression gallery.
- **Gallery changes** (`tests/gallery/`) require the user's approval of the rendered images.
- **Keep SKILL.md under ~150 lines.** Push detail into `references/`.
- **One rasterizer:** resvg via `acs_chem_svg.preview.rasterize`, with explicit font files only.
  Text metrics come from `acs_chem_svg.fontmetrics`, which applies the same pair kerning as resvg.

## Tools

- `uv run python scripts/measure_refs.py [--summary | --json out.json] [names]` measures ChemDraw exports.
- `uv run python scripts/compare_fonts.py` compares Helvetica with the Heros fallback (advances, kerning).
- `uv run python scripts/check_rasterizer.py` is the rasterizer/font acceptance test.
- `uv run python scripts/spike_rdkit.py` runs the M0 RDKit comparison (writes `spike/`).
- `uv run python scripts/measure_refs.py --fit-clip` refits the bond-to-label clip model (acs-style.md §2.1).
- `uv run python scripts/measure_refs.py --scheme` prints the scheme parameters per preset (acs-style.md §3).
- `uv run python scripts/measure_refs.py --mech` prints the mechanism parameters per preset (acs-style.md §4).
- `uv run python scripts/measure_refs.py --wrap` prints the multi-row scheme parameters per preset (acs-style.md §3.4).
- `uv run python scripts/measure_refs.py --resonance` prints the resonance arrow and bracket parameters per preset (acs-style.md §4.3).
- `uv run python scripts/measure_refs.py --abbrev` prints the abbreviation label parameters per preset, both forms (acs-style.md §3.5).
- `uv run python scripts/measure_refs.py --equilibrium` prints the equilibrium arrows and arrowhead variants per preset, incl. arrows_large_x1.0 (acs-style.md §3.1, §3.6).
- `uv run python scripts/cdxml_to_molfile.py` regenerates `templates/m1/` … `templates/m4/` and `templates/abbrev/` from the reference CDXMLs.
- `uv run python scripts/readme_images.py` re-renders the README's example images (`docs/images/`).
- `uv run python scripts/compare.py <ref>...` renders a reference's templates and prints ours vs
  ChemDraw, plus side-by-side PNGs in `build/compare/`. Scheme references also get our arrows and captions;
  mech_basic_x* also gets our curved arrows on ChemDraw's CurvePoints and our dots at its centers.
  `--spec examples/m2_scheme_multi_x1.0.yaml` adds our own layout's approximate items (arrow lengths, gaps);
  resonance_x* (with `examples/m4_resonance_x*.yaml`) and scheme_wrap_x* (`examples/m4_wrap_x*.yaml`) add
  bracket and row items. resonance_x* also gets our ↔ arrows and brackets on ChemDraw's. Arrows take the head
  variant of their CDXML HeadSize; equilibrium_x* and arrows_large_x* draw our arrows and captions only.
- `uv run python scripts/render.py spec.yaml [--out DIR] [--crops] [--maps] [--cdxml]` renders a spec: an M2–M4 scheme
  (`examples/m2_*.yaml` … `examples/m4_*.yaml`) or M1 species (`examples/m1_smiles_x*.yaml`). `--crops` adds
  16 px/pt crops of curved arrows, charged or dotted labels, ↔ arrows, equilibrium arrow ends and bracket corners; `--maps` adds
  `<stem>_maps.png` with map numbers beside the atoms (preview only). It prints each curved arrow's bow.
  `--cdxml` (M1 species only) adds `<stem>.cdxml` for ChemDraw (`cdxml_writer`, acs-style.md §13).
- `uv run python scripts/cdxml_headers.py` regenerates the CDXML header templates `presets/cdxml/` from m1_set_x*.
- `uv run python scripts/measure_refs.py --cdxml-symbols` prints the charge-symbol boxes per reference (acs-style.md §13).
- `uv run python scripts/cdxml_spike.py [--diff]` writes the CDXML gate files to `build/cdxml_spike/`; `--diff`
  compares ChemDraw's re-saves (`build/cdxml_spike/resaved/`) with them (acs-style.md §13, gate).
- `uv run python scripts/validate.py [--suggest] spec.yaml` checks schema, chemistry, electron bookkeeping, scheme
  overlaps and mechanism geometry, and prints paste-ready `suggested:` products on a mismatch (`--suggest`: for
  every step). `tests/broken/*.yaml` must each fail with their `# expect:` message.
- `examples/benchmarks/prompts.yaml` holds the M5 one-sentence benchmark prompts with reference specs; its
  README says how a trial is run (frozen copy of the skill, fresh agent, `build/m5_trials/<round>/`).
- The regression gallery is `tests/test_gallery.py` (part of `uv run pytest`), with CDXML text snapshots in `tests/gallery/cdxml/`. `scripts/gallery.py` reports which cases
  differ; `--review DIR` writes the new or changed renders for the user to look at. `--update` rewrites
  `tests/gallery/`; run it only after the user approves the new images.
- uv commands (including `uv run pytest`) run inside the sandbox; `~/.cache/uv` and `~/.local/share/uv` are allowed.
