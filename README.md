# acs-chem-svg

A [Claude Code](https://code.claude.com) skill that draws chemical structures, reaction
schemes, and arrow-pushing mechanisms in ACS style (ChemDraw's *ACS Document 1996*) as
editable **SVG**.

Ask in plain language:

> Draw the full arrow-pushing mechanism for the nitration of benzene with HNO3 and H2SO4:
> formation of the nitronium ion, attack by benzene, the resonance structures of the sigma
> complex, and deprotonation to nitrobenzene.

and Claude writes a short YAML spec, validates the chemistry and the electron bookkeeping,
renders it, looks at the rendered image, fixes what it sees, and hands back `.svg`, `.yaml`
and `.png` files:

![Nitration of benzene: mechanism with curved arrows, lone pairs and a bracketed set of resonance structures](docs/images/nitration.png)

![A three-step synthesis scheme with conditions and bold compound numbers](docs/images/synthesis.png)

## What it draws

- **Structures** from SMILES: rings, double and triple bonds, wedges and hashes, heteroatom
  labels with H counts and charges, abbreviations (OMe, OBn, CO₂Me, NO₂, Ph, …).
- **Schemes:** reaction arrows with conditions above and below (true subscripts and
  superscripts, italics), two-stage arrows, plus signs, bold compound numbers, products
  aligned to their reactants, rows wrapped to a page width.
- **Mechanisms:** two-electron and fishhook curved arrows, lone pairs, radicals, charges, and
  resonance structures in brackets joined by ↔.
- **Three sizes:** the ACS 1996 preset at 100 %, 80 % and 60 %. Bond lengths, line widths,
  double-bond spacing and label clearance are measured from ChemDraw exports and match them
  within 0.1 pt.

The SVG has real, editable text (Helvetica) and named groups, so it opens cleanly in a vector editor
(checked in Affinity) for final touches.

The output is a drawing, not ChemDraw-editable chemistry. Validation checks that every step's
curved arrows turn its reactants into its products exactly; it does not check that the
mechanism is the right one. Review research-level mechanisms yourself.

## Installation

The usual setup is one copy of the repository shared by any number of project folders, each
of which links to it. The skill is then available to every Claude Code session started in
those folders, and updating the repository updates all of them.

**Requirements:** [Claude Code](https://code.claude.com), [uv](https://docs.astral.sh/uv/)
(it installs Python 3.12+ and the dependencies, including RDKit), git. Developed on macOS; Linux should work with the bundled fallback font.

1. **Clone the repository** once, anywhere you like (here `~/tools`):

   ```bash
   git clone https://github.com/scotthartley/acs-chem-svg.git ~/tools/acs-chem-svg
   uv sync --project ~/tools/acs-chem-svg
   ```

   `uv sync` is optional (the first use does it too), but it gets the download out of the way.

2. **Link it into your project folder** (the folder you start Claude Code in):

   ```bash
   cd ~/path/to/my-project
   mkdir -p .claude/skills
   ln -sfn ~/tools/acs-chem-svg .claude/skills/acs-chem-svg
   ```

   Use `-sfn`: rerunning a plain `ln -s` on an existing link creates a stray link inside the
   repository. Repeat this step for every project folder that should have the skill.

3. **Start Claude Code** in the project folder and ask for a drawing. The skill loads on its
   own when the request is about structures, schemes or mechanisms; you can also invoke it
   with `/acs-chem-svg`.

To update, `git -C ~/tools/acs-chem-svg pull`.

Don't link the skill into the repository itself (a link back to the repo makes test
collection loop); development sessions there don't need it.

**Fonts.** On macOS, text is measured with the system Helvetica. Elsewhere the bundled TeX
Gyre Heros (a Helvetica clone, `fonts/`) is used, and Claude tells you so; the SVG still asks
for Helvetica, so it looks right on a machine that has it.

## Using it

Ask for what you want as you would ask a colleague. Name the preset if you want a smaller
one ("at 80 %", "at 60 %"), and give the conditions you want on the arrows; Claude uses only
the conditions you name. Some starting points:

- "Draw the SN2 mechanism for hydroxide reacting with bromomethane, with curved arrows."
- "Draw the radical bromination of isobutane with Br2 and light: initiation and both
  propagation steps, with fishhook arrows."
- "Draw a three-step synthesis scheme: methyl 4-hydroxybenzoate is benzylated with BnBr and
  K2CO3 in acetone, reduced with LiAlH4 in THF, and oxidized with PCC to the aldehyde; number
  the compounds."

Claude saves the files in your project folder. To change a drawing, ask for the change
("flip the benzene ring", "use the 80 % size", "make it fit 3.25 in"); Claude edits
the `.yaml` spec and re-renders, so keep the `.yaml` next to the `.svg`. The spec format is
documented in [`references/scheme-schema.md`](references/scheme-schema.md) if you want to
write or edit specs yourself; `examples/` has a spec for every kind of drawing.

When something goes wrong, Claude appends a note to `issues/inbox.md` in the repository
(local, not tracked by git). If it turns up a bug, please open a GitHub issue with the note.

## Limits

- Not yet supported: explicit CH₂/CH labels on internal carbons, charges beyond ±1,
  carbenes, more than two stages on one arrow, catalytic cycles.
- An abbreviation whose bond leaves to the right (MeO–, BnO–) isn't measured yet; Claude turns
  the structure or writes the group out.
- Curved-arrow placement is automatic with a few semantic overrides (bow side, curvature);
  crowded mechanisms can need a nudge in a vector editor.

## Documentation

| File | Contents |
|---|---|
| [`SKILL.md`](SKILL.md) | What Claude reads when it uses the skill |
| [`references/scheme-schema.md`](references/scheme-schema.md) | The YAML spec format |
| [`references/workflow.md`](references/workflow.md) | Worked example: from a request to a mechanism |
| [`references/acs-style.md`](references/acs-style.md) | The measured ACS style, with the ChemDraw reference behind each number |
| [`docs/design.md`](docs/design.md) | Design brief: architecture, validation, milestones |
| [`CLAUDE.md`](CLAUDE.md) | Rules for development sessions in this repository |

## Development

Open Claude Code in the repository itself to work on the renderer. `CLAUDE.md` holds the
rules (every style number traces to a ChemDraw reference; tests and the regression gallery
run before each commit). Tests: `uv run pytest`.

The ChemDraw reference exports that the style numbers were measured from (`reference_renders/`)
are not in the public repository; the tests that compare against them skip without them.

## License

MIT (`LICENSE`). The bundled TeX Gyre Heros fonts are under the GUST Font License
(`fonts/GUST-FONT-LICENSE.txt`).
