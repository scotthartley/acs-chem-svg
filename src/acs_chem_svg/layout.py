"""Scheme layout (design.md §10 M2, acs-style.md §3): one row of species, plus signs, reaction
arrows with conditions, and compound numbers.

Rules (values are presets' `scheme:` keys: round placement defaults chosen from
scheme_basic_x* and scheme_multi_x*, user decision 2026-10-24):
- Vertical: if every structure in the row has a ring, the ring nearest each structure's ink
  center is centered on the arrow axis (center_above_axis, 0 by default); otherwise the ink
  centers are. The plus sign's ink center sits there too.
- Horizontal: structure ink -> plus ink plus_gap; structure ink -> arrow tail and head -> next
  structure ink arrow_gap.
- Conditions: above, stacked upward with the last baseline conditions_above_baseline above
  the axis; below, stacked downward from conditions_below_baseline. Line pitch from
  richtext.line_pitch. Centered on tail -> notch; stages (numbered) left-justified with their
  ink at tail + multistep_indent, the first ⌈n/2⌉ above and the rest below.
- Arrow length = max(arrow_min_length, widest line + 2 conditions_pad + head_inset), or for
  stages max(min, indent + widest + pad + head_inset).
- Equilibrium arrows (a step's arrow_type, arrows.equilibrium_arrow): the same rules, with the
  above/below baselines measured from the outer shafts (equilibrium.shaft_spacing / 2 off the
  axis), the conditions centered on the arrow, and a head_inset at each end in the length
  (unmeasured rules, acs-style.md §3.6). Heads: spec arrow_head, else the preset's default.
- Compound numbers: bold, ink-centered under each structure, on one baseline: the tallest
  number's ink top is number_gap below the lowest structure ink in the row.

M3 mechanisms: the electron bookkeeping (electrons.check_spec) decides which atoms show their
lone pairs (electrons.lone_pair_atoms); the species are drawn with those dots, so the gaps
above count them. After placement, each step's curved arrows are laid out on its reactants
(mechanism.plan_step).

Groups (editor layers): species-<name> (with -bonds, -labels, -charges, -lone-pairs,
-radicals), plus-<i>, arrow-step-<n>, conditions-step-<n> (with -above / -below), numbers (with
number-<name>), curved-arrows-step-<n> (with curved-arrow-step-<n>-<index>).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from . import electrons, richtext
from .arrows import conditions_mid, heads_along, resonance_arrow, step_arrow
from .brackets import bracket_pair
from .chem import Molecule
from .draw_mol import DotSet, Drawer
from .geometry import Point, add, sub, union_box, unit
from .mechanism import AtomRef, MechArrow, Parts, plan_step, species_parts
from .preset import Preset, load_preset
from .scene import Group, Polyline, TextRun
from .spec import ArrowSet, ResonanceStep, Spec, load_spec, molecules, scheme_molecules


@dataclass
class Placed:
    name: str
    group: Group
    mol: Molecule
    ring_centers: list[Point]
    offset: Point = (0.0, 0.0)  # added to mol atom coordinates by the layout
    drawer: Drawer | None = None  # labels and dots, in mol coordinates

    def translate(self, dx: float, dy: float) -> None:
        self.group.translate(dx, dy)
        self.ring_centers = [(x + dx, y + dy) for x, y in self.ring_centers]
        self.offset = (self.offset[0] + dx, self.offset[1] + dy)


@dataclass
class StepLayout:
    n: int
    tail: Point
    head: Point
    above: list[richtext.Line]
    below: list[richtext.Line]
    staged: bool
    row: int = 0
    group: Group | None = None  # arrow-step-<n>
    kind: str = "forward"  # spec arrow_type
    variant: str | None = None  # head variant (spec arrow_head; None: the preset's default)

    def shaft_offset(self, preset: Preset) -> float:
        """How far the outer shafts sit off the axis (0 for a forward arrow)."""
        return 0.0 if self.kind == "forward" else preset.equilibrium.shaft_spacing / 2

    def mid_x(self, preset: Preset) -> float:
        """x the conditions (and reactants drawn above the arrow) are centered on."""
        return conditions_mid(self.tail, self.head, preset, self.kind, self.variant)


@dataclass
class ResonanceLayout:
    n: int
    contributors: list[str]
    row: int
    bar_x: tuple[float, float]  # stroke centers of the opening and closing brackets' bars
    arrows: list[tuple[Point, Point, Group]] = field(default_factory=list)  # double-headed: tail tip, head tip
    group: Group | None = None  # bracket-step-<n>, with -left and -right


@dataclass
class Unit:
    """An unbreakable piece of a row: the opening reactants, or a step's arrow with its
    conditions and products (a resonance bracket that follows included), plus the next step's
    extra reactants."""
    kind: str  # "reactants" | "step" | "products" (a mechanism step's products, which may start a row)
    step: int  # the step whose arrow leads it (reactants: the step they start, 0 for a bracket)
    row: int
    x0: float  # cursor before the unit (its leading gap starts here)
    lead: float  # leading gap, dropped when the unit starts a row
    x1: float = 0.0  # cursor after its last ink
    new_row: bool = False
    glue: bool = False  # stays on the row of the unit before it (a mechanism step's arrow)
    species: list[str] = field(default_factory=list)
    over: list[str] = field(default_factory=list)  # reactants drawn above a mechanism step's arrow
    pluses: list[TextRun] = field(default_factory=list)
    items: list = field(default_factory=list)

    @property
    def bracket(self) -> bool:
        return any(getattr(it, "id", "").startswith("bracket-step-") for it in self.items)


@dataclass
class SchemeLayout:
    root: Group
    preset: Preset
    species: dict[str, Placed] = field(default_factory=dict)
    steps: list[StepLayout] = field(default_factory=list)
    resonance: list[ResonanceLayout] = field(default_factory=list)
    pluses: list[TextRun] = field(default_factory=list)
    numbers: dict[str, list[TextRun]] = field(default_factory=dict)
    curved: list[MechArrow] = field(default_factory=list)
    checks: list[electrons.StepCheck] = field(default_factory=list)

    @property
    def width(self) -> float:
        b = self.root.bbox()
        return b[2] - b[0]


def m1_drawers(spec: Spec, base_dir: Path) -> list[Drawer]:
    """M1: one Drawer per species, in spec order."""
    preset = load_preset(spec.preset)
    return [Drawer(mol, preset, show_lone_pairs={a.idx for a in mol.atoms if a.map in sp.show_lone_pairs})
            for sp, mol in zip(spec.species.values(), molecules(spec, base_dir))]


def m1_row(spec: Spec, base_dir: Path, name: str) -> Group:
    """M1: species side by side, centered vertically, plus_gap between ink boxes."""
    preset = load_preset(spec.preset)
    groups, x = [], 0.0
    for d in m1_drawers(spec, base_dir):
        g = d.draw(d.mol.name)
        x0, y0, x1, y1 = g.bbox()
        g.translate(x - x0, -(y0 + y1) / 2)
        x += (x1 - x0) + preset.scheme.plus_gap
        groups.append(g)
    return Group(name, groups)


def _ring_centers(mol: Molecule) -> list[Point]:
    return [(sum(mol.atoms[i].x for i in r) / len(r), sum(mol.atoms[i].y for i in r) / len(r)) for r in mol.rings]


def _anchor_y(p: Placed, use_rings: bool) -> float:
    b = p.group.ink_bbox()
    cy = (b[1] + b[3]) / 2
    if use_rings and p.ring_centers:
        return min((c[1] for c in p.ring_centers), key=lambda y: abs(y - cy))
    return cy


def face_plus_width(preset: Preset) -> float:
    """Ink width of a '+' at caption size."""
    b = richtext.face_for("regular").ink_bounds("+", preset.text.caption_size)
    return b[2] - b[0]


def _plus(x: float, axis_y: float, preset: Preset, i: int) -> tuple[TextRun, float]:
    """A '+' (regular, caption size) with its ink left at x and ink center center_above_axis
    above the axis. Returns the run and its ink right edge."""
    size = preset.text.caption_size
    face = richtext.face_for("regular")
    b = face.ink_bounds("+", size)
    cy = axis_y - preset.scheme.center_above_axis
    run = TextRun(x - b[0], cy + (b[1] + b[3]) / 2, size, "+", face, source=f"plus:{i}")
    return run, x - b[0] + b[2]


def _conditions(step_lines: tuple[list[str], list[str]], staged: bool, preset: Preset, n: int):
    """Lay out above/below lines at x = 0 (left origin) relative to their own frame; returns
    (above, below, needed tail->notch length)."""
    size = preset.text.caption_size
    src = f"conditions:{n}"
    above = [richtext.layout_line(t, 0.0, 0.0, size, preset, src) for t in step_lines[0]]
    below = [richtext.layout_line(t, 0.0, 0.0, size, preset, src) for t in step_lines[1]]
    s = preset.scheme
    lines = above + below
    if not lines:
        return above, below, 0.0
    if staged:  # each line's origin sits at tail + indent - its ink offset
        need = max(s.multistep_indent - ln.ink()[0] + ln.advance for ln in lines) + s.conditions_pad
    else:
        need = max(ln.advance for ln in lines) + 2 * s.conditions_pad
    return above, below, need


def _place_conditions(st: StepLayout, preset: Preset) -> None:
    """Stack the step's lines around the arrow (module docstring)."""
    s = preset.scheme
    size = preset.text.caption_size
    axis, off = st.tail[1], st.shaft_offset(preset)
    mid = st.mid_x(preset)

    def x_for(ln):
        if st.staged:
            return st.tail[0] + s.multistep_indent - ln.ink()[0]
        return mid - ln.advance / 2

    base = axis - off - s.conditions_above_baseline
    for k in range(len(st.above) - 1, -1, -1):
        ln = st.above[k]
        if k < len(st.above) - 1:
            base -= richtext.line_pitch(ln, st.above[k + 1], size, preset)
        ln.shift(x_for(ln), base)
    base = axis + off + s.conditions_below_baseline
    for k, ln in enumerate(st.below):
        if k:
            base += richtext.line_pitch(st.below[k - 1], ln, size, preset)
        ln.shift(x_for(ln), base)


def _place(spec: Spec, out: SchemeLayout, breaks: set[tuple[str, int]], rings: dict[int, bool] | None) -> list[Unit]:
    """One placement pass in row-local coordinates (each row's axis at y = 0, its first unit at
    x = 0): units in order, with a new row at each new_row step and before each step in
    `breaks`. `rings` (row -> every structure has a ring) picks the vertical anchor; None in
    the measuring pass, where only x matters."""
    preset = out.preset
    s = preset.scheme
    out.steps, out.pluses, out.resonance = [], [], []
    units: list[Unit] = []
    x, row = 0.0, 0

    def unit(kind: str, step: int, lead: float, new_row: bool = False, glue: bool = False) -> None:
        nonlocal x, row
        if units:
            units[-1].x1 = x
            if new_row or (kind, step) in breaks:
                row += 1
                x = 0.0 if new_row else -lead  # a continuation row starts at 0 with its first ink
        units.append(Unit(kind, step, row, x, lead, new_row=new_row, glue=glue))

    def plus() -> None:
        nonlocal x
        run, right = _plus(x + s.plus_gap, 0.0, preset, len(out.pluses) + 1)
        out.pluses.append(run)
        units[-1].pluses.append(run)
        units[-1].items.append(Group(f"plus-{len(out.pluses)}", [run]))
        x = right + s.plus_gap

    def put(nm: str, own: ArrowSet | None = None) -> None:
        """Species nm from x. A resonance contributor's own curved arrows (`own`, planned on it
        alone) count as its ink here, so the double-headed arrows and brackets clear them."""
        nonlocal x
        p = out.species[nm]
        use_rings = rings[row] if rings is not None else False
        p.translate(0.0, -s.center_above_axis - _anchor_y(p, use_rings))
        b = p.group.ink_bbox()
        if own is not None:
            b = union_box([b] + [a.group.ink_bbox() for a in curved_arrows(out, spec, own, [nm], [])])
        p.translate(x - b[0], 0.0)
        units[-1].items.append(p.group)
        units[-1].species.append(nm)
        x += b[2] - b[0]

    def put_species(names: list[str]) -> None:
        for k, nm in enumerate(names):
            if k:
                plus()
            put(nm)

    def put_bracket(n: int) -> None:
        """Step n's contributors from x (the opening bracket's ink left); the brackets' heights
        come after the curved arrows (_finish_row)."""
        nonlocal x
        st, r, b = spec.steps[n - 1], preset.resonance, preset.bracket
        left = x + b.stroke_width / 2
        x = left + b.gap_x
        res = ResonanceLayout(n, list(st.resonance), row, (left, 0.0))
        own = {a.contributor: a for a in spec.arrow_sets() if a.step == n}
        for k, nm in enumerate(st.resonance):
            if k:
                tail = (x + r.arrow_gap, 0.0)
                head = (tail[0] + r.arrow_length, 0.0)
                g = resonance_arrow(tail, head, preset, f"resonance-arrow-step-{n}-{k}", spec.arrow_head)
                res.arrows.append((tail, head, g))
                units[-1].items.append(g)
                x = head[0] + r.arrow_gap
            put(nm, own.get(k))
        res.bar_x = (left, x + b.gap_x)
        x = res.bar_x[1] + b.stroke_width / 2
        res.group = Group(f"bracket-step-{n}")
        units[-1].items.append(res.group)
        out.resonance.append(res)

    def put_over(names: list[str], st: StepLayout) -> None:
        """Reactants a mechanism step adds, side by side (with '+') above its arrow: centered on
        tail -> notch, their ink bottom conditions_above_baseline above the arrow or above its
        conditions (unmeasured rule, acs-style.md §4.3)."""
        top = min([ln.ink()[1] for ln in st.above] + [-st.shaft_offset(preset)])
        base = top - s.conditions_above_baseline
        boxes = [out.species[n].group.ink_bbox() for n in names]
        width = sum(b[2] - b[0] for b in boxes) + (len(names) - 1) * (2 * s.plus_gap + face_plus_width(preset))
        cx = st.mid_x(preset) - width / 2
        mid_y = base - max(b[3] - b[1] for b in boxes) / 2
        for k, (n, b) in enumerate(zip(names, boxes)):
            if k:
                run, right = _plus(cx + s.plus_gap, mid_y + s.center_above_axis, preset, len(out.pluses) + 1)
                out.pluses.append(run)
                units[-1].pluses.append(run)
                units[-1].items.append(Group(f"plus-{len(out.pluses)}", [run]))
                cx = right + s.plus_gap
            p = out.species[n]
            p.translate(cx - b[0], base - b[3])
            units[-1].items.append(p.group)
            units[-1].over.append(n)
            cx += b[2] - b[0]

    after_bracket = False
    for i, st in enumerate(spec.steps, 1):
        if isinstance(st, ResonanceStep):
            if i == 1:
                unit("reactants", 0, 0.0)
                put_bracket(1)
                after_bracket = True
            continue  # otherwise drawn with the step before it
        # a mechanism step's arrow stays with its reactants (its curved arrows are drawn on them),
        # and the reactants it adds (or its `over`) sit above it
        over = spec.over_names(i)
        if i == 1 or st.new_row:
            unit("reactants", i, 0.0, new_row=i > 1)
            put_species([n for n in st.sources if n not in over])
        else:
            extra = [n for n in st.sources[spec.kept(i):] if n not in over]
            if extra:  # more reactants join the previous products: "+ X" (after a bracket, no plus)
                if after_bracket:
                    x += s.plus_gap
                else:
                    plus()
                put_species(extra)
        unit("step", i, s.arrow_gap, glue=bool(st.arrows))
        staged = bool(st.stages)
        above, below, need = _conditions(st.lines(), staged, preset, i)
        boxes = [out.species[n].group.ink_bbox() for n in over]
        if over:
            plus_w = face_plus_width(preset)
            need = max(need, sum(b[2] - b[0] for b in boxes) + (len(over) - 1) * (2 * s.plus_gap + plus_w)
                       + 2 * s.conditions_pad)
        length = max(s.arrow_min_length, need + heads_along(st.arrow_type) * preset.head(spec.arrow_head).head_inset)
        tail = (x + s.arrow_gap, 0.0)
        head = (tail[0] + length, 0.0)
        g = step_arrow(st.arrow_type, tail, head, preset, f"arrow-step-{i}", spec.arrow_head)
        step = StepLayout(i, tail, head, above, below, staged, row, g, st.arrow_type, spec.arrow_head)
        _place_conditions(step, preset)
        out.steps.append(step)
        units[-1].items.append(g)
        sides = [Group(f"conditions-step-{i}-{side}", [r for ln in lines for r in ln.runs])
                 for side, lines in (("above", above), ("below", below)) if lines]
        if sides:
            units[-1].items.append(Group(f"conditions-step-{i}", sides))
        if over:
            put_over(over, step)
        x = head[0]
        if st.arrows:
            unit("products", i, s.arrow_gap)
        x += s.arrow_gap
        if i < len(spec.steps) and isinstance(spec.steps[i], ResonanceStep):
            put_bracket(i + 1)  # opens with the first product; the other drawn products follow the bracket ("+ HX")
            after_bracket = len(st.drawn_products) == 1
            if not after_bracket:
                plus()
                put_species(st.drawn_products[1:])
        else:
            put_species(st.drawn_products)
            after_bracket = False
    units[-1].x1 = x
    return units


def _breaks(units: list[Unit], width: float) -> set[tuple[str, int]]:
    """Units that start a continuation row, as (kind, step): greedy, each row as full as
    `width` allows. A glued unit (a mechanism step's arrow) never starts a row."""
    blocks: list[list[Unit]] = []
    for u in units:
        if u.glue and blocks:
            blocks[-1].append(u)
        else:
            blocks.append([u])
    out: set[tuple[str, int]] = set()
    start = units[0]
    for b in blocks:
        first, last = b[0], b[-1]
        if first.new_row or first is units[0]:
            start = first
        elif last.x1 - (start.x0 + start.lead) > width + 1e-9:
            out.add((first.kind, first.step))
            start = first
        if first is start and last.x1 - (first.x0 + first.lead) > width + 1e-9:
            what = ("the resonance bracket of step 1" if first.kind == "reactants" and first.step == 0 else
                    f"the reactants of step {first.step}" if first.kind == "reactants" else
                    f"the products of step {first.step}" if first.kind == "products" else
                    f"step {first.step} (its arrow, conditions and products)")
            if len(b) > 1 or first.bracket:
                what += " with what must stay on its row (a bracket, a mechanism step's arrow)"
            raise ValueError(f"layout: {what} is "
                             f"{last.x1 - first.x0 - first.lead:.1f} pt wide, more than layout.width {width:g} pt, and a row "
                             "breaks only between steps; set a larger layout.width")
    return out


def build_scheme(spec: Spec, base_dir: Path, name: str = "scheme") -> SchemeLayout:
    preset = load_preset(spec.preset)
    s = preset.scheme
    mols = scheme_molecules(spec, base_dir)
    out = SchemeLayout(Group(name), preset)
    out.checks = electrons.check_spec(spec, base_dir)
    show = electrons.lone_pair_atoms(spec, out.checks)
    tails = tail_directions(spec, mols)
    for n, mol in mols.items():
        d = Drawer(mol, preset, show_lone_pairs={a.idx for a in mol.atoms if a.map in show[n]}, tails=tails[n])
        out.species[n] = Placed(n, d.draw(n), mol, _ring_centers(mol), drawer=d)

    # pass 1 measures the units on unbroken rows; pass 2 places them on the wrapped rows
    width = spec.layout.width or s.wrap_width
    breaks = _breaks(_place(spec, out, set(), None), width)
    units = _place(spec, out, breaks, None)
    rows: dict[int, list[str]] = {}
    for u in units:
        rows.setdefault(u.row, []).extend(u.species)
    rings = {r: all(out.species[n].ring_centers for n in names) for r, names in rows.items()}
    for u in units:
        rows[u.row].extend(u.over)
    units = _place(spec, out, breaks, rings)

    finished = [_finish_row(spec, out, [u for u in units if u.row == r], r, len(rows)) for r in sorted(rows)]
    groups = [g for g, _ in finished]
    if len(groups) == 1:
        out.root.children = groups[0].children
        return out
    bottom = groups[0].ink_bbox()[3]
    for r, g in enumerate(groups[1:], 1):  # rows stacked with row_gap between their full ink
        b = g.ink_bbox()
        dy = bottom + s.row_gap - b[1]
        g.translate(0.0, dy)
        _shift_row_records(out, r, rows[r], finished[r][1], dy)
        bottom = b[3] + dy
    out.root.children = groups
    return out


def _finish_row(spec: Spec, out: SchemeLayout, units: list[Unit], r: int, nrows: int) -> tuple[Group, list[MechArrow]]:
    """Curved arrows (on the arrow sets whose species are in this row: a step's reactants can
    end the row above its arrow), brackets and compound numbers of row r; its items as one
    group, and its curved arrows."""
    preset = out.preset
    names = [n for u in units for n in u.species + u.over]
    pluses = [p for u in units for p in u.pluses]
    items = [it for u in units for it in u.items]
    obstacles = {g.id: _arrow_parts(g) for g in items if g.id.startswith(("arrow-step-", "resonance-arrow-step-"))}
    curved: dict[int, list[MechArrow]] = {}
    for aset in spec.arrow_sets():
        if aset.sources[0] in names:
            curved.setdefault(aset.step, []).extend(curved_arrows(out, spec, aset, names, pluses, obstacles))
    for res in (x for x in out.resonance if x.row == r):
        box = union_box([out.species[n].group.ink_bbox() for n in res.contributors]
                        + [a.group.ink_bbox() for a in curved.get(res.n, [])])
        res.group.children = list(bracket_pair((res.bar_x[0] + preset.bracket.gap_x, box[1],
                                                res.bar_x[1] - preset.bracket.gap_x, box[3]),
                                               preset, f"bracket-step-{res.n}"))
    numbered = {n: str(spec.numbers[n]) for n in names if n in (spec.numbers or {})}
    if numbered:
        size = preset.text.caption_size
        lines = {n: richtext.layout_line(t, 0.0, 0.0, size, preset, f"number:{n}", style="bold")
                 for n, t in numbered.items()}
        bottom = max([out.species[n].group.ink_bbox()[3] for n in names]
                     + [res.group.ink_bbox()[3] for res in out.resonance if res.row == r])
        rise = max(-ln.ink()[1] for ln in lines.values())  # tallest ink above the baseline
        base = bottom + preset.scheme.number_gap + rise
        groups = []
        for n, ln in lines.items():
            b = out.species[n].group.ink_bbox()
            ink = ln.ink()
            ln.shift((b[0] + b[2]) / 2 - (ink[0] + ink[2]) / 2, base)
            out.numbers[n] = ln.runs
            groups.append(Group(f"number-{n}", ln.runs))
        items.append(Group("numbers" if nrows == 1 else f"numbers-row-{r + 1}", groups))
    for n, arrows in sorted(curved.items()):
        out.curved += arrows
        items.append(Group(f"curved-arrows-step-{n}", [a.group for a in arrows]))
    return Group(f"row-{r + 1}", items), [a for arrows in curved.values() for a in arrows]


def _arrow_parts(g: Group) -> Parts:
    """A reaction or resonance arrow as an obstacle for curved arrows: its shaft and the
    outlines of its heads."""
    out = Parts()
    for shaft in g.leaves():
        for leaf in [shaft, *shaft.head_leaves()]:
            pts = leaf.points
            closed = not isinstance(leaf, Polyline)
            hw = leaf.width / 2 if isinstance(leaf, Polyline) else 0.0
            out.bonds += [(a, b, hw, g.id) for a, b in zip(pts, pts[1:] + (pts[:1] if closed else []))]
    return out


def _shift_row_records(out: SchemeLayout, r: int, names: list[str], curved: list[MechArrow], dy: float) -> None:
    """After a row's group moved down by dy: move the records that hold coordinates."""
    for n in names:
        p = out.species[n]
        p.offset = (p.offset[0], p.offset[1] + dy)
        p.ring_centers = [(x, y + dy) for x, y in p.ring_centers]
    for st in out.steps:
        if st.row == r:
            st.tail, st.head = (st.tail[0], st.tail[1] + dy), (st.head[0], st.head[1] + dy)
    for res in out.resonance:
        if res.row == r:
            res.arrows = [((t[0], t[1] + dy), (h[0], h[1] + dy), g) for t, h, g in res.arrows]
    for a in curved:
        a.shift(0.0, dy)


def tail_directions(spec: Spec, mols: dict[str, Molecule]) -> dict[str, dict[int, dict[str, Point]]]:
    """species -> atom index -> {"pair" | "radical": direction to its arrow's target}, for
    lp(n) and electron(n) tails, so that dot set is drawn facing the target. A target in
    another species lies along the row: right if drawn later, left if earlier."""
    pos = {n: k for k, n in enumerate(spec.order())}
    out: dict[str, dict] = {n: {} for n in mols}
    for aset in spec.arrow_sets():
        where = {a.map: (n, a) for n in aset.sources for a in mols[n].atoms if a.map}
        for arr in aset.arrows:
            src, dst = arr.source, arr.target
            if src.kind not in ("lp", "electron") or src.atoms[0] not in where:
                continue
            n, a = where[src.atoms[0]]
            m = next((x for x in dst.atoms if x != src.atoms[0]), dst.atoms[0])
            if m not in where:
                continue
            tn, ta = where[m]
            d = unit(sub(ta.p, a.p)) if tn == n and ta.idx != a.idx else ((1.0, 0.0) if pos[tn] >= pos[n] else (-1.0, 0.0))
            out[n].setdefault(a.idx, {}).setdefault("pair" if src.kind == "lp" else "radical", d)
    return out


def _moved(ds: DotSet, d: Point) -> DotSet:
    return DotSet(ds.atom, ds.kind, add(ds.center, d), ds.direction, ds.side, [add(p, d) for p in ds.dots])


def curved_arrows(lay: SchemeLayout, spec: Spec, aset: ArrowSet, row: list[str] | None = None,
                  pluses: list[TextRun] | None = None, obstacles: dict[str, Parts] | None = None) -> list[MechArrow]:
    """One arrow set's curved arrows on its placed species (mechanism.plan_step), avoiding the
    species and plus signs of its row (default: all) and `obstacles` (its reaction and
    resonance arrows)."""
    atoms, bonds, centroids = {}, set(), {}
    for n in aset.sources:
        p = lay.species[n]
        mol, d, off = p.mol, p.drawer, p.offset
        sets = d.dot_sets()
        for a in mol.atoms:
            if a.map is None:
                continue
            box = None
            if a.idx in d.labels:
                b = d.labels[a.idx].elem_ink_box()
                box = (b[0] + off[0], b[1] + off[1], b[2] + off[0], b[3] + off[1])
            atoms[a.map] = AtomRef(n, a.idx, add(a.p, off), box, [_moved(ds, off) for ds in sets.get(a.idx, [])])
        for b in mol.bonds:
            ma, mb = mol.atoms[b.a].map, mol.atoms[b.b].map
            if ma and mb:
                bonds.add(frozenset((ma, mb)))
        centroids[n] = (sum(a.x for a in mol.atoms) / len(mol.atoms) + off[0],
                        sum(a.y for a in mol.atoms) / len(mol.atoms) + off[1])
    overrides = [spec.arrow_override(aset.step, aset.offset + k) for k in range(len(aset.arrows))]
    parts = {n: species_parts(lay.species[n].group) for n in (row if row is not None else lay.species)}
    parts |= obstacles or {}
    return plan_step(aset.step, list(aset.arrows), overrides, lay.preset, atoms, bonds, centroids, parts,
                     [r.bbox() for r in (pluses if pluses is not None else lay.pluses)], aset.offset)


def render_spec(path: Path) -> Group:
    """Scene for a spec file: M1 side-by-side species, or an M2 scheme when it has steps."""
    spec = load_spec(path)
    if not spec.steps:
        return m1_row(spec, path.parent, path.stem)
    return build_scheme(spec, path.parent, path.stem).root
