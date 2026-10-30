"""Catalytic cycles (acs-style.md §15): the circle solver, arc arrows, tails and the spec."""

import math
from pathlib import Path

import pytest
from pydantic import ValidationError

from acs_chem_svg.arrows import arc_arrow, arc_equilibrium_arrow
from acs_chem_svg.cycle import _box_at, _gap, arc_ends, build_cycle, check_cycle, solve, tail_items
from acs_chem_svg.geometry import arc_point, length, point_box_dist, sub
from acs_chem_svg.layout import render_spec
from acs_chem_svg.preset import load_preset
from acs_chem_svg.spec import CycleTail, Spec, load_spec

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = sorted((ROOT / "examples").glob("m6_*_cycle_x*.yaml"))
P = load_preset("acs1996")


def _cycle(steps, species=None, **kw):
    names = {s["from"] for s in steps} | {s["to"] for s in steps}
    sp = species or {n: {"smiles": "c1ccccc1"} for n in names}
    return Spec.model_validate({"species": sp, "cycle": {"steps": steps, **kw}})


def _ring(*names):
    return [{"from": a, "to": b} for a, b in zip(names, names[1:] + names[:1])]


@pytest.mark.parametrize("path", EXAMPLES, ids=lambda p: p.stem)
def test_examples_equal_free_arcs_start_on_top_and_clear(path):
    spec = load_spec(path)
    lay = build_cycle(spec, path.parent, path.stem)
    lengths = [lay.radius * abs(a.a1 - a.a0) for a in lay.arcs]
    assert max(lengths) - min(lengths) < 0.05
    assert lay.free == pytest.approx(lengths[0], abs=0.05)
    assert math.degrees(lay.theta[spec.cycle.members[0]]) == pytest.approx(90.0, abs=1e-6)
    c = P.cycle.clear
    for a in lay.arcs:  # no species within `clear` of any arc (its own ends sit exactly at `clear`)
        for j in range(101):
            q = arc_point(lay.center, lay.radius, a.a0 + (a.a1 - a.a0) * j / 100)
            for n, p in lay.species.items():
                if n in spec.cycle.members:
                    assert point_box_dist(q, p.group.ink_bbox()) > c - 0.05, (a.n, n)
    assert not [t for t in lay.tails if t.clash]
    assert check_cycle(lay) == []


def test_example_groups():
    root = render_spec(ROOT / "examples" / "m6_suzuki_cycle_x1.0.yaml")
    ids = {g.id for g in root.children}
    assert {"cycle-arc-1", "cycle-arc-2", "cycle-arc-3", "cycle-in-1-1", "cycle-in-2-1", "cycle-out-2-1",
            "cycle-out-3-1", "cycle-label-1", "cycle-center", "numbers", "species-pd0"} <= ids


def _radius(halves, need=40.0, labels=None):
    return solve(halves, need, P, -1.0, labels or [None] * len(halves))[0]


def test_radius_grows_with_count_and_size():
    small = [(10.0, 8.0)]
    rs = [_radius(small * n) for n in range(3, 9)]
    assert all(a < b for a, b in zip(rs, rs[1:]))
    assert _radius([(10.0, 8.0)] * 4) < _radius([(20.0, 8.0)] * 4) < _radius([(20.0, 16.0)] * 4)


def test_wide_species_takes_more_arc_at_the_top():
    """A wide box at 12 o'clock spans more of the circle than the same box at 3 o'clock would."""
    halves = [(40.0, 6.0), (6.0, 6.0), (6.0, 6.0), (6.0, 6.0)]
    R, (th, al, free) = solve(halves, 40.0, P, -1.0, [None] * 4)
    assert th[0] == pytest.approx(math.pi / 2)
    assert sum(al[0]) > 2 * sum(al[1])
    ends = arc_ends(th, al, -1.0)
    assert [R * abs(a1 - a0) for a0, a1 in ends] == pytest.approx([R * free] * 4, abs=0.01)


def test_labels_grow_the_radius_and_stay_clear():
    halves = [(12.0, 8.0)] * 3
    labels = ["oxidative addition", "transmetalation", "reductive elimination"]
    R0, R1 = _radius(halves), _radius(halves, labels=labels)
    assert R1 > R0
    R, (th, al, _) = solve(halves, 40.0, P, -1.0, labels)
    from acs_chem_svg.cycle import _label_line
    boxes = [_label_line(t, R, a0, a1, P, k).ink() for k, (t, (a0, a1)) in enumerate(zip(labels, arc_ends(th, al, -1.0)))]
    species = [_box_at(h, arc_point((0.0, 0.0), R, t)) for t, h in zip(th, halves)]
    assert all(_gap(a, b) >= P.cycle.label_clear - 1e-6 for a in boxes for b in species)
    for b in boxes:  # every corner label_gap inside the circle
        assert max(math.hypot(x, y) for x in (b[0], b[2]) for y in (b[1], b[3])) <= R - P.cycle.label_gap + 1e-3


def test_tail_text_does_not_change_the_radius(tmp_path):
    def radius(text):
        spec = _cycle([{"from": "a", "to": "b", "in": text}, {"from": "b", "to": "c"}, {"from": "c", "to": "a"}])
        return build_cycle(spec, tmp_path, "x").radius
    assert radius("X") == pytest.approx(radius("a much longer reagent, base, solvent"), abs=1e-6)


def test_arc_arrow_tip_notch_and_marker():
    c, r = (5.0, 7.0), 30.0
    g = arc_arrow(c, r, math.radians(200), math.radians(320), P, "x", "standard")
    shaft = g.children[0]
    head = shaft.head_leaves()[0]
    tip = arc_point(c, r, math.radians(320))
    assert min(length(sub(q, tip)) for q in head.points) < 1e-3  # the turn is rounded to 0.001° (marker id)
    notch = shaft.points[-1]
    assert length(sub(notch, c)) == pytest.approx(r, abs=1e-9)
    assert length(sub(tip, notch)) == pytest.approx(P.head("standard").head_inset, abs=1e-9)
    assert "-arc" in shaft.marker_end.id
    with pytest.raises(ValueError, match="shorter than its head"):
        arc_arrow(c, r, 0.0, 0.1, P, "x")


@pytest.mark.parametrize("a1, outer", [(math.radians(100), False), (math.radians(-100), True)])
def test_arc_equilibrium_forward_shaft_left_of_travel(a1, outer):
    """Counter-clockwise travel (a1 > a0): the forward shaft is inside; clockwise: outside."""
    c, r = (0.0, 0.0), 40.0
    fwd, back = arc_equilibrium_arrow(c, r, 0.0, a1, P, "eq").children
    rf, rb = length(sub(fwd.points[0], c)), length(sub(back.points[0], c))
    assert (rf > rb) is outer
    assert abs(rf - rb) == pytest.approx(P.equilibrium.shaft_spacing, abs=1e-9)


def test_tail_items_merge_text_keep_species():
    tails = [CycleTail(text="PhCHO"), CycleTail(text="H+"), CycleTail(species="s"), CycleTail(text="base")]
    assert tail_items(tails) == [("PhCHO, H+", None), (None, "s"), ("base", None)]


def test_start_and_direction(tmp_path):
    spec = _cycle(_ring("a", "b", "c", "d"), start="c", direction="counterclockwise")
    assert spec.cycle.members == ["c", "d", "a", "b"]
    assert [s.from_ for s in spec.cycle.ordered_steps()] == ["c", "d", "a", "b"]
    lay = build_cycle(spec, tmp_path, "x")
    assert math.degrees(lay.theta["c"]) == pytest.approx(90.0)
    assert lay.theta["d"] > lay.theta["c"]  # counter-clockwise: the angles grow


@pytest.mark.parametrize("cycle, species, message", [
    ({"steps": [{"from": "a", "to": "b"}, {"from": "b", "to": "c"}, {"from": "c", "to": "b"}]}, None,
     "must be the `from` of the first step"),
    ({"steps": [{"from": "a", "to": "b"}, {"from": "b", "to": "a"}]}, None, "at least 3"),
    ({"steps": _ring("a", "b", "c"), "start": "z"}, None, "start 'z' is not one of"),
    ({"steps": _ring("a", "b") + _ring("a", "c")}, None, "twice round the circle"),
    ({"steps": [{"from": "a", "to": "b", "in": {"species": "c"}}, {"from": "b", "to": "c"}, {"from": "c", "to": "a"}]},
     None, "on the circle"),
    ({"steps": _ring("a", "b", "c")}, {"a": {}, "b": {}, "c": {}, "d": {}}, "not used in the cycle: d"),
    ({"steps": [{"from": "a", "to": "b", "in": {"text": "x", "species": "d"}}, {"from": "b", "to": "c"},
                {"from": "c", "to": "a"}]}, None, "exactly one of text or species"),
])
def test_cycle_schema_errors(cycle, species, message):
    names = species or ({"a", "b", "c", "d"} if "'d'" in str(cycle) else {"a", "b", "c"})
    sp = {n: {"smiles": "C"} for n in names}
    with pytest.raises(ValidationError, match=message):
        Spec.model_validate({"species": sp, "cycle": cycle})


def test_cycle_and_steps_are_exclusive():
    with pytest.raises(ValidationError, match="either steps"):
        Spec.model_validate({"species": {n: {"smiles": "C"} for n in "abc"},
                             "steps": [{"from": "a", "to": "b"}], "cycle": {"steps": _ring("a", "b", "c")}})
