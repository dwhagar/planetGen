# tests/test_oddity_plausibility.py

"""TEST.490 (docs/design/systems/weird-worlds.md section 5): the oddity plausibility harness (`oddity_plausibility_cli.py`). It runs on stub rows and
passes them, fails on each deliberately broken row at the invariant that row breaks, holds every shipped row to the general invariants and its own
limits, and checks the isolation rule: the weirdness dial at 0 changes no draw but the oddity's own."""

import pathlib

import pytest

from planetgen.physics import oddities
from planetgen.util import draw
from tests import oddity_plausibility_cli as harness
from tests.oddity_plausibility_cli import BROKEN_CAUGHT_BY, BROKEN_RING_ROWS, BROKEN_ROWS, RING_CAUGHT_BY, STUB_ROWS

PLANETGEN = pathlib.Path(__file__).resolve().parent.parent / "planetgen"
BODIES = 25
"""int: Bodies forced per row in the test run (the harness's own default, 200, is for the command line and the build of a row)."""
SHIPPED_BODIES = 30


@pytest.mark.parametrize("row", STUB_ROWS, ids=lambda row: row.id)
def test_a_stub_row_keeps_every_invariant_and_its_limits(row):
    result = harness.run_row(row, BODIES)
    assert result["gated"] == BODIES, result
    assert not result["general"] and not result["limits"], result["examples"]


@pytest.mark.parametrize("row", BROKEN_ROWS, ids=lambda row: row.id)
def test_a_broken_row_is_caught_by_the_invariant_it_breaks(row):
    result = harness.run_row(row, 8)
    assert result["gated"] > 0, f"{row.id}: no body met the gate"
    assert BROKEN_CAUGHT_BY[row.id] in result["general"], (row.id, result["general"], result["examples"])


@pytest.mark.parametrize("row", BROKEN_RING_ROWS, ids=lambda row: row.id)
def test_a_broken_ring_row_is_caught_by_the_ring_invariant_it_breaks(row):
    """TEST.494: ring mass, the Roche limit, the optical depth as the formula's output, no ring inside 0.1 AU, the ring plane in the equator."""
    result = harness.run_row(row, 8)
    assert result["gated"] > 0, f"{row.id}: no body met the gate"
    assert RING_CAUGHT_BY[row.id] in result["general"], (row.id, result["general"], result["examples"])


def test_every_ring_invariant_is_broken_by_a_row():
    assert set(RING_CAUGHT_BY) == {row.id for row in BROKEN_RING_ROWS}
    assert len(set(RING_CAUGHT_BY.values())) == 5


def test_a_ring_beyond_the_roche_limit_passes_when_the_body_carries_a_resonance_row():
    """R4 (Quaoar, Haumea): a ring beyond the Roche limit is allowed for the body that names its resonance, and for no other."""
    row = [r for r in BROKEN_RING_ROWS if r.id == "R92"][0]
    body, star = harness.force(row, 0)
    assert any(name == "ring beyond the Roche limit" for name, _ in harness.ring_failures(body, star))
    assert not any(name == "ring beyond the Roche limit" for name, _ in harness.ring_failures(body, star, markers=("R4",)))


def test_a_body_with_no_ring_holds_no_ring_field_and_a_half_set_ring_is_refused():
    body, star = harness.draw_body(STUB_ROWS[0], 0)
    for name in harness.rings.BODY_FIELDS:
        setattr(body, name, None)
    assert harness.ring_failures(body, star) == []
    body.ring_outer_km = 5.0e4
    assert [name for name, _ in harness.ring_failures(body, star)] == ["ring fields"]


def test_every_broken_row_names_its_invariant_and_every_invariant_is_broken_somewhere():
    assert set(BROKEN_CAUGHT_BY) == {row.id for row in BROKEN_ROWS}
    assert set(BROKEN_CAUGHT_BY.values()) == {"energy", "structure", "stability", "air", "class", "no NaN"}


def test_a_row_limit_that_fails_is_reported_apart_from_the_general_invariants():
    broken = [row for row in BROKEN_ROWS if row.id == "B2"][0]      # a negative mass breaks the limit mass_earth >= 0
    result = harness.run_row(broken, 5)
    assert result["limits"] == result["gated"] > 0


def test_a_gate_that_is_never_met_is_recorded_not_forced():
    never = harness._stub("N1", "never", "stub_spin_retrograde", "L1", "obliquity_deg", [{"quantity": "obliquity_deg", "min": 0.0}],
                          requires=[{"quantity": "mass_earth", "min": 1.0e9}])
    assert harness.draw_body(never, 0) is None and harness.force(never, 0) is None
    assert harness.run_row(never, 3)["gated"] == 0


def test_a_forced_body_is_the_same_body_on_a_second_run():
    row = STUB_ROWS[0]
    first = harness.quantities_of(*harness.force(row, 3))
    second = harness.quantities_of(*harness.force(row, 3))
    assert first == second and first["obliquity_deg"] >= 90.0


def test_the_quantities_of_a_body_are_declared_in_the_table_vocabulary():
    body, star = harness.draw_body(STUB_ROWS[0], 0)
    assert set(harness.quantities_of(body, star)) <= set(oddities.QUANTITIES)


def test_every_shipped_row_has_an_effect_the_harness_can_force_and_keeps_the_invariants():
    for row in oddities.ROWS:
        assert row.effect in harness.EFFECTS, f"{row.id}: the module that implements it registers its effect in tests.oddity_plausibility_cli.EFFECTS"
        result = harness.run_row(row, SHIPPED_BODIES)
        assert not result["general"] and not result["limits"], result["examples"]


# --- the isolation rule (design section 4.4 and 5) ------------------------------------------------------------------------------

def _next_run_draws(roll_first, dial):
    """The next five draws of the run's stream after a body is drawn, with the row rolled (or not) first at `dial`."""
    with draw.bound("oddity-isolation"):
        if roll_first:
            oddities.roll(STUB_ROWS[0], "body-1", dial=dial)
            oddities.roll(STUB_ROWS[1], "body-1", dial=dial)
        return [draw.random() for _ in range(5)]


def test_rolling_a_row_reads_nothing_from_the_run_stream_at_any_dial():
    plain = _next_run_draws(False, 1.0)
    assert _next_run_draws(True, 1.0) == plain
    assert _next_run_draws(True, 0.0) == plain


def test_the_dial_at_zero_turns_every_row_off_and_leaves_the_rows_draw_where_it_was():
    row = STUB_ROWS[0]
    assert oddities.roll(row, "body-1", dial=1.0) is True
    assert oddities.roll(row, "body-1", dial=0.0) is False
    assert oddities.stream(row.id, "body-1").random() == oddities.stream(row.id, "body-1").random()


def test_one_rows_draw_does_not_move_when_another_row_is_added_or_rolled():
    alone = oddities.stream("S1", "body-9").random()
    oddities.roll(STUB_ROWS[1], "body-9", dial=1.0)
    assert oddities.stream("S1", "body-9").random() == alone
    assert oddities.stream("S2", "body-9").random() != alone


def test_a_body_is_the_same_body_with_the_dial_at_zero_and_at_one():
    """The generator itself (no row wired yet): nothing that draws a body reads the dial, so the body drawn with it at 0 or at 1 is the same."""
    row = STUB_ROWS[0]
    from planetgen import tuning
    original = tuning.WEIRDNESS_DIAL
    try:
        bodies = []
        for dial in (0.0, 1.0):
            tuning.WEIRDNESS_DIAL = dial
            body, star = harness.draw_body(row, 5)
            bodies.append(harness.quantities_of(body, star))
    finally:
        tuning.WEIRDNESS_DIAL = original
    assert bodies[0] == bodies[1]


def test_no_life_score_reads_an_oddity():
    """Design general invariant 5: the life scores do not read an oddity name. The modules that score life never mention the table."""
    files = [PLANETGEN / "generation" / "life.py"] + sorted((PLANETGEN / "physics").glob("life_*.py")) + [PLANETGEN / "physics" / "habitability.py"]
    reading = [path.name for path in files if "oddit" in path.read_text(encoding="utf-8").lower()]
    assert not reading, f"these read the oddity table: {reading}"
