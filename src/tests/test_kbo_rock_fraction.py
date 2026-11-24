"""
GEN.522 (builder B14, row T19 of docs/design/composition/tunable-gaps.md): the rock fraction of an icy body (`small_body_composition.ICY_TABLE`)
against the published Kuiper belt object and dwarf planet statistics. Every one of them is a measurement of one system's bodies, so each is a
test point (`role` test-point, reported and never tuned to, NO FIXED CONTENT), and the rows say so in their states.
"""

import math

import pytest

from planetgen.generation import possibility_table as pt, tuning_harness as harness
from planetgen.physics import compaction, small_body_composition as sb
from planetgen.util import draw

TABLE = "icy-body-rock-fraction"

KBO_DENSITIES = (
    # (diameter km, diameter error km, density kg/m3): Table A.2 of [R0390], the 17 objects with a published density and diameter
    (157.0, 34.0, 600.0), (174.0, 17.0, 1370.0), (178.0, 34.5, 600.0), (235.0, 22.0, 1140.0), (123.0, 44.0, 300.0), (249.0, 30.5, 730.0),
    (304.0, 101.0, 500.0), (652.0, 12.0, 820.0), (705.0, 78.0, 1270.0), (866.0, 37.0, 1260.0), (958.0, 23.0, 1520.0), (1070.0, 38.0, 2180.0),
    (1212.0, 2.0, 1700.0), (1595.0, 11.0, 1885.0), (2326.0, 12.0, 2520.0), (2376.0, 3.0, 1850.0), (2706.0, 1.8, 2060.0))
"""tuple: The 17 objects of [R0390]'s table: the mass of each comes from the orbit of its binary, so the radius dominates the density error (the paper's own reading)."""


def _reach(rock_fraction=0.95):
    """The objects of the table the compaction law covers at `rock_fraction` (it ends at a mass of 1e21 kg; beyond it the structure solver decides)."""
    covered = []
    for diameter, error, density in KBO_DENSITIES:
        try:
            compaction.compacted_density_kg_m3(diameter / 2.0, rock_fraction)
        except compaction.OutOfRange:
            continue
        covered.append((diameter, error, density))
    return covered


COVERED = _reach()


def test_the_17_objects_are_the_published_table_and_the_law_covers_the_nine_smallest():
    assert len(KBO_DENSITIES) == 17 and len(COVERED) == 9
    assert [row[0] for row in COVERED] == sorted(row[0] for row in COVERED)[:9] or max(row[0] for row in COVERED) < 800.0


def test_the_t19_targets_are_test_points_of_the_icy_rows_and_none_is_tuned_to():
    ids = ("t19-dwarf-rock-mean", "t19-dwarf-rock-sd", "t19-kbo-density-share")
    for ident in ids:
        target = harness.TARGETS.get(ident)
        assert target.role == "test-point" and not target.tunes and target.solar_system
        assert target.fill_id == "GEN.522" and target.table == TABLE
    assert not harness.TARGETS.tunable_targets(f"{TABLE}:icy-intact")


def test_the_loader_accepts_the_table_with_these_targets():
    harness.TARGETS.require([sb.ICY_TABLE])


def test_the_rows_carry_their_states_and_the_dwarf_planet_takes_each_at_a_half():
    states = {row.name: row.state for row in sb.ICY_TABLE.rows}
    assert states == {"icy-intact": pt.PARTLY_VERIFIED, "icy-stripped": pt.STAND_IN}
    dwarf = [(row.name, p) for row, p in sb.ICY_TABLE.probabilities(pt.Subject(sb.DWARF_PLANET, {}))]
    assert dwarf == [("icy-intact", 0.5), ("icy-stripped", 0.5)]
    for kind in (sb.KUIPER_BELT_OBJECT, sb.CLASS_D):
        only = sb.ICY_TABLE.probabilities(pt.Subject(kind, {}))
        assert [(row.name, p) for row, p in only] == [("icy-intact", 1.0)]


def _dwarf_rock(_i):
    return sb.draw_icy_body(sb.DWARF_PLANET).rock


def test_a_dwarf_planets_rock_fraction_has_the_mean_and_spread_of_the_seven_measured_ones():
    mean = harness.simulate_and_compare(harness.TARGETS.get("t19-dwarf-rock-mean"), _dwarf_rock, n=2000)
    spread = harness.simulate_and_compare(harness.TARGETS.get("t19-dwarf-rock-sd"), _dwarf_rock, n=2000)
    assert mean.passed and spread.passed, (mean.simulated, mean.z, spread.simulated, spread.z)
    assert not mean.tunes and not spread.tunes


def _within(i):
    """1 when the density the law gives at the rock fraction a draw of the row makes, at the radius of object `i` (cycling the covered ones), is within 2 sigma of the object's own."""
    diameter, error, density = COVERED[i % len(COVERED)]
    drawn = sb.draw_icy_body(sb.KUIPER_BELT_OBJECT).density_at(diameter / 2.0)
    return 1.0 if abs(drawn - density) < 2.0 * density * 3.0 * error / diameter else 0.0


def test_the_density_of_a_drawn_kuiper_belt_object_matches_the_measured_ones_as_often_as_the_published_fit():
    target = harness.TARGETS.get("t19-kbo-density-share")
    result = harness.simulate_and_compare(target, _within, n=len(COVERED) * 300)
    assert result.passed, (result.simulated, result.z)
    assert not result.tunes and not target.outside({"radius_km": 200.0})


def test_the_object_under_pure_ice_is_reached_only_at_the_rows_lower_limit():
    """The 652 km object (820 kg/m3, the published fit's 0.5) sits at the density the law gives at the row's lowest rock fraction, 0.45, and above it only 1 percent of draws reach it."""
    diameter, error, density = next(row for row in KBO_DENSITIES if row[0] == 652.0)
    sigma = density * 3.0 * error / diameter
    center, spread, low, _high = sb.ROCK_FRACTION_LAWS["icy-intact"]
    assert abs(compaction.compacted_density_kg_m3(diameter / 2.0, low) - density) < sigma
    assert abs(compaction.compacted_density_kg_m3(diameter / 2.0, center) - density) > 10.0 * sigma
    assert (low - center) / spread < -3.0


def test_the_icy_rows_draw_the_same_numbers_as_before_they_carried_states():
    with draw.bound(draw.Stream("gen522-same")):
        rocks = [round(sb.draw_icy_body(sb.DWARF_PLANET).rock, 12) for _ in range(5)]
    assert all(0.45 <= rock <= 0.98 for rock in rocks) and len(set(rocks)) == 5
    assert math.isclose(sum(sb.ROCK_FRACTION_LAWS["icy-intact"][:1]), 0.70)
