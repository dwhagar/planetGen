# tests/test_prevalence_oddities.py

"""TEST.492: the oddity census harness (`tests/audit/prevalence_oddities.py`): each flag trips on a deliberately broken input. The live census runs in `test_prevalence_systems_census.py`."""

import types

import pytest

from planetgen.physics import oddities
from tests.audit import prevalence_oddities as odd

ROW = {row.id: row for row in oddities.ROWS}


def _planet(tilt, carried=()):
    star = types.SimpleNamespace(mass=1.98847e30, luminosity=3.828e26, age=4.5, temperature=5772.0, radius=6.957e5)
    return types.SimpleNamespace(body_type="t", mass=5.972e24, radius=6371.0, distance=1.0, eccentricity=0.0167, star=star, is_moon=False, moons=[],
                                 axial_tilt_deg=tilt, oddities=tuple(carried))


def test_a_body_is_counted_gated_when_its_gate_holds_and_carried_when_it_has_the_id():
    system = types.SimpleNamespace(all_planets=[_planet(95.0, ("W1",)), _planet(95.0), _planet(10.0), types.SimpleNamespace(body_type="a")])
    assert odd.record_of(system)["W1"] == [2, 1]


def test_a_body_that_carries_a_row_its_gate_no_longer_holds_still_counts_as_gated():
    system = types.SimpleNamespace(all_planets=[_planet(10.0, ("W1",))])
    assert odd.record_of(system)["W1"] == [1, 1]


def test_the_summary_adds_the_records_that_finished_and_skips_the_abandoned_and_the_raised():
    records = [{"oddities": {"W1": [4, 1]}}, {"oddities": {"W1": [6, 3]}}, {"abandoned": True}, {"error": "x", "oddities": {"W1": [100, 100]}}]
    assert odd.summarize(records)["W1"] == {"gated": 10, "carried": 4, "share": 0.4}


def test_a_share_at_the_rate_raises_no_flag_and_one_far_from_it_does():
    rate = oddities.probability(ROW["W5"])
    inside = {"W5": {"gated": 400, "carried": round(400 * rate), "share": round(400 * rate) / 400}}
    assert odd.compare(inside, rows=[ROW["W5"]]) == []
    for share in (0.0, min(1.0, rate + 0.4)):
        far = {"W5": {"gated": 400, "carried": round(400 * share), "share": share}}
        assert [f.kind for f in odd.compare(far, rows=[ROW["W5"]])] in (["missing"], ["out-of-range"])


def test_a_row_expected_to_appear_that_never_does_is_missing():
    never = {"W5": {"gated": 400, "carried": 0, "share": 0.0}}
    assert [f.kind for f in odd.compare(never, rows=[ROW["W5"]])] == ["missing"]


def test_a_row_with_too_few_gated_bodies_is_reported_and_not_judged():
    few = {"W5": {"gated": odd.MIN_GATED - 1, "carried": 0, "share": 0.0}}
    assert odd.compare(few, rows=[ROW["W5"]]) == []
    assert "too few gated bodies" in odd.format_report(few, [], rows=[ROW["W5"]])


def test_a_row_that_can_be_declined_has_no_lower_bound_and_still_has_an_upper_one():
    low, high = odd.expected_range(ROW["W2"], 200)
    assert low == 0.0 and high >= oddities.probability(ROW["W2"]) and odd.SHARE_CAN_FALL_SHORT == {"W2", "W19", "W20"}
    none = {"W2": {"gated": 200, "carried": 0, "share": 0.0}}
    assert odd.compare(none, rows=[ROW["W2"]]) == []


def test_a_known_anomaly_names_the_item_that_explains_the_flag():
    off = {"W5": {"gated": 400, "carried": 400, "share": 1.0}}
    assert [f.known for f in odd.compare(off, rows=[ROW["W5"]], known={"W5": "GEN.9999"})] == ["GEN.9999"]
