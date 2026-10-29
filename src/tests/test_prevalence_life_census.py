# tests/test_prevalence_life_census.py

"""
TEST.435: the life census (part G of the prevalence suite). The life worlds of whole seeded systems are put on the highest rung they
reached and the shares are checked against Boss's prevalence targets for the life class system (`tests/data/prevalence_life_expectations.json`).
The live census fails on a flag no known anomaly explains.
"""

import copy
import re
import types

import pytest

from tests.audit import prevalence as census
from tests.audit import prevalence_life as life

ITEM = re.compile(r"^[A-Z]+\.\d+$")


def _body(**fields):
    base = dict(life_sentient=False, life_fauna=0, life_flora=0, life_stage="abiogenesis")
    base.update(fields)
    return types.SimpleNamespace(**base)


def _result(**rungs):
    life_worlds = sum(rungs.values())
    return {"systems": 10, "worlds": 100, "life": life_worlds, "rungs": rungs, "classes": {}, "abandoned": 0, "errors": [], "error_count": 0}


def _expectations(**changes):
    bands = {"stops_before_animals": (0.5, 1.0), "reaches_animals": (0.15, 0.35), "reaches_sentient": (0.0, 0.02),
             "sentient_of_animal_worlds": (0.0, 0.05)}
    measures = {name: {"min": low, "max": high, "kind": "stand-in", "source": "test"} for name, (low, high) in bands.items()}
    for name, extra in changes.items():
        measures[name].update(extra)
    return {"measures": measures, "anomalies": {"GEN.1": "test"}}


# --- The rungs --------------------------------------------------------------------------------------------------------------


def test_a_life_world_is_put_on_its_highest_rung():
    assert life.rung(_body()) == "origin"
    assert life.rung(_body(life_stage="photosynthesis", life_flora=2)) == "producers"
    assert life.rung(_body(life_stage="complex_cells")) == "complex"
    assert life.rung(_body(life_stage="multicellularity", life_flora=3)) == "multicellular"
    assert life.rung(_body(life_stage="multicellularity", life_flora=3, life_fauna=2)) == "animals"
    assert life.rung(_body(life_stage="technological_civilization", life_fauna=4, life_sentient=True)) == "sentient"


def test_the_shares_count_sentient_worlds_among_the_animal_ones():
    values, counts = life.summarize(_result(origin=30, producers=30, animals=30, sentient=10))
    assert values["stops_before_animals"] == 0.6 and values["reaches_animals"] == 0.4 and values["reaches_sentient"] == 0.1
    assert values["sentient_of_animal_worlds"] == 0.25 and counts["sentient_of_animal_worlds"] == 40


# --- The harness: each flag trips on a deliberately broken input ------------------------------------------------------------


def test_shares_inside_their_bands_raise_no_flag():
    assert life.compare(_result(origin=60, producers=10, animals=29, sentient=1), _expectations()) == []


def test_a_share_outside_its_band_is_flagged_both_ways():
    flags = life.compare(_result(origin=100, animals=800, sentient=100), _expectations())
    assert {f.measure for f in flags} >= {"stops_before_animals", "reaches_animals", "reaches_sentient", "sentient_of_animal_worlds"}
    assert [f.kind for f in life.compare(_result(origin=400, animals=600), _expectations()) if f.measure == "stops_before_animals"] == [
        "out-of-range"]


def test_no_animals_at_all_is_flagged_missing():
    flags = life.compare(_result(origin=100), _expectations(reaches_animals={"min": 0.5}))
    assert [f.kind for f in flags] == ["missing"]


def test_a_measure_that_rests_on_no_world_is_not_judged():
    assert life.compare(_result(origin=100), _expectations(reaches_animals={"min": 0.0})) == []


def test_a_run_with_no_life_world_is_flagged():
    assert [f.kind for f in life.compare(_result(), _expectations())] == ["missing"]


def test_the_band_is_widened_by_the_sampling_of_the_run():
    wide = life.allowed({"min": 0.15, "max": 0.35}, 25)
    narrow = life.allowed({"min": 0.15, "max": 0.35}, 2500)
    assert wide[0] < narrow[0] < 0.15 < 0.35 < narrow[1] < wide[1]


def test_a_known_anomaly_is_named_and_a_stale_tag_is_listed():
    spec = _expectations(reaches_animals={"anomaly": "GEN.1"})
    assert [f.known for f in life.compare(_result(origin=100), spec) if f.measure == "reaches_animals"] == ["GEN.1"]
    ok = life.compare(_result(origin=60, producers=10, animals=29, sentient=1), spec)
    assert life.stale_anomalies(spec, ok) == ["reaches_animals (GEN.1)"]


def test_a_system_that_raised_or_hung_is_flagged():
    result = _result(origin=60, animals=39, sentient=1)
    result.update(error_count=2, errors=["ValueError: x"], abandoned=1)
    assert {f.kind for f in life.compare(result, _expectations())} == {"error", "hang"}


def test_the_report_names_the_flags_the_shares_and_the_rungs():
    result = _result(origin=100)
    result["classes"] = {"None": 60, "3a": 40}
    report = life.format_report(result, life.compare(result, _expectations()))
    assert "reaches_animals" in report and "Highest rung: origin 100" in report and "3a 40" in report


# --- The stored bands --------------------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def expectations():
    return life.load_expectations()


def test_every_measure_has_a_band_a_kind_a_source_and_a_known_item(expectations):
    assert set(expectations["measures"]) == set(life.MEASURES)
    for measure, bounds in expectations["measures"].items():
        assert 0.0 <= bounds["min"] <= bounds["max"] <= 1.0 and bounds["kind"] == "stand-in" and bounds["source"], measure
        if "anomaly" in bounds:
            assert ITEM.match(bounds["anomaly"]) and bounds["anomaly"] in expectations["anomalies"], measure


def test_the_bands_are_boss_targets(expectations):
    bands = expectations["measures"]
    assert (bands["reaches_animals"]["min"], bands["reaches_animals"]["max"]) == (0.15, 0.35)
    assert bands["reaches_sentient"]["max"] == 0.02 and bands["sentient_of_animal_worlds"]["max"] == 0.05
    assert bands["stops_before_animals"]["min"] == 0.5


# --- The census ------------------------------------------------------------------------------------------------------------


def test_the_draw_repeats_and_counts_life_worlds_by_rung():
    first = life.draw_life([0, 1, 2])
    assert first == life.draw_life([0, 1, 2])
    assert first["life"] == sum(first["rungs"].values()) <= first["worlds"] and first["systems"] == 3


def test_a_broken_band_trips_the_live_flag(expectations):
    result = life.draw_life([0, 1, 2, 3, 4, 5])
    broken = copy.deepcopy(expectations)
    broken["measures"]["stops_before_animals"].update(min=0.0, max=0.0)
    for bounds in broken["measures"].values():
        bounds.pop("anomaly", None)
    assert any(f.measure == "stops_before_animals" for f in life.compare(result, broken))


@pytest.mark.slow
def test_the_life_census_raises_no_new_flag(expectations):
    """Generates the seeded systems, writes the latest-run report and fails on a flag that no known anomaly explains."""
    result = life.run()
    findings = life.compare(result, expectations)
    header = f"Targets: {expectations['reference']['kind']}."
    stale = life.stale_anomalies(expectations, findings)
    if stale:
        header += "\n\nStale anomaly tags (the band holds again; remove the tag): " + ", ".join(stale)
    path = census.write_report(life.format_report(result, findings, header), name=life.REPORT_NAME)
    new = [finding for finding in findings if not finding.known]
    assert not new, f"new life census flags (report: {path}):\n" + "\n".join(f.line() for f in new)
