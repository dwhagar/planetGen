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


def _result(oxygenic=None, multi=None, living=None, traits=None, trait_worlds=0, **rungs):
    """A run whose life worlds stand on `rungs`; the oxygenic worlds of the host systems stand on `oxygenic` (default: the same rungs)."""
    life_worlds = sum(rungs.values())
    host = rungs if oxygenic is None else oxygenic
    return {"systems": 10, "worlds": 100, "life": life_worlds, "rungs": rungs, "classes": {}, "host_systems": 10, "host_oxygenic": host,
            "multi": int(0.8 * life_worlds) if multi is None else multi, "living": living or {}, "traits": traits or {}, "trait_worlds": trait_worlds, "minor_sentient": 0, "abandoned": 0, "errors": [], "error_count": 0}


def _expectations(**changes):
    bands = {"stops_before_animals": (0.5, 1.0), "reaches_animals": (0.15, 0.35), "reaches_sentient": (0.0, 0.02),
             "sentient_of_animal_worlds": (0.0, 0.05), "two_or_more_classes": (0.6, 1.0), "any_trait": (0.01, 0.25), "extreme_of_traits": (0.0, 0.75)}
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


def test_the_share_of_worlds_with_several_living_classes_is_over_every_life_world():
    values, counts = life.summarize(_result(origin=30, producers=30, animals=30, sentient=10, multi=60))
    assert values["two_or_more_classes"] == 0.6 and counts["two_or_more_classes"] == 100


def test_a_list_that_collapses_to_one_class_is_flagged():
    flags = life.compare(_result(origin=60, producers=10, animals=29, sentient=1, multi=0), _expectations())
    assert [f.measure for f in flags] == ["two_or_more_classes"] and flags[0].kind == "missing"


def test_the_report_names_the_share_of_worlds_each_class_lives_on():
    result = _result(origin=60, producers=10, animals=29, sentient=1, living={"3a": 80, "4b": 60})
    report = life.format_report(result, life.compare(result, _expectations()))
    assert "3a 0.800" in report and "4b 0.600" in report and "beside the dominant one reached a civilization" in report


def test_the_report_names_the_life_classes_no_sampled_world_lived_on():
    from planetgen.physics import life_classes

    result = _result(origin=60, producers=10, animals=29, sentient=1, living={"3a": 80, "4b": 60})
    report = life.format_report(result, life.compare(result, _expectations()))
    line = next(line for line in report.splitlines() if line.startswith("Life classes that lived on no sampled world"))
    assert "3a" not in line.split(": ", 1)[1] and "4b" not in line.split(": ", 1)[1]
    assert all(class_id in line for class_id in life_classes.IDS if class_id not in ("3a", "4b"))
    everything = _result(origin=60, producers=10, animals=29, sentient=1, living={class_id: 5 for class_id in life_classes.IDS})
    assert next(line for line in life.format_report(everything, []).splitlines() if line.startswith("Life classes that lived on no")).endswith("none")


# --- The harness: each flag trips on a deliberately broken input ------------------------------------------------------------


def test_shares_inside_their_bands_raise_no_flag():
    assert life.compare(_result(origin=60, producers=10, animals=29, sentient=1), _expectations()) == []


def test_a_share_outside_its_band_is_flagged_both_ways():
    flags = life.compare(_result(origin=100, animals=800, sentient=100), _expectations())
    assert {f.measure for f in flags} >= {"stops_before_animals", "reaches_animals", "reaches_sentient", "sentient_of_animal_worlds"}
    assert [f.kind for f in life.compare(_result(origin=400, animals=600), _expectations()) if f.measure == "stops_before_animals"] == [
        "out-of-range"]


def test_the_animal_shares_are_taken_over_the_oxygenic_worlds_and_not_over_all_life_worlds():
    result = _result(origin=950, producers=50, oxygenic={"producers": 3, "animals": 6, "sentient": 1})
    values, counts = life.summarize(result)
    assert values["stops_before_animals"] == 1.0 and counts["stops_before_animals"] == 1000
    assert values["reaches_animals"] == 0.7 and counts["reaches_animals"] == 10
    assert values["reaches_sentient"] == 0.1 and values["sentient_of_animal_worlds"] == 1 / 7 and counts["sentient_of_animal_worlds"] == 7


def test_a_host_run_that_found_no_oxygenic_world_is_flagged_missing():
    result = _result(origin=100, oxygenic={})
    assert any(f.measure == "oxygenic life worlds" and f.kind == "missing" for f in life.compare(result, _expectations()))


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
    assert "Oxygenic life worlds of the default systems: 0 of 100" in report


def test_the_report_counts_the_oxygenic_life_worlds_and_the_rung_they_reached():
    result = _result(origin=96, animals=4)
    result["oxygenic"] = {"producers": 3, "animals": 1}
    report = life.format_report(result, life.compare(result, _expectations()))
    assert "Oxygenic life worlds of the default systems: 4 of 100, highest rung: origin 0, producers 3" in report
    assert "around G2V, K5V" in report and "animals 1" in report


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
    result = life.draw_life([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])      # GEN.1334: a planet under an envelope holds no life, so 6 systems held one living world
    broken = copy.deepcopy(expectations)
    broken["measures"]["stops_before_animals"].update(min=0.0, max=0.0)
    for bounds in broken["measures"].values():
        bounds.pop("anomaly", None)
    assert any(f.measure == "stops_before_animals" for f in life.compare(result, broken))


@pytest.mark.slow
@pytest.mark.timeout(3600)    # PERF.241: the census draws every cell as its reference did, on every core: 886 to 1,036 s on 4 vCPUs against the default 900 s
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


# --- TEST.436: traits, chemistries per star letter, animal grades ---------------------------------------------------------------


def test_the_trait_shares_are_over_the_life_worlds_and_over_the_traits_drawn():
    result = _result(origin=30, producers=30, animals=30, sentient=10, traits={"radiation:2": 6, "radiation:3": 2, "sulfide:2": 2}, trait_worlds=8)
    values, counts = life.summarize(result)
    assert values["any_trait"] == 0.08 and counts["any_trait"] == 100
    assert values["extreme_of_traits"] == 0.2 and counts["extreme_of_traits"] == 10


def test_a_census_in_which_the_trait_draw_went_silent_is_flagged():
    result = _result(origin=60000, producers=10000, animals=29000, sentient=1000, multi=80000)
    flags = life.compare(result, _expectations())
    assert [f.measure for f in flags] == ["any_trait"] and flags[0].kind == "missing"


def test_a_census_in_which_every_life_world_carries_a_trait_is_flagged():
    result = _result(origin=30, producers=30, animals=30, sentient=10, traits={"radiation:2": 90}, trait_worlds=90)
    assert "any_trait" in [f.measure for f in life.compare(result, _expectations())]


def test_the_star_letter_of_a_body_is_the_first_letter_of_its_spectral_type():
    assert life.star_letter(types.SimpleNamespace(star=types.SimpleNamespace(star_type="g2V"))) == "G"
    assert life.star_letter(types.SimpleNamespace(star=types.SimpleNamespace(star_type=None))) == "?"
    assert life.star_letter(types.SimpleNamespace()) == "?"


def test_the_draw_counts_traits_star_letters_and_animal_grades_of_the_life_worlds():
    found = life.draw_life([0, 1, 2, 3])
    assert sum(found["by_letter"].values()) == found["life"] == sum(found["grades"].values())
    assert found["trait_worlds"] <= found["life"] and sum(found["traits"].values()) >= found["trait_worlds"]
    assert all(key.split(":")[1].isalpha() or key.endswith(":?") for key in found["by_letter"])
    assert all(name.split(":")[0] in ("radiation", "sulfide", "carbon_dioxide") and name.endswith((":2", ":3")) for name in found["traits"])


def test_the_report_lists_the_traits_the_classes_by_star_letter_and_the_animal_grades():
    result = _result(origin=60, producers=10, animals=29, sentient=1, traits={"radiation:2": 3}, trait_worlds=3)
    result.update(classes={"3a": 40, "4b": 60}, by_letter={"3a:K": 30, "3a:G": 10, "4b:G": 60}, grades={"3a:0": 40, "4b:2": 20, "4b:0": 40})
    report = life.format_report(result, life.compare(result, _expectations()))
    assert "radiation:2 3" in report and "- 3a: G 10, K 30" in report and "- 4b: grade 0 40, grade 2 20" in report
