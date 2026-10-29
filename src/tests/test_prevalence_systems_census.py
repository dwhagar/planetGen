# tests/test_prevalence_systems_census.py

"""
TEST.431: the system census (part B of the prevalence suite). Seeded full systems of each census star type, summed into the figures
a whole system decides (planets, moons, belts and comets per system, the share with no planet, the share that are pairs, the
share with a habitable-class or a Blue or Green world, the share where life reaches each milestone), checked against the stored
ranges of `tests/data/prevalence_systems_expectations.json`. The live census fails on a flag with no known anomaly.
"""

import copy
import re

import pytest

from tests.audit import prevalence as census
from tests.audit import prevalence_occurrence as occurrence
from tests.audit import prevalence_systems as systems

ITEM = re.compile(r"^[A-Z]+\.\d+$")


def _summary(star="G2V", n=100, abandoned=0, errors=(), **values):
    return {star: {"n": n, "abandoned": abandoned, "errors": list(errors), "error_count": len(errors), "values": values, "sd": {}}}


def _spec(**measures):
    return {"types": {"G2V": {"measures": {name: {"min": lo, "max": hi, "kind": "stand-in", "source": "test"}
                                            for name, (lo, hi) in measures.items()}}}}


# --- The harness: each flag trips on a deliberately broken input ------------------------------------------------------


def test_a_census_inside_its_ranges_raises_no_flag():
    assert systems.compare(_summary(planets=8.0, binary=0.4), _spec(planets=(5.0, 11.0), binary=(0.3, 0.6))) == []


def test_a_figure_outside_its_range_is_flagged():
    findings = systems.compare(_summary(planets=2.0, binary=0.9), _spec(planets=(5.0, 11.0), binary=(0.3, 0.6)))
    assert {(f.kind, f.measure) for f in findings} == {("out-of-range", "planets"), ("out-of-range", "binary")}


def test_a_milestone_nobody_reaches_is_flagged_missing():
    findings = systems.compare(_summary(life_technological_civilization=0.0), _spec(life_technological_civilization=(0.01, 0.3)))
    assert [(f.kind, f.measure) for f in findings] == [("missing", "life_technological_civilization")]


def test_a_system_over_the_time_limit_and_one_that_raised_are_flagged():
    findings = systems.compare(_summary(abandoned=2, errors=["ZeroDivisionError: x"]), _spec())
    assert {f.kind for f in findings} == {"hang", "error"}


def test_a_known_anomaly_is_named_and_no_type_drawn_is_flagged():
    spec = _spec(planets=(5.0, 11.0))
    spec["types"]["G2V"]["measures"]["planets"]["anomaly"] = "GEN.855"
    assert [f.known for f in systems.compare(_summary(planets=1.0), spec)] == ["GEN.855"]
    assert [f.kind for f in systems.compare({}, spec)] == ["missing"]


def test_a_share_expected_five_times_must_appear_whatever_the_band():
    low, high = systems.measure_range("life_complex_cells", 0.14, 150, 100)
    assert low >= 1.0 / 100 and high > 0.14
    assert systems.measure_range("life_complex_cells", 0.02, 150, 100)[0] == 0.0
    assert systems.measure_range("planets", 8.0, 150, 100, sd=3.0)[0] > 5.0


def test_a_record_summary_gives_shares_and_means():
    records = [{"star": "G2V", "planets": 4, "belts": 1, "moons": 4, "comets": 2, "binary": True, "habitable_class": True,
                "blue_green": False, "reached": dict.fromkeys(systems.MILESTONES, False)},
               {"star": "G2V", "planets": 0, "belts": 0, "moons": 0, "comets": 0, "binary": False, "habitable_class": False,
                "blue_green": False, "reached": dict.fromkeys(systems.MILESTONES, False)},
               {"star": "G2V", "abandoned": True}, {"star": "G2V", "error": "ValueError: x"}]
    summary = systems.summarize(records)["G2V"]
    assert (summary["n"], summary["abandoned"], summary["error_count"]) == (2, 1, 1)
    assert summary["values"]["planets"] == 2.0 and summary["values"]["zero_planets"] == 0.5
    assert summary["values"]["binary"] == 0.5 and summary["values"]["moons_per_planet"] == 1.0


def test_a_type_draw_repeats():
    first = systems.draw_types(2, ("M2V",))
    assert first == systems.draw_types(2, ("M2V",)) and len(first) == 2 and not any(r.get("error") for r in first)


# --- The stored ranges -----------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def expectations():
    return systems.load_expectations()


def test_every_star_type_has_every_measure_with_a_kind_and_a_source(expectations):
    assert set(expectations["types"]) == set(census.STAR_TYPES)
    for star, spec in expectations["types"].items():
        assert set(spec["measures"]) == {*systems.MEAN_MEASURES, *systems.SHARE_MEASURES}, star
        for name, bounds in spec["measures"].items():
            where = f"{star} {name}"
            assert bounds["kind"] in {"observed", "norm", "stand-in"} and bounds["source"], where
            assert 0.0 <= bounds["min"] <= bounds["max"], where
            if name in systems.SHARE_MEASURES:
                assert bounds["max"] <= 1.0, where
            if "anomaly" in bounds:
                assert ITEM.match(bounds["anomaly"]) and bounds["anomaly"] in expectations["anomalies"], where


def test_the_reference_records_which_commit_and_seeds_it_was_measured_on(expectations):
    reference = expectations["reference"]
    assert re.fullmatch(r"[0-9a-f]{40}", reference["commit"]) and reference["systems_per_type"] >= 50
    assert reference["seeds"]


# --- The census ------------------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def occurrence_expectations():
    return occurrence.load_expectations()


@pytest.mark.slow
def test_the_system_and_occurrence_census_raise_no_new_flag(expectations, occurrence_expectations):
    """Draws the systems as the reference did (the same count per star type, once: the occurrence census reads the same systems),
    writes both latest-run reports and fails on a flag no known anomaly explains, in either."""
    records = systems.draw_types(expectations["reference"]["systems_per_type"])
    summary = systems.summarize(records)
    findings = systems.compare(summary, expectations)
    header = (f"Reference: {expectations['reference']['label']} at {expectations['reference']['commit'][:9]}, "
              f"{expectations['reference']['systems_per_type']} systems per star type.")
    stale = systems.stale_anomalies(expectations, findings)
    if stale:
        header += "\n\nStale anomaly tags (the range holds again; remove the tag): " + ", ".join(stale)
    path = census.write_report(systems.format_report(summary, findings, header), name=systems.REPORT_NAME)

    grouped = occurrence.summarize(records)
    gaps = occurrence.compare(grouped, occurrence_expectations)
    gap_header = "Observed rates recalled by Aux 16 (tier T2/T3, unverified until the papers are read)."
    stale_gaps = systems.stale_anomalies(occurrence_expectations, gaps)
    if stale_gaps:
        gap_header += "\n\nStale anomaly tags (the rate is met now; remove the tag): " + ", ".join(stale_gaps)
    gap_path = census.write_report(occurrence.format_report(grouped, gaps, gap_header), name=occurrence.REPORT_NAME)

    new = [finding for finding in findings if not finding.known]
    new_gaps = [finding for finding in gaps if not finding.known]
    assert not new and not new_gaps, (
        f"new system census flags (report: {path}):\n" + "\n".join(f.line() for f in new)
        + f"\nnew gaps against observation (report: {gap_path}):\n" + "\n".join(f.line() for f in new_gaps))


def test_a_broken_range_trips_the_live_flag(expectations):
    summary = systems.summarize(systems.draw_types(3, ("G2V",)))
    broken = copy.deepcopy(expectations)
    broken["types"]["G2V"]["measures"]["planets"].update(min=100.0, max=200.0)
    assert any(f.measure == "planets" and f.kind == "out-of-range" for f in systems.compare(summary, broken))
