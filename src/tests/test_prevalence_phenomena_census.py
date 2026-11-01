# tests/test_prevalence_phenomena_census.py

"""
TEST.433: the phenomena census (part D of the prevalence suite). The galaxy-wide totals of the classes the phenomenon scatter places
(stellar-mass and intermediate-mass black holes, neutron stars, magnetars, planetary nebulae, supernova remnants) over the Milky Way-shaped
reference galaxy, checked against the researched ranges of `tests/data/prevalence_phenomena_expectations.json`. The live census fails on a
flag no known anomaly explains.
"""

import copy
import re

import pytest

from tests.audit import prevalence as census
from tests.audit import prevalence_phenomena as phenomena

ITEM = re.compile(r"^[A-Z]+\.\d+$")
KINDS = {"observed", "norm", "stand-in"}


def _spec(measure="supernova_remnants", low=294, high=3000):
    return {"measures": {measure: {"min": low, "max": high, "kind": "observed", "source": "test"}}}


# --- The harness: each flag trips on a deliberately broken input ------------------------------------------------------


def test_a_total_inside_its_range_raises_no_flag():
    assert phenomena.compare({"supernova_remnants": 1000.0}, _spec()) == []


def test_a_total_outside_its_range_is_flagged_both_ways():
    assert [f.kind for f in phenomena.compare({"supernova_remnants": 260.0}, _spec())] == ["out-of-range"]
    assert [f.kind for f in phenomena.compare({"supernova_remnants": 9000.0}, _spec())] == ["out-of-range"]


def test_a_class_with_no_objects_is_flagged_missing():
    assert [f.kind for f in phenomena.compare({"supernova_remnants": 0.0}, _spec())] == ["missing"]
    assert [f.kind for f in phenomena.compare({}, _spec())] == ["missing"]


def test_a_known_anomaly_is_named_and_a_stale_tag_is_listed():
    spec = _spec()
    spec["measures"]["supernova_remnants"]["anomaly"] = "GEN.900"
    assert [f.known for f in phenomena.compare({"supernova_remnants": 260.0}, spec)] == ["GEN.900"]
    assert phenomena.stale_anomalies(spec, phenomena.compare({"supernova_remnants": 1000.0}, spec)) == [
        "supernova_remnants (GEN.900)"]
    assert phenomena.stale_anomalies(spec, phenomena.compare({"supernova_remnants": 260.0}, spec)) == []


def test_a_total_that_rests_on_draws_is_widened_by_their_poisson_error():
    bounds = {"min": 1000, "max": 1500}
    assert phenomena.allowed("x", bounds) == (1000, 1500)
    low, high = phenomena.allowed("x", bounds, hits=400)
    assert low == pytest.approx(800) and high == pytest.approx(1800)
    spec = _spec("wolf_rayet_stars", 1000, 1500)
    assert phenomena.compare({"wolf_rayet_stars": 850.0}, spec, {"wolf_rayet_stars": 400}) == []
    assert [f.kind for f in phenomena.compare({"wolf_rayet_stars": 850.0}, spec)] == ["out-of-range"]


def test_the_neutron_stars_total_includes_the_magnetars():
    totals = {("neutron-star", None): 1000.0, ("neutron-star", "magnetar"): 5.0, ("black-hole", "stellar"): 7.0}
    values = phenomena.summarize(totals)
    assert values["neutron_stars"] == 1005.0 and values["magnetars"] == 5.0 and values["stellar_black_holes"] == 7.0
    assert values["planetary_nebulae"] == 0.0


def test_the_report_names_the_flags_and_the_totals():
    values = {"supernova_remnants": 260.0}
    report = phenomena.format_report(values, phenomena.compare(values, _spec()))
    assert "1 new, 0 known anomalies" in report and "OUT-OF-RANGE supernova_remnants: 260 not in [294, 3,000]" in report
    assert "- supernova_remnants: 260" in report


def test_the_first_shells_total_is_repeatable_and_has_a_total_for_the_classes_that_reach_it():
    first = phenomena.galaxy_totals(workers=1, shells=3)
    assert first == phenomena.galaxy_totals(workers=1, shells=3)
    assert first and all(value >= 0.0 for value in first.values())
    assert set(first) <= set(phenomena.CLASS_MEASURES)


# --- The stored ranges -----------------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def expectations():
    return phenomena.load_expectations()


def test_every_measure_of_the_census_has_a_range(expectations):
    assert set(expectations["measures"]) == set(phenomena.CLASS_MEASURES.values()) | {"wolf_rayet_stars"} | set(phenomena.RATE_MEASURES)


def test_every_range_has_a_kind_a_source_and_a_sane_width(expectations):
    for measure, bounds in expectations["measures"].items():
        assert bounds["kind"] in KINDS and bounds["source"], measure
        assert 0 < bounds["min"] < bounds["max"], measure
        if "anomaly" in bounds:
            assert ITEM.match(bounds["anomaly"]) and bounds["anomaly"] in expectations["anomalies"], measure


def test_the_reference_says_which_galaxy_and_which_papers(expectations):
    reference = expectations["reference"]
    assert reference["kind"].startswith("norm") and "milky_way_like_shape" in reference["galaxy"]
    assert set(reference["sources"]) >= {"Green2019", "McGill", "FrewParker", "Rosslowe", "Mroz2017", "Kirkpatrick2021", "CantatGaudin2020"}


def test_the_reference_galaxy_holds_the_stars_its_history_forms_and_forms_the_milky_ways_rate():
    """GEN.900 to GEN.902: the totals of the remnants, magnetars and Wolf-Rayet stars follow the young stars, so the reference galaxy's
    density, disk history and mass must agree, and its star formation must be the Milky Way's (1.65 +- 0.19 solar masses a year,
    [R2680, R2681], recalled): the old fixture held a bulge as heavy as its thin disk and ran its history on a larger galaxy."""
    import numpy as np

    from planetgen.galaxy import density, disk_history

    shape = phenomena.milky_way_like_shape()
    masses = density.component_masses(shape)
    assert sum(masses.values()) == pytest.approx(phenomena.MILKY_WAY_STELLAR_MASS_SOLAR, rel=1.0e-6)
    assert masses["bulge"] / sum(masses.values()) == pytest.approx(phenomena.MILKY_WAY_BULGE_TO_TOTAL, abs=0.02)
    history = disk_history.history_of(shape)
    radii = history.radius_scale_lengths
    area = 2.0 * np.pi * radii * shape.disk_scale_length_pc * (radii[1] - radii[0]) * shape.disk_scale_length_pc
    assert float((history.stars_now * area).sum()) == pytest.approx(masses["thin_disk"], rel=0.05)
    rate = float((history.rate[:, -1] * area).sum()) / 1.0e9
    assert 1.0 < rate < 3.0, rate


# --- The census ------------------------------------------------------------------------------------------------------------


@pytest.mark.slow
def test_the_phenomena_census_raises_no_new_flag(expectations):
    """Totals the scatter's expectation over every shell of the reference galaxy, writes the latest-run report and fails on a flag
    that no known anomaly explains. Known anomalies are in the report; one that no longer trips is listed as stale there."""
    values = phenomena.summarize(phenomena.galaxy_totals())
    values["wolf_rayet_stars"], hits = phenomena.wolf_rayet_total()
    values.update(phenomena.rate_values())
    findings = phenomena.compare(values, expectations, {"wolf_rayet_stars": hits})
    header = f"Norms: {expectations['reference']['kind']}. Galaxy: {expectations['reference']['galaxy']}."
    stale = phenomena.stale_anomalies(expectations, findings)
    if stale:
        header += "\n\nStale anomaly tags (the range holds again; remove the tag): " + ", ".join(stale)
    path = census.write_report(phenomena.format_report(values, findings, header), name=phenomena.REPORT_NAME)
    new = [finding for finding in findings if not finding.known]
    assert not new, f"new phenomena census flags (report: {path}):\n" + "\n".join(phenomena.line(f) for f in new)


def test_the_wolf_rayet_count_repeats_and_is_positive_where_the_young_stars_are():
    total, hits = phenomena.wolf_rayet_total(workers=2, draws=2000, shells=40)
    assert (total, hits) == phenomena.wolf_rayet_total(workers=2, draws=2000, shells=40)
    assert total >= 0.0 and hits >= 0


def test_a_broken_range_trips_the_live_flag(expectations):
    values = phenomena.summarize(phenomena.galaxy_totals(workers=1, shells=3))
    broken = copy.deepcopy(expectations)
    broken["measures"]["stellar_black_holes"].update(min=1.0e12, max=1.0e13)
    assert any(f.measure == "stellar_black_holes" and f.kind in ("out-of-range", "missing")
               for f in phenomena.compare(values, broken))


def test_the_rates_are_read_where_the_fill_reads_them_and_repeat(expectations):
    values = phenomena.rate_values()
    assert values == phenomena.rate_values()
    assert set(values) == set(phenomena.RATE_MEASURES)
    assert phenomena.compare(values, {"measures": {m: b for m, b in expectations["measures"].items() if m in values}}) == []


def test_a_rate_scale_that_dials_a_kind_down_trips_its_flag(monkeypatch, expectations):
    from planetgen import tuning

    monkeypatch.setitem(tuning.PHENOMENON_RATE_SCALE, "rogue-planet", 0.1)
    monkeypatch.setitem(tuning.PHENOMENON_RATE_SCALE, "brown-dwarf", 2.0)
    values = phenomena.rate_values()
    flagged = {f.measure for f in phenomena.compare(values, expectations) if f.measure in phenomena.RATE_MEASURES}
    assert {"rogue_terrestrial_per_star", "rogue_terrestrial_per_pc3", "brown_dwarfs_per_pc3"} <= flagged


def test_the_open_clusters_of_the_reference_galaxy_are_alive_young_and_many(expectations):
    from planetgen import tuning

    count = phenomena.rate_values()["open_clusters"]
    bounds = expectations["measures"]["open_clusters"]
    assert bounds["min"] <= count <= bounds["max"] <= tuning.OPEN_CLUSTER_MAX_COUNT
