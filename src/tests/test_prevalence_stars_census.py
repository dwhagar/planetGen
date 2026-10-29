# tests/test_prevalence_stars_census.py

"""
TEST.432: the star census (part C of the prevalence suite). The star type mix at the default galaxy's reference solar position
and the multiplicity of each primary type, checked against the researched norms of `tests/data/prevalence_stars_expectations.json`
(the RECONS 10 parsec census, Mamajek's number densities, Duchene and Kraus, Raghavan, De Rosa, Moe and Di Stefano). The live census
fails on a flag no known anomaly explains.
"""

import copy
import re

import pytest

from tests.audit import prevalence as census
from tests.audit import prevalence_stars as stars

ITEM = re.compile(r"^[A-Z]+\.\d+$")
KINDS = {"observed", "norm"}


def _mix(n=40000, abandoned=0, errors=(), **values):
    return {"n": n, "values": values, "errors": list(errors), "error_count": len(errors), "abandoned": abandoned}


def _summary(mix=None, **multiplicity):
    return {"mix": mix or _mix(), "multiplicity": multiplicity}


def _spec(measure="K", low=0.10, high=0.14, star="G2V", companion=None):
    spec = {"mix": {measure: {"norm_min": low, "norm_max": high, "kind": "observed", "source": "test"}},
            "multiplicity": {}}
    if companion:
        spec["multiplicity"][star] = {"companion": {"norm_min": companion[0], "norm_max": companion[1], "kind": "norm",
                                                    "source": "test"}}
    return spec


# --- The harness: each flag trips on a deliberately broken input ------------------------------------------------------


def test_a_census_inside_its_norms_raises_no_flag():
    summary = _summary(_mix(K=0.12), G2V={"n": 2500, "values": {"companion": 0.46}, "errors": [], "error_count": 0, "abandoned": 0})
    assert stars.compare(summary, _spec(companion=(0.42, 0.50))) == []


def test_a_share_outside_its_norm_is_flagged():
    findings = stars.compare(_summary(_mix(K=0.19)), _spec())
    assert [(f.kind, f.cell, f.measure) for f in findings] == [("out-of-range", "mix", "K")]


def test_a_kind_of_star_that_never_appears_is_flagged_missing():
    findings = stars.compare(_summary(_mix(A=0.0)), _spec("A", 0.004, 0.016))
    assert [(f.kind, f.measure) for f in findings] == [("missing", "A")]


def test_a_multiplicity_outside_its_norm_is_flagged_for_its_star_type():
    summary = _summary(G2V={"n": 2500, "values": {"companion": 0.15}, "errors": [], "error_count": 0, "abandoned": 0})
    findings = stars.compare(summary, {"mix": {}, "multiplicity": _spec(companion=(0.42, 0.50))["multiplicity"]})
    assert [(f.kind, f.cell, f.measure) for f in findings] == [("out-of-range", "G2V", "companion")]


def test_a_draw_that_raised_or_ran_over_the_time_limit_is_flagged():
    findings = stars.compare(_summary(_mix(K=0.12, errors=["ValueError: x"], abandoned=1)), _spec())
    assert {f.kind for f in findings} == {"error", "hang"}


def test_a_known_anomaly_is_named_and_a_type_with_no_systems_is_flagged():
    spec = _spec()
    spec["mix"]["K"]["anomaly"] = "GEN.855"
    assert [f.known for f in stars.compare(_summary(_mix(K=0.2)), spec)] == ["GEN.855"]
    no_mix = {"mix": {}, "multiplicity": _spec(companion=(0.4, 0.5))["multiplicity"]}
    assert [(f.kind, f.cell) for f in stars.compare(_summary(), no_mix)] == [("missing", "G2V")]


def test_a_stale_tag_is_listed_when_its_range_holds_again():
    spec = _spec()
    spec["mix"]["K"]["anomaly"] = "GEN.855"
    assert stars.stale_anomalies(spec, stars.compare(_summary(_mix(K=0.12)), spec)) == ["mix K (GEN.855)"]
    assert stars.stale_anomalies(spec, stars.compare(_summary(_mix(K=0.2)), spec)) == []


def test_the_sampling_slack_shrinks_as_the_draw_grows():
    bounds = {"norm_min": 0.10, "norm_max": 0.14}
    small, large = stars.allowed(bounds, 1000), stars.allowed(bounds, 100000)
    assert small[0] < large[0] < 0.10 < 0.14 < large[1] < small[1]
    assert stars.allowed({"norm_min": 0.0, "norm_max": 0.0001}, 100)[0] == 0.0


def test_the_report_names_the_flags_and_the_figures():
    summary = _summary(_mix(K=0.19, M=0.68, G=0.03), G2V={"n": 100, "values": {"companion": 0.4}, "errors": [], "error_count": 0,
                                                          "abandoned": 0})
    report = stars.format_report(summary, stars.compare(summary, _spec()))
    assert "1 new, 0 known anomalies" in report and "OUT-OF-RANGE mix K" in report and "K 19.00%" in report
    assert "G2V (n=100" in report


def test_the_mix_draw_repeats_and_the_multiplicity_draw_counts_stars():
    first = stars.draw_mix(60)
    assert first == stars.draw_mix(60) and first["n"] == 60 and sum(first["counts"].values()) == 60
    assert set(first["counts"]) <= set(stars.MIX_MEASURES)
    one = stars.draw_multiplicity("G2V", 6)
    assert one == stars.draw_multiplicity("G2V", 6) and sum(one["stars"].values()) == 6
    assert all(int(count) >= 1 for count in one["stars"])


def test_the_solar_position_is_a_disk_place_with_every_population_density_given():
    populations = stars.solar_populations()
    assert populations["old"] > 0 and {"young", "intermediate", "old"} <= set(populations)


# --- The stored norms ------------------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def expectations():
    return stars.load_expectations()


def _all_bounds(expectations):
    for measure, bounds in expectations["mix"].items():
        yield f"mix {measure}", bounds
    for star, measures in expectations["multiplicity"].items():
        for measure, bounds in measures.items():
            yield f"{star} {measure}", bounds


def test_every_measure_the_census_counts_has_a_norm(expectations):
    assert set(expectations["mix"]) == set(stars.MIX_MEASURES)
    assert set(expectations["multiplicity"]) == set(stars.MULTIPLICITY_TYPES)
    assert all("companion" in measures for measures in expectations["multiplicity"].values())


def test_every_norm_has_a_kind_a_source_and_a_sane_range(expectations):
    for where, bounds in _all_bounds(expectations):
        assert bounds["kind"] in KINDS and bounds["source"], where
        assert 0.0 <= bounds["norm_min"] <= bounds["norm_max"] <= 1.0, where
        if "anomaly" in bounds:
            assert ITEM.match(bounds["anomaly"]) and bounds["anomaly"] in expectations["anomalies"], where


def test_the_reference_says_what_the_norms_are_read_from(expectations):
    reference = expectations["reference"]
    assert reference["kind"].startswith("norm") and reference["mix_draws"] == stars.MIX_DRAWS
    assert reference["systems_per_type"] == stars.SYSTEMS_PER_TYPE and set(reference["sources"]) >= {"RECONS", "Mamajek"}


def test_the_observed_mix_norms_cover_what_a_census_of_stars_adds_up_to(expectations):
    """The lowest shares can total at most 1 and the highest at least 1, else no mix could meet them all."""
    mix = expectations["mix"]
    assert sum(b["norm_min"] for b in mix.values()) <= 1.0 <= sum(b["norm_max"] for b in mix.values())


# --- The census ------------------------------------------------------------------------------------------------------------


@pytest.mark.slow
def test_the_star_census_raises_no_new_flag(expectations):
    """Draws the stars and the systems as the norms file says, writes the latest-run report and fails on a flag that no known
    anomaly explains. Known anomalies are in the report; one that no longer trips is listed as stale there."""
    reference = expectations["reference"]
    summary = stars.summarize(stars.run(reference["mix_draws"], reference["systems_per_type"]))
    findings = stars.compare(summary, expectations)
    header = (f"Norms: {reference['kind']}. {reference['mix_draws']} stars at {reference['point']}; "
              f"{reference['systems_per_type']} systems per primary type.")
    stale = stars.stale_anomalies(expectations, findings)
    if stale:
        header += "\n\nStale anomaly tags (the range holds again; remove the tag): " + ", ".join(stale)
    path = census.write_report(stars.format_report(summary, findings, header), name=stars.REPORT_NAME)
    new = [finding for finding in findings if not finding.known]
    assert not new, f"new star census flags (report: {path}):\n" + "\n".join(f.line() for f in new)


def test_a_broken_norm_trips_the_live_flag(expectations):
    summary = stars.summarize({"mix": stars.draw_mix(200), "multiplicity": {"G2V": stars.draw_multiplicity("G2V", 20)}})
    broken = copy.deepcopy(expectations)
    broken["mix"]["M"].update(norm_min=0.99, norm_max=1.0)
    assert any(f.cell == "mix" and f.measure == "M" and f.kind == "out-of-range" for f in stars.compare(summary, broken))
