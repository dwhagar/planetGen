# tests/test_mass_function_curve.py

"""
TEST.454: the star mass function census (item 1 of the star mass function chain, GEN.922). The generator's stars, bin by bin in mass,
against the published curves held as data in `tests/data/mass_function_curve.json` [R0705, R2498]. The live census
fails on a flag no known anomaly explains; a bin the chain's later items fix carries the item's number until it does.
"""

import copy
import math
import re

import pytest

from planetgen import tuning
from planetgen.physics import imf
from tests.audit import mass_curve as curve
from tests.audit import prevalence as census

ITEM = re.compile(r"^[A-Z]+\.\d+$")


@pytest.fixture(scope="module")
def spec():
    return curve.load_curve()


def _draws(primaries=None, stars=None, **kwargs):
    """A `run`-shaped result with `primaries` and `stars` as the six bins' counts."""
    primaries = primaries or [100, 100, 100, 100, 100, 100]
    stars = stars or primaries
    return {"n": sum(primaries), "primary_bins": primaries, "star_bins": stars, "systems": sum(primaries), "stars": sum(stars),
            "brown_dwarfs": 0, "phase_bins": [{} for _ in range(6)], "wd_mass_sum": 0.0, "wd_count": 0, "errors": [],
            "error_count": 0, **kwargs}


def _no_anomalies(spec):
    quiet = copy.deepcopy(spec)
    quiet["anomalies"] = {}
    return quiet


# --- The curves as data ----------------------------------------------------------------------------------------------------


def test_each_curve_gives_bin_shares_that_add_to_one(spec):
    for name, shape in spec["curves"].items():
        shares = curve.bin_shares(shape, spec["bins"]["edges_sol"])
        assert math.isclose(sum(shares), 1.0, rel_tol=1e-9), name
        assert all(share > 0 for share in shares), name


def test_the_system_function_peaks_higher_in_mass_than_the_single_star_function(spec):
    """Chabrier's systems peak near 0.22 solar masses, single stars near 0.08: the system curve has the lighter share of its
    low bins and the heavier share of its upper bins."""
    edges = spec["bins"]["edges_sol"]
    system = curve.bin_shares(spec["curves"]["chabrier2003_system"], edges)
    single = curve.bin_shares(spec["curves"]["chabrier2003_individual"], edges)
    assert system[0] < single[0] and system[-1] > single[-1]


def test_the_chabrier_individual_curve_is_the_shape_the_generator_rows_carry_today(spec):
    """The curve file's copy of [R0705] (individual stars) is the function `tuning` holds (GEN.927), so a change of the
    rows is seen against the published curve and not against itself: the generator's bins are the curve's, to the
    chord error of the segments."""
    shape = spec["curves"]["chabrier2003_individual"]
    assert (shape["m_c"], shape["sigma"], shape["break_sol"]) == (
        tuning.IMF_BELL_CENTRE_SOL, tuning.IMF_BELL_WIDTH_DEX, tuning.IMF_BELL_TOP_SOL)
    assert shape["power_law_x"] + 1.0 == tuning.IMF_HIGH_SLOPE
    edges = spec["bins"]["edges_sol"]
    published = curve.bin_shares(shape, edges)
    held = [imf.number(a, b) for a, b in zip(edges, edges[1:])]
    held = [h / sum(held) for h in held]
    assert held == pytest.approx(published, rel=0.01)


def test_the_power_law_curve_matches_a_hand_integral(spec):
    shape = {"breaks_sol": [0.1, 1.0, 10.0], "slopes": [1.0, 2.0]}
    shares = curve.bin_shares(shape, [0.1, 1.0, 10.0])
    # dN/dM = M^-1 up to 1 (ln 10 stars), then M^-2 (continuous at 1): integral from 1 to 10 of M^-2 = 0.9
    assert shares[0] == pytest.approx(math.log(10) / (math.log(10) + 0.9), rel=1e-3)


def test_the_bins_are_ascending_and_inside_the_function(spec):
    edges = spec["bins"]["edges_sol"]
    assert edges == sorted(edges) and edges[0] >= tuning.IMF_BREAKS_SOL[0] and edges[-1] <= 1.0
    above = spec["bins"]["evolved_report_edges_sol"]
    assert above == sorted(above) and above[0] == edges[-1] and above[-1] <= tuning.IMF_BREAKS_SOL[-1]


def test_every_sample_names_curves_that_exist_and_a_sane_tolerance(spec):
    for name, sample in spec["samples"].items():
        assert sample["curves"] and set(sample["curves"]) <= set(spec["curves"]), name
        assert 0.0 <= sample["tolerance_relative"] < 0.5, name


def test_every_curve_and_norm_says_where_it_comes_from(spec):
    for name, shape in spec["curves"].items():
        assert "T2" in shape["source"] or "recalled" in shape["source"], name
    assert spec["norms"]["white_dwarf_mean_mass_sol"]["source"]
    assert any("[R0705]" in text for text in spec["reference"]["sources"])


def test_every_known_anomaly_is_an_item_with_a_note(spec):
    for cell, item in spec["anomalies"].items():
        assert ITEM.match(item) and item in spec["anomaly_notes"], cell


def test_the_ranges_widen_by_the_tolerance_and_stay_inside_zero_and_one(spec):
    bare = copy.deepcopy(spec)
    for sample in spec["samples"]:
        bare["samples"][sample]["tolerance_relative"] = 0.0
        tolerance = spec["samples"][sample]["tolerance_relative"]
        for (low, high, share_min, share_max), (_l, _h, raw_min, raw_max) in zip(
                curve.expected_ranges(spec, sample), curve.expected_ranges(bare, sample)):
            assert 0.0 < share_min < share_max <= 1.0
            assert share_min == pytest.approx(raw_min * (1.0 - tolerance))
            assert share_max == pytest.approx(min(1.0, raw_max * (1.0 + tolerance)))


# --- The harness: each flag trips on a deliberately broken input -----------------------------------------------------------


def test_a_census_that_follows_the_curves_raises_no_flag(spec):
    quiet = _no_anomalies(spec)
    system = curve.bin_shares(spec["curves"]["chabrier2003_system"], spec["bins"]["edges_sol"])
    single = curve.bin_shares(spec["curves"]["chabrier2003_individual"], spec["bins"]["edges_sol"])
    draws = _draws([round(share * 50000) for share in system], [round(share * 50000) for share in single])
    assert curve.compare(draws, quiet) == []


def test_a_bin_with_too_many_stars_is_flagged_for_its_sample(spec):
    quiet = _no_anomalies(spec)
    draws = _draws([5000, 100, 100, 100, 100, 100])
    findings = curve.compare(draws, quiet)
    assert ("out-of-range", "primaries") in {(f.kind, f.cell) for f in findings}
    assert all(f.measure.endswith("solar masses") for f in findings)


def test_an_empty_bin_is_flagged_missing(spec):
    quiet = _no_anomalies(spec)
    findings = curve.compare(_draws([0, 100, 100, 100, 100, 100]), quiet)
    assert ("missing", "primaries") in {(f.kind, f.cell) for f in findings}


def test_a_sample_with_no_stars_is_flagged(spec):
    findings = curve.compare(_draws([0] * 6, [0] * 6), _no_anomalies(spec))
    assert {f.kind for f in findings} == {"missing"} and {f.cell for f in findings} == {"primaries", "all_stars"}


def test_a_known_anomaly_is_named_and_a_stale_tag_is_listed(spec):
    tagged = copy.deepcopy(spec)
    tagged["anomalies"] = {"primaries 0.08-0.12": "GEN.928"}
    findings = curve.compare(_draws([5000, 100, 100, 100, 100, 100]), tagged)
    assert [f.known for f in findings if f.cell == "primaries" and f.measure.startswith("0.08")] == ["GEN.928"]
    assert curve.stale_anomalies(tagged, []) == ["primaries 0.08-0.12 (GEN.928)"]
    assert curve.stale_anomalies(tagged, findings) == []


def test_white_dwarfs_outside_their_mean_mass_are_flagged(spec):
    quiet = _no_anomalies(spec)
    system = curve.bin_shares(spec["curves"]["chabrier2003_system"], spec["bins"]["edges_sol"])
    single = curve.bin_shares(spec["curves"]["chabrier2003_individual"], spec["bins"]["edges_sol"])
    even = dict(primaries=[round(share * 50000) for share in system], stars=[round(share * 50000) for share in single])
    flagged = curve.compare(_draws(wd_count=10, wd_mass_sum=9.0, **even), quiet)
    assert [(f.cell, f.measure) for f in flagged] == [("white dwarfs", "mean mass (solar masses)")]
    assert curve.compare(_draws(wd_count=10, wd_mass_sum=6.2, **even), quiet) == []


def test_a_draw_that_raised_is_flagged(spec):
    system = curve.bin_shares(spec["curves"]["chabrier2003_system"], spec["bins"]["edges_sol"])
    draws = _draws([round(share * 50000) for share in system], error_count=2, errors=["ValueError: x"])
    kinds = {f.kind for f in curve.compare(draws, _no_anomalies(spec))}
    assert "error" in kinds


def test_the_report_names_the_flags_the_bins_and_the_phases(spec):
    draws = _draws([5000, 100, 100, 100, 100, 100], phase_bins=[{"main sequence": 3, "white dwarf": 1}] + [{} for _ in range(5)],
                   wd_count=2, wd_mass_sum=1.3)
    draws.update(systems=7, stars=9)
    report = curve.format_report(draws, curve.compare(draws, _no_anomalies(spec)), spec)
    assert "OUT-OF-RANGE primaries" in report and "| 0.08 to 0.12 |" in report
    assert "75.0%" in report and "mean mass 0.650" in report and "| 0.8 to 1.2 (4) |" in report


def test_the_draw_repeats_and_counts_every_system():
    first = curve.draw_stars(40, label="test")
    assert first == curve.draw_stars(40, label="test")
    assert first["systems"] + first["error_count"] == 40
    assert sum(first["primary_bins"]) <= 40 and sum(first["star_bins"]) <= first["stars"]


# --- The census ------------------------------------------------------------------------------------------------------------


@pytest.mark.slow
def test_the_mass_function_census_raises_no_new_flag(spec):
    """Draws the primaries and their systems as the curve file says, writes the latest-run report and fails on a flag no known
    anomaly explains. A tag whose bin holds again is listed as stale in the report."""
    draws = curve.run(spec["reference"]["primaries"])
    findings = curve.compare(draws, spec)
    header = f"Curves: {spec['reference']['kind']}. {spec['reference']['primaries']} systems at {spec['reference']['point']}."
    stale = curve.stale_anomalies(spec, findings)
    if stale:
        header += "\n\nStale anomaly tags (the bin holds again; remove the tag): " + ", ".join(stale)
    path = census.write_report(curve.format_report(draws, findings, spec, header), name=curve.REPORT_NAME)
    new = [finding for finding in findings if not finding.known]
    assert not new, f"new mass function flags (report: {path}):\n" + "\n".join(f.line() for f in new)
