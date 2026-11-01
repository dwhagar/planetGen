# tests/test_wd_pollution.py

"""
GEN.508: the white-dwarf pollution detection model (`planetgen.physics.wd_pollution`): the sample of 24 polluted white dwarfs as
the papers have it, the loader's refusals, the absolute abundance error fitted to it, one body observed with the papers' rule
for what is reported, the sample a population would show (quadrature, no random numbers), and the registration as the selection model
the tuning harness reads. The Done test: a mock population applied to the model recovers the published recovery classes and mean.
"""

import copy
import json
import math
import statistics
import subprocess
import sys

import pytest

from planetgen.generation import tuning_harness as th
from planetgen.physics import wd_pollution as wd
from planetgen.physics.elements import DataError
from planetgen.util import draw


def written(tmp_path, mutate=None):
    with open(wd.DATA_FILE, encoding="utf-8") as handle:
        document = json.load(handle)
    if mutate:
        mutate(document)
    path = tmp_path / "sample.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    return str(path)


def refused(tmp_path, match, mutate):
    with pytest.raises(DataError, match=match):
        wd.load(written(tmp_path, mutate))


# --- the sample ------------------------------------------------------------------------------------------------------------


def test_the_sample_is_the_24_stars_with_the_published_statistics():
    sample = wd.SAMPLE
    assert len(sample.stars) == 24 and len({s.name for s in sample.stars}) == 24
    assert sample.counts == {"both": 16, "upper": 6, "limit": 2}
    values = sample.values
    assert len(values) == 22
    assert statistics.mean(values) == pytest.approx(-1.248, abs=0.001)
    assert statistics.median(values) == pytest.approx(-1.256, abs=0.001)
    assert statistics.stdev(values) == pytest.approx(0.447, abs=0.001)
    assert sample.teff_range == (10120.0, 21200.0) and wd.VALIDITY["teff_k"] == (10120.0, 21200.0)
    assert {s.source for s in sample.stars} == {"doyle20", "doyle23"}
    assert sum(1 for s in sample.stars if s.table == "1") == 16 and sum(1 for s in sample.stars if s.table == "4") == 8


def test_a_star_with_no_median_is_an_upper_limit_and_a_star_with_no_lower_bound_has_a_median_and_an_upper_bound():
    by_name = {s.name: s for s in wd.SAMPLE.stars}
    assert by_name["G 29-38"].recovery == "limit" and by_name["WD 1415+234"].recovery == "limit"
    assert by_name["Ton 345"].recovery == "upper" and by_name["GD 40"].recovery == "upper"
    assert by_name["WD 1350-162"].recovery == "both"
    for star in wd.SAMPLE.stars:
        assert (star.median_diw is None) == (star.recovery == "limit")
        assert (star.minus is None) == (star.recovery != "both")


def test_the_loader_refuses_a_broken_sample(tmp_path):
    refused(tmp_path, "listed twice", lambda d: d["stars"][1].update(name=d["stars"][0]["name"]))
    refused(tmp_path, "not in the sources", lambda d: d["stars"][0].update(source="nobody"))
    refused(tmp_path, "not positive", lambda d: d["stars"][1].update(plus=0))
    refused(tmp_path, "not positive", lambda d: d["stars"][1].update(minus=-0.1))
    refused(tmp_path, "an upper limit has no lower bound", lambda d: d["stars"][8].update(minus=0.3))
    refused(tmp_path, "not a finite number", lambda d: d["stars"][1].update(teff_k="hot"))
    refused(tmp_path, "out of range", lambda d: d.update(x_metal=1.5))
    refused(tmp_path, "out of range", lambda d: d.update(negative_fraction_limit_only=0.1))
    refused(tmp_path, "no citation", lambda d: d["sources"]["doyle20"].pop("citation"))
    refused(tmp_path, "no stars", lambda d: d.update(stars=[]))
    with pytest.raises(DataError, match="cannot read"):
        wd.load(str(tmp_path / "missing.json"))


def test_the_fe_o_mole_fraction_and_the_oxygen_fugacity_convert_back_and_forth():
    assert wd.x_of(0.0) == pytest.approx(0.85)
    assert wd.x_of(-1.25) == pytest.approx(0.85 * 10 ** -0.625)
    for diw in (-3.0, -1.25, 0.4):
        assert wd.diw_of(wd.x_of(diw)) == pytest.approx(diw)


# --- the abundance error ---------------------------------------------------------------------------------------------------


def test_the_error_data_read_the_16_stars_with_both_bounds_and_bound_the_six_and_two_others():
    exact, censored = wd.error_data()
    assert len(exact) == 16 and all(0.02 < v < 0.2 for v in exact)
    assert len(censored) == 8 and censored[-2:] == [max(censored[:6])] * 2
    z = -statistics.NormalDist().inv_cdf(0.165)
    ton = wd.x_of(-2.098) / z
    assert any(c == pytest.approx(ton) for c in censored)


def test_the_error_model_is_a_deterministic_maximum_likelihood_fit_and_the_bounded_one_is_the_optimistic_bracket():
    model = wd.fit_error_model()
    assert model == wd.ERROR_MODEL and wd.fit_error_model() == model
    assert model.exact == 16 and model.censored == 8 and model.flagged
    assert model.mu == pytest.approx(-2.18, abs=0.1) and model.sigma == pytest.approx(1.0, abs=0.1)
    bounded = wd.BOUNDED_ERROR_MODEL
    assert not bounded.flagged and bounded.censored == 0
    assert bounded.mu < model.mu and bounded.sigma < model.sigma
    exact, censored = wd.error_data()
    assert (wd._log_likelihood(model.mu, model.sigma, exact, censored)
            >= wd._log_likelihood(bounded.mu, bounded.sigma, exact, censored))


def test_the_fit_needs_stars_with_both_bounds(tmp_path):
    def only_limits(document):
        for star in document["stars"]:
            star["minus"] = None
    sample = wd.load(written(tmp_path, only_limits))
    with pytest.raises(DataError, match="at least 3 stars"):
        wd.fit_error_model(sample)


# --- one body observed -----------------------------------------------------------------------------------------------------


def test_the_papers_rule_classifies_a_measurement_by_the_share_of_negative_draws():
    assert wd.classify(0.0, 0.1)[0] == "limit" and wd.classify(-0.05, 0.1)[0] == "limit"
    recovery, negative = wd.classify(0.05, 0.1)             # 0.5 sigma above zero: 31 percent negative
    assert recovery == "upper" and negative == pytest.approx(0.3085, abs=0.001)
    assert wd.classify(0.2, 0.1)[0] == "both"               # 2 sigma: 2.3 percent
    z = -statistics.NormalDist().inv_cdf(0.165)
    assert wd.classify(z * 0.1 * 1.001, 0.1)[0] == "both" and wd.classify(z * 0.1 * 0.999, 0.1)[0] == "upper"


def test_an_observation_is_two_draws_repeats_on_a_seed_and_hides_the_median_of_a_limit():
    with draw.bound(draw.Stream("wd:a")):
        first = wd.observe(-1.25)
    with draw.bound(draw.Stream("wd:a")):
        again = wd.observe(-1.25)
    assert first == again and first.diw_true == -1.25
    results = []
    for i in range(400):
        with draw.bound(draw.Stream(f"wd:b{i}")):
            results.append(wd.observe(-2.5))
    classes = {r.recovery for r in results}
    assert classes == {"both", "upper", "limit"}
    for r in results:
        assert (r.diw is None) == (r.recovery == "limit")
        item = r.as_item()
        assert item["recovery"] == r.recovery and item["diw"] == r.diw and 0.0 <= item["negative_fraction"] <= 1.0
        if r.recovery == "limit":
            assert r.negative_fraction >= 0.5


def test_an_oxidised_body_is_reported_with_both_bounds_and_a_reduced_one_often_is_not():
    oxidised = wd.class_probabilities(-0.3)
    reduced = wd.class_probabilities(-3.0)
    assert sum(oxidised.values()) == pytest.approx(1.0) and sum(reduced.values()) == pytest.approx(1.0)
    assert oxidised["both"] > 0.8 and reduced["limit"] > 0.3
    assert reduced["limit"] > oxidised["limit"] and reduced["both"] < oxidised["both"]
    ladder = [wd.class_probabilities(d)["limit"] for d in (-0.5, -1.0, -1.5, -2.0, -3.0)]
    assert ladder == sorted(ladder)


# --- the sample a population shows -----------------------------------------------------------------------------------------


def test_the_quadrature_agrees_with_many_observations():
    counts = {r: 0 for r in wd.RECOVERIES}
    reported = []
    n = 3000
    for i in range(n):
        with draw.bound(draw.Stream(f"wd:mc{i}")):
            body = draw.gauss(-1.25, 0.2)
            obs = wd.observe(body)
        counts[obs.recovery] += 1
        if obs.diw is not None:
            reported.append(obs.diw)
    expected = wd.expected_sample(-1.25, 0.2)
    for recovery in wd.RECOVERIES:
        assert counts[recovery] / n == pytest.approx(expected.fractions[recovery], abs=0.03), recovery
    assert statistics.mean(reported) == pytest.approx(expected.reported_mean, abs=0.06)


def test_the_expected_sample_is_pure_and_its_fractions_sum_to_one():
    first = wd.expected_sample(-1.25, 0.2)
    assert first == wd.expected_sample(-1.25, 0.2)
    assert sum(first.fractions.values()) == pytest.approx(1.0)
    assert wd.expected_sample(-1.25, 0.0).reported_sd > 0


def test_a_mock_population_recovers_the_published_recovery_classes_and_mean_within_their_intervals():
    """The Done test: the model applied to bodies centred where the sample is gives the classes and the mean the papers report."""
    expected = wd.expected_sample(-1.25, 0.2)
    counts = wd.SAMPLE.counts
    n = len(wd.SAMPLE.stars)
    for recovery in wd.RECOVERIES:
        low, high = th.wilson_interval(counts[recovery], n, z=2.0)
        assert low <= expected.fractions[recovery] <= high, (recovery, expected.fractions[recovery], low, high)
    standard_error = statistics.stdev(wd.SAMPLE.values) / math.sqrt(len(wd.SAMPLE.values))
    assert abs(expected.reported_mean - statistics.mean(wd.SAMPLE.values)) <= 2.0 * standard_error


def test_the_published_spread_lies_between_the_two_error_models_at_zero_width():
    """The spread of the reported medians is explained by the measurement alone: the intrinsic width is not identified by the sample."""
    observed = statistics.stdev(wd.SAMPLE.values)
    optimistic = wd.expected_sample(-1.25, 0.0, wd.BOUNDED_ERROR_MODEL).reported_sd
    pessimistic = wd.expected_sample(-1.25, 0.0, wd.ERROR_MODEL).reported_sd
    assert optimistic < observed < pessimistic
    assert wd.expected_sample(-1.25, 0.45, wd.BOUNDED_ERROR_MODEL).reported_sd > observed, "a width of 0.45 is too wide even there"


# --- the harness -----------------------------------------------------------------------------------------------------------


def target(**overrides):
    raw = {"id": "wd-mean", "row": "tab:r", "statistic": "mean of dIW over the polluted white dwarfs", "kind": "mean",
           "value": -1.25, "sigma": 0.1, "n": 22, "source": "[R1152: Table 1] and [R1153: Table 4]", "label": "COMPUTED",
           "selection_model": wd.MODEL_ID, "field": "diw", "validity": {"teff_k": list(wd.VALIDITY["teff_k"])}}
    raw.update(overrides)
    return th.TargetTable({"version": 1, "targets": [raw]}).get("wd-mean")


def test_the_model_is_registered_under_its_id_and_reads_the_recovery_of_the_item():
    assert th.SELECTION_MODELS[wd.MODEL_ID] is wd.reported_median
    assert th.selection_model(wd.MODEL_ID) is wd.reported_median
    assert wd.reported_median({"recovery": "both"}) == 1.0 and wd.reported_median({"recovery": "upper"}) == 1.0
    assert wd.reported_median({"recovery": "limit"}) == 0.0
    with pytest.raises(DataError, match="none of"):
        wd.reported_median({"recovery": "maybe"})


def test_a_target_naming_the_model_imports_it_on_first_use():
    code = ("from planetgen.generation import tuning_harness as t; "
            "assert not t.SELECTION_MODELS; "
            "assert t.selection_model('wd-pollution-detection') is not None; "
            "assert set(t.SELECTION_MODELS) == {'wd-pollution-detection'}, 'only the named model is imported'; print('ok')")
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert result.stdout.strip() == "ok"


def test_the_harness_cuts_the_observed_sample_with_the_model_and_compares_the_reported_medians():
    t = target()
    sampler = wd.observed_sampler(lambda i: draw.gauss(-1.25, 0.2))
    items = th.simulate(t, sampler, 600)
    assert 0.80 < len(items) / 600 < 0.98, "about 10 percent of bodies give a limit only and are dropped"
    assert all(item["recovery"] != "limit" and item["diw"] is not None for item in items)
    result = th.compare(t, items)
    assert result.passed and abs(result.simulated + 1.25) < 0.2


def test_a_tuned_width_is_compared_with_the_reported_sample_and_not_the_bodies():
    """Bodies at one dIW still give a spread of reported medians, wider than the bodies, and the sd target sees it."""
    sd_target = target(id="wd-mean", kind="sd", statistic="sample standard deviation of dIW over the same stars",
                       value=0.447, sigma=0.069)
    sampler = wd.observed_sampler(lambda i: -1.25)
    items = th.simulate(sd_target, sampler, 800)
    values = th.values_of(sd_target, items)
    assert statistics.stdev(values) > 0.3, "the measurement alone spreads the reported medians"
    assert statistics.stdev(values) == pytest.approx(wd.expected_sample(-1.25, 0.0).reported_sd, abs=0.08)


def test_the_error_model_changes_nothing_in_the_shipped_recipes_or_targets():
    redox = [t for t in th.TARGETS if t.row_id in ("a3", "a4")]
    assert {t.id for t in redox} == {"a3-mean", "a3-sd", "a4-mean", "a4-sd"}
    assert all(t.selection_model == "wd-pollution-detection" and t.field == "diw" for t in redox)
    document = copy.deepcopy(th.TARGETS.document)
    assert th.TargetTable(document).digest() == th.TARGETS.digest()


# --- the redox rows (GEN.513) ------------------------------------------------------------------------------------------------


def _a4_recipes(sigma):
    import copy
    from planetgen.physics.atmosphere_recipes import RECIPES, Recipes
    document = copy.deepcopy(RECIPES.document)
    next(e for e in document["entries"] if e["id"] == "a4")["dist"]["sigma"] = sigma
    return Recipes(document, register=False)


def _compare(target_id, sigma, model, n=3000):
    from planetgen.generation import tuning_harness as th
    from planetgen.generation.possibility_table import Subject
    subject = Subject("terrestrial", {"host_type": "WD"})
    sampler = wd.observed_sampler(th.recipe_sampler("a4", subject, _a4_recipes(sigma)), model)
    return th.simulate_and_compare(th.TARGETS.get(target_id), sampler, n=n, quantities=subject)


def test_the_shipped_width_is_the_intrinsic_spread_of_the_well_constrained_stars_and_inside_the_accepted_range():
    from planetgen.physics.atmosphere_recipes import RECIPES
    for row in ("a3", "a4"):
        dist = RECIPES.entry(row).spec["dist"]
        assert dist["mu"] == pytest.approx(-1.25) and dist["sigma"] == pytest.approx(0.22)
    # the deterministic expected sample: the reported spread at the shipped width is inside 2 sigma of the published 0.447
    bounded = wd.fit_error_model(flagged=False)
    assert abs(wd.expected_sample(-1.25, 0.22, bounded).reported_sd - 0.447) < 2 * 0.069
    assert wd.expected_sample(-1.25, 0.0, bounded).reported_sd > 0.447 - 2 * 0.069, "zero width is accepted too: the width is not identified"
    assert wd.expected_sample(-1.25, 0.45, bounded).reported_sd > 0.447 + 2 * 0.069, "the earlier width 0.45 is too wide"


def test_the_spread_is_tuned_with_the_error_of_the_stars_with_both_bounds_and_the_whole_sample_model_over_disperses():
    bounded = wd.fit_error_model(flagged=False)
    assert _compare("a4-sd", 0.22, bounded).passed
    assert not _compare("a4-sd", 0.22, wd.ERROR_MODEL).passed, "the whole-sample error model fails the spread at every width"
    assert not _compare("a4-sd", 0.01, wd.ERROR_MODEL).passed
    assert not _compare("a4-sd", 0.45, bounded).passed, "the earlier width"
    assert _compare("a4-mean", 0.22, wd.ERROR_MODEL).passed and _compare("a4-mean", 0.22, bounded).passed


# --- the power of the redox targets (GEN.665) --------------------------------------------------------------------------------


def test_the_redox_targets_refuse_a_centre_wrong_by_0_3_dex_and_a_spread_wrong_by_a_factor_of_two_but_not_a_zero_width():
    bounded = wd.fit_error_model(flagged=False)
    mean = th.TARGETS.get("a4-mean")
    assert _compare("a4-mean", 0.22, bounded).passed
    for shift in (-0.3, 0.3):
        assert abs((mean.value + shift) - mean.value) > th.SIGMAS * mean.sigma, "a centre 0.3 dex off is more than 2 sigma of the 22 values"
    sd = th.TARGETS.get("a4-sd")
    for wrong in (sd.value * 2, sd.value / 2):
        assert abs(wrong - sd.value) > th.SIGMAS * sd.sigma
    assert _compare("a4-sd", 0.01, bounded).passed, "the spread does not tell a near-zero intrinsic width from 0.22: the sample of 7 does not identify it"
    assert not _compare("a4-sd", 0.45, bounded).passed
