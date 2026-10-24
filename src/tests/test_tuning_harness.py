# tests/test_tuning_harness.py

"""
GEN.503: the shared tuning harness (`planetgen.generation.tuning_harness`): the target table and its refusals, the loader
check that a tunable row has a target, the comparison within 2 sigma (a mean, a standard deviation, a median, a proportion by
Wilson intervals), the validity range and the flagged extrapolation, the Solar System flag, the selection-model hook, and the
worked example (the white-dwarf redox rows a3 and a4 drawn at many seeds against their published statistics).
"""

import copy
import math

import pytest

from planetgen.generation import possibility_table as pt
from planetgen.generation import tuning_harness as th
from planetgen.generation.possibility_table import Row, Subject, Table
from planetgen.generation.tuning_harness import TargetError, TargetTable
from planetgen.physics import atmosphere_recipes as ar
from planetgen.physics.atmosphere_recipes import RECIPES, Recipes
from planetgen.util import draw


def target(**overrides):
    raw = {"id": "t1", "row": "tab:r", "statistic": "mean of x over stars", "kind": "mean", "value": 1.0, "sigma": 0.1,
           "n": 20, "source": "Test et al. 2020 Table 1", "label": "READ", "validity": {"host_type": "WD"}}
    raw.update(overrides)
    return raw


def table(*targets, version=1):
    return TargetTable({"version": version, "targets": list(targets)})


def one(**overrides):
    return table(target(**overrides)).get("t1")


def refused(match, **overrides):
    with pytest.raises(TargetError, match=match):
        table(target(**overrides))


# --- the shipped table ------------------------------------------------------------------------------------------------


def test_the_shipped_table_holds_the_white_dwarf_redox_targets_and_each_names_its_fill_and_falsifier():
    targets = th.TARGETS
    assert {t.id for t in targets} == {"a3-mean", "a3-sd", "a4-mean", "a4-sd"}
    for t in targets:
        assert t.tunes and t.label in th.USABLE_LABELS and t.fill_id.startswith("GEN.") and t.falsifier
        assert t.outside({"host_type": "WD"}) == []
        assert t.outside({"host_type": "G"}) == ["host_type"]


def test_the_shipped_digest_is_stable_and_changes_with_a_value():
    document = copy.deepcopy(th.TARGETS.document)
    assert TargetTable(document).digest() == th.TARGETS.digest()
    document["targets"][0]["value"] += 0.01
    assert TargetTable(document).digest() != th.TARGETS.digest()


# --- the table refuses -----------------------------------------------------------------------------------------------


def test_a_table_needs_a_version_and_targets_and_no_duplicate_id():
    with pytest.raises(TargetError, match="version"):
        TargetTable({"targets": [target()]})
    with pytest.raises(TargetError, match="no targets"):
        TargetTable({"version": 1, "targets": []})
    with pytest.raises(TargetError, match="twice"):
        table(target(), target())


@pytest.mark.parametrize("overrides, match", [
    ({"source": ""}, "no source"),
    ({"statistic": " "}, "no statistic"),
    ({"row": "r"}, "not '<table>:<row>'"),
    ({"kind": "mode"}, "kind"),
    ({"label": "HEARSAY"}, "label"),
    ({"sigma": 0.0}, "sigma must be positive"),
    ({"sigma": float("nan")}, "not a finite number"),
    ({"value": "1"}, "not a finite number"),
    ({"n": 1}, "whole number of 2"),
    ({"validity": {}}, "validity range"),
    ({"validity": {"fe_h": [0.5, -0.5]}}, "not \\[minimum, maximum\\]"),
    ({"selection_model": "magic"}, "selection model"),
    ({"role": "guess"}, "role"),
    ({"colour": "red"}, "unknown fields"),
    ({"counts": [1, 4]}, "only a proportion"),
])
def test_a_target_with_a_gap_is_refused_at_load(overrides, match):
    refused(match, **overrides)


def test_a_proportion_states_its_count_and_the_value_follows_from_it():
    refused("count", kind="proportion", value=0.5, sigma=0.1)
    refused("0 <= k <= n", kind="proportion", counts=[5, 4])
    refused("is not the count", kind="proportion", counts=[1, 4], value=0.5)
    raw = target(kind="proportion", counts=[13, 25])
    raw.pop("value"), raw.pop("sigma")
    t = table(raw).get("t1")
    assert t.value == pytest.approx(0.52) and 0.05 < t.sigma < 0.15


def test_a_solar_system_value_is_a_test_point_with_a_flag_and_never_tuned_to():
    refused("Solar System body", statistic="dIW of the Earth upper mantle")
    refused("Solar System body", source="Halley's comet table")
    refused("never a target", solar_system="Earth only")
    refused("states why", role="test-point")
    t = one(statistic="dIW of the Earth upper mantle", role="test-point", solar_system="measured on Earth only")
    assert not t.tunes
    assert one().tunes


def test_a_snippet_or_a_recalled_value_is_listed_and_never_tuned_to():
    for label in ("SNIPPET", "RECALLED", "NOT REACHED"):
        assert not one(label=label).tunes
    for label in ("READ", "DATA", "COMPUTED", "2ND-HAND"):
        assert one(label=label).tunes


# --- the loader check ------------------------------------------------------------------------------------------------


def tab(*rows):
    return Table("tab", 1, rows)


def trow(name, tunable):
    return Row(name, ("rocky",), 1.0, "Test et al. 2020", tunable=tunable)


def test_a_tunable_row_without_a_target_is_refused_naming_the_row():
    targets = table(target())
    targets.require([tab(trow("r", True))])
    with pytest.raises(TargetError, match="tab:s"):
        targets.require([tab(trow("r", True), trow("s", True))])


def test_a_row_whose_only_targets_are_test_points_or_untunable_labels_still_has_no_target():
    flagged = table(target(role="test-point", statistic="Mars value", solar_system="Mars only"))
    with pytest.raises(TargetError, match="tab:r"):
        flagged.require([tab(trow("r", True))])
    snippet = table(target(label="SNIPPET"))
    with pytest.raises(TargetError, match="tab:r"):
        snippet.require([tab(trow("r", True))])


def test_a_target_for_a_row_that_is_not_tunable_or_does_not_exist_is_refused():
    targets = table(target())
    with pytest.raises(TargetError, match="not marked tunable"):
        targets.require([tab(trow("r", False))])
    with pytest.raises(TargetError, match="no such row"):
        targets.require([tab(trow("q", False))])
    targets.require([Table("other", 1, [Row("r", ("rocky",), 1.0, "Test et al. 2020")])])


def test_a_test_point_may_name_a_row_that_is_not_tunable():
    table(target(role="test-point", statistic="Mars value", solar_system="Mars only")).require([tab(trow("r", False))])


def test_the_registry_check_passes_for_every_registered_table_and_refuses_a_new_unbacked_tunable_row():
    assert RECIPES  # the atmosphere tables are registered
    pt.REGISTRY.validate(targets=th.TARGETS)
    registry = pt.Registry()
    registry.register(tab(trow("r", True), trow("s", True)))
    with pytest.raises(TargetError, match="tab:s"):
        registry.validate(targets=table(target()))


def test_the_recipes_loader_refuses_a_tunable_row_with_no_target():
    document = copy.deepcopy(RECIPES.document)
    next(e for e in document["entries"] if e["id"] == "a5")["tunable"] = True
    recipes = Recipes(document, register=False)
    with pytest.raises(TargetError, match="atmosphere-a:a5"):
        th.TARGETS.require(recipes.tables)
    assert {e.id for e in RECIPES.entries if e.tunable} == {"a3", "a4"}


# --- the comparison --------------------------------------------------------------------------------------------------


def normal_values(mean, sd, n, label):
    stream = draw.Stream(label)
    return [stream.gauss(mean, sd) for _ in range(n)]


def test_a_mean_within_two_sigma_passes_and_one_outside_fails():
    t = one(value=0.0, sigma=0.1)
    ok = th.compare(t, normal_values(0.05, 1.0, 4000, "ok"))
    assert ok.passed and abs(ok.z) <= 2 and ok.draws == 4000 and ok.tunes
    far = th.compare(t, normal_values(0.6, 1.0, 4000, "far"))
    assert not far.passed and far.z > 2
    assert far.chi_square == pytest.approx(far.z ** 2)


def test_the_sigma_combines_the_published_uncertainty_and_the_sampling_error():
    t = one(value=0.0, sigma=0.3)
    values = normal_values(0.0, 2.0, 100, "combine")
    result = th.compare(t, values)
    sampling = result.sigma ** 2 - 0.3 ** 2
    assert math.sqrt(sampling) == pytest.approx(2.0 / math.sqrt(100), rel=0.2)


def test_a_standard_deviation_and_a_median_are_compared_the_same_way():
    sd = one(kind="sd", value=0.45, sigma=0.07)
    assert th.compare(sd, normal_values(0.0, 0.46, 3000, "sd-ok")).passed
    assert not th.compare(sd, normal_values(0.0, 0.8, 3000, "sd-far")).passed
    median = one(kind="median", value=10.0, sigma=0.5)
    assert th.compare(median, normal_values(10.2, 3.0, 3000, "med-ok")).passed
    assert not th.compare(median, normal_values(13.0, 3.0, 3000, "med-far")).passed


def test_a_proportion_is_compared_by_wilson_intervals_on_both_sides():
    raw = target(kind="proportion", counts=[13, 25])
    raw.pop("value"), raw.pop("sigma")
    t = table(raw).get("t1")
    assert th.compare(t, [], count=(520, 1000)).passed
    assert th.compare(t, [True] * 55 + [False] * 45).passed
    assert not th.compare(t, [], count=(900, 1000)).passed
    assert not th.compare(t, [], count=(50, 1000)).passed
    with pytest.raises(TargetError, match="no simulated draws"):
        th.compare(t, [])


def test_the_wilson_interval_has_the_published_values_and_keeps_inside_0_and_1():
    low, high = th.wilson_interval(13, 25, z=1.96)
    assert low == pytest.approx(0.3355, abs=1e-3) and high == pytest.approx(0.6989, abs=1e-3)
    assert th.wilson_interval(0, 10)[0] == 0.0 and th.wilson_interval(10, 10)[1] == pytest.approx(1.0)
    with pytest.raises(ValueError):
        th.wilson_interval(5, 4)


def test_a_comparison_needs_two_values_and_flags_an_extrapolation():
    t = one()
    with pytest.raises(TargetError, match="at least 2"):
        th.compare(t, [1.0])
    inside = th.compare(t, [1.0, 1.1], quantities={"host_type": "WD"})
    outside = th.compare(t, [1.0, 1.1], quantities=Subject("E", {"host_type": "G"}))
    assert inside.extrapolated == () and outside.extrapolated == ("host_type",)
    missing = th.compare(t, [1.0, 1.1], quantities={})
    assert missing.extrapolated == ("host_type",)


def test_a_numeric_validity_range_has_open_ends():
    t = one(validity={"fe_h": [-0.5, None], "teff_k": [3000.0, 4000.0]})
    assert t.outside({"fe_h": 0.2, "teff_k": 3500.0}) == []
    assert t.outside({"fe_h": -0.6, "teff_k": 4500.0}) == ["fe_h", "teff_k"]


# --- the selection model ---------------------------------------------------------------------------------------------


def test_a_target_that_names_a_selection_model_nobody_built_refuses_to_simulate():
    t = one(selection_model="kepler-recovery")
    with pytest.raises(TargetError, match="not registered"):
        th.simulate(t, lambda i: 1.0, 5)


def test_a_registered_selection_model_cuts_the_sample_by_its_probability(monkeypatch):
    monkeypatch.setitem(th.SELECTION_MODELS, "kepler-recovery", lambda value: 1.0 if value > 0 else 0.0)
    t = one(selection_model="kepler-recovery")
    values = th.simulate(t, lambda i: float(i % 2) - 0.5, 200)
    assert len(values) == 100 and all(v > 0 for v in values)
    with pytest.raises(TargetError):
        th.register_selection_model("magic", lambda v: 1.0)


def test_a_sampler_that_returns_mappings_is_read_at_the_targets_field_and_cut_by_the_model_on_the_whole_item(monkeypatch):
    monkeypatch.setitem(th.SELECTION_MODELS, "kepler-recovery", lambda item: 1.0 if item["period_d"] < 100 else 0.0)
    t = one(selection_model="kepler-recovery", field="radius_re", value=2.0, sigma=0.1)
    items = th.simulate(t, lambda i: {"period_d": 10.0 + 20.0 * (i % 10), "radius_re": 2.0}, 100)
    assert len(items) == 50 and all(item["period_d"] < 100 for item in items)
    result = th.compare(t, items)
    assert result.simulated == 2.0 and result.passed
    with pytest.raises(TargetError, match="names no field"):
        th.compare(one(), items)


# --- simulate --------------------------------------------------------------------------------------------------------


def test_simulate_gives_each_draw_its_own_seeded_stream_and_repeats():
    t = one()
    first = th.simulate(t, lambda i: draw.random(), 50)
    again = th.simulate(t, lambda i: draw.random(), 50)
    other = th.simulate(one(id="t1", statistic="mean of y over stars"), lambda i: draw.random(), 50)
    assert first == again and len(set(first)) == 50 and first == other
    assert th.simulate(table(target(id="t2")).get("t2"), lambda i: draw.random(), 50) != first


def test_simulate_skips_a_draw_whose_gate_did_not_pass():
    values = th.simulate(one(), lambda i: None if i % 2 else 1.0, 10)
    assert values == [1.0] * 5


# --- the worked example ----------------------------------------------------------------------------------------------


@pytest.mark.parametrize("target_id, body_class", [("a3-mean", "C"), ("a3-sd", "C"), ("a4-mean", "E"), ("a4-sd", "E")])
def test_the_tunable_redox_rows_drawn_at_many_seeds_match_their_published_statistics(target_id, body_class):
    t = th.TARGETS.get(target_id)
    subject = Subject(body_class, {"host_type": "WD"})
    result = th.simulate_and_compare(t, th.recipe_sampler(t.row_id, subject), n=2000, quantities=subject)
    assert result.draws == 2000 and result.extrapolated == ()
    assert result.passed, (target_id, result.simulated, result.z)


def test_a_knob_retuned_away_from_the_target_fails_the_comparison():
    document = copy.deepcopy(RECIPES.document)
    entry = next(e for e in document["entries"] if e["id"] == "a4")
    entry["dist"]["mu"] = -0.7
    recipes = Recipes(document, register=False)
    t = th.TARGETS.get("a4-mean")
    subject = Subject("E", {"host_type": "WD"})
    result = th.simulate_and_compare(t, th.recipe_sampler("a4", subject, recipes), n=2000)
    assert not result.passed and result.z > 2


def test_the_harness_holds_no_solar_system_value_as_a_target_without_the_flag():
    for t in th.TARGETS:
        if t.role == "target":
            assert not th._solar_system_hits(t.statistic) and not th._solar_system_hits(t.source)
        else:
            assert t.solar_system
