# tests/test_tuned_draws.py

"""
The tuned draws (GEN.509 to GEN.525): every target of the tuning-target table that names the body it is drawn for and needs no selection
model is checked against the shipped recipe (`tuning_harness.check_target`: the row drawn at 2000 to 3000 seeds, then compared within 2
sigma), the knobs that were retuned are pinned, and a retune away from the data fails.
"""

import copy

import pytest

from planetgen.generation import tuning_harness as th
from planetgen.physics.atmosphere_recipes import RECIPES, Recipes

CHECKED = [t.id for t in th.TARGETS if t.role == "target" and t.subject is not None and t.selection_model is None]


def retuned(row_id, mutate):
    document = copy.deepcopy(RECIPES.document)
    mutate(next(e for e in document["entries"] if e["id"] == row_id))
    return Recipes(document, register=False)


@pytest.mark.parametrize("target_id", CHECKED)
def test_the_shipped_recipe_reproduces_the_published_statistic(target_id):
    target = th.TARGETS.get(target_id)
    result = th.check_target(target, n=3000)
    assert result.extrapolated == (), "the subject lies inside the validity range of the data"
    assert result.passed, (target_id, result.simulated, target.value, result.z)


def test_every_target_with_a_subject_is_checked_here_or_through_its_selection_model():
    assert CHECKED, "no target is checked"
    for target in th.TARGETS:
        if target.subject is not None:
            assert target.selection_model is None or target.selection_model in th.SELECTION_MODEL_ITEMS


# --- giants (GEN.516) --------------------------------------------------------------------------------------------------------


def test_the_giant_metallicity_follows_the_later_larger_fit_and_not_the_first_relation():
    spec = RECIPES.entry("f5").spec["dist"]
    assert spec["mu"] == "log10(2.6/mass_jupiter + 3.3)" and spec["sigma"] == 0.26
    old = retuned("f5", lambda e: e["dist"].update(mu="log10(9.7*mass_jupiter**-0.45)"))
    assert not th.check_target(th.TARGETS.get("f5-z-1p0"), recipes=old, n=3000).passed, "the first relation is 1.6 times the fit at one M_J"
    assert th.check_target(th.TARGETS.get("f5-z-0p1"), recipes=old, n=3000).passed, "and the two meet at 0.1 M_J"


def test_the_giant_core_width_comes_from_the_population_and_not_from_one_planet():
    assert RECIPES.entry("f7").spec["dist"]["core"]["sigma"] == 0.21
    narrow = retuned("f7", lambda e: e["dist"]["core"].update(sigma=0.08))
    assert not th.check_target(th.TARGETS.get("f7-low"), recipes=narrow, n=3000).passed, "12 of 69 planets are below 0.35"
    wide_tail = retuned("f7", lambda e: e["constants"].update(f_tail=[0.3, 0.34, "log"]))
    assert not th.check_target(th.TARGETS.get("f7-tail"), recipes=wide_tail, n=3000).passed, "the old tail bracket is outside the data"


# --- the harness's statistic fields ------------------------------------------------------------------------------------------


def test_a_log_scale_a_threshold_and_a_subject_are_read_from_the_target():
    log_sd = th.TARGETS.get("f5-width")
    assert log_sd.scale == "log10" and log_sd.subject == {"body_class": "I", "quantities": {"mass_jupiter": 1.0}}
    assert th.compare(log_sd, [1.0, 10.0, 100.0, 10.0]).simulated == pytest.approx(0.8165, abs=1e-3), "sample sd of 0, 1, 2, 1"
    below = th.TARGETS.get("f7-low")
    assert below.below and below.threshold == 0.35
    result = th.compare(below, [0.1, 0.2, 0.5, 0.6], count=None)
    assert result.simulated == 0.5
    above = th.TARGETS.get("f7-tail")
    assert not above.below and th.compare(above, [0.5, 0.9, 1.2, 0.3]).simulated == 0.5
    with pytest.raises(th.TargetError, match="no logarithm"):
        th.compare(log_sd, [1.0, -2.0])
    assert th.subject_of(log_sd).body_class == "I" and th.subject_of(th.TARGETS.get("a4-mean")) is not None


def test_the_table_refuses_a_malformed_scale_threshold_or_subject():
    base = dict(th.TARGETS.document["targets"][0])
    for key, value, match in (("scale", "ln", "scale"), ("threshold", 0.5, "threshold"), ("below", True, "below"),
                              ("subject", {"quantities": {}}, "subject"), ("subject", "x", "subject")):
        document = {"version": 1, "targets": [dict(base, **{key: value})]}
        with pytest.raises(th.TargetError, match=match):
            th.TargetTable(document)
    proportion = {"id": "p", "row": "t:r", "statistic": "share", "kind": "proportion", "counts": [1, 10], "source": "s",
                  "label": "DATA", "validity": {"q": [0, 1]}}
    th.TargetTable({"version": 1, "targets": [dict(proportion, threshold=0.5, below=True)]})
    with pytest.raises(th.TargetError, match="needs a threshold"):
        th.TargetTable({"version": 1, "targets": [dict(proportion, below=True)]})
    with pytest.raises(th.TargetError, match="names no subject"):
        th.check_target(th.TargetTable({"version": 1, "targets": [proportion]}).get("p"))
