# tests/test_radius_valley.py

"""GEN.509: the radius valley of a simulated Kepler survey and the tuning of rows e1 and e2 (`planetgen.physics.radius_valley`)."""

import math

import pytest

from planetgen.generation import tuning_harness as th
from planetgen.generation.possibility_table import Subject
from planetgen.physics import radius_valley as rv
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.util import draw

CENTRE = 0.045
SURVEYS = 24


def _compare(target_id, f_env0, eta, field):
    target = th.TARGETS.get(target_id)
    return th.simulate_and_compare(target, rv.valley_sampler(f_env0, eta, field), n=SURVEYS)


def test_a_survey_repeats_on_a_seed_and_finds_a_valley_between_one_and_three_earth_radii():
    with draw.bound(draw.Stream("valley:1")):
        first = rv.valley_of(rv.survey(CENTRE, rv.LOSS_PER_ENVELOPE * CENTRE))
    with draw.bound(draw.Stream("valley:1")):
        second = rv.valley_of(rv.survey(CENTRE, rv.LOSS_PER_ENVELOPE * CENTRE))
    assert first == second and first is not None
    slope, radius_10d = first
    assert slope < 0 and 1.4 < radius_10d < 2.4


def test_a_sparse_survey_has_no_valley():
    assert rv.valley_of([(5.0, 2.0, 1.0)] * 10) is None


def test_the_valley_slope_and_location_follow_the_published_values_along_the_tuned_ridge():
    for f_env0 in (0.035, CENTRE, 0.08):
        assert _compare("e1-valley-slope", f_env0, rv.LOSS_PER_ENVELOPE * f_env0, "slope").passed, f_env0
        assert _compare("e1-valley-radius", f_env0, rv.LOSS_PER_ENVELOPE * f_env0, "radius_10d").passed, f_env0


def test_a_weak_loss_puts_the_valley_too_low_and_a_strong_one_too_high():
    assert not _compare("e1-valley-radius", CENTRE, 0.1 * CENTRE, "radius_10d").passed, "no photoevaporation: the valley sits near 1.6"
    assert not _compare("e1-valley-radius", CENTRE, 40 * CENTRE, "radius_10d").passed


def test_the_recipe_brackets_the_ridge_and_both_rows_are_marked_as_tuned():
    assert RECIPES.document["groups"]["e"]["constants"]["f_env0"][:2] == [0.035, 0.08]
    entry = RECIPES.entry("e1")
    assert entry.tunable and entry.state == "partly verified" and "degenerate ridge" in entry.limit
    assert {t.row for t in th.TARGETS.for_row("atmosphere-e:e1")} >= {"atmosphere-e:e1"}


def test_the_row_formula_is_the_survey_formula_for_the_same_planet():
    expression = RECIPES.entry("e1").lets["f_lost"].text
    assert "10*f_env0" in expression and "(15/v_esc_kms)**2" in expression and "3.1e16" in expression
    core_mass, insolation, f_env0 = 4.0, 60.0, 0.04
    core_radius = core_mass ** 0.25
    v_esc = math.sqrt(2 * rv.G * core_mass * rv.M_EARTH / (core_radius * rv.R_EARTH)) / 1000.0
    row = 10 * f_env0 * (15 / v_esc) ** 2 * math.pi * (core_radius * rv.R_EARTH) ** 3 * 3.1e16 * insolation / (rv.G * (core_mass * rv.M_EARTH) ** 2)
    survey = (rv.LOSS_PER_ENVELOPE * f_env0 * (v_esc / 15.0) ** -2 * math.pi * (core_radius * rv.R_EARTH) ** 3
              * rv.XUV_ENERGY_J_M2 * insolation / (rv.G * (core_mass * rv.M_EARTH) ** 2))
    assert row == pytest.approx(survey, rel=1e-9)
