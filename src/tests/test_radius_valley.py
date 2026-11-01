# tests/test_radius_valley.py

"""GEN.509: the radius valley of a simulated Kepler survey and the tuning of rows e1 and e2 (`planetgen.physics.radius_valley`)."""

import math

import pytest

from planetgen.generation import tuning_harness as th
from planetgen.generation.possibility_table import Subject
from planetgen.physics import envelope_escape, radius_valley as rv
from planetgen.physics.atmosphere_recipes import RECIPES
from planetgen.util import draw

CENTRE = 0.1
SURVEYS = 24


MASS_SURVEYS = 8
"""int: Surveys per comparison of a stellar-mass slope (each survey is 12000 planets, three times the period survey)."""


def _compare(target_id, f_env0, field):
    target = th.TARGETS.get(target_id)
    return th.simulate_and_compare(target, rv.valley_sampler(f_env0, field), n=SURVEYS)


def _compare_mass(target_id, f_env0, field):
    target = th.TARGETS.get(target_id)
    return th.simulate_and_compare(target, rv.mass_slope_sampler(f_env0, field), n=MASS_SURVEYS)


def test_a_survey_repeats_on_a_seed_and_finds_a_valley_between_one_and_three_earth_radii():
    with draw.bound(draw.Stream("valley:1")):
        first = rv.valley_of(rv.triples(rv.survey(CENTRE)))
    with draw.bound(draw.Stream("valley:1")):
        second = rv.valley_of(rv.triples(rv.survey(CENTRE)))
    assert first == second and first is not None
    slope, radius_10d = first
    assert slope < 0 and 1.4 < radius_10d < 2.4


def test_a_sparse_survey_has_no_valley():
    assert rv.valley_of([(5.0, 2.0, 1.0)] * 10) is None


def test_the_valley_slope_and_location_follow_the_published_values_across_the_recipe_bracket():
    for f_env0 in (0.05, CENTRE, 0.2):
        assert _compare("e1-valley-slope", f_env0, "slope").passed, f_env0
        assert _compare("e1-valley-radius", f_env0, "radius_10d").passed, f_env0


def test_a_weaker_or_stronger_loss_moves_the_valley_out_of_the_published_range():
    """No fit to a second knob here: the onset age is the one that moves the loss (module constant, patched for the test)."""
    for onset in (0.3, 3.0):
        old = envelope_escape.ONSET_GY
        envelope_escape.ONSET_GY = onset
        try:
            assert not _compare("e1-valley-radius", CENTRE, "radius_10d").passed, onset
        finally:
            envelope_escape.ONSET_GY = old


def test_the_recipe_brackets_the_tuned_range_and_both_rows_are_marked_as_tuned():
    assert RECIPES.document["groups"]["e"]["constants"]["f_env0"][:2] == [0.05, 0.2]
    entry = RECIPES.entry("e1")
    assert entry.tunable and entry.state == "partly verified" and "onset age" in entry.limit
    assert {t.row for t in th.TARGETS.for_row("atmosphere-e:e1")} >= {"atmosphere-e:e1"}


def test_the_row_formula_is_the_function_of_the_module_for_the_same_planet():
    expression = RECIPES.entry("e1").lets["f_lost"].text
    assert expression.startswith("envelope_lost_fraction(") and "3.1e16" not in expression
    row = RECIPES.entry("e2").lets["f_lost"].text
    assert row == expression, "the two rows weigh the same share, one by it and one by its complement"


def test_the_valley_against_insolation_and_stellar_mass_follows_berger_2023_across_the_bracket():
    """GEN.722 and GEN.822: the hydro-based loss, with no stellar-mass dependence of its own beyond the star's XUV history, still puts alpha against
    insolation and beta against stellar mass inside the measurements of [R0354] (the paper's two mechanisms are both within 1.1 sigma of
    beta, so the survey cannot choose one; the design note, section 4)."""
    for f_env0 in (0.05, CENTRE, 0.2):
        assert _compare_mass("e1-valley-flux-slope", f_env0, "flux_slope").passed, f_env0
        assert _compare_mass("e1-valley-mass-slope", f_env0, "mass_slope").passed, f_env0


def test_a_survey_over_a_wider_range_of_hosts_draws_them_from_the_range():
    with draw.bound(draw.Stream("valley:hosts")):
        found = rv.survey(CENTRE, 400, rv.WIDE_HOST_MASS_RANGE_SOLAR)
    masses = [c.host_mass_solar for c in found]
    assert min(masses) >= 0.5 and max(masses) <= 1.4 and min(masses) < 0.6 and max(masses) > 1.3
    assert all(c.insolation_earth > 0 and c.radius_earth > 0 for c in found)
    assert rv.mass_slopes_of(found[:20]) is None
