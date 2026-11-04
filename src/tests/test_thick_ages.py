# tests/test_thick_ages.py

"""GEN.1025: the old population's ages reach the thick disc's 13 Gy (`planetgen/galaxy/thick_ages.py`)."""

import random

import pytest

from planetgen import tuning
from planetgen.galaxy import thick_ages
from planetgen.galaxy.age_law import AgeLaw

THIN = AgeLaw((3.0, 6.5, 10.0), (0.5, 1.0, 0.8))


def test_with_no_thick_disc_the_law_is_the_thin_laws_own():
    assert thick_ages.mixed_old_law(THIN, 0.0) is THIN


def test_the_mixture_reaches_the_oldest_thick_disc_stars_and_keeps_the_thin_discs_start():
    law = thick_ages.mixed_old_law(THIN, 0.2)
    assert law.window == (3.0, tuning.THICK_DISC_AGE_RANGE_GY[1])


def test_the_share_of_stars_older_than_the_thin_disc_is_the_thick_discs_share_of_its_window():
    low, high = tuning.THICK_DISC_AGE_RANGE_GY
    for share in (0.05, 0.2, 0.5):
        law = thick_ages.mixed_old_law(THIN, share)
        beyond = law.probability(10.0 + thick_ages.EDGE_GY, high)
        assert beyond == pytest.approx(share * (high - 10.0 - thick_ages.EDGE_GY) / (high - low), abs=0.003), share


def test_the_thick_discs_stars_start_at_its_branch_and_none_are_older_than_thin_before_it():
    law = thick_ages.mixed_old_law(THIN, 0.3)
    thin_only = AgeLaw(THIN.ages, THIN.rates)
    # below the thick disc's start the shape is the thin law's, scaled by its weight
    ratio = law.probability(3.0, 5.0) / thin_only.probability(3.0, 5.0)
    assert ratio == pytest.approx(0.7, rel=0.01)


def test_samples_fall_in_the_window_and_the_old_tail_is_populated():
    law = thick_ages.mixed_old_law(THIN, 0.25)
    rng = random.Random(7)
    ages = [law.sample(3.0, 13.0, rng) for _ in range(20000)]
    assert min(ages) >= 3.0 and max(ages) <= 13.0
    tail = sum(age > 10.0 for age in ages) / len(ages)
    assert tail == pytest.approx(0.25 * 3.0 / 5.0, abs=0.01)


def test_the_thick_share_is_the_thick_units_over_the_old_units_and_held_to_zero_and_one():
    assert thick_ages.thick_share_of_old(10.0, 2.0) == pytest.approx(0.2)
    assert thick_ages.thick_share_of_old(0.0, 2.0) == 0.0
    assert thick_ages.thick_share_of_old(1.0, 5.0) == 1.0


def test_a_milky_way_galaxy_has_old_stars_from_10_to_13_gy():
    from tests.test_disk_history import SHAPE
    from planetgen.galaxy import density, disk_history
    share = thick_ages.thick_share_of_old(density.population_units(SHAPE)["old"], density.thick_disk_units(SHAPE))
    assert 0.0 < share < 0.6
    law = disk_history.population_age_laws(SHAPE)["old"]
    assert thick_ages.old_age_law(SHAPE, disk_history.history_of(SHAPE).age_law(*tuning.THIN_DISC_OLD_WINDOW_GY)).window == law.window
    assert law.window == (3.0, tuning.THICK_DISC_AGE_RANGE_GY[1])
    assert law.probability(10.0 + thick_ages.EDGE_GY, 13.0) == pytest.approx(share * 0.6, abs=0.01)
    # the thin disc part is unchanged where the thick disc is absent
    assert law.probability(3.0, 7.0) > 0.4
