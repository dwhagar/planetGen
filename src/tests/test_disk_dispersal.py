# src/tests/test_disk_dispersal.py

"""GEN.1551 slice 0: a star's own disk dispersal time (`physics/disk_dispersal.py`), the law and its draw, not yet wired."""

import math

import pytest

from planetgen.physics import disk_dispersal as dd

MYR = 1.0e6
UNIFORM = [(i + 0.5) / 20000 for i in range(20000)]      # an even grid of draws, so no check can flake


def _share_with_disk(mass, age_yr, close=False):
    return sum(dd.dispersal_time_yr(mass, u, close) > age_yr for u in UNIFORM) / len(UNIFORM)


def test_the_time_of_a_star_is_an_exponential_of_the_characteristic_time_of_its_mass():
    assert dd.characteristic_time_yr(0.5) == pytest.approx(2.5 * MYR)
    assert dd.characteristic_time_yr(0.01) == pytest.approx(3.0 * MYR), "a brown dwarf: the long end"
    assert dd.characteristic_time_yr(5.0) == pytest.approx(1.0 * MYR), "over 1.3 solar masses: the short end"
    times = [dd.characteristic_time_yr(m) for m in (0.08, 0.1, 0.3, 0.5, 0.8, 1.0, 1.3)]
    assert times == sorted(times, reverse=True), "heavier stars lose their disks sooner"
    assert dd.dispersal_time_yr(0.5, math.exp(-1.0)) == pytest.approx(2.5 * MYR)
    assert dd.dispersal_time_yr(0.5, 1.0) == 0.0


def test_the_share_of_stars_with_a_disk_follows_the_survey_exponential():
    for age_myr in (1.0, 2.5, 3.0, 5.5, 10.0):
        assert _share_with_disk(0.5, age_myr * MYR) == pytest.approx(math.exp(-age_myr / 2.5), abs=0.002)
        assert dd.disk_fraction(0.5, age_myr * MYR) == pytest.approx(math.exp(-age_myr / 2.5))
    assert dd.disk_fraction(0.5, 1.7 * MYR) == pytest.approx(0.5, abs=0.02), "the half-life of the survey, 1.7 Myr"
    assert dd.disk_fraction(0.5, 10.0 * MYR) < 0.03, "near none by 10 Myr"
    assert dd.disk_fraction(0.5, 0.0) == 1.0


def test_two_thirds_of_the_stars_with_a_close_companion_lose_the_disk_within_1_myr_and_none_keep_it_longer():
    assert _share_with_disk(0.5, 1.0 * MYR, close=True) == pytest.approx(1.0 / 3.0, abs=0.002)
    assert _share_with_disk(0.5, 1.0 * MYR, close=False) == pytest.approx(math.exp(-1.0 / 2.5), abs=0.002)
    for mass in (0.1, 0.5, 1.0, 2.0):
        for u in (0.001, 0.1, 0.5, 0.9):
            assert dd.dispersal_time_yr(mass, u, close_companion=True) <= dd.dispersal_time_yr(mass, u)
    assert dd.dispersal_time_yr(0.5, 0.3, True) / dd.dispersal_time_yr(0.5, 0.3) < 1.0, "the same draw, a shorter time: not a second draw"


def test_the_time_is_the_stars_own_and_does_not_move_with_anything_but_its_mass_and_companion():
    first = dd.draw_dispersal_time_yr("system-7:star-0", 0.8)
    assert dd.draw_dispersal_time_yr("system-7:star-0", 0.8) == first
    assert dd.draw_dispersal_time_yr("system-7:star-1", 0.8) != first
    assert dd.draw_dispersal_time_yr("system-7:star-0", 0.8, close_companion=True) <= first
    times = [dd.draw_dispersal_time_yr(f"star-{i}", 0.5) for i in range(4000)]
    assert sum(t > 3.0 * MYR for t in times) / len(times) == pytest.approx(math.exp(-3.0 / 2.5), abs=0.03)
    assert all(t >= 0.0 for t in times)


def test_a_draw_outside_the_open_interval_or_a_mass_of_nothing_is_an_error():
    with pytest.raises(ValueError):
        dd.dispersal_time_yr(0.5, 0.0)
    with pytest.raises(ValueError):
        dd.dispersal_time_yr(0.5, 1.5)
    with pytest.raises(ValueError):
        dd.characteristic_time_yr(0.0)
