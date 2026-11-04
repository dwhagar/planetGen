# tests/test_field_rates.py

"""The physical birth rates of interstellar asteroid fields and the host rule of a supernova remnant (GEN.1101, `physics/field_rates.py`): the design's standing numbers, the half-mass rule, and
the rule that nothing is drawn at a made-up rate."""

import pytest

from planetgen.physics import constants, field_dispersal, field_rates as fr

EARTH = constants.EARTH_MASS_TO_KG
MSUN = constants.SOLAR_MASS_TO_KG
AU = constants.AU_M
R_EARTH = 6.371e6


def test_the_unbound_belts_of_the_galaxy_stand_at_about_two_hundred():
    speed = field_dispersal.unbound_belt_speed_ms(10 * MSUN, 1.4 * MSUN, 50 * AU)
    rate = fr.unbound_belt_birth_rate_per_year(1.0e7)
    assert rate == pytest.approx(0.02 * 0.41)
    assert 90 < fr.standing_number(rate, field_dispersal.retirement_age_yr(speed)) < 250          # the design: 90 to 220
    assert fr.standing_per_sector(200.0, 9.0, 2.7e11) == pytest.approx(6.7e-9, rel=0.05)           # 3e-9 to 8e-9 per sector


def test_the_belt_share_falls_with_the_age_of_the_progenitor():
    shares = [fr.belt_share(age) for age in (1e7, 3e7, 1e8, 3e8, 1e9)]
    assert shares == [0.41, 0.41, 0.25, 0.15, 0.02]
    with pytest.raises(fr.FieldRatesError):
        fr.belt_share(-1.0)


def test_a_belt_is_unbound_only_when_over_half_the_mass_leaves():
    assert fr.belt_is_unbound(10 * MSUN, 1.4 * MSUN)
    assert not fr.belt_is_unbound(25 * MSUN, 13 * MSUN)                                         # a heavy black hole keeps it
    assert fr.unbound_belt_birth_rate_per_year(1e7, unbound_share=0.0) == 0.0
    with pytest.raises(fr.FieldRatesError):
        fr.unbound_belt_birth_rate_per_year(1e7, unbound_share=1.5)


def test_a_young_remnant_hosts_a_belt_and_an_old_one_does_not():
    args = (0.2, 1.0e7, 10 * MSUN, 1.4 * MSUN, 50 * AU)
    belt = fr.hosted_belt(*args, 5000.0)
    assert belt.age_yr == 5000.0 and belt.belt_radius_m == 50 * AU
    assert belt.speed_ms / 1000 == pytest.approx(11.3, abs=0.1)
    assert fr.hosted_belt(*args, 20000.0) is not None
    assert fr.hosted_belt(*args, 40000.0) is None                                              # past a light year: retired, not hosted
    assert fr.hosted_belt(0.5, 1.0e7, 10 * MSUN, 1.4 * MSUN, 50 * AU, 5000.0) is None          # 0.5 is over the belt share of 0.41
    assert fr.hosted_belt(0.2, 1.0e7, 25 * MSUN, 13 * MSUN, 50 * AU, 5000.0) is None          # the remnant keeps the belt
    with pytest.raises(fr.FieldRatesError):
        fr.hosted_belt(1.0, 1.0e7, 10 * MSUN, 1.4 * MSUN, 50 * AU, 5000.0)


def test_the_share_of_remnants_that_host_a_belt_is_about_a_quarter():
    assert fr.hosted_age_share((50.0, 1.0e5), 11.3e3) == pytest.approx(0.265, abs=0.01)
    assert fr.hosted_age_share((50.0, 1.0e5), 0.5e3) == 1.0                                     # slow enough to keep for the whole life of a remnant
    assert fr.hosted_age_share((50.0, 1.0e5), 1.0e6) < 0.01
    with pytest.raises(fr.FieldRatesError):
        fr.hosted_age_share((10.0, 5.0), 1e4)


def test_a_collision_field_in_a_sector_is_astronomically_rare():
    rate = fr.collision_field_rate_per_sector_per_year(0.78, 64.0, EARTH, EARTH, R_EARTH, R_EARTH, 40e3, 25e3)
    assert 1e-23 < rate < 1e-21
    assert fr.standing_number(rate, 6e4) < 1e-16                                               # a field for every 1e16 sectors at the least
    faster = fr.collision_field_rate_per_sector_per_year(0.78, 64.0, EARTH, EARTH, R_EARTH, R_EARTH, 40e3, 100e3)
    assert faster < rate                                                                      # a hot population is past the rubble window
    assert fr.collision_field_rate_per_sector_per_year(0.78, 128.0, EARTH, EARTH, R_EARTH, R_EARTH, 40e3, 25e3) == pytest.approx(2 * rate)


def test_a_young_cluster_strips_some_belts_and_the_swarm_lasts_tens_of_kyr():
    rate = fr.flyby_stripping_rate_per_star_per_year(1000.0, 2e3, 50 * AU, MSUN)
    assert rate * 1e6 == pytest.approx(0.02, rel=0.3)                                          # about 20 per Myr per 1000 stars, 0.02 per Myr per star
    assert fr.flyby_stripping_rate_per_star_per_year(2000.0, 2e3, 50 * AU, MSUN) == pytest.approx(2 * rate)
    assert rate < fr.flyby_stripping_rate_per_star_per_year(1000.0, 8e3, 50 * AU, MSUN) < 2 * rate    # faster flies more, but focusing carries most of the slow rate
    assert fr.flyby_stripping_rate_per_star_per_year(1000.0, 2e3, 100 * AU, MSUN) > 2 * rate       # a wider belt is stripped by a wider flyby (focusing: the area grows with q, not q^2)
    assert fr.cluster_swarm_lifetime_yr(50 * AU, MSUN) == pytest.approx(7.1e4, rel=0.1)
    assert fr.cluster_swarm_standing(1000, 1000.0, 2e3, 50 * AU, MSUN) == pytest.approx(0.5, rel=0.3)
    assert fr.cluster_swarm_standing(0, 1000.0, 2e3, 50 * AU, MSUN) == 0.0


def test_bad_arguments_are_refused():
    with pytest.raises(fr.FieldRatesError):
        fr.standing_number(-1.0, 1.0)
    with pytest.raises(fr.FieldRatesError):
        fr.standing_per_sector(1.0, 1.0, 0.0)
    with pytest.raises(fr.FieldRatesError):
        fr.flyby_stripping_rate_per_star_per_year(0.0, 2e3, AU, MSUN)
    with pytest.raises(fr.FieldRatesError):
        fr.collision_field_rate_per_sector_per_year(0.78, 0.0, EARTH, EARTH, R_EARTH, R_EARTH, 40e3, 25e3)
