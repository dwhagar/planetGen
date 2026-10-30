# tests/test_rogue_moons.py

"""GEN.1112 RP4: whether a moon leaves with the planet that is thrown out and what it does there (`planetgen.physics.rogue_moons`)."""

import math

import pytest

from planetgen.physics import constants, rogue_moons as rm, tidal_heating

EARTH = constants.EARTH_MASS_TO_KG
JUPITER = constants.JUPITER_MASS_TO_KG
SUN = constants.SOLAR_MASS_TO_KG
MOON = 7.346e22
R_EARTH = 6.371e6


def test_the_survival_curve_passes_through_the_two_points_the_paper_prints():
    assert rm.closest_approach_survival(0.01) == pytest.approx(0.82)
    assert rm.closest_approach_survival(0.1) == pytest.approx(0.48)
    assert rm.closest_approach_survival(0.089) == pytest.approx(0.5, abs=0.02)      # their median
    assert rm.closest_approach_survival(1e-6) == 1.0 and rm.closest_approach_survival(100.0) == 0.0
    shares = [rm.closest_approach_survival(10.0 ** k) for k in range(-4, 3)]
    assert shares == sorted(shares, reverse=True)


def test_the_draw_inverts_the_survival_curve():
    for u in (0.1, 0.4, 0.5, 0.7):
        assert rm.closest_approach_survival(rm.closest_approach_hill(u)) == pytest.approx(1.0 - u)
    assert rm.closest_approach_hill(0.5) == pytest.approx(0.089, rel=0.1)
    assert rm.closest_approach_hill(0.0) < rm.closest_approach_hill(1.0)


def test_a_moon_follows_its_planet_when_its_orbit_is_inside_the_encounter_scale():
    r_min = 1.0e9
    reach = rm.survival_distance_m(r_min, EARTH, JUPITER)
    assert reach == pytest.approx(r_min * (EARTH / (EARTH + JUPITER)) ** (1.0 / 3.0))
    assert rm.moon_survives(0.99 * reach, r_min, EARTH, JUPITER)
    assert not rm.moon_survives(1.01 * reach, r_min, EARTH, JUPITER)
    # the orbital time against the encounter time of eq. 4: equal at the limit
    t_moon = math.sqrt(reach ** 3 / (constants.G * EARTH))
    t_encounter = math.sqrt(r_min ** 3 / (constants.G * (EARTH + JUPITER)))
    assert t_moon == pytest.approx(t_encounter)


def test_a_heavier_scatterer_makes_it_harder_for_the_moon():
    assert rm.survival_distance_m(1e9, EARTH, 10 * JUPITER) < rm.survival_distance_m(1e9, EARTH, JUPITER)


def test_the_earth_moon_orbit_survives_most_encounters_and_a_wide_orbit_few():
    hill = rm.mutual_hill_radius_m(5.2 * constants.AU_M, 4.0 * constants.AU_M, JUPITER, EARTH, SUN)
    assert hill == pytest.approx(0.5 * 9.2 * constants.AU_M * ((JUPITER + EARTH) / (3 * SUN)) ** (1.0 / 3.0))
    near = rm.survival_probability(3.84e8, EARTH, JUPITER, hill)
    wide = rm.survival_probability(0.5 * hill, EARTH, JUPITER, hill)
    assert 0.4 < near < rm.SURVIVAL_AT_0P01_HILL and wide < 0.2             # under the paper's 0.82 at 0.01 R_H: the light planet's reach is cut by (m / (m + M))^(1/3)


def test_the_circularisation_heat_of_the_debes_and_sigurdsson_median_pair():
    """Their eq. 1 for the median pair (e 0.13, a 30 Earth radii, a lunar moon) with their own rheology (k2 0.299, mu 4, Q 12)."""
    a = 30.0 * R_EARTH
    n = math.sqrt(constants.G * (EARTH + MOON) / a ** 3)
    debes = 63.0 * 0.13 ** 2 * n * R_EARTH ** 5 * constants.G * MOON ** 2 / (4.0 * 4.0 * 12.0 * a ** 6)
    flux = rm.circularisation_flux_w_m2(EARTH, R_EARTH, MOON, a, 0.13, k2_over_q=0.299 / 12.0)
    assert flux * 4.0 * math.pi * R_EARTH ** 2 == pytest.approx(debes, rel=0.25)
    assert 1e-3 < flux < 1e-2                                          # inside the old default's range of 1e-4 to 1e-2 W/m2


def test_the_default_rheology_is_the_rock_of_tidal_heating_and_the_flux_scales_with_the_square_of_the_eccentricity():
    a = 30.0 * R_EARTH
    flux = rm.circularisation_flux_w_m2(EARTH, R_EARTH, MOON, a)
    assert flux == pytest.approx(rm.circularisation_flux_w_m2(EARTH, R_EARTH, MOON, a, k2_over_q=tidal_heating.rigid_k2(EARTH, R_EARTH) / tidal_heating.ROCK_TIDAL_Q))
    assert rm.circularisation_flux_w_m2(EARTH, R_EARTH, MOON, a, 0.26) == pytest.approx(4.0 * flux)
    assert rm.circularisation_flux_w_m2(EARTH, R_EARTH, MOON, a, 0.0) == 0.0
    assert rm.circularisation_flux_w_m2(EARTH, R_EARTH, MOON, 2.0 * a) < flux


def test_bad_inputs_are_refused():
    with pytest.raises(rm.RogueMoonError):
        rm.closest_approach_survival(0.0)
    with pytest.raises(rm.RogueMoonError):
        rm.closest_approach_hill(1.5)
    with pytest.raises(rm.RogueMoonError):
        rm.moon_survives(0.0, 1e9, EARTH, JUPITER)
    with pytest.raises(rm.RogueMoonError):
        rm.circularisation_flux_w_m2(EARTH, R_EARTH, MOON, 1e8, 1.0)
