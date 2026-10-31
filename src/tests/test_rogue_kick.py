# tests/test_rogue_kick.py

"""
GEN.1112 RP6: the ejection speed of a rogue planet and its velocity through the Galaxy (`planetgen.physics.rogue_kick`): escape and orbital speeds, the Safronov number, the speed at
infinity from the scatterer, the isotropic direction, the velocity as the host's plus the drift's, and the speed window of an interstellar object.
"""

import math
import random

import pytest

from planetgen.physics import constants, rogue_drift
from planetgen.physics import rogue_kick as rk

SUN = constants.SOLAR_MASS_TO_KG
JUPITER = (1.8982e27, 7.1492e7)
NEPTUNE = (1.0243e26, 2.4622e7)
EARTH = (5.972e24, 6.371e6)


def test_the_suns_escape_and_circular_speeds_at_one_au():
    assert rk.escape_speed_km_s(SUN, 1.0) == pytest.approx(42.1, abs=0.1)
    assert rk.circular_speed_km_s(SUN, 1.0) == pytest.approx(29.78, abs=0.05)
    assert rk.planet_escape_speed_km_s(*JUPITER) == pytest.approx(59.5, abs=0.5)
    assert rk.planet_escape_speed_km_s(*EARTH) == pytest.approx(11.19, abs=0.05)


def test_the_safronov_number_splits_the_scatterers_that_eject_from_those_that_stir():
    assert rk.safronov_number(*JUPITER, SUN, 5.2) == pytest.approx(10.4, rel=0.03)
    assert rk.safronov_number(*EARTH, SUN, 1.0) == pytest.approx(0.071, rel=0.05)
    assert rk.ejects(*JUPITER, SUN, 5.2) and rk.ejects(*NEPTUNE, SUN, 5.2)
    assert not rk.ejects(*EARTH, SUN, 1.0)
    # Theta = (M_p / M_*) (a / R_p)
    assert rk.safronov_number(*JUPITER, SUN, 5.2) == pytest.approx(JUPITER[0] / SUN * 5.2 * constants.AU_M / JUPITER[1], rel=1e-3)


def _mean_speed(planet_mass, distance_au, ejected_mass=0.0, n=20000, seed=3):
    rng = random.Random(seed)
    return sum(rk.ejection_speed_km_s(rng.random(), planet_mass, distance_au, ejected_mass) for _ in range(n)) / n


def test_the_scale_is_huang_and_lais_three_km_a_second_for_a_jupiter_mass_at_one_au():
    assert rk.scatter_speed_km_s(constants.JUPITER_MASS_TO_KG, 1.0) == pytest.approx(3.0)
    assert rk.scatter_speed_km_s(constants.JUPITER_MASS_TO_KG, 5.2) == pytest.approx(1.316, abs=0.01)
    assert rk.scatter_speed_km_s(10.0 * constants.JUPITER_MASS_TO_KG, 5.2) == pytest.approx(4.16, abs=0.02)
    # the ejected mass lowers it by m_p / (m_p + m_ej), and the star does not enter at all
    assert rk.scatter_speed_km_s(1e27, 2.0, 1e27) == pytest.approx(0.5 * rk.scatter_speed_km_s(1e27, 2.0))


def test_the_mean_speeds_are_kramer_and_profumos_for_one_and_ten_jupiter_masses_at_5_2_au():
    mj = constants.JUPITER_MASS_TO_KG
    assert _mean_speed(mj, 5.2) == pytest.approx(1.15, rel=0.15)                   # Sec. 4.2: 1.15 km/s
    assert _mean_speed(10.0 * mj, 5.2) == pytest.approx(3.5, rel=0.12)             # Sec. 4.1: 3.4 to 3.6 km/s
    exponent = math.log(_mean_speed(10.0 * mj, 5.2) / _mean_speed(mj, 5.2)) / math.log(10.0)
    assert 0.4 < exponent < 0.55                                                   # their fit over 1 to 75 Jupiter masses: 0.42


def test_the_ejected_mass_makes_no_difference_over_four_decades():
    mj = constants.JUPITER_MASS_TO_KG
    light, heavy = (_mean_speed(10.0 * mj, 5.2, ejected) for ejected in (1e-5 * mj, 1e-1 * mj))
    assert heavy == pytest.approx(light, rel=0.1)                                  # Kramer and Profumo Fig. 1: 3.4 to 3.6 km/s throughout


def test_the_spread_tops_out_near_two_and_a_half_times_the_scale():
    scale = rk.scatter_speed_km_s(constants.JUPITER_MASS_TO_KG, 1.0)
    assert rk.ejection_speed_km_s(0.99, constants.JUPITER_MASS_TO_KG, 1.0) == pytest.approx(2.5 * scale, rel=0.05)
    assert rk.ejection_speed_km_s(0.5, constants.JUPITER_MASS_TO_KG, 1.0) < scale


def test_a_giants_ejections_are_about_a_kilometre_a_second_and_slower_further_out_inside_sayeds_range():
    mj = constants.JUPITER_MASS_TO_KG
    medians = {a: rk.ejection_speed_km_s(0.5, mj, a) for a in (1.0, 5.2, 30.0)}
    assert medians[1.0] > medians[5.2] > medians[30.0]
    rng = random.Random(6)
    speeds = sorted(rk.ejection_speed_km_s(rng.random(), mj * 10 ** rng.uniform(-0.5, 1.0), 10 ** rng.uniform(0.0, 1.5)) for _ in range(5000))
    assert 0.7 < speeds[len(speeds) // 2] < 10.7                                   # Sayed, Profumo and Smyth Table 10: medians 0.708 to 10.716 km/s
    assert speeds[int(0.9 * len(speeds))] < 10.7


def test_the_speed_rises_with_the_random_number_and_is_small_at_the_threshold():
    speeds = [rk.ejection_speed_km_s(i / 100.0, JUPITER[0], 5.2) for i in range(100)]
    assert speeds == sorted(speeds)
    assert rk.ejection_speed_km_s(0.0, JUPITER[0], 5.2) < 0.05 * rk.scatter_speed_km_s(JUPITER[0], 5.2)


def test_a_moon_survives_when_it_is_inside_half_the_hill_sphere_the_closest_approach_leaves_it():
    m_ej, m_p = 1.0 * constants.JUPITER_MASS_TO_KG, 10.0 * constants.JUPITER_MASS_TO_KG
    reach = rk.MOON_LOSS_FRACTION * 1.0e10 * (m_ej / (m_ej + m_p)) ** (1.0 / 3.0)
    assert rk.moon_survives(0.99 * reach, 1.0e10, m_p, m_ej) and not rk.moon_survives(1.01 * reach, 1.0e10, m_p, m_ej)


def test_the_closest_approach_has_huang_and_lais_median_and_spread():
    mj = constants.JUPITER_MASS_TO_KG
    hill = rk.mutual_hill_radius_m(SUN, 5.2, 10.0 * mj, mj)
    assert hill == pytest.approx(5.2 * constants.AU_M * (11.0 * mj / (3.0 * SUN)) ** (1.0 / 3.0))
    rng = random.Random(1)
    draws = sorted(rk.closest_approach_m(rng.random(), SUN, 5.2, 10.0 * mj, mj) / hill for _ in range(20000))
    assert draws[len(draws) // 2] == pytest.approx(0.089, rel=0.07)
    assert sum(1 for d in draws if d > 0.1) / len(draws) == pytest.approx(0.48, abs=0.03)       # Sec. 7: 48 percent above 0.1
    assert sum(1 for d in draws if d > 0.01) / len(draws) == pytest.approx(0.82, abs=0.03)      # and 82 percent above 0.01
    light = rk.closest_approach_m(0.5, SUN, 5.2, 1.0 * mj, 0.1 * mj) / rk.mutual_hill_radius_m(SUN, 5.2, 1.0 * mj, 0.1 * mj)
    assert light == pytest.approx(0.0066, rel=0.05)                                            # Fig. 2: 1 Jupiter mass scatterer


def test_the_moon_survival_probability_falls_with_the_moons_distance_and_matches_the_draw():
    mj = constants.JUPITER_MASS_TO_KG
    args = (SUN, 5.2, 10.0 * mj, mj)
    hill = rk.mutual_hill_radius_m(*args[:2], *args[2:])
    near, far = (rk.moon_survival_probability(f * hill, *args) for f in (0.001, 0.1))
    assert 0.85 < near < 1.0 and 0.0 < far < near
    assert 0.5 < rk.moon_survival_probability(0.01 * hill, *args) < 0.85                       # a Galilean-moon radius: the paper's 70 to 82 percent
    rng = random.Random(5)
    r_moon = 0.02 * hill
    kept = sum(rk.moon_survives(r_moon, rk.closest_approach_m(rng.random(), *args), args[2], args[3]) for _ in range(20000)) / 20000
    assert kept == pytest.approx(rk.moon_survival_probability(r_moon, *args), abs=0.02)


def test_the_isotropic_direction_is_a_unit_vector_whose_components_average_to_nothing():
    rng = random.Random(4)
    vectors = [rk.isotropic_direction(rng.random(), rng.random()) for _ in range(6000)]
    assert all(math.sqrt(sum(c * c for c in v)) == pytest.approx(1.0) for v in vectors)
    for axis in range(3):
        assert abs(sum(v[axis] for v in vectors) / len(vectors)) < 0.03
        assert sum(v[axis] ** 2 for v in vectors) / len(vectors) == pytest.approx(1.0 / 3.0, abs=0.02)
    assert rk.isotropic_direction(1.0, 0.0)[2] == pytest.approx(1.0)


def test_the_rogues_velocity_is_the_hosts_plus_the_drifts_relative_velocity():
    host = (10.0, -5.0, 3.0)
    kick = rk.kick(3.0, (0.0, 1.0, 0.0))
    assert rk.rogue_velocity(host, kick, 0.0) == pytest.approx((10.0, -2.0, 3.0))
    later = rk.rogue_velocity(host, kick, 60.0)
    drift = rogue_drift.drift(kick, 60.0)
    assert later == pytest.approx((host[0] + drift.vx_km_s, host[1] + drift.vy_km_s, host[2] + drift.vz_km_s))
    assert rk.relative_speed_km_s(later, host) == pytest.approx(math.sqrt(drift.vx_km_s ** 2 + drift.vy_km_s ** 2 + drift.vz_km_s ** 2))


def test_the_speed_window_of_a_rogue_passing_the_sun_holds_the_known_interstellar_objects():
    # hosts with the thin disc's velocity dispersion (about 35, 25, 20 km/s), a giant's kick, a random age: the rogue's speed relative to the Sun
    rng = random.Random(8)
    speeds = []
    for _ in range(4000):
        host = (rng.gauss(0, 35), rng.gauss(0, 25), rng.gauss(0, 20))
        kick = rk.kick(rk.ejection_speed_km_s(rng.random(), JUPITER[0], rng.choice((3.0, 5.2, 10.0, 30.0))),
                       rk.isotropic_direction(rng.random(), rng.random()))
        speeds.append(rk.relative_speed_km_s(rk.rogue_velocity(host, kick, rng.uniform(10.0, 8000.0))))
    speeds.sort()
    p10, p50, p90 = (speeds[int(f * len(speeds))] for f in (0.1, 0.5, 0.9))
    assert 20.0 < p50 < 60.0
    assert p10 < 26.0 < p90 and p10 < 32.0 < p90                          # 'Oumuamua 26 and Borisov 32 km/s relative to the Sun, both inside the window


def test_inputs_the_model_cannot_use_are_refused():
    with pytest.raises(rk.RogueKickError):
        rk.escape_speed_km_s(0.0, 1.0)
    with pytest.raises(rk.RogueKickError):
        rk.planet_escape_speed_km_s(1e27, 0.0)
    with pytest.raises(rk.RogueKickError):
        rk.ejection_speed_km_s(1.0, JUPITER[0], 5.2)
    with pytest.raises(rk.RogueKickError):
        rk.isotropic_direction(1.5, 0.0)
    with pytest.raises(rk.RogueKickError):
        rk.rogue_velocity((1.0, 2.0), (0.0, 0.0, 0.0), 1.0)
