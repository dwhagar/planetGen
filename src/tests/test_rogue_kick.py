# tests/test_rogue_kick.py

"""
GEN.1112 RP6: the ejection speed of a rogue planet and its velocity through the Galaxy (`planetgen.physics.rogue_kick`): escape and orbital speeds, the Safronov number, the speed at
infinity from the energy step, the isotropic direction, the velocity as the host's plus the drift's, and the speed window of an interstellar object.
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


def test_the_energy_step_is_capped_at_the_scale_for_a_scatterer_that_ejects():
    assert rk.energy_step(*JUPITER, SUN, 5.2) == pytest.approx(rk.ENERGY_STEP_SCALE)
    assert rk.energy_step(*EARTH, SUN, 1.0) < rk.ENERGY_STEP_SCALE


def test_a_giants_ejections_are_a_few_to_ten_kilometres_a_second_and_slower_further_out():
    medians = {a: rk.ejection_speed_km_s(0.5, SUN, a, *JUPITER) for a in (1.0, 5.2, 30.0)}
    assert medians[1.0] > medians[5.2] > medians[30.0]
    assert 3.0 < medians[5.2] < 10.0 and 1.0 < medians[30.0] < 5.0                    # the note's 3 km/s and Sayed et al.'s 0.7 to 10.7
    tail = rk.ejection_speed_km_s(0.999, SUN, 5.2, *JUPITER)
    assert 10.0 < tail < 30.0                                                         # Kramer and Profumo: a ceiling near 23 km/s, rare kicks past it


def test_the_speed_rises_with_the_random_number_and_is_zero_at_the_threshold():
    assert rk.ejection_speed_km_s(0.0, SUN, 5.2, *JUPITER) == 0.0
    speeds = [rk.ejection_speed_km_s(i / 100.0, SUN, 5.2, *JUPITER) for i in range(100)]
    assert speeds == sorted(speeds)


def test_the_energy_beyond_ejection_is_exponentially_distributed():
    rng = random.Random(2)
    step = rk.energy_step(*JUPITER, SUN, 5.2)
    v_esc = rk.escape_speed_km_s(SUN, 5.2)
    excess = [0.5 * rk.ejection_speed_km_s(rng.random(), SUN, 5.2, *JUPITER) ** 2 for _ in range(20000)]
    mean_expected = step * v_esc ** 2 / 4.0
    assert sum(excess) / len(excess) == pytest.approx(mean_expected, rel=0.05)


def test_a_heavier_star_gives_a_faster_ejection_at_the_same_distance():
    assert rk.ejection_speed_km_s(0.5, 2.0 * SUN, 5.2, *JUPITER) > rk.ejection_speed_km_s(0.5, SUN, 5.2, *JUPITER)


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
        kick = rk.kick(rk.ejection_speed_km_s(rng.random(), SUN, rng.choice((3.0, 5.2, 10.0, 30.0)), *JUPITER),
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
        rk.ejection_speed_km_s(1.0, SUN, 5.2, *JUPITER)
    with pytest.raises(rk.RogueKickError):
        rk.isotropic_direction(1.5, 0.0)
    with pytest.raises(rk.RogueKickError):
        rk.rogue_velocity((1.0, 2.0), (0.0, 0.0, 0.0), 1.0)
