# tests/test_rprocess_enrichment.py

"""The thorium-to-uranium enrichment factor from the disc's history (GEN.663)."""

import pytest

from planetgen.galaxy import disk_history
from planetgen.physics import rprocess_enrichment as rp


def _history(stellar_mass, halo_mass):
    return disk_history.history(0.8 * stellar_mass, 2600.0, halo_mass)


def test_the_decay_constants_and_the_constant_input_factor_are_the_rows_former_numbers():
    l232, l238 = rp.decay_constants_per_gy()
    assert l238 - l232 == pytest.approx(0.105802, rel=1e-3)
    assert rp.constant_input_factor(3.0) == pytest.approx(1.16254, rel=1e-3)
    assert rp.constant_input_factor(13.0) == pytest.approx(1.71724, rel=1e-3)
    assert rp.effective_duration_gy(rp.constant_input_factor(6.5)) == pytest.approx(6.5, abs=1e-6)


def test_a_constant_star_formation_and_prompt_sources_reproduce_the_constant_input_factor():
    times = [i * 0.05 for i in range(201)]
    factor = rp.th_u_factor(times, [1.0] * len(times), 10.0, tau_min_gy=0.0501)
    # a 1/tau law over a flat history is flat apart from its start-up: within 7 percent (the start-up) of the 10 Gyr constant input
    assert factor == pytest.approx(rp.constant_input_factor(10.0), rel=0.07)
    assert rp.th_u_factor(times, [1.0] * len(times), 0.0) == 1.0


def test_an_older_star_has_less_earlier_gas_and_a_younger_star_more_thorium():
    history = _history(6.1e10, 3.1e12)
    factors = [rp.th_u_factor_of_history(history, age) for age in (9.0, 7.0, 4.5, 2.0, 0.0)]
    assert factors == sorted(factors), "younger stars formed from gas whose earlier ejecta had longer to decay"
    assert 1.0 <= factors[0] and factors[-1] < 2.9, "below the ratio of a single burst at the start of the disc, exp(0.1058 x 10) = 2.88"


def test_a_galaxy_that_formed_its_stars_early_gives_the_gas_of_a_late_star_more_thorium():
    early = rp.th_u_factor_of_history(_history(6.1e10, 3.1e12), 1.0)       # heavy: early history
    late = rp.th_u_factor_of_history(_history(1.0e10, 3.2e11), 1.0)        # light: rising history
    assert early > late


def test_the_solar_system_datum_is_a_test_point_the_disc_alone_stays_under():
    # CI chondrites need an effective constant input of 6.5 Gyr; the thin disc alone (no halo or thick disc gas) gives less (a lower bound)
    found = rp.effective_duration_gy(rp.th_u_factor_of_history(_history(6.1e10, 3.1e12), 4.57))
    assert 3.5 < found < 6.5


def test_the_merger_delay_floor_moves_the_factor_little_and_prompt_sources_least():
    history = _history(6.1e10, 3.1e12)
    prompt = rp.th_u_factor_of_history(history, 4.5, tau_min_gy=0.005)
    slow = rp.th_u_factor_of_history(history, 4.5, tau_min_gy=0.3)
    assert abs(prompt - slow) / prompt < 0.05
