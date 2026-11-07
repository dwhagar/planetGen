# tests/test_surface_air_steam_limit.py

"""
GEN.1329: a body under its runaway greenhouse limit cannot hold its water as steam.

A hot steam air radiates about 280 W/m2 whatever it absorbs [R3236], so a body that gets less than the runaway flux of its star and mass
(`orbits.runaway_greenhouse_flux`) cools and its water returns to the surface. The fixed point had no such limit: a dense, hydrogen-rich first step (28 bar, which the
table warms to 1,000 K) carried the air past the critical temperature of water (647 K), `freeze_out` then let the whole ocean join the air (a gas above its critical
point has no saturation) and the steam air of 1e10 to 1e11 Pa stood at 2 W/m2 absorbed: 44 of the 63 steam airs in 450 generated systems were under the limit, one at 0.008
of it, and 8 of 14 Earth-mass bodies in the habitable zone of a Sun-like star ended at 1,100 to 1,500 K. `Setting.steam_possible` (`surface_air._steam_possible`) now
holds the surface of such a body with water under `surface_air.MOIST_LIMIT_K`.
"""

import math
import types

import pytest

from planetgen import tuning
from planetgen.physics import orbits, surface_air as sa
from tests.test_surface_air import _earth_like_planet


def _body(mass_earth=1.0, temperature_k=3400.0, own_w_m2=0.0, tidal_w_m2=0.0):
    return types.SimpleNamespace(star=types.SimpleNamespace(temperature=temperature_k), mass=mass_earth * 5.972e24,
                                 internal_heat_flux_w_m2=own_w_m2, tidal_flux_w_m2=tidal_w_m2)


def test_steam_is_possible_only_past_the_runaway_limit_of_the_star_and_the_mass():
    limit = orbits.runaway_greenhouse_flux(3400.0, 6.0)
    s = tuning.RUNAWAY_FLUX_REFERENCE_W_M2
    body = _body(mass_earth=6.0)
    assert sa._steam_possible(body, 0.3, 1.2 * limit * s)
    assert not sa._steam_possible(body, 0.3, 0.8 * limit * s)
    assert not sa._steam_possible(body, 0.3, 0.008 * limit * s)


def test_the_heat_of_the_body_counts_toward_the_limit_as_the_starlight_that_would_absorb_the_same():
    limit = orbits.runaway_greenhouse_flux(3400.0, 6.0)
    s = tuning.RUNAWAY_FLUX_REFERENCE_W_M2
    at_limit_w_m2 = limit * s * (1.0 - 0.3) / 4.0                  # what the body absorbs at the limit
    assert sa._steam_possible(_body(mass_earth=6.0, tidal_w_m2=1.1 * at_limit_w_m2), 0.3, 0.0)
    assert not sa._steam_possible(_body(mass_earth=6.0, tidal_w_m2=0.1), 0.3, 0.0)


def test_a_body_with_no_star_takes_the_limit_of_the_sun():
    body = types.SimpleNamespace(star=types.SimpleNamespace(), mass=5.972e24, internal_heat_flux_w_m2=0.0, tidal_flux_w_m2=0.0)
    limit = orbits.runaway_greenhouse_flux(sa.SUN_TEMPERATURE_K, 1.0)
    s = tuning.RUNAWAY_FLUX_REFERENCE_W_M2
    assert sa._steam_possible(body, 0.3, 1.1 * limit * s) and not sa._steam_possible(body, 0.3, 0.9 * limit * s)


def test_the_moist_limit_is_under_the_critical_temperature_of_water():
    assert 340.0 < sa.MOIST_LIMIT_K < 647.0


@pytest.mark.parametrize("seed", range(1, 10))
def test_an_earth_mass_body_in_the_habitable_zone_is_never_a_steam_air(seed):
    """Seeds 1, 3, 4, 5, 6, 7, 8, 9 of this body ended at 1,100 to 1,500 K with 1e9 to 2e10 Pa of steam before GEN.1329."""
    planet = _earth_like_planet(seed)
    assert planet.surface_temperature <= sa.MOIST_LIMIT_K * (1.0 + 1e-6)
    assert planet.atmospheric_pressure < 1.0e7
