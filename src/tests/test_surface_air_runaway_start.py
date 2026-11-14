# tests/test_surface_air_runaway_start.py

"""
GEN.1238: a rocky body with surface water under starlight past its runaway greenhouse limit holds its water as steam, not as a cool ocean.

The air's fixed point has two answers for such a body: a cool one (a few bar of vapor at 350 to 420 K, which the greenhouse table, fitted to a radiative-convective
model, gives no outgoing-flux limit to) and the steam one (the water in the air). Which it finds depended on where it started, and it started from the airless
temperature, so a Venus twin (Venus's mass at Venus's 1.9 times the Sun's flux at the Earth) kept a liquid ocean in more than half of its draws although the water-loss
history (`runaway_history`) already called it past the limit (about 1.1 S for a body of its mass around a Sun-like star). `air_final.past_runaway_with_water` now
starts such a body hot (`tuning.RUNAWAY_START_TEMPERATURE_K`), so the fixed point finds the steam air.
"""

import math
import types

import pytest

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import air_final, orbits
from planetgen.util import draw
from tests import class_support

VENUS_MASS_KG = 4.867e24
VENUS_EFFECTIVE_FLUX = 1.91
"""float: Venus's flux, in the Sun's at the Earth."""
SOLAR_L_W = 3.828e26


def _planet(temperature_k=5772.0, mass_kg=VENUS_MASS_KG, body_type="t"):
    return types.SimpleNamespace(star=types.SimpleNamespace(temperature=temperature_k), mass=mass_kg, body_type=body_type)


def test_a_body_with_water_past_its_runaway_limit_starts_hot_and_one_below_it_or_without_water_does_not():
    limit = orbits.runaway_greenhouse_flux(5772.0, VENUS_MASS_KG / 5.972e24)
    s = tuning.RUNAWAY_FLUX_REFERENCE_W_M2
    assert air_final.past_runaway_with_water(_planet(), 1.2 * limit * s, 1e-3)
    assert not air_final.past_runaway_with_water(_planet(), 0.8 * limit * s, 1e-3)
    assert not air_final.past_runaway_with_water(_planet(), 1.2 * limit * s, 0.0)
    assert not air_final.past_runaway_with_water(_planet(), 1.2 * limit * s, None)
    assert not air_final.past_runaway_with_water(_planet(body_type="g"), 1.2 * limit * s, 1e-3)


def test_the_limit_is_the_one_the_water_loss_history_follows():
    """The Sun's flux at the Earth is the unit of `orbits.runaway_greenhouse_flux`, so the reference is the solar constant."""
    assert tuning.RUNAWAY_FLUX_REFERENCE_W_M2 == pytest.approx(3.828e26 / (4.0 * math.pi * 1.495978707e11 ** 2), rel=2e-3)


@pytest.fixture(scope="module")
def venus_twins():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    with draw.bound("gen1238-star"):
        star = Star(config)
    distance_au = math.sqrt(star.luminosity / SOLAR_L_W / VENUS_EFFECTIVE_FLUX)
    return [class_support.draw_planet(star, distance_au, seed, config, mass=VENUS_MASS_KG, moon_count=0) for seed in range(12)]


def test_a_wet_body_past_the_runaway_limit_holds_its_water_as_steam(venus_twins):
    # A body whose recipe draws a retention regime with a pressure cap (`surface_air._retention_cap_pa`, rows d1 to d4 of the recipe, GEN.1340) keeps an air under that cap whatever its water,
    # so it cannot hold its water as steam; which regime a seed draws moves with every change of the draws before it (the nitrogen of GEN.1584 made seed 8 draw d4, 0.1 bar), so the
    # bodies the test speaks of are the ones no cap holds, and the capped ones are held to their cap.
    wet = [p for p in venus_twins if p.water_mass_fraction > 1.0e-3 and p._retention_cap_pa is None]
    assert wet, "the sample holds no wet body, so it says nothing"
    assert [p.hydrosphere for p in venus_twins if p.hydrosphere == "surface ocean"] == []
    for planet in (p for p in venus_twins if p._retention_cap_pa is not None):
        assert planet.atmospheric_pressure <= planet._retention_cap_pa * (1.0 + 1e-9)
    for planet in wet:
        assert planet.surface_temperature > 550.0, (planet.water_mass_fraction, planet.surface_temperature, planet.atmospheric_pressure)
        assert planet.atmospheric_pressure > 20.0e5
