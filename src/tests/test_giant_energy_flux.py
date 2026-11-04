"""
GEN.784: a giant's `energy_flux_w_m2` is the energy it radiates through its 1 bar level (the starlight it absorbs and its own heat,
`sigma T_eff^4`), and its `flux_variation` is the swing of its starlight, as every planet's.
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import climate_draws, constants, giant_draw, habitability_world
from planetgen.util import draw


def _giants():
    found = []
    for seed in range(12):
        with draw.bound(seed):
            system = StarSystem(SystemConfig())
        found += [planet for planet in system.planets if planet.body_type == "g"]
    if not found:
        pytest.skip("no giant in the sample systems")
    return found


def test_a_giants_energy_flux_is_what_it_absorbs_and_what_its_interior_gives():
    for planet in _giants():
        absorbed = (1.0 - giant_draw.BOND_ALBEDO) * giant_draw.starlight_w_m2(planet) / 4.0
        assert planet.energy_flux_w_m2 == pytest.approx(absorbed + giant_draw.own_flux_w_m2(planet), rel=1e-6)
        assert planet.energy_flux_w_m2 == pytest.approx(habitability_world.giant_energy_flux_w_m2(planet), rel=1e-12)


def test_a_giants_flux_variation_is_filled_and_inside_its_bounds():
    for planet in _giants():
        assert planet.flux_variation is not None and 0.0 <= planet.flux_variation <= 1.0


def test_a_giant_with_no_starlight_radiates_its_own_heat():
    planet = _giants()[0]
    stripped = type("Moonlike", (), {})()
    stripped.is_moon = True
    stripped.mass, stripped.radius = planet.mass, planet.radius
    stripped.star = planet.star
    flux = giant_draw.own_flux_w_m2(stripped)
    assert habitability_world.giant_energy_flux_w_m2(stripped) == pytest.approx(
        flux + constants.STEFAN_BOLTZMANN_CONSTANT * climate_draws.equilibrium_temperature_k(giant_draw.BOND_ALBEDO, 0.0) ** 4, rel=1e-9)
