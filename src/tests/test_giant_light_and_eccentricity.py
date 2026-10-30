# tests/test_giant_light_and_eccentricity.py

"""
GEN.977: the table radius of a giant is read at the light of its semi-major axis (`giant_draw`), not at the mean light of its orbit
(`orbit_elements.mean_flux_factor`), because the eccentricity is drawn after the radius. That is within 1% of the mean-light radius
unless the giant gets more than 300 times Earth's light and an eccentricity over 0.3 together (GEN.727, GEN.842). This guards the
corner: a seeded sweep fails if the generator ever makes a table-sized giant in it, and the corner itself is pinned so the bound in
the design note stays the one where the error passes 1%.
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, giant_draw, giant_radius, orbit_elements
from planetgen.util import draw

LIGHT_BOUND_EARTH = 300.0
ECCENTRICITY_BOUND = 0.3
SWEEP_SEED = 977
SWEEP_SYSTEMS = 80
CORNER_MASSES_EARTH = (50.0, 300.0, 1000.0)
CORE_MASS_EARTH = 10.0
AGE_GY = 4.5


def _table_sized_giants(systems, seed):
    """`(planet, light)` for every planet of `systems` seeded systems that is a giant the radius tables cover, `light` being its
    starlight at the semi-major axis in Earth units."""
    draw.set_run_seed(seed)
    found = []
    for _ in range(systems):
        for planet in StarSystem(system_config=SystemConfig()).planets:
            if getattr(planet, "body_type", None) != "g" or planet.is_moon:
                continue
            if not giant_radius.in_range(planet.mass / constants.EARTH_MASS_TO_KG):
                continue
            light = giant_draw.starlight_w_m2(planet) / giant_radius.flux_w_m2(1.0)
            found.append((planet, light))
    return found


def _mean_light_error(mass_earth, light_earth, eccentricity):
    """The share the table radius moves by when read at the mean light of the orbit instead of the light at its semi-major axis."""
    at_axis = giant_radius.radius_km(mass_earth, CORE_MASS_EARTH, light_earth, AGE_GY)
    at_mean = giant_radius.radius_km(mass_earth, CORE_MASS_EARTH, light_earth * orbit_elements.mean_flux_factor(eccentricity), AGE_GY)
    return abs(at_mean - at_axis) / at_axis


@pytest.mark.slow
def test_no_generated_table_giant_gets_more_than_300_times_earths_light_with_an_eccentricity_over_0_3():
    giants = _table_sized_giants(SWEEP_SYSTEMS, SWEEP_SEED)
    assert giants, "the sweep made no giant the tables cover"
    inside = [(planet.name, light, planet.eccentricity) for planet, light in giants
              if light > LIGHT_BOUND_EARTH and (planet.eccentricity or 0.0) > ECCENTRICITY_BOUND]
    assert not inside, f"table giants in the corner where the radius is read at the wrong light: {inside}"


def test_below_the_bound_the_mean_light_moves_a_table_radius_by_under_a_percent():
    for mass in CORNER_MASSES_EARTH:
        for light in (1.0, 30.0, LIGHT_BOUND_EARTH):
            assert _mean_light_error(mass, light, ECCENTRICITY_BOUND * 0.5) < 0.01, (mass, light)
        assert _mean_light_error(mass, 1.0, 0.9) < 0.01, mass


def test_at_the_bound_the_most_sensitive_giant_passes_a_percent():
    worst = max(_mean_light_error(mass, light, ECCENTRICITY_BOUND) for mass in CORNER_MASSES_EARTH for light in (300.0, 1000.0))
    assert worst > 0.01
