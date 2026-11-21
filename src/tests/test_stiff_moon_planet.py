"""
A planet whose moon history is stiff (GEN.1815): seed 9101 of `tests/test_aquifers` has a rocky planet of 8.9e24 kg and a moon of 5.0e23 kg born at 3.8 planet radii, which the three BDF runs of
`moon_orbits.STRATEGIES` could not integrate. The moons then stayed as born over a spin that had not felt them, and their tide put 2.5e7 W/m2 (4e8 times the Earth's heat flow) into the planet. Radau finishes the history.
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import moon_orbits
from tests.fuzz_support import deterministic_entropy

IO_FLUX_W_M2 = 2.0
"""float: The surface heat flow of Io, the most tidally heated body there is (about 2 W/m2)."""


@pytest.fixture(scope="module")
def system():
    config = SystemConfig()
    config.MOONS = True
    with deterministic_entropy(9101):
        return StarSystem(system_config=config)


def _bodies(system):
    for planet in system.planets + system.secondary_planets:
        if planet.body_type == "a":
            continue
        yield planet
        yield from planet.moons


def test_the_last_integrator_of_the_ladder_is_the_implicit_one_that_finishes_a_stiff_history():
    assert moon_orbits.STRATEGIES[-1][0] == "Radau"


def test_no_body_of_the_system_is_heated_by_a_tide_more_than_io_is(system):
    flux = {body.name: (getattr(body, "tidal_flux_w_m2", None) or 0.0) for body in _bodies(system)}
    assert max(flux.values()) < 10.0 * IO_FLUX_W_M2, {name: value for name, value in flux.items() if value >= 10.0 * IO_FLUX_W_M2}


def test_the_planet_with_the_stiff_history_lost_its_close_moon_to_the_tide(system):
    planet = next(p for p in system.planets if getattr(p, "name", "").endswith("Miorutor III"))
    assert not planet.moons and planet.fallen_moons, "the 5e23 kg moon born at 3.8 radii falls to its planet under the tide (1.4 Gyr)"
    assert abs(planet.fallen_moons[0].mass_kg - 5.0e23) < 1.0e22
