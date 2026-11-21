"""
The orbit a planet's position and velocity describe is the orbit its eccentricity names (GEN.1817). A chain planet's forced eccentricity (oddity W5) can move by less than
`planets.ECCENTRICITY_SETTLE_TOLERANCE`, which leaves the air and the class as they are; the position, the moons' host eccentricity and the light go with it all the same.
"""

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.system import StarSystem
from planetgen.physics import planets
from planetgen.util import draw

SMALL_MOVE = planets.ECCENTRICITY_SETTLE_TOLERANCE / 2.0
"""float: An eccentricity change under the settle tolerance."""


@pytest.fixture(scope="module")
def planet():
    draw.set_run_seed(2117)
    for _ in range(30):
        system = StarSystem(system_config=SystemConfig())
        for body in system.primary.orbiters:
            if isinstance(body, Planet) and body.body_type in ("t", "g") and body.eccentricity > 0.0:
                return body
    pytest.fail("no planet on an eccentric orbit in 30 systems")


def test_a_move_under_the_tolerance_keeps_the_vector_on_the_eccentricity(planet, monkeypatch):
    target = planet.eccentricity + SMALL_MOVE

    def force(body):
        body.eccentricity = target

    monkeypatch.setattr(planets, "apply_orbit_oddities", force)
    assert planets.finish_orbit_oddities(planet) is False      # the air is not made again
    assert planet.orbit_from_vector()["eccentricity"] == pytest.approx(planet.eccentricity, abs=1e-9)
    for moon in planet.moons:
        assert moon.host_eccentricity == planet.eccentricity


def test_no_move_leaves_the_planet_as_it_is(planet, monkeypatch):
    before = (planet.position_x, planet.position_y, planet.position_z)
    monkeypatch.setattr(planets, "apply_orbit_oddities", lambda body: None)
    assert planets.finish_orbit_oddities(planet) is False
    assert (planet.position_x, planet.position_y, planet.position_z) == before
