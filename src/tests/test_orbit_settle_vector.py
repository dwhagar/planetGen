"""
The orbit a planet's position and velocity describe is the orbit its eccentricity names (GEN.1817). A chain planet's forced eccentricity (oddity W5) can move by less than
`planets.ECCENTRICITY_SETTLE_TOLERANCE`, which leaves the air and the class as they are; the position, the moons' host eccentricity and the light go with it all the same.
"""

import math

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.system import StarSystem
from planetgen.physics import planets, spin
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


@pytest.fixture(scope="module")
def planet_with_moons():
    draw.set_run_seed(2118)
    for _ in range(40):
        system = StarSystem(system_config=SystemConfig())
        for body in system.primary.orbiters:
            if isinstance(body, Planet) and body.body_type in ("t", "g") and body.moons and spin.spin_axis_of(body) is not None:
                return body
    pytest.fail("no planet with moons and a spin in 40 systems")


def _offset_from_plane(body):
    """How far, AU, the position of a planet or moon lies from the plane its inclination and node name (a moon's in the frame of its planet's equator, GEN.1352)."""
    position = (body.position_x, body.position_y, body.position_z)
    axis = body.build.frame_axis
    if axis is not None:
        position = spin.turn_from_axis(position, axis)
    inclination, node = math.radians(body.orbital_inclination_deg), math.radians(body.orbital_ascending_node_deg)
    normal = (math.sin(node) * math.sin(inclination), -math.cos(node) * math.sin(inclination), math.cos(inclination))
    return abs(sum(c * n for c, n in zip(position, normal)))


def test_a_row_that_turns_the_plane_after_the_exchange_takes_the_planet_with_it(planet_with_moons, monkeypatch):
    """The final roll of the orbit rows (`finish_orbit_oddities`) gives a polar orbit (W4) a fresh draw of the inclination, after the secular exchange or the Kozai pump has moved the plane, and with
    the eccentricity unmoved nothing else moves the position: the giant of M50 at 92 AU lay 40 AU off the plane it named."""
    body = planet_with_moons
    old_inclination = body.orbital_inclination_deg
    polar = 100.0 if abs(old_inclination - 100.0) > 5.0 else 80.0
    obliquity_before = _obliquity_deg(body)
    radius_before = math.sqrt(body.position_x ** 2 + body.position_y ** 2 + body.position_z ** 2)
    eccentricity = body.eccentricity

    def turn(planet):
        planet.orbital_inclination_deg = polar

    monkeypatch.setattr(planets, "apply_orbit_oddities", turn)
    assert planets.finish_orbit_oddities(body) is False
    assert body.orbital_inclination_deg == polar
    assert body.eccentricity == eccentricity
    assert _offset_from_plane(body) < 1.0e-9 * radius_before
    assert math.sqrt(body.position_x ** 2 + body.position_y ** 2 + body.position_z ** 2) == pytest.approx(radius_before, rel=1e-9)      # it keeps its place along an orbit of the same size
    assert _obliquity_deg(body) == pytest.approx(obliquity_before, abs=1.0e-6)      # the spin axis turned with the plane
    for moon in body.moons:
        assert _offset_from_plane(moon) < 1.0e-9 * math.sqrt(moon.position_x ** 2 + moon.position_y ** 2 + moon.position_z ** 2)


def _obliquity_deg(body):
    """The angle between the spin axis of a planet and the normal of its orbit, degrees."""
    axis, normal = spin.spin_axis_of(body), spin.orbit_normal(body.orbital_inclination_deg, body.orbital_ascending_node_deg)
    return math.degrees(math.acos(max(-1.0, min(1.0, sum(a * n for a, n in zip(axis, normal))))))
