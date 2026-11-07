# tests/test_moon_frame.py

"""
GEN.1352: a moon's orbit is drawn against its planet's equator and set in the system's frame by a turn onto the planet's spin axis
(`planets.orient_moons`), so the moons of a planet lying on its side go round in the tipped equator.
"""

import math

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import planets, spin


def _angle_deg(a, b):
    return spin.obliquity_between(a, b)


def _axis(planet):
    return (planet.spin_axis_x, planet.spin_axis_y, planet.spin_axis_z)


def _orbit_normal(body):
    return spin.orbit_normal_from_state((body.position_x, body.position_y, body.position_z),
                                        (body.velocity_x_kms, body.velocity_y_kms, body.velocity_z_kms))


def _planets_with_moons(count, giants_only=False):
    found = []
    for _ in range(count * 40):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.MOONS = True
        cfg.MAX_PLANETS = True
        for planet in StarSystem(system_config=cfg).planets:
            if planet.body_type != "a" and planet.moons and (planet.body_type == "g" or not giants_only):
                found.append(planet)
        if len(found) >= count:
            return found
    pytest.fail("could not generate planets with moons")


@pytest.fixture(scope="module")
def giants():
    return _planets_with_moons(4, giants_only=True)


def test_a_turn_onto_an_axis_and_back_is_the_identity_and_an_upright_pole_leaves_a_vector_alone():
    vector = (0.3, -0.5, 0.8)
    assert spin.turn_to_axis(vector, (0.0, 0.0, 1.0)) == vector
    for axis in ((0.6, 0.0, 0.8), (0.0, 1.0, 0.0), (-0.48, 0.64, -0.6), (0.0, 0.0, -1.0)):
        there = spin.turn_to_axis(vector, axis)
        assert math.isclose(sum(c * c for c in there), sum(c * c for c in vector))
        assert spin.turn_from_axis(there, axis) == pytest.approx(vector)
    assert spin.turn_to_axis((0.0, 0.0, 1.0), (0.6, 0.0, 0.8)) == pytest.approx((0.6, 0.0, 0.8))


def test_every_moon_goes_round_in_its_planets_equator_by_the_inclination_it_keeps(giants):
    for planet in giants:
        for moon in planet.moons:
            assert _angle_deg(_orbit_normal(moon), _axis(planet)) == pytest.approx(moon.orbital_inclination_deg, abs=1.0e-6)
            assert moon.spin_axis_x is not None


def test_the_moons_of_a_giant_lying_on_its_side_lie_in_its_tipped_equator(giants):
    planet = giants[0]
    before = {id(m): (m.orbital_inclination_deg, m.orbital_ascending_node_deg, m.distance) for m in planet.moons}
    spin.set_spin(planet, spin.orbit_normal(planet.orbital_inclination_deg, planet.orbital_ascending_node_deg), 98.0)
    assert planet.axial_tilt_deg == 98.0
    planets.orient_moons(planet)
    planet_normal = spin.orbit_normal(planet.orbital_inclination_deg, planet.orbital_ascending_node_deg)
    for moon in planet.moons:
        assert (moon.orbital_inclination_deg, moon.orbital_ascending_node_deg, moon.distance) == before[id(moon)]
        normal = _orbit_normal(moon)
        assert _angle_deg(normal, _axis(planet)) == pytest.approx(moon.orbital_inclination_deg, abs=1.0e-6)
        assert abs(_angle_deg(normal, planet_normal) - 98.0) <= moon.orbital_inclination_deg + 1.0e-6     # on its side, not in the plane of the planet's orbit
        assert _angle_deg(_axis(moon), normal) == pytest.approx(moon.axial_tilt_deg, abs=1.0e-6)          # the moon's own tilt is against its own orbit
        assert math.sqrt(moon.position_x ** 2 + moon.position_y ** 2 + moon.position_z ** 2) == pytest.approx(
            moon.distance, rel=max(1.0e-9, (moon.eccentricity or 0.0) + 1.0e-9))


def test_orienting_the_moons_again_does_not_turn_them_twice(giants):
    planet = giants[-1]
    planets.orient_moons(planet)
    state = [(m.position_x, m.position_y, m.position_z, m.spin_axis_x) for m in planet.moons]
    planets.orient_moons(planet)
    for moon, kept in zip(planet.moons, state):
        assert (moon.position_x, moon.position_y, moon.position_z, moon.spin_axis_x) == pytest.approx(kept)


def test_the_moons_of_an_upright_planet_lie_in_the_plane_of_its_orbit_to_its_inclination():
    for planet in _planets_with_moons(3):
        planet_normal = spin.orbit_normal(planet.orbital_inclination_deg, planet.orbital_ascending_node_deg)
        for moon in planet.moons:
            assert _angle_deg(_orbit_normal(moon), planet_normal) <= planet.axial_tilt_deg + moon.orbital_inclination_deg + 1.0e-6
