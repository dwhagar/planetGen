# tests/test_orbit_from_vector.py

"""GEN.122: a planet's, moon's or comet's orbit is worked out from its current state vector, never stored, so it
follows the vector; the closed ellipse is its projected course."""

import math

import pytest

from planetgen.generation.comet import Comet
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.system import StarSystem
from planetgen.physics import state_vectors as sv
from planetgen.physics.planets import update_orbital_position

TWO_PI = 2.0 * math.pi


def _system(predicate):
    for _ in range(400):
        system = StarSystem(system_config=SystemConfig())
        if predicate(system):
            return system
    pytest.fail("no system matching")


def _angle_close(a, b, tol=1e-6):
    return abs((a - b + math.pi) % TWO_PI - math.pi) < tol


def test_a_planets_orbit_from_its_vector_is_the_circle_it_was_generated_on():
    system = _system(lambda s: any(isinstance(p, Planet) for p in s.planets))
    planet = next(p for p in system.planets if isinstance(p, Planet))
    orbit = planet.orbit_from_vector()
    assert orbit["kind"] == "elliptical"
    assert orbit["eccentricity"] == pytest.approx(0.0, abs=1e-9)
    assert orbit["semi_major_axis"] == pytest.approx(planet.distance, rel=1e-9)
    assert orbit["period"] == pytest.approx(planet.period, rel=1e-9)
    assert orbit["inclination"] == pytest.approx(math.radians(planet.orbital_inclination_deg), abs=1e-9)
    if planet.orbital_inclination_deg > 1e-6:
        assert _angle_close(orbit["ascending_node"], math.radians(planet.orbital_ascending_node_deg))
    assert _angle_close(orbit["true_anomaly"], math.radians(planet.orbital_phase_deg))


def test_a_moons_orbit_from_its_vector_is_about_its_planet():
    system = _system(lambda s: any(getattr(p, "moons", None) for p in s.planets))
    planet = next(p for p in system.planets if getattr(p, "moons", None))
    moon = planet.moons[0]
    orbit = moon.orbit_from_vector()
    assert orbit["semi_major_axis"] == pytest.approx(moon.distance, rel=1e-9)
    assert orbit["period"] == pytest.approx(moon.period, rel=1e-9)


def test_the_projected_course_is_the_closed_ellipse_round_the_primary():
    system = _system(lambda s: any(isinstance(p, Planet) and p.distance > 0 for p in s.planets))
    planet = next(p for p in system.planets if isinstance(p, Planet))
    points = planet.projected_orbit_au(64)
    assert len(points) == 64
    assert all(math.dist(point, (0.0, 0.0, 0.0)) == pytest.approx(planet.distance, rel=1e-9) for point in points)
    # The body is on its course: its position is a point of the ellipse it is on.
    orbit = planet.orbit_from_vector()
    position, _velocity = sv.state_from_elements(
        planet.orbit_mu_au3_per_year2(), orbit["periapsis_distance"], orbit["eccentricity"], orbit["inclination"],
        orbit["ascending_node"], orbit["arg_periapsis"], orbit["true_anomaly"])
    assert position == pytest.approx((planet.position_x, planet.position_y, planet.position_z), rel=1e-6, abs=1e-9)


def test_an_elliptical_comets_orbit_from_its_vector_matches_its_elements():
    system = _system(lambda s: any(c.orbit_type == "elliptical" for c in s.comets))
    comet = next(c for c in system.comets if c.orbit_type == "elliptical")
    orbit = comet.orbit_from_vector()
    assert orbit["kind"] == "elliptical"
    assert orbit["periapsis_distance"] == pytest.approx(comet.perihelion_distance_au, rel=1e-6)
    assert orbit["eccentricity"] == pytest.approx(comet.eccentricity, rel=1e-6, abs=1e-9)
    assert orbit["inclination"] == pytest.approx(math.radians(comet.inclination_deg), rel=1e-6, abs=1e-9)
    assert sv.mean_anomaly_from_true(orbit["true_anomaly"], orbit["eccentricity"]) == pytest.approx(
        math.radians(comet.mean_anomaly_deg) % TWO_PI, abs=1e-5)
    assert len(comet.projected_orbit_au(32)) == 32


def test_a_comet_that_does_not_close_has_no_projected_ellipse():
    system = _system(lambda s: any(c.orbit_type != "elliptical" for c in s.comets))
    comet = next(c for c in system.comets if c.orbit_type != "elliptical")
    assert comet.orbit_from_vector()["kind"] in ("parabolic", "hyperbolic")
    with pytest.raises(ValueError):
        comet.projected_orbit_au()


def test_the_orbit_follows_the_vector_when_it_changes_not_the_stored_elements():
    system = _system(lambda s: any(isinstance(p, Planet) for p in s.planets))
    planet = next(p for p in system.planets if isinstance(p, Planet))
    circle = planet.orbit_from_vector()
    # A kick along the motion (a flyby's pull) makes the orbit an ellipse with periapsis where it was.
    planet.set_velocity_kms(*(1.1 * v for v in (planet.velocity_x_kms, planet.velocity_y_kms, planet.velocity_z_kms)))
    kicked = planet.orbit_from_vector()
    assert kicked["eccentricity"] == pytest.approx(0.21, abs=1e-6)  # 1.1**2 - 1
    assert kicked["periapsis_distance"] == pytest.approx(planet.distance, rel=1e-9)
    assert kicked["semi_major_axis"] > circle["semi_major_axis"]
    assert planet.distance == pytest.approx(circle["semi_major_axis"], rel=1e-12)  # the generated elements are untouched
    # The orbital update puts the circle back, and the derived orbit with it.
    update_orbital_position(planet)
    assert planet.orbit_from_vector()["eccentricity"] == pytest.approx(0.0, abs=1e-9)


def test_a_body_with_no_position_has_no_orbit_yet():
    comet = Comet.__new__(Comet)
    comet.spatial = None
    with pytest.raises(ValueError):
        comet.orbit_from_vector()
