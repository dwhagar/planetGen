# tests/test_mass_object.py

"""GEN.565: the shared core of every body with a mass (`mass_object.MassObject`)."""

import pytest

from planetgen.galaxy.galactic_orbit import generate_galactic_orbit_fields
from planetgen.generation.comet import Comet
from planetgen.generation.mass_object import GALACTIC_MOTION_FIELDS, MassObject
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.util import draw

FAMILIES = [Star, BlackHole, NeutronStar, Planet, Comet, RoguePlanet, InterstellarComet, AsteroidField]


@pytest.mark.parametrize("family", FAMILIES)
def test_every_family_is_a_mass_object_with_unique_field_names(family):
    assert issubclass(family, MassObject)
    assert len(family.SERIALIZABLE_FIELDS) == len(set(family.SERIALIZABLE_FIELDS))


@pytest.mark.parametrize("family", [Star, BlackHole, NeutronStar, RoguePlanet, InterstellarComet, AsteroidField])
def test_a_body_in_a_galactic_orbit_persists_the_four_motion_fields(family):
    assert set(GALACTIC_MOTION_FIELDS) <= set(family.SERIALIZABLE_FIELDS)


def test_set_galactic_motion_draws_what_the_orbit_function_draws():
    expected = None
    with draw.bound(draw.Stream("mass-object-test")):
        expected = generate_galactic_orbit_fields(30000.0, None)
    body = MassObject()
    with draw.bound(draw.Stream("mass-object-test")):
        body.set_galactic_motion(30000.0)
    assert tuple(getattr(body, field) for field in GALACTIC_MOTION_FIELDS) == expected


def test_to_dict_applies_encode_and_from_dict_applies_decode():
    class Body(MassObject):
        SERIALIZABLE_FIELDS = ["tags"]

        def _encode(self, data):
            data["tags"] = list(self.tags)

        def _decode(self, data):
            self.tags = tuple(data["tags"])

    body = object.__new__(Body)
    body.tags = ("a", "b")
    data = body.to_dict()
    assert data == {"tags": ["a", "b"]}
    again = Body.from_dict(data, "config")
    assert again.tags == ("a", "b") and again.system_config == "config"
    assert Body._rebuild(data, "config", star="s").star == "s"


# --- Motion (GEN.565b): one advance() for a planet, a moon, a comet, a star and a wanderer -----------------------------

import math

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import body_state, constants


def _first(predicate):
    for _ in range(400):
        system = StarSystem(system_config=SystemConfig())
        found = predicate(system)
        if found:
            return found
    pytest.fail("no system matching")


def test_a_planets_advance_is_the_state_step_its_batch_update_takes():
    planet = _first(lambda s: next((p for p in s.planets if isinstance(p, Planet)), None))
    start = (planet.position_x, planet.position_y, planet.position_z)
    velocity = (planet.velocity_x_kms, planet.velocity_y_kms, planet.velocity_z_kms)
    seconds = 0.3 * planet.period * constants.SECONDS_PER_YEAR
    km = constants.AU_M / 1000.0
    expected, _ = body_state.advance(tuple(v * km for v in start), velocity, planet.distance * km, planet.period, seconds)
    planet.advance(seconds)
    assert (planet.position_x * km, planet.position_y * km, planet.position_z * km) == pytest.approx(expected, rel=1e-9)


def test_a_planet_comes_back_after_one_period_and_a_step_back_undoes_a_step():
    planet = _first(lambda s: next((p for p in s.planets if isinstance(p, Planet)), None))
    start = (planet.position_x, planet.position_y, planet.position_z)
    period_s = planet.period * constants.SECONDS_PER_YEAR
    planet.advance(period_s)
    assert (planet.position_x, planet.position_y, planet.position_z) == pytest.approx(start, rel=1e-6, abs=1e-9)
    planet.advance(0.37 * period_s)
    planet.advance(-0.37 * period_s)
    assert (planet.position_x, planet.position_y, planet.position_z) == pytest.approx(start, rel=1e-6, abs=1e-9)


def test_a_moon_advances_about_its_planet_by_the_same_function():
    planet = _first(lambda s: next((p for p in s.planets if getattr(p, "moons", None)), None))
    moon = planet.moons[0]
    start = (moon.position_x, moon.position_y, moon.position_z)
    moon.advance(0.25 * moon.period * constants.SECONDS_PER_YEAR)
    assert (moon.position_x, moon.position_y, moon.position_z) != pytest.approx(start)
    distance = math.sqrt(moon.position_x ** 2 + moon.position_y ** 2 + moon.position_z ** 2)
    assert distance == pytest.approx(moon.distance, rel=1e-6)


def test_a_comet_advances_on_its_hosts_pull():
    comet = _first(lambda s: next((c for c in getattr(s, "comets", []) if isinstance(c, Comet)
                                   and c.orbit_type == "elliptical"), None))
    before = comet.orbit_from_vector()
    comet.advance(0.1 * before["period"] * constants.SECONDS_PER_YEAR)
    after = comet.orbit_from_vector()
    assert after["semi_major_axis"] == pytest.approx(before["semi_major_axis"], rel=1e-6)
    assert after["eccentricity"] == pytest.approx(before["eccentricity"], abs=1e-6)


@pytest.mark.parametrize("family", [Star, NeutronStar, RoguePlanet])
def test_a_body_that_goes_round_the_galaxy_advances_its_phase_by_the_orbits_period(family):
    with draw.bound(565):
        body = family(SystemConfig())
    start, period_years = body.galactic_orbital_phase_deg, body.galactic_orbital_period_gy * 1.0e9
    body.advance(0.25 * period_years * constants.SECONDS_PER_YEAR)
    assert body.galactic_orbital_phase_deg == pytest.approx((start + 90.0) % 360.0)
    body.advance(period_years * constants.SECONDS_PER_YEAR)
    assert body.galactic_orbital_phase_deg == pytest.approx((start + 90.0) % 360.0)


def test_a_body_with_no_trajectory_cannot_advance():
    with pytest.raises(ValueError):
        MassObject().advance(1.0)
