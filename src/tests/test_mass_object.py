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
