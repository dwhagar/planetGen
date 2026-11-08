"""
GEN.1334: a planet with a hydrogen-helium envelope of 0.5 percent of its mass or more has no ground (design: docs/design/composition/radius-valley-tuning.md section 8).

`planet_type.has_ground` is the one test. A planet under a deep envelope (a sub-Neptune, a gas dwarf) has its air solved at the 1 bar level as a giant's is (a hydrogen and helium
mixture with its metals), holds no ocean, land, ground, surface dose or liquid, and is scored for no life.
"""

import collections

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import (atmosphere, habitability_world, hydrosphere, life_classes, planet_type, radiation, regional_temperature, surface_composition,
                               surface_liquid_world)
from planetgen.util import draw

SYSTEMS = 12
SEED = 20261108

GROUND_FIELDS = tuple(dict.fromkeys((
    *hydrosphere.HYDROSPHERE_FIELDS, *surface_composition.COLUMN_FIELDS, *radiation.BODY_FIELDS, *habitability_world.BODY_FIELDS, *habitability_world.LIKELIHOOD_FIELDS,
    *surface_liquid_world.COLUMN_FIELDS, *life_classes.BODY_FIELDS, *regional_temperature.BAND_COLUMNS)))
"""tuple: What a giant leaves empty because it has no surface; a planet under a deep envelope leaves the same."""

# the fields a body keeps whatever it stands on: the bulk water of its composition, the swing of its starlight and the energy it receives
KEPT = {"water_mass_fraction", "flux_variation", "energy_flux_w_m2"}


@pytest.fixture(scope="module")
def sample():
    """Every planet and moon of `SYSTEMS` G2V systems."""
    draw.set_run_seed(SEED)
    bodies = []
    for _ in range(SYSTEMS):
        config = SystemConfig()
        config.STAR_TYPE = "G2V"
        for planet in StarSystem(system_config=config).all_planets:
            if getattr(planet, "body_type", None) in ("t", "g"):
                bodies.append(planet)
                bodies.extend(getattr(planet, "moons", None) or [])
    return bodies


def test_a_body_has_ground_unless_it_is_a_giant_or_holds_half_a_percent_of_its_mass_as_envelope():
    assert planet_type.has_ground({"body_type": "t", "envelope_mass_fraction": 0.0})
    assert planet_type.has_ground({"body_type": "t", "envelope_mass_fraction": 0.004999})
    assert not planet_type.has_ground({"body_type": "t", "envelope_mass_fraction": 0.005})
    assert not planet_type.has_ground({"body_type": "t", "envelope_mass_fraction": 0.2})
    assert not planet_type.has_ground({"body_type": "g", "envelope_mass_fraction": 0.9})
    assert planet_type.has_ground({"body_type": "t"})            # an envelope share that is not stored: the body type's
    assert not planet_type.has_ground({"body_type": "g"})


def test_the_sample_holds_planets_under_an_envelope_and_planets_with_ground(sample):
    covered = [body for body in sample if body.body_type == "t" and not planet_type.has_ground(body)]
    grounded = [body for body in sample if body.body_type == "t" and planet_type.has_ground(body)]
    assert len(covered) >= 5 and len(grounded) >= 5, (len(covered), len(grounded))


def test_no_planet_under_an_envelope_carries_a_surface_field(sample):
    carried = collections.Counter(field for body in sample if body.body_type == "t" and not planet_type.has_ground(body)
                                  for field in GROUND_FIELDS if field not in KEPT and getattr(body, field, None) not in (None, {}, [], ()))
    assert not carried, dict(carried)


def test_a_planet_with_ground_still_carries_its_surface_fields(sample):
    grounded = [body for body in sample if body.body_type == "t" and planet_type.has_ground(body) and body.atmospheric_pressure]
    assert grounded and all(body.ocean_fraction is not None and body.ground_type is not None for body in grounded)


def test_the_air_of_a_planet_under_an_envelope_is_the_hydrogen_and_helium_of_its_1_bar_level(sample):
    covered = [body for body in sample if body.body_type == "t" and not planet_type.has_ground(body)]
    for body in covered:
        partials = atmosphere.gas_partials_pa(body)
        total = sum(partials.values())
        assert body.atmospheric_pressure == pytest.approx(1.0e5, rel=1e-3), body.atmospheric_pressure
        assert partials.get("H2", 0.0) / total > 0.5 and partials.get("He", 0.0) > 0.0, partials


def test_a_planet_under_an_envelope_is_scored_for_no_life_and_has_none(sample):
    for body in sample:
        if body.body_type == "t" and not planet_type.has_ground(body):
            assert body.life_class is None and body.life_living in (None, {}, []) and body.phi4 is None
