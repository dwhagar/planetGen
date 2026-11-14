"""GEN.582 (light): the birth and present-state classification of the stored fields of stars, planets and moons (`generation/body_state_fields.py`)."""

import copy

import pytest

from planetgen.generation import body_state_fields as fields
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import envelope_share, planets as planet_physics
from planetgen.util import draw


@pytest.mark.parametrize("kind", [Star, Planet])
def test_every_birth_name_is_a_stored_field_once(kind):
    birth = fields.birth_fields(kind)
    assert len(birth) == len(set(birth))
    assert set(birth) <= set(kind.SERIALIZABLE_FIELDS), sorted(set(birth) - set(kind.SERIALIZABLE_FIELDS))


@pytest.mark.parametrize("kind", [Star, Planet])
def test_every_stored_field_is_birth_or_present_and_never_both(kind):
    birth, present = set(fields.birth_fields(kind)), set(fields.present_fields(kind))
    assert not birth & present
    assert birth | present == set(kind.SERIALIZABLE_FIELDS)
    assert present and birth


def test_the_star_follows_its_age_for_the_state_and_keeps_the_draws_of_its_making():
    present = set(fields.present_fields(Star))
    for field in ("luminosity", "radius", "temperature", "yerkes_class", "habitable_zone", "snow_line_au", "rotation_period_hours", "l_xuv_w", "disk_kind"):
        assert field in present, field
    for field in ("age", "metallicity_fe_h", "mass", "lifespan"):
        assert fields.is_birth(Star, field), field


def test_the_planet_keeps_its_composition_and_orbit_draws_and_remakes_the_air():
    present = set(fields.present_fields(Planet))
    for field in ("surface_temperature", "atmospheric_pressure", "planet_class", "hydrosphere", "ocean_fraction", "radius", "ring_mass_kg", "life_class", "tidal_flux_w_m2"):
        assert field in present, field
    for field in ("mass", "distance", "rock_mg_si_ratio", "core_s_fraction", "c_u0_ppm", "moon_origin", "interior_kind", "born_beyond_snow_line"):
        assert fields.is_birth(Planet, field), field


AIR_SHARE_FIELDS = (*envelope_share.SOLIDS, "envelope_mass_fraction", "envelope_helium_fraction", "envelope_z_fraction")
"""tuple: The birth fields that the share of the air moves (`envelope_share.fill`)."""
ORBIT_DRAWS = ("distance", "eccentricity", "orbital_inclination_deg", "orbital_ascending_node_deg", "host_eccentricity")


def _snapshot(body):
    return {field: copy.deepcopy(getattr(body, field, None)) for field in Planet.SERIALIZABLE_FIELDS}


def test_a_remake_after_an_orbit_change_keeps_every_birth_value_and_works_the_present_state_out_again():
    # a rocky planet moved 12 percent outward and made again (`remake_after_orbit_change`): its birth fields (but the orbit draw it was moved off) are untouched, and the present fields
    # that depend on the orbit move
    moved, tested = set(), 0
    for seed in range(6):
        config = SystemConfig()
        config.PLANETS = True
        with draw.bound(seed):
            system = StarSystem(system_config=config)
        for planet in system.all_planets:
            if planet.body_type != "t":
                continue
            bodies = [planet, *planet.moons]
            before = {id(body): _snapshot(body) for body in bodies}
            with draw.bound(seed + 1000):
                planet.distance *= 1.12
                planet_physics.remake_after_orbit_change(planet)
            for body in bodies:
                tested += 1
                for field in fields.birth_fields(Planet):
                    if field in ORBIT_DRAWS:
                        continue            # the orbit the test moved; a moon's own orbit follows its parent's change a little (`follow_orbit_change`): the orbit draws keep their first values in GEN.582's draw record
                    if field in AIR_SHARE_FIELDS:       # the air's share comes out of the solid fractions (`envelope_share.fill`), so an air that appears or goes at the new distance moves them a little
                        assert getattr(body, field, None) == pytest.approx(before[id(body)][field], abs=1e-6), (seed, body.name, field)
                        continue
                    assert getattr(body, field, None) == before[id(body)][field], (seed, body.name, field)
                moved |= {field for field in fields.present_fields(Planet) if getattr(body, field, None) != before[id(body)][field]}
    assert tested >= 3
    for field in ("period", "surface_temperature", "orbital_speed_kms"):
        assert field in moved, field
