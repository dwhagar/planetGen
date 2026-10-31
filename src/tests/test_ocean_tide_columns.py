# tests/test_ocean_tide_columns.py

"""
GEN.1034: the stored ocean-tide columns (`ocean_tides.compute` / `generate`): an Earth with its Moon and Sun comes out near the published tide, a sea with no
moon shows only the star's, a world locked to its star shows none from it, an ice-covered ocean is zero, a body without a free sea is empty, and a moon's
sea is raised by its planet.
"""

import math
from types import SimpleNamespace

import pytest

from planetgen.physics import constants, ocean_tides, spin

EARTH_KG, EARTH_KM = 5.972e24, 6371.0
MOON_KG, MOON_AU = 7.342e22, 384400.0 / 149597870.7
SUN = SimpleNamespace(mass=1.989e30, radius=695700.0, name="star")


def _body(**given):
    values = dict(mass=EARTH_KG, radius=EARTH_KM, hydrosphere="surface ocean", rotation_period_hours=24.0, distance=1.0, eccentricity=0.0167,
                  orbital_inclination_deg=0.0, orbital_ascending_node_deg=0.0, spin_axis_x=0.0, spin_axis_y=0.0, spin_axis_z=1.0, axial_tilt_deg=23.4,
                  body_type="t", planet_class="H", star=SUN, moons=[], water_mass_fraction=0.0003, is_moon=False, period=1.0)
    values.update(given)
    return SimpleNamespace(**values)


def _moon(**given):
    values = dict(mass=MOON_KG, distance=MOON_AU, eccentricity=0.055, orbital_inclination_deg=5.0, orbital_ascending_node_deg=0.0)
    values.update(given)
    return SimpleNamespace(**values)


def test_the_fields_are_the_two_ranges_and_the_liquid():
    assert ocean_tides.FIELDS == ("tide_range_spring_m", "tide_range_neap_m", "ocean_liquid")


def test_an_earth_with_its_moon_and_sun_has_a_spring_range_of_the_right_order_and_a_smaller_neap_range():
    earth = _body(moons=[_moon()])
    found = ocean_tides.compute(earth)
    assert 0.3 < found["tide_range_spring_m"] < 1.2, found       # the published equilibrium spring range is about 0.8 m
    assert 0.0 < found["tide_range_neap_m"] < found["tide_range_spring_m"]
    assert found["ocean_liquid"] == "water"


def test_a_sea_with_no_moon_shows_only_the_stars_tide_and_it_has_one_range():
    found = ocean_tides.compute(_body())
    assert 0.1 < found["tide_range_spring_m"] < 0.35            # the Sun's 0.25 m less the floor's share
    assert found["tide_range_neap_m"] == pytest.approx(found["tide_range_spring_m"])


def test_a_world_locked_to_its_star_has_no_tide_from_it_and_a_moon_only_the_moons():
    year_h = 365.25 * 24.0
    locked = ocean_tides.compute(_body(rotation_period_hours=year_h))
    assert locked["tide_range_spring_m"] == 0.0 and locked["ocean_liquid"] == "water"
    with_moon = ocean_tides.compute(_body(rotation_period_hours=year_h, moons=[_moon()]))
    alone = ocean_tides.compute(_body(moons=[_moon()]))
    assert 0.0 < with_moon["tide_range_spring_m"] < alone["tide_range_spring_m"]


def test_a_retrograde_orbit_cannot_lock_the_body_it_goes_round_the_other_way_from():
    year_h = 365.25 * 24.0
    against = _body(rotation_period_hours=year_h, spin_axis_z=-1.0)       # spins opposite to the planet's orbit
    assert ocean_tides.compute(against)["tide_range_spring_m"] > 0.0


def test_an_ice_covered_ocean_is_zero_and_a_body_without_a_free_sea_is_empty():
    ice = ocean_tides.compute(_body(hydrosphere="ice-covered ocean", moons=[_moon()]))
    assert ice == {"tide_range_spring_m": 0.0, "tide_range_neap_m": 0.0, "ocean_liquid": "water"}
    for state in ("dry", "vapour", "ice", "brine seep", "hycean", None):
        assert ocean_tides.compute(_body(hydrosphere=state, moons=[_moon()])) == dict.fromkeys(ocean_tides.FIELDS)
    assert ocean_tides.compute(_body(rotation_period_hours=None)) == dict.fromkeys(ocean_tides.FIELDS)
    assert ocean_tides.compute(_body(spin_axis_x=None)) == dict.fromkeys(ocean_tides.FIELDS)


def test_a_moons_sea_is_raised_by_its_planet_and_the_star():
    giant = _body(mass=1.898e27, radius=69911.0, hydrosphere=None, distance=5.2, body_type="g")
    moon = _body(mass=1.5e22, radius=1560.0, hydrosphere="surface ocean", distance=0.0044, rotation_period_hours=3.55 * 24.0 * 0.0 + 85.2,
                 eccentricity=0.01, is_moon=True, planet_class="H")
    free = ocean_tides.compute(moon, giant)
    assert free["tide_range_spring_m"] > 0.0
    locked_moon = ocean_tides.compute(_body(mass=1.5e22, radius=1560.0, distance=0.0044, eccentricity=0.01, is_moon=True,
                                            rotation_period_hours=2.0 * math.pi / math.sqrt(constants.G * (1.898e27 + 1.5e22) / (0.0044 * constants.AU_TO_M) ** 3) / 3600.0),
                                      giant)
    assert 0.0 <= locked_moon["tide_range_spring_m"] < free["tide_range_spring_m"]


def test_generate_sets_the_attributes_on_the_body():
    earth = _body(moons=[_moon()])
    ocean_tides.generate(earth)
    assert earth.tide_range_spring_m > earth.tide_range_neap_m > 0.0 and earth.ocean_liquid == "water"
    dry = _body(hydrosphere="dry")
    ocean_tides.generate(dry)
    assert dry.tide_range_spring_m is None and dry.ocean_liquid is None


def test_the_unit_registry_knows_the_columns_in_metres():
    from planetgen.db import column_units
    for column in ("tide_range_spring_m", "tide_range_neap_m"):
        assert column_units.unit_of("planets", column) == "m" or column_units.unit_of("moons", column) == "m", column


def test_a_finished_system_has_the_tide_of_its_final_water():
    """The final-age pass makes the water again (a brine seep can end as an ice-covered ocean), so the tides are made after it: every body
    of a finished system holds what `ocean_tides.compute` gives for it now (seed 7100 had an ice-covered ocean on a class E planet with none)."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    draw.set_run_seed(7100)
    config = SystemConfig()
    config.PLANETS = True
    config.MAX_PLANETS = True
    system = StarSystem(system_config=config, galactic_center_dist_ly=26000.0)
    checked = 0
    for planet in system.all_planets:
        if planet.body_type == "a":
            continue
        for body, parent in [(planet, None)] + [(moon, planet) for moon in planet.moons]:
            assert {name: getattr(body, name) for name in ocean_tides.FIELDS} == ocean_tides.compute(body, parent)
            checked += body.hydrosphere == "ice-covered ocean"
    assert checked, "the seed made no ice-covered ocean"
