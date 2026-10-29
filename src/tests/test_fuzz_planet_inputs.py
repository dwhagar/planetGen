# tests/test_fuzz_planet_inputs.py

"""
Every way a caller can pin a planet's radius and mass before generation (`planets.generate_planet_properties`'s four branches: neither,
radius, mass, both), in every zone, with values drawn from inside, on the edge of and outside what bodies have -- plus the scalar
physics helpers at their boundaries.

The contract under attack: a planet either comes out whole (a class derived from what it is, radius, mass, distance and period finite
and positive, a given radius or mass kept) or construction fails with a plain `ValueError`. Nothing else may escape, and nothing
half-built may come back.
"""

import math

import pytest
from hypothesis import assume, example, given, settings
from hypothesis import strategies as st

from planetgen.physics import planets
from planetgen import tuning as prog_c
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from tests.fuzz_support import any_float

RADII_KM = [1.0, prog_c.HOSTING_SMALL_BODY_MAX_RADIUS_KM, 2000.0, 6371.0, 15000.0, prog_c.GIANT_MIN_RADIUS_KM, 70000.0, 200000.0]
MASSES_KG = [prog_c.PLANET_MIN_MASS_KG * 0.5, prog_c.PLANET_MIN_MASS_KG, 1e23, 5.97e24, 1.2e26, 1.9e27, 1e29]


@pytest.fixture(scope="module")
def star():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    return Star(config)


def _distance(star, zone):
    inner, outer = star.habitable_zone
    return {"h": inner * 0.5, "e": (inner + outer) / 2.0, "c": outer * 2.0}[zone]


def _near(values):
    low, high = min(values), max(values)
    return st.one_of(st.sampled_from(values), st.floats(low * 0.5, high * 1.5))


@st.composite
def pinned_inputs(draw):
    zone = draw(st.sampled_from("hec"))
    give_radius, give_mass = draw(st.tuples(st.booleans(), st.booleans()))
    radius = draw(_near(RADII_KM)) if give_radius else None
    mass = draw(_near(MASSES_KG)) if give_mass else None
    habitable = draw(st.sampled_from([None, True, False]))
    return zone, radius, mass, habitable


def _check_whole(planet):
    assert planet.planet_class in prog_c.PLANET_CLASSES
    assert planet.planet_class == planets.planet_classes.classify(planet)
    for value in (planet.radius, planet.mass, planet.distance, planet.period):
        assert math.isfinite(value) and value > 0


@settings(max_examples=30, deadline=None)
@given(inputs=pinned_inputs())
def test_pinned_planets_come_out_whole_or_raise_value_error(star, inputs):
    zone, radius, mass, habitable = inputs
    config = SystemConfig()
    config.HABITABLE_WORLD = habitable
    try:
        planet = Planet(config, star, star.habitable_zone, _distance(star, zone), radius=radius, mass=mass, moon_count=0)
    except ValueError:
        return
    _check_whole(planet)
    if radius is not None:
        assert planet.radius == radius
    if mass is not None:
        assert planet.mass == mass
    if habitable is False and planet.zone == "e":
        assert planet.planet_class not in prog_c.HABITABLE_PLANET_CLASSES


def test_a_planet_class_can_no_longer_be_pinned(star):
    with pytest.raises(TypeError):
        Planet(SystemConfig(), star, star.habitable_zone, _distance(star, "e"), planet_class="M")


def test_every_class_has_a_description_composition_and_type():
    for cls, data in prog_c.PLANET_CLASSES.items():
        assert data["description"] and data["composition"] and data["type"] in ("t", "g"), cls


# --- orbital period ------------------------------------------------------------

@given(st.floats(1e-6, 1e6), st.floats(1e20, 1e33))
def test_orbital_period_obeys_keplers_third_law(distance, mass):
    period = planets.calculate_orbital_period_years(distance, mass)
    assert math.isfinite(period) and period > 0
    doubled = planets.calculate_orbital_period_years(distance * 4, mass)
    assert doubled == pytest.approx(period * 8, rel=1e-9)


@given(any_float, any_float)
@example(math.nan, 1e30)
@example(1.0, math.nan)
@example(math.inf, 1e30)
def test_orbital_period_rejects_non_positive_inputs_cleanly(distance, mass):
    """Non-positive, NaN or infinite inputs are all a ValueError (NaN used
    to slip past the `<= 0` check and come back as NaN)."""
    assume(not (distance > 0 and mass > 0 and math.isfinite(distance) and math.isfinite(mass)))
    with pytest.raises(ValueError):
        planets.calculate_orbital_period_years(distance, mass)
