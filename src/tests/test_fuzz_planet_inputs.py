# tests/test_fuzz_planet_inputs.py

"""
Every way a caller can pin a planet's class, radius and mass before
generation (`planets.generate_planet_properties`'s eight branches:
none, class, radius, mass, class+radius, class+mass, radius+mass, all
three), in every zone, with values drawn from inside, on the edge of and
outside each class's declared ranges -- plus the scalar physics helpers
at their boundaries.

The contract under attack: a planet either comes out whole (a class
valid for its zone, radius and mass finite, positive and inside that
class's ranges) or construction fails with a plain `ValueError`. Nothing
else may escape, and nothing half-built may come back.
"""

import math

import pytest
from hypothesis import assume, example, given
from hypothesis import strategies as st

from planetgen.physics import planets
from planetgen import tuning as prog_c
from stellarObjects.config import SystemConfig
from stellarObjects.planetData import Planet
from stellarObjects.starData import Star
from tests.fuzz_support import any_float

CLASSES = sorted(prog_c.PLANET_CLASSES)
MASS_RANGES = planets.get_planet_mass_ranges()


@pytest.fixture(scope="module")
def star():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    return Star(config)


def _distance(star, zone):
    inner, outer = star.habitable_zone
    return {"h": inner * 0.5, "e": (inner + outer) / 2.0, "c": outer * 2.0}[zone]


def _radius_near(cls):
    low, high = prog_c.PLANET_CLASSES[cls]["radius_range"]
    return st.one_of(
        st.sampled_from([low, high, low - 1e-6, high + 1e-6, (low + high) / 2]),
        st.floats(low * 0.5, high * 1.5),
    )


def _mass_near(cls):
    low, high = MASS_RANGES[cls]
    return st.one_of(
        st.sampled_from([low, high, low * (1 - 1e-9), high * (1 + 1e-9), (low + high) / 2]),
        st.floats(low * 0.5, high * 1.5),
    )


@st.composite
def pinned_inputs(draw):
    cls = draw(st.sampled_from(CLASSES))
    other = draw(st.sampled_from(CLASSES))
    zone = draw(st.sampled_from("hec"))
    give_class, give_radius, give_mass = draw(st.tuples(st.booleans(), st.booleans(), st.booleans()))
    # Radius/mass drawn around a (possibly different) class's range, so
    # mismatched combinations get tried as often as matching ones.
    radius = draw(_radius_near(other if draw(st.booleans()) else cls)) if give_radius else None
    mass = draw(_mass_near(other if draw(st.booleans()) else cls)) if give_mass else None
    habitable = draw(st.sampled_from([None, True, False]))
    return zone, (cls if give_class else None), radius, mass, habitable


def _check_whole(planet):
    assert planet.planet_class in prog_c.PLANET_CLASSES
    for value in (planet.radius, planet.mass, planet.distance, planet.period):
        assert math.isfinite(value) and value > 0
    low, high = prog_c.PLANET_CLASSES[planet.planet_class]["radius_range"]
    assert low <= planet.radius <= high


@given(inputs=pinned_inputs())
def test_pinned_planets_come_out_whole_or_raise_value_error(star, inputs):
    zone, cls, radius, mass, habitable = inputs
    config = SystemConfig()
    config.HABITABLE_WORLD = habitable
    try:
        planet = Planet(config, star, star.habitable_zone, _distance(star, zone),
                        planet_class=cls, radius=radius, mass=mass)
    except ValueError:
        return
    _check_whole(planet)
    if cls is not None and radius is not None:
        assert planet.radius == radius
    if habitable is False and planet.zone == "e":
        assert planet.planet_class not in prog_c.HABITABLE_PLANET_CLASSES


@given(st.sampled_from(CLASSES), st.sampled_from("hec"))
def test_a_class_pinned_in_a_zone_it_does_not_allow_is_rejected(star, cls, zone):
    assume(not prog_c.PLANET_CLASSES[cls][zone])
    with pytest.raises(ValueError):
        Planet(SystemConfig(), star, star.habitable_zone, _distance(star, zone), planet_class=cls)


@given(st.sampled_from(CLASSES), st.sampled_from("hec"), st.data())
def test_a_radius_outside_the_pinned_class_is_rejected(star, cls, zone, data):
    assume(prog_c.PLANET_CLASSES[cls][zone])
    low, high = prog_c.PLANET_CLASSES[cls]["radius_range"]
    radius = data.draw(st.one_of(st.floats(1, low, exclude_max=True), st.floats(high, high * 3, exclude_min=True)))
    with pytest.raises(ValueError):
        Planet(SystemConfig(), star, star.habitable_zone, _distance(star, zone), planet_class=cls, radius=radius)


def test_unknown_planet_class_is_rejected(star):
    for cls in ("", "Z", "a", "EE", "?", "\x00"):
        if cls in prog_c.PLANET_CLASSES:
            continue
        with pytest.raises((ValueError, KeyError)):
            Planet(SystemConfig(), star, star.habitable_zone, _distance(star, "e"), planet_class=cls)


def test_every_class_has_a_sane_mass_and_radius_range():
    for cls, data in prog_c.PLANET_CLASSES.items():
        low, high = data["radius_range"]
        assert 0 < low < high, cls
        m_low, m_high = MASS_RANGES[cls]
        assert 0 < m_low < m_high and math.isfinite(m_high), cls
        assert any(data[z] for z in "hec"), f"class {cls} is allowed in no zone"


@given(st.sampled_from(CLASSES))
def test_sampled_class_radius_stays_in_range(cls):
    low, high = prog_c.PLANET_CLASSES[cls]["radius_range"]
    for _ in range(25):
        assert low <= planets._sample_class_radius(cls, low, high) <= high


@given(st.lists(st.sampled_from(CLASSES), min_size=1, unique=True))
def test_weighted_class_choice_only_returns_offered_classes(classes):
    assert planets._choose_weighted_planet_class(classes) in classes


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
