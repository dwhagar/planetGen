# tests/test_envelope_capture.py

"""
GEN.1323: a rocky planet of 1 to 20 Earth masses captures a hydrogen-helium envelope and keeps it unless its star strips it
(`physics/envelope_capture.py`, `planets.settle_envelope`, the envelope pass of `StarSystem`): sub-Neptunes and the radius valley. The loss
itself is tested in test_envelope_escape.py and the survey that tuned it in test_radius_valley.py.
"""

import math
import statistics
from types import SimpleNamespace

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import bulk_structure, constants, envelope_capture, envelope_escape, envelope_share, planet_type, planets, radius_valley
from planetgen.util import draw
from tests import class_support

EARTH = constants.EARTH_MASS_TO_KG


def _body(mass_earth, distance_au, star, eccentricity=0.0):
    return SimpleNamespace(star=star, mass=mass_earth * EARTH, distance=distance_au, eccentricity=eccentricity)


@pytest.fixture()
def sun():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    with draw.bound("envelope-capture:star"):
        star = Star(cfg)
    star.age = 4.5
    return star


def test_the_start_is_centred_where_the_valley_was_tuned_and_falls_below_it():
    assert envelope_capture.REFERENCE_MASS_EARTH == radius_valley.CORE_MASS_MEDIAN_EARTH
    assert envelope_capture.START_SPREAD == radius_valley.ENVELOPE_SIGMA
    low, high = 0.05, 0.2           # the recipe bracket of f_env0 (scripts/atmosphere_recipes_data.py)
    assert low <= envelope_capture.MEDIAN_START <= high
    assert envelope_capture.median_start(10.0) == envelope_capture.median_start(envelope_capture.REFERENCE_MASS_EARTH) == envelope_capture.MEDIAN_START
    assert envelope_capture.median_start(1.0) < 0.005 < envelope_capture.median_start(2.0) < envelope_capture.median_start(3.0)


def test_only_planets_of_one_to_twenty_earth_masses_draw_a_start_and_the_start_is_capped():
    with draw.bound("envelope-capture:start"):
        assert envelope_capture.draw_start(0.5) == 0.0 and envelope_capture.draw_start(25.0) == 0.0
        starts = [envelope_capture.draw_start(8.0) for _ in range(400)]
    assert all(0.0 < value <= envelope_capture.MAX_START <= envelope_escape.MAX_ENVELOPE for value in starts)
    assert statistics.median(starts) == pytest.approx(envelope_capture.MEDIAN_START, rel=0.2)
    assert max(starts) == envelope_capture.MAX_START      # the tail of the log-normal is held at the top of the radius law's range


def test_a_planet_close_to_its_star_loses_what_one_far_from_it_keeps(sun):
    near = envelope_capture.kept(_body(5.0, 0.05, sun), 0.1)
    far = envelope_capture.kept(_body(5.0, 1.0, sun), 0.1)
    assert near is None
    fraction, radius_earth = far
    assert 0.05 < fraction < 0.12
    assert radius_earth > 2.2 > bulk_structure.core_radius_earth(5.0)       # a sub-Neptune, not the rocky 1.5
    assert envelope_capture.kept(_body(5.0, 1.0, sun), 0.0) is None


def test_the_orbit_is_read_at_its_mean_light(sun):
    circular = envelope_capture.kept(_body(8.0, 0.1, sun), 0.1)
    eccentric = envelope_capture.kept(_body(8.0, 0.1, sun, eccentricity=0.6), 0.1)
    assert (eccentric is None) or (circular is not None and eccentric[0] <= circular[0])      # more light, no more envelope


def _planet(star, distance, mass_earth, seed=0):
    return class_support.draw_planet(star, distance, seed, SystemConfig(), mass=mass_earth * EARTH, zone_override="c", moon_count=0)


def _kept_planet(star, distance, mass_earth):
    """A planet that was drawn with an envelope (the draw gives it when the planet is made), and the radius of its drawn core."""
    for seed in range(60):
        planet = _planet(star, distance, mass_earth, seed)
        if planet._captured_envelope:
            return planet, planet._bare_core["radius"]
    pytest.fail("no planet kept an envelope in 60 draws")


def test_a_planet_is_made_with_the_envelope_it_kept_its_radius_density_gravity_and_fractions(sun):
    planet, bare_radius = _kept_planet(sun, 1.5, 6.0)
    assert planet.envelope_mass_fraction == pytest.approx(planet._captured_envelope) and 0.0 < planet.envelope_mass_fraction < 0.5
    assert planet.radius > bare_radius * 1.15
    assert planet.volume == pytest.approx(4 / 3 * math.pi * planet.radius ** 3, rel=1e-12)
    assert planet.gravity == pytest.approx(constants.G * planet.mass / (planet.radius * 1000) ** 2 / constants.EARTH_GRAVITY, rel=1e-6)
    assert planet.density == pytest.approx(planet.mass / (4 / 3 * math.pi * (planet.radius * 1000) ** 3) / 1000, rel=1e-9)
    assert planet.density < planet.mass / (4 / 3 * math.pi * (bare_radius * 1000) ** 3) / 1000
    total = sum(getattr(planet, name) or 0.0 for name in (*envelope_share.SOLIDS, "envelope_mass_fraction"))
    assert total == pytest.approx(1.0, abs=1e-9)
    assert planet.planet_class is not None and planet.scale_height is not None
    assert 0.2 < planet.envelope_helium_fraction < 0.35 and 0.0 < planet.envelope_z_fraction < 0.05       # the star's gas
    expected = planet_type.TRANSITIONAL if planet.envelope_mass_fraction >= planet_type.TERRESTRIAL_MAX_ENVELOPE else planet_type.TERRESTRIAL
    assert planet_type.planet_type(planet) == expected


def _stripped_by_age(star):
    """A planet close in whose envelope the star's X-rays thin as it ages: drawn at 0.8 Gyr (before the loss starts), settled at 6."""
    star.age = 0.8
    for seed in range(60):
        planet = _planet(star, 0.12, 12.0, seed)
        if planet._captured_envelope:
            star.age = 6.0
            return planet
    pytest.fail("no planet kept an envelope at 0.8 Gyr in 60 draws")


def test_settling_at_the_final_age_takes_the_loss_again_and_makes_the_air_and_class_again(sun):
    planet = _stripped_by_age(sun)
    young = (planet.envelope_mass_fraction, planet.radius)
    assert planets.settle_envelope(planet) is True
    assert planet.envelope_mass_fraction < young[0] and planet.radius < young[1]
    assert planet.gravity == pytest.approx(constants.G * planet.mass / (planet.radius * 1000) ** 2 / constants.EARTH_GRAVITY, rel=1e-6)
    total = sum(getattr(planet, name) or 0.0 for name in (*envelope_share.SOLIDS, "envelope_mass_fraction"))
    assert total == pytest.approx(1.0, abs=1e-9)


def test_a_settled_planet_is_settled_for_good_and_a_second_pass_does_not_compound(sun):
    planet = _stripped_by_age(sun)
    planets.settle_envelope(planet)
    stored = (planet.radius, planet.gravity, planet.envelope_mass_fraction, planet.metal_mass_fraction, planet.composition)
    assert planets.settle_envelope(planet) is False
    assert (planet.radius, planet.gravity, planet.envelope_mass_fraction, planet.metal_mass_fraction, planet.composition) == stored


def test_the_star_stripping_all_of_it_gives_the_drawn_core_back(sun):
    planet = _stripped_by_age(sun)
    bare_radius = planet._bare_core["radius"]
    sun.age = 10.0
    sun.luminosity *= 30.0       # a star that has turned on: the loss takes it all
    assert planets.settle_envelope(planet) is True
    assert planet._captured_envelope == 0.0 and planet.radius == pytest.approx(bare_radius) and not planet.envelope_mass_fraction
    assert planet.envelope_helium_fraction is None and planet.envelope_z_fraction is None


def test_a_bare_core_and_a_planet_out_of_range_are_left_as_drawn(sun):
    close = _planet(sun, 0.04, 8.0)
    before = (close.radius, close.gravity, close.envelope_mass_fraction, close.planet_class)
    assert planets.settle_envelope(close) is False
    assert (close.radius, close.gravity, close.envelope_mass_fraction, close.planet_class) == before
    light = _planet(sun, 1.5, 0.5)
    assert planets.settle_envelope(light) is False and not light._captured_envelope


def test_settling_draws_nothing_from_the_run_stream(sun):
    planet = _stripped_by_age(sun)
    draw.set_run_seed(5)
    planets.settle_envelope(planet)
    after = draw.random()
    draw.set_run_seed(5)
    assert after == draw.random()


def test_the_envelope_is_not_counted_twice_as_a_hydrogen_air(sun):
    planet, _bare = _kept_planet(sun, 1.5, 6.0)
    assert envelope_share.air_share(planet) == 0.0


def test_a_planet_with_no_star_keeps_the_whole_envelope_it_started_with():
    """A rogue planet: no light, no loss (`envelope_escape`), the radius read at the least light the law is fitted for."""
    star = SimpleNamespace(luminosity=0.0, mass=1.989e30, age=4.5)
    body = SimpleNamespace(star=star, mass=8.0 * EARTH, distance=1.0e5, eccentricity=0.0)
    fraction, radius_earth = envelope_capture.kept(body, 0.1)
    assert fraction == pytest.approx(0.1) and radius_earth > 2.5


@pytest.mark.slow
def test_the_generated_planets_show_the_radius_valley_a_system_sweep_finds():
    """The sweep of GEN.1323: inside 100 days, about half the planets of 1 to 4 Earth radii are above 2.0 (Fulton et al. 2017) and the
    planets of 1 to 2 Earth masses are mostly rocky."""
    draw.set_run_seed(20261107)
    inner, big, light, light_rocky = 0, 0, 0, 0
    for _ in range(24):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        from planetgen.generation.system import StarSystem
        for planet in StarSystem(system_config=cfg).all_planets:
            if getattr(planet, "body_type", None) != "t" or planet.is_moon:
                continue
            radius_earth = planet.radius / constants.EARTH_RADIUS_KM
            mass_earth = planet.mass / EARTH
            if planet.period * 365.25 < 100.0 and 1.0 <= radius_earth <= 4.0:
                inner += 1
                big += radius_earth > 2.0
            if 1.0 <= mass_earth < 1.5:
                light += 1
                light_rocky += (planet.envelope_mass_fraction or 0.0) < 0.01
    assert inner >= 8 and 0.2 <= big / inner <= 0.75
    assert light >= 10 and light_rocky / light >= 0.6
