"""
Planet generation regression tests.

Builds a `Planet` directly (bypassing `StarSystem`'s orbit-placement loop) at a distance in each zone, from several seeded
streams, and checks the resulting physical properties and life-chemical assignment. The class is derived from what the body
comes out as (GEN.646), so every check reads it from the planet. One shared G2V star is used as the host, since planet physics
only reads a handful of scalar properties off it (luminosity, mass, radius, temperature, habitable zone) -- spectral-type
coverage of the star itself belongs to test_star_matrix.py.

Run with: pytest tests/test_planets.py
"""
import math

import pytest

from planetgen.describe import build
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen import tuning as prog_c
from planetgen.physics import constants, gases, planet_classes, surface_air
from planetgen.generation import habitable_request, life as planetLife
from planetgen.util import draw
from tests import class_support

ZONE_CHARS = "hec"
SEEDS = range(4)
ZONE_SEEDS = [(zone, seed) for zone in ZONE_CHARS for seed in SEEDS]
ZONE_SEED_IDS = [f"{zone}-{seed}" for zone, seed in ZONE_SEEDS]

# GEN.830: the mass share of a body that its air can hold beyond its water: the carbon, nitrogen and sulfur a mantle degasses, bounded by
# the 9 percent that CI chondrites carry of them in all (Lodders); the degassed share of that is far smaller.
VOLATILE_ALLOWANCE = 0.09

# TEST.485: the least pressure a stored air has. A thin regime caps the air at `surface_air.VACUUM_PRESSURE_PA` (1 Pa) and the solve then drops every
# gas whose partial is under `STORE_MIN_FRACTION` of the total, so the stored sum sits below the cap by that share of each gas dropped, at most
# all the gases the registry knows (a draw of 1 Pa stored as 0.99999999993 Pa).
THIN_AIR_FLOOR_PA = surface_air.VACUUM_PRESSURE_PA * (1.0 - surface_air.STORE_MIN_FRACTION * len(gases.GASES))

@pytest.fixture(scope="module")
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    with draw.bound("test-planets-host-star"):      # TEST.497: a star drawn from the run's entropy gave a different air to the checks each run
        return Star(cfg)


def distance_for_zone(star, zone):
    inner, outer = star.habitable_zone
    if zone == 'h':
        return inner * 0.5
    if zone == 'c':
        return outer * 2.0
    return (inner + outer) / 2.0


def make_planet(host_star, zone, seed=0, **kwargs):
    return class_support.draw_planet(host_star, distance_for_zone(host_star, zone), seed, **kwargs)


@pytest.mark.parametrize("zone,seed", ZONE_SEEDS, ids=ZONE_SEED_IDS)
def test_planet_is_in_its_zone_with_a_known_class_and_the_type_of_that_class(host_star, zone, seed):
    planet = make_planet(host_star, zone, seed)
    assert planet.radius > 0
    assert planet.planet_class in prog_c.PLANET_CLASSES
    assert planet.zone == zone
    shown = prog_c.PLANET_CLASSES[planet.planet_class]
    assert planet.body_type in (shown["type"], shown.get("also_type"))
    assert planet.planet_class == planet_classes.classify(planet)


def test_the_thin_air_floor_allows_the_trace_trim_and_no_more():
    assert THIN_AIR_FLOOR_PA < 0.99999999993 < surface_air.VACUUM_PRESSURE_PA
    assert THIN_AIR_FLOOR_PA > surface_air.VACUUM_PRESSURE_PA * (1.0 - 1.0e-6)


@pytest.mark.parametrize("zone,seed", ZONE_SEEDS, ids=ZONE_SEED_IDS)
def test_planet_physical_properties_are_finite_and_sane(host_star, zone, seed):
    planet = make_planet(host_star, zone, seed)
    cls = planet.planet_class

    for attr in ("mass", "density", "gravity", "surface_temperature", "volume", "period", "hill_radius", "min_orbit_distance"):
        value = getattr(planet, attr)
        assert math.isfinite(value), f"{cls}-{zone}: {attr}={value} is not finite"
        assert value > 0, f"{cls}-{zone}: {attr}={value} should be positive"

    if planet.atmospheric_pressure == 0.0:
        # GEN.411: the air is solved from the body's own interior, sinks and temperature, so a body can be airless where nothing
        # stays (a hot small body); the stored fields then agree that there is none.
        assert planet.atmosphere == "None"
        assert planet.atm_density is None
        assert planet.atm_molar_density is None
    else:
        # GEN.85: the text is written from the stored gases.
        assert planet.atmosphere.startswith("a ") and " mix of " in planet.atmosphere
        # A thin, Mars-like air up to a thick, Venus-like one; a unit-conversion bug would show as near-zero pressure.
        # GEN.830: the weight of the air on the surface (g * m / area) is its mass, which cannot exceed the water the body has plus the
        # degassed share of its other volatiles. This replaces the flat ceiling: a sealed ocean world that crosses the critical point of
        # water holds all its surface water as steam (several GPa), which is its own inventory and no more.
        # The floor is the 1 Pa a thin regime caps the air at (`surface_air.VACUUM_PRESSURE_PA`); the solve then drops every gas whose partial is
        # under `STORE_MIN_FRACTION` of the total, so the stored sum can sit below the cap by that share of each gas dropped (TEST.485:
        # 0.99999999993 Pa), at most all the gases the registry knows.
        assert math.isfinite(planet.atmospheric_pressure) and THIN_AIR_FLOOR_PA <= planet.atmospheric_pressure
        air_kg = planet.atmospheric_pressure * 4.0 * math.pi * (planet.radius * 1000.0) ** 2 / (planet.gravity * constants.EARTH_GRAVITY)
        assert air_kg <= planet.mass * ((planet.water_mass_fraction or 0.0) + VOLATILE_ALLOWANCE)
        assert math.isfinite(planet.atm_density) and planet.atm_density > 0
        assert math.isfinite(planet.atm_molar_density) and planet.atm_molar_density > 0
        assert math.isfinite(planet.scale_height) and planet.scale_height > 0


@pytest.mark.parametrize("zone,seed", ZONE_SEEDS, ids=ZONE_SEED_IDS)
def test_planet_life_class_is_one_its_conditions_support(host_star, zone, seed):
    """
    A planet's life class (once planetLife.apply_life_data runs) is one its own stored conditions support: it
    is among the classes the body scores 1 or more for (GEN.889). A gas giant carries none.
    """
    planet = make_planet(host_star, zone, seed)
    cls = planet.planet_class
    for _ in range(5):
        planetLife.apply_life_data(planet)
        if planet.body_type == "g":
            assert planet.life_class is None
        elif planet.life_class_id is not None:
            assert planet.life_class_scores.get(planet.life_class_id, 0) > 0, f"{cls}-{zone}: {planet.life_class_id!r} is not supported"


def test_the_class_that_lives_on_a_world_is_stored_among_the_classes_it_supports_even_when_its_life_used_up_its_fuel():
    """TEST.486: a cold world whose life drew down what it lived on is scored 0 for its class once the air is built; the class still lives there
    and is still stored among the supported ones, at the score it was chosen at. The host star is seeded, so the same stars come every run."""
    kept = 0
    for trial in range(40):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        with draw.bound(trial):
            star = Star(cfg)
        for seed in range(4):
            planet = make_planet(star, "c", seed)
            planetLife.apply_life_data(planet)
            if planet.body_type != "g" and planet.life_class_id is not None:
                kept += 1
                assert planet.life_class_scores.get(planet.life_class_id, 0) > 0, (trial, seed, planet.life_class_id, planet.life_class_scores)
    assert kept


def test_decide_flavor_text_is_generation_time_only(host_star, monkeypatch):
    """
    Regression test for the Phase 0 bug (see TODO.md): flavor text used to be
    rolled -- and system_config.system_flavor_count/recent_flavor_texts
    mutated -- at render time, so describing a planet twice double-rolled and
    double-mutated shared state. planetLife.decide_flavor_text now makes that decision
    once, at generation time, and the Document is a pure read of the
    already-decided self.flavor_text.
    """
    planet = make_planet(host_star, "e")
    planet.life_class, planet.life_class_variant, planet.life_stage = 3, "b", None  # Life, but not multicellular.
    monkeypatch.setattr(prog_c, "FLAVOR_CHANCE_PLANET", 1.0)

    planetLife.decide_flavor_text(planet)

    assert planet.flavor_text is not None
    assert planet.flavor_text_count == 1
    assert planet.system_config.system_flavor_count == 1

    first_render = build._planet_section(planet)
    second_render = build._planet_section(planet)

    assert first_render == second_render
    assert planet.flavor_text_count == 1
    assert planet.system_config.system_flavor_count == 1


def test_the_flavour_line_is_stored_as_a_key_and_read_back_from_it(host_star, monkeypatch):
    """DB.51: the line is kept as `category:index` into the `tuning` list it was drawn from."""
    monkeypatch.setattr(prog_c, "FLAVOR_CHANCE_PLANET", 1.0)
    planet = make_planet(host_star, "e")
    planet.life_class, planet.life_class_variant, planet.life_stage = 3, "b", "multicellularity"
    planetLife.decide_flavor_text(planet)
    category, index = planet.flavor_key.split(":")
    assert category == "habitable"
    assert planet.flavor_text == prog_c.HABITABLE_FLAVOR[int(index)]
    assert planetLife.flavor_key("planet", prog_c.PLANET_FLAVOR[3]) == "planet:3"
    assert planetLife.flavor_line(None) is None


def test_habitable_world_false_forbids_habitable_classes_in_ecosphere(host_star):
    cfg = SystemConfig()
    cfg.HABITABLE_WORLD = False
    distance = distance_for_zone(host_star, 'e')
    for seed in range(15):
        planet = class_support.draw_planet(host_star, distance, seed, cfg)
        assert not habitable_request.meets(planet)         # GEN.916: the bar is on the conditions of a habitable world, not a class letter


@pytest.mark.parametrize("zone", list(ZONE_CHARS))
@pytest.mark.parametrize("habitable_world", [True, False, None])
def test_fully_random_planet_respects_zone_and_habitable_world(host_star, zone, habitable_world):
    cfg = SystemConfig()
    cfg.HABITABLE_WORLD = habitable_world
    distance = distance_for_zone(host_star, zone)
    for seed in range(6):
        planet = class_support.draw_planet(host_star, distance, seed, cfg)
        assert planet.zone == zone
        if habitable_world is False and zone == 'e':  # the bar is on the ecosphere; a hot or cold body's conditions are its own
            assert not habitable_request.meets(planet)


def test_a_moon_keeps_its_orbit_around_its_parent(host_star):
    """
    A planet may be redrawn to a new distance only for the star's own zones -- a moon's `distance` is its orbit around its
    *parent planet*, an entirely different (much smaller) AU-scale value, so building a moon must not move it.
    """
    cfg = SystemConfig()
    parent = Planet(cfg, host_star, host_star.habitable_zone, distance_for_zone(host_star, 'e'), moon_count=0)
    moon_distance = 384400 / 149_597_870.7  # real Moon's distance, in AU

    moon = Planet(
        cfg, host_star, host_star.habitable_zone, moon_distance,
        mass=7.3e22, zone_override=parent.zone, distance_override=parent.distance, is_moon=True,
        primary_mass_kg=parent.mass,
    )

    assert moon.distance == pytest.approx(moon_distance)


def test_no_planet_class_names_a_life_class():
    """GEN.647: life follows a body's stored conditions, not its class letter."""
    assert not any("life_class" in data or "life_chemical" in data for data in prog_c.PLANET_CLASSES.values())
    assert not hasattr(prog_c, "PLANET_CLASS_MAX_LIFE_STAGE")
    assert not hasattr(prog_c, "LIFE_CHEMICALS")


def test_a_world_no_class_supports_carries_no_class(host_star):
    """A body that supports no life class carries none whatever its air (GEN.206)."""
    seen_none = False
    for seed in range(12):
        planet = make_planet(host_star, "e", seed)
        planetLife.apply_life_data(planet)
        if not planet.life_class_scores:  # GEN.411: the class no longer fixes the air, so not every body supports a class
            seen_none = True
            assert planet.life_class is None and planet.life_stage is None
    if not seen_none:
        pytest.skip("no draw in zone e supported no life class")
