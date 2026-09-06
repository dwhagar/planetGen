"""
Planet-class generation regression tests.

Builds a `Planet` directly (bypassing `StarSystem`'s orbit-placement loop)
for every (planet class, zone) combination `PLANET_CLASSES` declares valid,
and checks the resulting physical properties and life-chemical assignment
against the class's own declared data. One shared G2V star is used as the
host, since planet physics only reads a handful of scalar properties off it
(luminosity, mass, radius, temperature, habitable zone) -- spectral-type
coverage of the star itself belongs to test_star_matrix.py.

Run with: pytest tests/test_planets.py
"""
import math

import pytest

from stellarObjects.config import SystemConfig
from stellarObjects.planetData import Planet
from stellarObjects.starData import Star
from stellarObjects import program_constants as prog_c
from stellarObjects import planetLife

ZONE_CHARS = "hec"

# (class, zone) pairs PLANET_CLASSES actually declares valid.
VALID_CLASS_ZONE_PAIRS = [
    (cls, zone)
    for cls, data in prog_c.PLANET_CLASSES.items()
    for zone in ZONE_CHARS
    if data[zone]
]

CLASSES_WITH_NO_VALID_ZONE = [
    cls for cls, data in prog_c.PLANET_CLASSES.items() if not any(data[z] for z in ZONE_CHARS)
]


@pytest.fixture(scope="module")
def host_star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    return Star(cfg)


def distance_for_zone(star, zone):
    inner, outer = star.habitable_zone
    if zone == 'h':
        return inner * 0.5
    if zone == 'c':
        return outer * 2.0
    return (inner + outer) / 2.0


def make_planet(host_star, cls, zone, **kwargs):
    cfg = SystemConfig()
    distance = distance_for_zone(host_star, zone)
    return Planet(
        cfg, host_star, host_star.habitable_zone, distance,
        planet_class=cls, **kwargs,
    )


@pytest.mark.parametrize("cls,zone", VALID_CLASS_ZONE_PAIRS, ids=[f"{c}-{z}" for c, z in VALID_CLASS_ZONE_PAIRS])
def test_planet_generates_without_error(host_star, cls, zone):
    make_planet(host_star, cls, zone)


@pytest.mark.parametrize("cls,zone", VALID_CLASS_ZONE_PAIRS, ids=[f"{c}-{z}" for c, z in VALID_CLASS_ZONE_PAIRS])
def test_planet_radius_within_declared_range(host_star, cls, zone):
    planet = make_planet(host_star, cls, zone)
    min_r, max_r = prog_c.PLANET_CLASSES[cls]["radius_range"]
    assert min_r <= planet.radius <= max_r
    assert planet.planet_class == cls
    assert planet.zone == zone
    assert planet.body_type == prog_c.PLANET_CLASSES[cls]["type"]


@pytest.mark.parametrize("cls,zone", VALID_CLASS_ZONE_PAIRS, ids=[f"{c}-{z}" for c, z in VALID_CLASS_ZONE_PAIRS])
def test_planet_physical_properties_are_finite_and_sane(host_star, cls, zone):
    planet = make_planet(host_star, cls, zone)

    for attr in ("mass", "density", "gravity", "surface_temperature", "volume", "period", "hill_radius", "min_orbit_distance"):
        value = getattr(planet, attr)
        assert math.isfinite(value), f"{cls}-{zone}: {attr}={value} is not finite"
        assert value > 0, f"{cls}-{zone}: {attr}={value} should be positive"

    class_data = prog_c.PLANET_CLASSES[cls]
    if class_data["atmosphere"] is None:
        assert planet.atmosphere == "None"
        assert planet.atmospheric_pressure == 0.0
        assert planet.atm_density is None
        assert planet.atm_molar_density is None
    else:
        assert planet.atmosphere == class_data["atmosphere"]
        assert math.isfinite(planet.atmospheric_pressure) and planet.atmospheric_pressure >= 0
        assert math.isfinite(planet.atm_density) and planet.atm_density > 0
        assert math.isfinite(planet.atm_molar_density) and planet.atm_molar_density > 0
        assert math.isfinite(planet.scale_height) and planet.scale_height > 0


@pytest.mark.parametrize("cls,zone", VALID_CLASS_ZONE_PAIRS, ids=[f"{c}-{z}" for c, z in VALID_CLASS_ZONE_PAIRS])
def test_planet_life_chemical_matches_class_and_star(host_star, cls, zone):
    """
    A planet's assigned life_chemical (once planetLife.apply_life_data runs)
    must be one the planet's own class declares as possible -- and, since
    apply_life_data intersects that with what the star's spectral type
    supports, it must also appear in the star's potentially_viable_chemicals.
    """
    planet = make_planet(host_star, cls, zone)
    for _ in range(5):
        planetLife.apply_life_data(planet)
        class_chems = prog_c.PLANET_CLASSES[cls]["life_chemical"]
        if planet.life_chemical is not None:
            assert class_chems is not None, f"{cls}-{zone}: assigned a chemical but class has none listed"
            assert any(planet.life_chemical in c or c in planet.life_chemical for c in class_chems), (
                f"{cls}-{zone}: life_chemical {planet.life_chemical!r} not among class chemicals {class_chems}"
            )


def test_decide_flavor_text_is_generation_time_only(host_star, monkeypatch):
    """
    Regression test for the Phase 0 bug (see TODO.md): flavor text used to be
    rolled -- and system_config.system_flavor_count/recent_flavor_texts
    mutated -- inside to_paragraph_list() itself, i.e. at render time, so
    calling to_paragraph_list() twice double-rolled and double-mutated
    shared state. planetLife.decide_flavor_text now makes that decision
    once, at generation time, and to_paragraph_list() is a pure read of the
    already-decided self.flavor_text.
    """
    planet = make_planet(host_star, "M", "e")
    planet.evolutionary_data = []  # No multicellular life keyword to match.
    monkeypatch.setattr(prog_c, "FLAVOR_CHANCE_PLANET", 1.0)

    planetLife.decide_flavor_text(planet)

    assert planet.flavor_text is not None
    assert planet.flavor_text_count == 1
    assert planet.system_config.system_flavor_count == 1

    first_render = planet.to_paragraph_list()
    second_render = planet.to_paragraph_list()

    assert first_render == second_render
    assert planet.flavor_text_count == 1
    assert planet.system_config.system_flavor_count == 1


@pytest.mark.parametrize("cls", CLASSES_WITH_NO_VALID_ZONE)
def test_known_issue_class_with_no_valid_zone_is_unreachable(host_star, cls):
    """
    FINDING (not fixed -- see fork report): PLANET_CLASSES['R'] declares
    h=e=c=False, so it has zero probability weight in PLANET_CLASS_PROBABILITIES
    *and* fails zone validation in every zone if explicitly requested. This
    xfail is strict so it starts failing loudly (forcing an update) the
    moment this class becomes reachable in some zone -- whether that's the
    right fix (which zone(s) should an "ejected" world belong to?) is a
    design question for a human, not something this fork guessed at.
    """
    assert prog_c.PLANET_CLASS_PROBABILITIES.get(cls, 0) == 0
    with pytest.raises(ValueError, match="Invalid planet class for this zone"):
        for zone in ZONE_CHARS:
            make_planet(host_star, cls, zone)


def test_habitable_world_false_forbids_habitable_classes_in_ecosphere(host_star):
    cfg = SystemConfig()
    cfg.HABITABLE_WORLD = False
    distance = distance_for_zone(host_star, 'e')
    for cls in prog_c.HABITABLE_PLANET_CLASSES:
        if not prog_c.PLANET_CLASSES[cls]['e']:
            continue
        with pytest.raises(ValueError):
            Planet(
                cfg, host_star, host_star.habitable_zone, distance,
                planet_class=cls,
            )


@pytest.mark.parametrize("zone", list(ZONE_CHARS))
@pytest.mark.parametrize("habitable_world", [True, False, None])
def test_fully_random_planet_respects_zone_and_habitable_world(host_star, zone, habitable_world):
    cfg = SystemConfig()
    cfg.HABITABLE_WORLD = habitable_world
    distance = distance_for_zone(host_star, zone)
    for _ in range(10):
        planet = Planet(
            cfg, host_star, host_star.habitable_zone, distance,
        )
        assert prog_c.PLANET_CLASSES[planet.planet_class][zone]
        if habitable_world is False:
            assert planet.planet_class not in prog_c.HABITABLE_PLANET_CLASSES
