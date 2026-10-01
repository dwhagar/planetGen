"""
Tests for the central validate module, `stellarObjects/validation.py`
(TODO ADM.5): its checks find nothing wrong with generated systems and
find each kind of fault, and the stabilize pass fixes what an edit broke.
"""
import pytest

from stellarObjects import physical_constants, planetPhysics, validation
from stellarObjects.config import SystemConfig
from stellarObjects.planetData import Planet
from stellarObjects.starData import Star
from stellarObjects.systemData import StarSystem


def make_system(star_type="G2V", **overrides):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.MOONS = False
    for key, value in overrides.items():
        setattr(cfg, key, value)
    return StarSystem(cfg)


@pytest.fixture
def star():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    return Star(cfg)


def planet_at(star, distance_au, planet_class, moon_count=0):
    return Planet(star.system_config, star, star.habitable_zone, distance_au,
                  planet_class=planet_class, moon_count=moon_count)


@pytest.mark.parametrize("star_type", ["G2V", "M3V", "K1V", "F5V", "A0V"])
def test_generated_systems_validate(star_type):
    for _ in range(5):
        system = make_system(star_type, PLANETS=True)
        assert validation.check_star_system(system) == []


def test_generated_wide_binaries_validate():
    for _ in range(5):
        system = make_system("K1V", BINARY_SYSTEM=True, WIDE_BINARY=True)
        assert validation.check_star_system(system) == []


def test_check_planet_reports_a_class_out_of_its_zone(star):
    planet = planet_at(star, star.habitable_zone[1] * 3, "J")
    assert validation.check_planet(planet) == []
    planet.planet_class = "M"  # Earth-like, far out in the cold zone
    messages = [p.message for p in validation.check_planet(planet)]
    assert any("can't exist in zone 'c'" in m for m in messages)


def test_check_planet_reports_a_stale_zone(star):
    planet = planet_at(star, star.habitable_zone[1] * 3, "J")
    planet.zone = "h"
    assert any("its zone is 'h'" in p.message for p in validation.check_planet(planet))


def test_check_orbit_spacing_and_space_orbits(star):
    inner = planet_at(star, 5.0, "J")
    outer = planet_at(star, 5.01, "J")
    planets = [inner, outer]
    assert any("too close" in p.message for p in validation.check_orbit_spacing(planets))
    moved = validation.space_orbits(planets)
    assert moved == [outer]
    assert outer.distance > 5.01
    assert validation.check_orbit_spacing(planets) == []


def test_space_orbits_keeps_a_pinned_class(star):
    inner = planet_at(star, star.habitable_zone[0] * 1.01, "J")
    # An ecosphere-only class right next to a gas giant: the push carries
    # it out of the ecosphere.
    outer = planet_at(star, star.habitable_zone[0] * 1.02, "M")
    outer.distance = star.habitable_zone[1] * 0.999
    inner.distance = star.habitable_zone[1] * 0.99
    validation.space_orbits([inner, outer], pinned=[outer])
    assert outer.planet_class == "M"
    assert outer.zone == "c"


def test_trim_to_orbit_ceiling_returns_what_it_removed(star):
    planets = [planet_at(star, 6.0, "J"), planet_at(star, 50.0, "J")]
    removed = validation.trim_to_orbit_ceiling(planets, 10.0)
    assert [p.distance for p in removed] == [50.0]
    assert len(planets) == 1


def test_stabilize_lunar_system_respaces_crowded_moons(star):
    planet = planet_at(star, 5.0, "J", moon_count=3)
    assert len(planet.moons) >= 2
    planet.moons[1].distance = planet.moons[0].distance  # stack two moons
    assert any("too close" in p.message for p in validation.check_lunar_system(planet))
    moved = validation.stabilize_lunar_system(planet)
    assert moved
    assert validation.check_lunar_system(planet) == []


def test_check_lunar_system_reports_a_moon_outside_the_stable_range(star):
    planet = planet_at(star, 5.0, "J", moon_count=1)
    _low, high = planetPhysics.moon_orbit_bounds_km(planet)
    planet.moons[0].distance = high * 2 / physical_constants.AU_TO_KM
    assert any("outside" in p.message for p in validation.check_lunar_system(planet))


def test_stabilize_star_system_fixes_a_crowded_system():
    system = make_system("G2V", PLANETS=True, MAX_PLANETS=True)
    planets = [p for p in system.planets if p.body_type != 'a']
    if len(planets) < 2:
        pytest.skip("needs two planets")
    planets[1].distance = planets[0].distance * 1.0001
    report = validation.stabilize_star_system(system)
    assert report.moved
    assert report.removed == []
    assert [p for p in report.problems if "too close" in p.message] == []


def test_stabilize_reports_rather_than_removes_by_default():
    system = make_system("G2V", PLANETS=True)
    if not system.planets:
        pytest.skip("needs a body")
    outer = system.planets[-1]
    ceiling = validation.orbit_ceiling_au(system.star)
    if outer.body_type == 'a':
        span = outer.upper_limit - outer.lower_limit
        outer.lower_limit = ceiling * 2
        outer.upper_limit = outer.lower_limit + span
        outer.distance = outer.lower_limit + span / 2
    else:
        outer.distance = ceiling * 2
    report = validation.stabilize_star_system(system)
    assert report.removed == []
    assert any("beyond the farthest stable orbit" in p.message for p in report.problems)

    report = validation.stabilize_star_system(system, allow_removal=True)
    assert report.removed
    assert not any("beyond" in p.message for p in report.problems)


def test_system_methods_still_delegate():
    system = make_system("G2V")
    assert system._orbit_ceiling_au(system.star) == validation.orbit_ceiling_au(system.star)


def test_stabilizing_an_untouched_system_moves_nothing():
    # A body generation put exactly on its spacing limit must not be
    # pushed out again by a float rounding error: an edit elsewhere in
    # the system would otherwise move it (and a belt by a whole Hill
    # sphere).
    for _ in range(40):
        system = make_system("G2V", PLANETS=True, ASTEROID_BELT=True)
        if validation.check_star_system(system):
            continue
        report = validation.stabilize_star_system(system)
        assert report.moved == [] and report.reclassified == []
