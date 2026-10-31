"""
Statistical checks on what generation produces over many systems (each
system draws fresh random numbers, so every system is its own seed):

- GEN.37: the zone mix of planets (hot, ecosphere, cold).
- GEN.34: gas and ice giants have real bulk densities, and super-Jupiters
  occur.
- GEN.35: rocky planets get more than Class D moons.
- GEN.36 / GEN.25: a moon made again after its planet moves is a rocky body
  its planet can hold, and every system passes the validator with no
  moon-size problems.
- GEN.45: the rogue planet split of terrestrial and gas giants.
"""
import collections
import statistics

import pytest

from planetgen.physics import constants as pc
from planetgen.generation import validation
from planetgen.physics import planet_classes, planets as planetPhysics
from planetgen import tuning as prog_c
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem

N_SYSTEMS = 250   # every body solves its own air (GEN.411), so a system takes seconds; the bands below hold at this size


@pytest.fixture(scope="module")
def systems():
    return [StarSystem(system_config=SystemConfig()) for _ in range(N_SYSTEMS)]


def _planets(systems):
    for system in systems:
        for _star, planets in validation.star_lists(system):
            for body in planets:
                if body.body_type != 'a':
                    yield body


def test_zone_mix_is_spread(systems):
    """About 30% hot, 12% ecosphere and 50 to 55% cold (it was 96% cold; the ecosphere was 6% before GEN.146 widened it to
    the Kopparapu edges, 0.95 to 1.68 AU for the Sun; the mass law of GEN.646 puts a few more planets close in)."""
    zones = collections.Counter(planet.zone for planet in _planets(systems))
    total = sum(zones.values())
    assert total > 3 * N_SYSTEMS
    assert 0.20 < zones['h'] / total < 0.42
    assert 0.06 < zones['e'] / total < 0.20
    assert 0.40 < zones['c'] / total < 0.75


def test_giants_have_real_densities_and_super_jupiters_occur(systems):
    giants = [p for p in _planets(systems) if p.body_type == 'g']
    densities = [p.density for p in giants]
    assert len(giants) > N_SYSTEMS // 6
    # Real gas and ice giants run 0.69-1.64 g/cm^3 (Saturn to Neptune);
    # the median was 0.25 before GEN.34.
    assert 0.69 <= statistics.median(densities) <= 1.64
    # GEN.460: a light giant of the radius tables, or a hot one swollen by its starlight, can be very puffy (measured
    # super-puffs run down to 0.03 g/cm^3, Masuda 2014).
    assert all(0.05 < d < 40 for d in densities)
    super_jupiters = [p for p in giants if p.mass > pc.JUPITER_MASS_TO_KG]
    assert len(super_jupiters) / len(giants) > 0.05


def test_rocky_planets_get_more_than_class_d_moons(systems):
    classes = collections.Counter(moon.planet_class for planet in _planets(systems)
                                  if planet.body_type == 't' for moon in planet.moons)
    assert classes['C'] > 0 and classes['D'] > 0


def test_every_moon_is_a_rocky_body_that_fits_its_planet(systems):
    for planet in _planets(systems):
        max_radius, max_mass = planetPhysics.moon_size_limits(planet)
        for moon in planet.moons:
            assert moon.body_type == 't'
            assert moon.radius <= max_radius * (1 + 1e-9)
            assert moon.mass <= max_mass * (1 + 1e-9)


def test_generated_systems_have_no_moon_size_problems(systems):
    """GEN.25's goal: no moon too large for its planet in the sample."""
    problems = [problem for system in systems for problem in validation.check_star_system(system)
                if "too large" in problem.message]
    assert problems == []


def test_moved_planet_makes_its_moons_again_as_moons():
    """GEN.36/GEN.25: a giant pushed from the ecosphere into the cold zone has its moons made again for the new place, and
    every one is a rocky body of its planet's zone, of the class its conditions give, that its planet can hold."""
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    star = Star(cfg)
    inner, outer = star.habitable_zone
    giant_kg = 300 * pc.EARTH_MASS_TO_KG
    moons_seen = 0
    for _ in range(25):
        planet = Planet(cfg, star, star.habitable_zone, (inner + outer) / 2, mass=giant_kg, zone_override='e')
        planet.distance = outer * 3
        validation.reconcile_moved_planet(planet)
        max_radius, max_mass = planetPhysics.moon_size_limits(planet)
        for moon in planet.moons:
            moons_seen += 1
            assert moon.zone == planet.zone == 'c' and moon.body_type == 't'
            assert moon.planet_class == planet_classes.classify(moon)
            assert moon.radius <= max_radius * (1 + 1e-9)
            assert moon.mass <= max_mass * (1 + 1e-9)
        assert validation.check_lunar_system(planet) == []
    assert moons_seen > 0


def test_rogue_planet_split_follows_the_mass_function():
    """GEN.45, GEN.1133: about 95% terrestrial and 5% giants (the line is the 10 Earth masses a core needs to hold a
    gas envelope), as dN/dlogM ~ M^-0.65 from 0.1 Earth masses to 13 Jupiter masses gives."""
    slope = prog_c.ROGUE_PLANET_MASS_FUNCTION_SLOPE
    low = min(b[0] for b in prog_c.ROGUE_PLANET_MASS_BINS.values())
    high = max(b[1] for b in prog_c.ROGUE_PLANET_MASS_BINS.values())
    threshold = prog_c.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER * pc.JUPITER_MASS_TO_KG / pc.EARTH_MASS_TO_KG
    expected = (threshold ** -slope - high ** -slope) / (low ** -slope - high ** -slope)
    assert expected == pytest.approx(0.049, abs=0.002)
    assert sum(rate for _lo, _hi, rate in prog_c.ROGUE_PLANET_MASS_BINS.values()) == pytest.approx(5.8)

    cfg = SystemConfig()
    n = 20000
    gas = sum(1 for _ in range(n) if RoguePlanet(cfg).planet_type == 'g')
    assert abs(gas / n - expected) < 0.008
