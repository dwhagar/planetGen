"""
Statistical checks on what generation produces over many systems (each
system draws fresh random numbers, so every system is its own seed):

- GEN.37: the zone mix of planets (hot, ecosphere, cold).
- GEN.34: gas and ice giants have real bulk densities, and super-Jupiters
  occur.
- GEN.35: rocky planets get more than Class D moons.
- GEN.36 / GEN.25: a moon regenerated after its planet moves only gets a
  class a moon of that planet may have, and every system passes the
  validator with no moon-size problems.
- GEN.45: the rogue planet split of terrestrial and gas giants.
"""
import collections
import math
import statistics

import pytest

from stellarObjects import physical_constants as pc
from stellarObjects import planetPhysics, validation
from stellarObjects import program_constants as prog_c
from stellarObjects.config import SystemConfig
from stellarObjects.planetData import Planet
from stellarObjects.roguePlanetData import RoguePlanet
from stellarObjects.starData import Star
from stellarObjects.systemData import StarSystem

N_SYSTEMS = 1000


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
    """About 30% hot, 6% ecosphere and 64% cold (it was 96% cold)."""
    zones = collections.Counter(planet.zone for planet in _planets(systems))
    total = sum(zones.values())
    assert total > 3000
    assert 0.20 < zones['h'] / total < 0.42
    assert 0.03 < zones['e'] / total < 0.12
    assert 0.50 < zones['c'] / total < 0.75


def test_giants_have_real_densities_and_super_jupiters_occur(systems):
    giants = [p for p in _planets(systems) if p.body_type == 'g']
    densities = [p.density for p in giants]
    assert len(giants) > 500
    # Real gas and ice giants run 0.69-1.64 g/cm^3 (Saturn to Neptune);
    # the median was 0.25 before GEN.34.
    assert 0.69 <= statistics.median(densities) <= 1.64
    assert all(0.2 < d < 40 for d in densities)
    super_jupiters = [p for p in giants if p.mass > pc.JUPITER_MASS_TO_KG]
    assert len(super_jupiters) / len(giants) > 0.05


def test_rocky_planets_get_more_than_class_d_moons(systems):
    classes = collections.Counter(moon.planet_class for planet in _planets(systems)
                                  if planet.body_type == 't' for moon in planet.moons)
    assert classes['C'] > 0 and classes['D'] > 0


def test_every_moon_fits_its_planet(systems):
    for planet in _planets(systems):
        max_radius, max_mass = planetPhysics.moon_size_limits(planet)
        for moon in planet.moons:
            data = prog_c.PLANET_CLASSES[moon.planet_class]
            assert data["type"] == 't' and moon.planet_class not in prog_c.MOON_BLACKLIST
            assert moon.radius <= max_radius * (1 + 1e-9)
            assert moon.mass <= max_mass * (1 + 1e-9)


def test_generated_systems_have_no_moon_size_problems(systems):
    """GEN.25's goal: no moon too large for its planet in 1,000 systems."""
    problems = [problem for system in systems for problem in validation.check_star_system(system)
                if "too large" in problem.message]
    assert problems == []


def test_moon_class_options_size_rules():
    """An Earth-sized rocky planet can hold a Class C moon (up to the
    Moon's scale) as well as Class D, never a gas giant or a blacklisted
    class, and each ceiling keeps the moon under a tenth of its mass."""
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    star = Star(cfg)
    earth = Planet(cfg, star, star.habitable_zone, 2.0, planet_class="C", radius=6371, moon_count=0)
    options = planetPhysics.moon_class_options(earth, earth.zone)
    assert {'C', 'D'} <= set(options)
    assert all(prog_c.PLANET_CLASSES[c]["type"] == 't' and c not in prog_c.MOON_BLACKLIST for c in options)
    max_radius, max_mass = planetPhysics.moon_size_limits(earth)
    densest = pc.PLANET_DENSITY['t'][1] * 1000
    for ceiling in options.values():
        assert ceiling <= max_radius
        assert (4 / 3) * math.pi * (ceiling * 1000) ** 3 * densest <= max_mass * (1 + 1e-9)


def test_moved_planet_regenerates_moons_as_moons():
    """GEN.36/GEN.25: a giant pushed from the ecosphere into the cold zone
    regenerates its habitable-class moons, and every regenerated moon is a
    class and size a moon of it may have."""
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    star = Star(cfg)
    inner, outer = star.habitable_zone
    regenerated = 0
    for _ in range(60):
        planet = Planet(cfg, star, star.habitable_zone, (inner + outer) / 2, planet_class="J", zone_override='e')
        before = [moon.planet_class for moon in planet.moons]
        planet.distance = outer * 3
        validation.reconcile_moved_planet(planet, keep_class=True)
        regenerated += sum(1 for old in before if not prog_c.PLANET_CLASSES[old]['c'])
        max_radius, max_mass = planetPhysics.moon_size_limits(planet)
        for moon in planet.moons:
            data = prog_c.PLANET_CLASSES[moon.planet_class]
            assert data['c'] and data["type"] == 't' and moon.planet_class not in prog_c.MOON_BLACKLIST
            assert moon.radius <= max_radius * (1 + 1e-9)
            assert moon.mass <= max_mass * (1 + 1e-9)
        assert validation.check_lunar_system(planet) == []
    assert regenerated > 0


def test_rogue_planet_split_follows_the_mass_function():
    """GEN.45: about 96% terrestrial and 4% gas giants, as dN/dlogM ~
    M^-0.65 from 0.1 Earth masses to 13 Jupiter masses gives."""
    slope = prog_c.ROGUE_PLANET_MASS_FUNCTION_SLOPE
    low = min(b[0] for b in prog_c.ROGUE_PLANET_MASS_BINS.values())
    high = max(b[1] for b in prog_c.ROGUE_PLANET_MASS_BINS.values())
    threshold = prog_c.ROGUE_PLANET_GAS_GIANT_MASS_THRESHOLD_JUPITER * pc.JUPITER_MASS_TO_KG / pc.EARTH_MASS_TO_KG
    expected = (threshold ** -slope - high ** -slope) / (low ** -slope - high ** -slope)
    assert expected == pytest.approx(0.036, abs=0.002)
    assert sum(rate for _lo, _hi, rate in prog_c.ROGUE_PLANET_MASS_BINS.values()) == pytest.approx(6.5)

    cfg = SystemConfig()
    n = 20000
    gas = sum(1 for _ in range(n) if RoguePlanet(cfg).planet_type == 'g')
    assert abs(gas / n - expected) < 0.008
