# tests/test_belt_planet_moon_tides.py

"""
GEN.1541 (GEN.1464 part C): a belt member the planet generator made a planet (GEN.1510) raises tides on its own moons like any planet, and is raised on by the star and its sibling planets (the belts' tides are the effect ledger's, GEN.1464 part B).

There is no second way: the promoted body is in the star's list of planets, so `tide_sources` already finds its sources and the sources of its moons. These tests pin that, and that the belt it came
from is not counted twice (the member left the belt's mass when it was promoted).
"""

import pytest

from planetgen.generation import irregular_moons
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import tide_sources
from planetgen.util import draw

SEEDS = (0, 6)


def _config():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = True
    cfg.MOONS = True
    cfg.BINARY_SYSTEM = False
    return cfg


def _siblings(system):
    """The planets `tide_sources.register` hands each body: every body of the system's list that is not a belt (a belt is not a source there, GEN.1464)."""
    return [p for p in system.primary.orbiters if getattr(p, "body_type", None) != "a"]


def _with_moons(seed):
    draw.set_run_seed(seed)
    system = StarSystem(system_config=_config())
    return system, [p for p in system.primary.orbiters if getattr(p, "belt_rank", None) is not None and getattr(p, "moons", None)]


@pytest.fixture(scope="module")
def cases():
    found = []
    for seed in SEEDS:
        system, planets = _with_moons(seed)
        found.extend((system, planet) for planet in planets[:2])
    assert found, "no promoted planet with a moon in the test seeds"
    return found


def test_a_promoted_planet_is_raised_on_by_the_star_and_its_siblings(cases):
    for system, planet in cases:
        kinds = [source.kind for source in tide_sources.sources_of(planet, None, _siblings(system))]
        assert kinds.count("star") == 1
        assert kinds.count("planet") >= 1


def test_its_moon_is_raised_on_by_the_planet_and_the_star(cases):
    for system, planet in cases:
        moon = planet.moons[0]
        sources = tide_sources.sources_of(moon, planet, _siblings(system))
        kinds = [source.kind for source in sources]
        assert sources[0].body is planet
        assert kinds.count("star") == 1


def test_the_planet_raises_the_tide_the_ocean_tide_model_reads_on_its_moon(cases):
    for system, planet in cases:
        moon = planet.moons[0]
        own = tide_sources.own_sources(moon, planet)
        assert [source.body for source in own] == [planet]
        assert own[0].g3 > 0.0 and own[0].mass_kg == float(planet.mass)


def test_the_promoted_planet_is_not_counted_again_in_its_belt(cases):
    for system, planet in cases:
        listed = [source.body for source in tide_sources.sources_of(planet, None, _siblings(system))]
        assert planet not in listed
        own = irregular_moons.belts_of(system.primary.orbiters)[planet.belt_order]        # the planet's own belt: a rank is the member's place in ITS belt, so another belt may hold a member of the same rank
        assert planet.belt_rank not in [member["rank"] for member in own.members]
