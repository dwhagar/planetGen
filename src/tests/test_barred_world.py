"""
The bar on habitable worlds (GEN.1250): with `HABITABLE_WORLD` False no finished planet or moon has the conditions a request for a habitable
world asks for (`habitable_request.meets`: temperate, an air of at least the minimum pressure, a surface ocean).

A planet is never drawn again for the bar (GEN.1397): the finished bodies are judged on their conditions, after the star's final age, the light of every
star and the oxygen of a life are in, and a body that meets the conditions or is rated habitable is dropped from the system
(`StarSystem._drop_barred_bodies`, counted in `barred_world_dropped`).
"""

import pytest

from planetgen.generation import habitable_request
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

# Systems the bar let through before GEN.1250 (found by a sample of 200 systems): the air remade at the star's final age (104), a planet outside
# the ecosphere zone (160), a thin air that the oxygen of its life lifted over the minimum pressure (270, 277).
SEEDS_THAT_BROKE_THE_BAR = (104, 160, 270, 277)
GALACTIC_DISTANCE_LY = 26000.0


def barred_system(seed):
    draw.set_run_seed(seed)
    config = SystemConfig()
    config.HABITABLE_WORLD = False
    config.PLANETS = True
    return StarSystem(system_config=config, galactic_center_dist_ly=GALACTIC_DISTANCE_LY)


def bodies_that_meet(system):
    return [(body.name, body.zone, body.surface_temperature, body.atmospheric_pressure, body.hydrosphere)
            for planet in system.all_planets if planet.body_type != "a" for body in (planet, *planet.moons) if habitable_request.meets(body)]


@pytest.mark.parametrize("seed", SEEDS_THAT_BROKE_THE_BAR)
def test_a_system_that_bars_habitable_worlds_has_none_that_met_the_request_before(seed):
    system = barred_system(seed)
    assert bodies_that_meet(system) == []
    assert system.barred_world_dropped >= 0


def test_the_bar_does_not_touch_a_system_that_does_not_bar_habitable_worlds():
    draw.set_run_seed(104)
    config = SystemConfig()
    config.PLANETS = True
    system = StarSystem(system_config=config, galactic_center_dist_ly=GALACTIC_DISTANCE_LY)
    assert not hasattr(system, "barred_world_dropped")


SAMPLE_CHUNKS = 6
SAMPLE_FIRST, SAMPLE_COUNT = 100, 300


@pytest.mark.slow
@pytest.mark.parametrize("chunk", range(SAMPLE_CHUNKS))
def test_over_300_sampled_systems_none_breaks_the_bar(chunk):
    """The census guard (GEN.1250): 300 systems of seeds 100 to 399, in six runs of 50 so that each stays inside the 900 s test limit (PERF.255)."""
    size = SAMPLE_COUNT // SAMPLE_CHUNKS
    first = SAMPLE_FIRST + chunk * size
    broken = {}
    for seed in range(first, first + size):
        found = bodies_that_meet(barred_system(seed))
        if found:
            broken[seed] = found
    assert not broken, f"{len(broken)} of {size} systems (seeds {first} to {first + size - 1}) hold a world that meets the request: {broken}"


def test_no_tide_or_heating_row_counts_a_body_the_bar_dropped():
    """GEN.1604: the bar drops a body before the tide pass, so the ledger of the survivors names no source that is not in the system (seed 433 dropped one and, before, left 40 rows naming it)."""
    system = barred_system(433)
    assert system.barred_world_dropped >= 1
    alive = {id(body) for planet in system.all_planets if planet.body_type != "a" for body in (planet, *planet.moons)}
    stars = {id(star) for star in (system.star, getattr(system, "secondary_star", None)) if star is not None}
    for effect in ("tide", "heating"):
        for _receiver, sources in system.effect_ledger.contributions(effect):
            for source, _value in sources:
                if hasattr(source, "is_moon"):
                    assert id(source) in alive, f"{effect}: a dropped body is still a source"
