"""GEN.458: each star draws its own snow line, 0.2 dex around the 182 K blackbody position (Lodders 2003), and stores it once; the fixed 2.7 AU
constant is gone (`physics/formation.py`, `Star.snow_line_au`, `stars.snow_line_km`, `planets/moons.born_beyond_snow_line`)."""

import math
import concurrent.futures
import multiprocessing
import os
import statistics

import pytest

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import condensate_shares, constants, formation
from planetgen.util import draw

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")


def _config(star_type="G2V", **overrides):
    cfg = SystemConfig()
    cfg.STAR_TYPE = star_type
    cfg.PLANETS = True
    cfg.BINARY_SYSTEM = False
    for name, value in overrides.items():
        setattr(cfg, name, value)
    return cfg


def _system(seed, star_type="G2V", **overrides):
    with draw.bound(seed):
        return StarSystem(system_config=_config(star_type, **overrides))


# --- the centre and the scatter ---------------------------------------------------------------------------------

def test_the_centre_is_the_182_k_blackbody_position_and_scales_with_the_root_of_the_luminosity():
    sun = formation.snow_line_au(constants.SOLAR_LUMINOSITY)
    assert sun == pytest.approx(2.339, abs=0.005)
    # an independent route: the temperature of a blackbody at that distance is 182 K
    distance_m = sun * constants.AU_M
    temperature = (constants.SOLAR_LUMINOSITY / (16.0 * math.pi * constants.STEFAN_BOLTZMANN_CONSTANT * distance_m ** 2)) ** 0.25
    assert temperature == pytest.approx(condensate_shares.ICE_CONDENSATION_K["H2O"], rel=1e-9)
    assert formation.snow_line_au(4 * constants.SOLAR_LUMINOSITY) == pytest.approx(2 * sun)


def test_the_drawn_line_has_the_median_and_scatter_stated():
    luminosity = constants.SOLAR_LUMINOSITY
    centre = formation.snow_line_au(luminosity)
    with draw.bound(458):
        logs = sorted(math.log10(formation.draw_snow_line_au(luminosity) / centre) for _ in range(20000))
    assert statistics.median(logs) == pytest.approx(0.0, abs=0.01)
    # 0.2 dex, a little under it for the cut at 3 sigma
    assert statistics.pstdev(logs) == pytest.approx(tuning.SNOW_LINE_SCATTER_DEX, abs=0.005)
    limit = tuning.SNOW_LINE_SCATTER_DEX * tuning.SNOW_LINE_SCATTER_MAX_SIGMA
    assert -limit <= logs[0] and logs[-1] <= limit
    assert logs[0] < -0.5 and logs[-1] > 0.5, "the tails are reached"


def test_the_same_star_draws_the_same_line_from_its_own_stream():
    stream = lambda: draw.Stream("snow_line:test")
    first = formation.draw_snow_line_au(constants.SOLAR_LUMINOSITY, stream())
    assert first == formation.draw_snow_line_au(constants.SOLAR_LUMINOSITY, stream())


# --- in a built system ------------------------------------------------------------------------------------------

def _facts(seed):
    system = _system(seed)
    star = system.star
    return (star.snow_line_au, [(p.distance, p.born_beyond_snow_line, [m.born_beyond_snow_line for m in p.moons])
                                for p in system.planets if p.body_type != 'a'])


def test_a_system_gives_the_same_line_and_flags_from_its_seed_at_one_and_four_workers():
    seeds = list(range(6))
    serial = [_facts(s) for s in seeds]
    assert serial == [_facts(s) for s in seeds]
    assert len({line for line, _ in serial}) == len(seeds), "stars differ"
    with concurrent.futures.ProcessPoolExecutor(max_workers=4, mp_context=multiprocessing.get_context("spawn")) as pool:
        assert list(pool.map(_facts, seeds)) == serial


def test_a_star_stores_a_line_inside_the_scatter_of_its_blackbody_position():
    for seed in range(12):
        star = _system(seed).star
        centre = formation.snow_line_au(star.luminosity)
        assert abs(math.log10(star.snow_line_au / centre)) <= tuning.SNOW_LINE_SCATTER_DEX * tuning.SNOW_LINE_SCATTER_MAX_SIGMA + 1e-12


def test_a_planet_is_born_beyond_the_line_of_its_star_when_it_orbits_at_or_beyond_it_and_a_moon_takes_its_planets():
    seen = set()
    for seed in range(40):
        system = _system(seed, MAX_PLANETS=True)
        for planet in system.planets:
            if planet.body_type == 'a':
                continue
            assert planet.born_beyond_snow_line == (planet.distance >= planet.star.snow_line_au)
            seen.add(planet.born_beyond_snow_line)
            assert all(moon.born_beyond_snow_line == planet.born_beyond_snow_line for moon in planet.moons)
    assert seen == {True, False}


def test_a_close_pairs_line_is_the_root_of_the_sum_of_the_squares_of_its_stars():
    found = 0
    for seed in range(60):
        system = _system(seed, BINARY_SYSTEM=True)
        if getattr(system, "binary_type", None) != "close":
            continue
        pair = system.star
        assert pair.snow_line_au == pytest.approx(math.hypot(*(s.snow_line_au for s in pair.stars)))
        for planet in system.planets:
            if planet.body_type != 'a':
                assert planet.born_beyond_snow_line == (planet.distance >= pair.snow_line_au)
        found += 1
    assert found


def test_the_pair_survives_a_round_trip_through_a_dict():
    for seed in range(60):
        system = _system(seed, BINARY_SYSTEM=True)
        if getattr(system, "binary_type", None) != "close":
            continue
        again = type(system.star).from_dict(system.star.to_dict(), system.system_config)
        assert again.snow_line_au == pytest.approx(system.star.snow_line_au)
        return
    pytest.fail("no close pair in 60 seeds")


# --- nothing reads the old constant -----------------------------------------------------------------------------

