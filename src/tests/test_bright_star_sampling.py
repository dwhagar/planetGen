"""
Bright/dim star sampling (`planetgen.generation.star_population`), stellar
populations by position (`density.population_densities`) and
pre-placed systems (`SpaceSector.add_preplaced_system`): GEN.21
and GEN.22, the Physics part of bright-star pre-placement.

Run with: pytest src/tests/test_bright_star_sampling.py
"""
import collections
import math
import random

import pytest

from planetgen.galaxy import density
from planetgen.physics import constants as phys_c
from planetgen import tuning as prog_c
from planetgen.physics import stellar_evolution as se
from planetgen.generation import star_population as sp
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.sector import SpaceSector, distance_between, required_separation_ly
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem

THRESHOLD = 500.0

SHAPE = density.build_galaxy_shape(
    disk_scale_length_pc=2800.0, disk_scale_height_pc=350.0, bulge_scale_radius_pc=200.0,
    bulge_amplitude=1.0, arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4,
)


def _brute_force(population, draws, seed):
    """Whole-model draws and the bright ones among them."""
    rng = random.Random(seed)
    bright = []
    for _ in range(draws):
        mass, age, state = se.sample_living_star(rng=rng, population=population)
        if state["luminosity_sol"] >= THRESHOLD:
            bright.append(state)
    return bright


@pytest.fixture(scope="module")
def young_brute_force():
    draws = 150_000
    return draws, _brute_force("young", draws, seed=11)


def test_threshold_must_be_at_least_every_white_dwarf():
    with pytest.raises(ValueError):
        sp.bright_star_fraction(prog_c.WD_LUMINOSITY_RANGE_SOL[1] * 0.99)
    assert sp.bright_star_fraction(prog_c.WD_LUMINOSITY_RANGE_SOL[1]) > 0.0
    with pytest.raises(ValueError):
        sp.sample_bright_stars(1, float("nan"))
    with pytest.raises(ValueError):
        sp.bright_star_fraction(THRESHOLD, "halo")


def test_bright_fraction_matches_the_whole_model(young_brute_force):
    draws, bright = young_brute_force
    expected = sp.bright_star_fraction(THRESHOLD, "young")
    sigma = math.sqrt(expected * draws) / draws
    assert abs(len(bright) / draws - expected) < 4 * sigma


def test_bright_fraction_of_the_whole_disk_is_well_under_a_thousandth():
    # The galaxy studies' estimate: ~7e-4 of all stars at >= 500 Lsun.
    assert 4e-4 < sp.bright_star_fraction(THRESHOLD) < 1e-3
    assert sp.bright_star_fraction(THRESHOLD, "young") > 5 * sp.bright_star_fraction(THRESHOLD, "old")


def test_bright_stars_are_all_bright_and_match_the_whole_models_mix(young_brute_force):
    _, brute = young_brute_force
    stars = sp.sample_bright_stars(20_000, THRESHOLD, "young", random.Random(3))
    assert all(s["luminosity_w"] >= THRESHOLD * phys_c.SOLAR_LUMINOSITY * (1 - 1e-12) for s in stars)
    sampled = collections.Counter(s["yerkes_class"] == "V" for s in stars)
    reference = collections.Counter(state["yerkes_class"] == "V" for state in brute)
    share = sampled[True] / len(stars)
    reference_share = reference[True] / len(brute)
    sigma = math.sqrt(reference_share * (1 - reference_share) / len(brute))
    assert abs(share - reference_share) < 4 * sigma + 0.01


@pytest.mark.parametrize("population", ["old", "bulge"])
def test_old_populations_have_only_bright_giants(population):
    stars = sp.sample_bright_stars(2000, THRESHOLD, population, random.Random(4))
    assert {s["yerkes_class"] for s in stars} <= {"III", "II"}
    low = se.population_age_range_gy(population)[0]
    assert all(s["age_gy"] >= low for s in stars)


def test_seeded_draws_repeat():
    first = sp.sample_bright_stars(50, THRESHOLD, None, random.Random(99))
    second = sp.sample_bright_stars(50, THRESHOLD, None, random.Random(99))
    assert first == second


def test_dim_stars_are_all_dim():
    rng = random.Random(5)
    for _ in range(2000):
        params = sp.sample_dim_star(THRESHOLD, "young", rng)
        assert params["luminosity_w"] < THRESHOLD * phys_c.SOLAR_LUMINOSITY
        assert params["age_gy"] < se.population_age_range_gy("young")[1]


def test_bright_star_params_rebuild_a_full_system():
    params = sp.sample_bright_stars(1, THRESHOLD, "intermediate", random.Random(6))[0]
    system = StarSystem(system_config=SystemConfig(), primary_star_params=params)
    stars = [system.primary_star, getattr(system, "secondary_star", None)]
    assert any(star is not None and star.luminosity == params["luminosity_w"] for star in stars)


@pytest.mark.parametrize("position", [(8000, 0, 0), (8000, 0, 300), (500, 0, 0), (0, 0, 0), (8000, 0, 20000)])
def test_population_densities_sum_to_the_total(position):
    densities = density.population_densities(position, SHAPE)
    assert set(densities) == set(sp.POPULATIONS)
    assert all(value >= 0 for value in densities.values())
    assert sum(densities.values()) == pytest.approx(density.relative_density(position, SHAPE), rel=1e-12)


def test_young_stars_hug_the_plane_and_the_bulge_is_old():
    def shares(position):
        densities = density.population_densities(position, SHAPE)
        total = sum(densities.values())
        return {name: value / total for name, value in densities.items()}

    plane, high = shares((8000, 0, 0)), shares((8000, 0, 1000))
    assert plane["young"] > 100 * high["young"]
    assert high["old"] > 0.9
    assert shares((0, 0, 2000))["bulge"] > 0.5
    assert shares((8000, 0, 0))["bulge"] < 1e-6


def test_pick_population_follows_the_densities():
    rng = random.Random(8)
    counts = collections.Counter(sp.pick_population({"young": 1, "old": 3, "bulge": 0}, rng) for _ in range(8000))
    assert counts["bulge"] == 0
    assert counts["old"] / 8000 == pytest.approx(0.75, abs=0.03)
    assert sp.pick_population({"young": 0.0}) is None


def test_config_population_and_luminosity_cap_steer_random_stars():
    cfg = SystemConfig()
    cfg.POPULATION = "young"
    cfg.MAX_STAR_LUMINOSITY_SOL = 1.0
    for _ in range(200):
        star = Star(cfg)
        assert star.age < se.population_age_range_gy("young")[1]
        # White dwarfs are exempt from the cap (`sample_living_star`): a
        # young massive star's remnant can still be brighter than the Sun.
        if star.yerkes_class != "VII":
            assert star.luminosity < phys_c.SOLAR_LUMINOSITY


@pytest.mark.parametrize("population, flag", [("bulge", "LARGE_STAR"), ("young", "HABITABLE_WORLD")])
def test_impossible_population_requests_fall_back_to_the_whole_disk(population, flag):
    cfg = SystemConfig()
    cfg.POPULATION = population
    setattr(cfg, flag, True)
    for _ in range(20):
        Star(cfg)  # would never find a star in the population alone


def _system():
    cfg = SystemConfig()
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    return StarSystem(system_config=cfg)


def test_preplaced_systems_keep_their_position_and_clear_later_systems():
    sector = SpaceSector("preplaced", edge_ly=40.0)
    bright = StarSystem(system_config=SystemConfig(),
                        primary_star_params=sp.sample_bright_stars(1, THRESHOLD, "young", random.Random(9))[0])
    entry = sector.add_preplaced_system(bright, (1.0, -2.0, 3.0))
    assert entry.preplaced and entry.position == (1.0, -2.0, 3.0)
    for _ in range(5):
        other = sector.add_system(_system())
        assert not other.preplaced
        assert distance_between(other, entry) >= required_separation_ly(other.star_system, bright)

    reloaded = SpaceSector.from_dict(sector.to_dict())
    assert [e.preplaced for e in reloaded.entries] == [e.preplaced for e in sector.entries]


def test_preplaced_position_must_be_inside_the_sector():
    sector = SpaceSector("preplaced", edge_ly=10.0)
    with pytest.raises(ValueError):
        sector.add_preplaced_system(_system(), (6.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        sector.add_preplaced_system(_system(), (float("nan"), 0.0, 0.0))


def test_the_dim_cap_never_removes_a_white_dwarf(monkeypatch):
    # The hottest white dwarfs are clamped to exactly the lowest allowed
    # cap; they are never pre-placed, so the dim draw must keep them.
    from planetgen.physics import stellar_evolution
    cap = prog_c.WD_LUMINOSITY_RANGE_SOL[1]
    white_dwarf = {"yerkes_class": "VII", "luminosity_sol": cap}
    monkeypatch.setattr(stellar_evolution, "evolve_star", lambda mass, age, rng: dict(white_dwarf))
    assert stellar_evolution.sample_living_star(max_luminosity_sol=cap)[2] == white_dwarf
    giant = {"yerkes_class": "III", "luminosity_sol": cap}
    monkeypatch.setattr(stellar_evolution, "evolve_star", lambda mass, age, rng: dict(giant))
    with pytest.raises(ValueError):
        stellar_evolution.sample_living_star(max_luminosity_sol=cap)
