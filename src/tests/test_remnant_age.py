# tests/test_remnant_age.py

"""GEN.553: a neutron star's or stellar black hole's age is a remnant's in its population: the time since the star that
made it collapsed, which is the star's birth age (the population's window and age law, as any star's) less the time the
progenitor took to reach its end."""

import pytest

from planetgen import tuning
from planetgen.galaxy.age_law import AgeLaw
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.physics import stellar_evolution as evolution
from planetgen.util import draw


def _ages(kind, population=None, age_law=None, age_bias=None, count=300):
    with draw.bound(553):
        return [evolution.collapse_age_gy(kind, age_bias, draw, population, age_law) for _ in range(count)]


@pytest.mark.parametrize("kind", ["neutron-star", "black-hole"])
@pytest.mark.parametrize("population", [None, "young", "intermediate", "old", "bulge", "halo"])
def test_a_remnant_is_never_older_than_its_populations_oldest_star(kind, population):
    low, high = tuning.STELLAR_POPULATION_AGE_RANGES_GY[population] if population else tuning.STAR_FORMATION_AGE_RANGE_GY
    ages = _ages(kind, population)
    assert all(tuning.COMPACT_REMNANT_AGE_RANGE_GY[0] <= age <= high for age in ages)


@pytest.mark.parametrize("kind", ["neutron-star", "black-hole"])
def test_an_old_populations_remnants_are_old_and_a_young_populations_are_young(kind):
    old, young = _ages(kind, "old"), _ages(kind, "young")
    assert min(old) > tuning.STELLAR_POPULATION_AGE_RANGES_GY["old"][0] - 0.1
    assert max(young) <= tuning.STELLAR_POPULATION_AGE_RANGES_GY["young"][1]
    assert sum(old) / len(old) > 20 * sum(young) / len(young)


def test_the_age_is_the_birth_age_less_the_time_the_progenitor_took_to_collapse():
    """Over the whole disk's window a remnant's age is uniform from nothing to the window's end less the progenitor's
    collapse time, and the progenitors of both kinds (8 solar masses and heavier) collapse inside 0.1 Gy."""
    high = tuning.STAR_FORMATION_AGE_RANGE_GY[1]
    for kind in ("neutron-star", "black-hole"):
        ages = _ages(kind, count=3000)
        assert max(ages) < high
        assert sum(ages) / len(ages) == pytest.approx(0.5 * high, abs=0.25)


def test_a_remnant_follows_the_populations_age_law():
    window = tuning.STELLAR_POPULATION_AGE_RANGES_GY["old"]
    front = AgeLaw(ages=(window[0], window[1]), rates=(1.0, 0.0))
    early = _ages("neutron-star", "old", front)
    flat = _ages("neutron-star", "old", None)
    assert sum(early) / len(early) < window[0] + 0.4 * (window[1] - window[0])
    assert sum(flat) / len(flat) > window[0] + 0.4 * (window[1] - window[0])


def test_a_window_no_progenitor_can_have_collapsed_in_gives_the_youngest_age():
    """A black hole's progenitors (20 solar masses and up) take millions of years: none can have collapsed in a window
    that closes sooner, so the remnant is as young as a remnant is."""
    floor = tuning.COMPACT_REMNANT_AGE_RANGE_GY[0]
    assert evolution.youngest_progenitor_sol(1.0) > tuning.BLACK_HOLE_PROGENITOR_SOL[1]
    assert evolution.collapse_age_gy("black-hole", population="cluster", age_law=None) >= floor
    assert evolution.collapse_age_gy("black-hole", None, draw, "open_cluster_00") == floor


def test_the_formation_time_falls_with_the_progenitors_mass():
    times = [evolution.remnant_formation_time_gy(mass) for mass in (8.0, 15.0, 30.0, 100.0)]
    assert times == sorted(times, reverse=True)
    assert 0.01 < times[0] < 0.1


def test_a_neutron_star_and_a_stellar_black_hole_take_the_config_population_but_other_black_holes_do_not():
    config = SystemConfig()
    config.POPULATION = "young"
    with draw.bound(553):
        neutron = NeutronStar(config)
        stellar = BlackHole(config, mass_class="stellar")
        intermediate = BlackHole(config, mass_class="intermediate")
    assert neutron.age <= 0.1 and stellar.age <= 0.1
    with draw.bound(553):
        ages = [BlackHole(SystemConfig(), mass_class="intermediate").age for _ in range(40)]
    assert max(ages) > 1.0 and 1.0e-5 <= min(ages)
    assert intermediate.lifespan == float("inf")


def test_a_given_age_is_kept():
    config = SystemConfig()
    config.POPULATION = "young"
    with draw.bound(553):
        assert NeutronStar(config, age_gy=7.5).age == 7.5


def test_a_sectors_ejected_remnants_take_the_age_of_the_population_its_fill_gives_them():
    from types import SimpleNamespace
    from planetgen.generation import fill_stages, run_galaxy, run_sector
    from planetgen.galaxy.sector import SpaceSector

    def young(config, rng=draw):
        config.POPULATION = "young"
        return config

    fill = SimpleNamespace(profile=SimpleNamespace(ejecta={"neutron-star": 25.0, "black-hole": 0.0}),
                           apply_population=young)
    args = run_galaxy._default_generation_args()
    args.num_systems = 0
    args.density = None
    with draw.bound(5):
        sector = SpaceSector("Remnant Age Test", edge_ly=100.0)
        fill_stages.begin(sector)
        entries = run_sector.add_ejecta(sector, args, fill)
    assert len(entries) > 10
    assert all(entry.phenomenon.age <= tuning.STELLAR_POPULATION_AGE_RANGES_GY["young"][1] for entry in entries)
