# tests/test_age_law.py

"""GEN.546: a disk star's age inside its population's window follows the disk's star formation rate there, not an even
spread. The law is the history's own rate (`galaxy.disk_history`, GEN.469) as a piecewise-linear density; the samplers
(`stellar_evolution`, `star_population`) draw, weigh and count with the same law, so the scatter and a sector's own
draw agree. A population without a law draws evenly, as every population did before."""

import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import density, disk_history, mass_chain
from planetgen.galaxy.age_law import AgeLaw
from planetgen.generation import bright_stars, star_population
from planetgen.generation.config import SystemConfig
from planetgen.physics import stellar_evolution
from planetgen.util import draw
from tests.toy_shape import toy_shape


def _shape(stellar_mass):
    return density.build_galaxy_shape(
        disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
        arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4,
        masses=(stellar_mass, mass_chain.median_halo_mass(stellar_mass)))


SHAPE = _shape(5.0e10)
LAWS = disk_history.population_age_laws(SHAPE)
FALLING = AgeLaw((3.0, 10.0), (1.0, 0.0))
"""AgeLaw: A rate that starts at its highest and falls to nothing at the window's far end."""


def test_a_uniform_law_weighs_the_overlap_and_draws_like_rng_uniform():
    law = AgeLaw.uniform((0.1, 3.0))
    assert law.window == (0.1, 3.0) and law.span == pytest.approx(2.9)
    assert law.weight(0.0, 1.0) == pytest.approx(0.9)
    assert law.weight(2.0, 99.0) == pytest.approx(1.0)
    assert law.weight(4.0, 5.0) == 0.0
    first, second = draw.Stream(7), draw.Stream(7)
    assert law.sample(0.5, 2.5, first) == second.uniform(0.5, 2.5)


def test_a_law_needs_increasing_ages_and_a_positive_rate_somewhere():
    for ages, rates in (((1.0,), (1.0,)), ((1.0, 1.0), (1.0, 1.0)), ((1.0, 2.0), (1.0,)), ((1.0, 2.0), (0.0, 0.0)),
                        ((1.0, 2.0), (-1.0, 1.0))):
        with pytest.raises(ValueError):
            AgeLaw(ages, rates)


def test_the_weights_of_neighbouring_intervals_add_up_to_the_windows_length():
    for law in (FALLING, AgeLaw((0.1, 1.0, 3.0), (0.4, 1.0, 0.7)), LAWS["intermediate"], LAWS["old"]):
        low, high = law.window
        edges = np.linspace(low, high, 7)
        assert sum(law.weight(a, b) for a, b in zip(edges, edges[1:])) == pytest.approx(law.span)
        assert law.probability(low, high) == pytest.approx(1.0)


def test_a_falling_rate_puts_more_weight_where_it_is_high():
    assert FALLING.weight(3.0, 5.0) > FALLING.weight(8.0, 10.0)
    assert FALLING.probability(3.0, 6.5) == pytest.approx(0.75)    # the area of a triangle: 1 - (1 - 0.5) ** 2


def test_drawn_ages_follow_the_rate_across_the_window():
    """The age histogram of many draws is the law's own probability in each bin."""
    rng = draw.Stream(546)
    for law in (FALLING, LAWS["intermediate"], LAWS["old"]):
        ages = np.array([law.sample(*law.window, rng) for _ in range(40_000)])
        low, high = law.window
        assert ages.min() >= low and ages.max() <= high
        edges = np.linspace(low, high, 9)
        counts, _ = np.histogram(ages, bins=edges)
        for count, a, b in zip(counts, edges, edges[1:]):
            expected = law.probability(a, b)
            assert count / len(ages) == pytest.approx(expected, abs=4.0 * math.sqrt(expected / len(ages)) + 1e-3)


def test_a_draw_inside_a_subinterval_stays_in_it_and_follows_the_law_there():
    rng = draw.Stream(5)
    ages = np.array([FALLING.sample(4.0, 6.0, rng) for _ in range(20_000)])
    assert ages.min() >= 4.0 and ages.max() <= 6.0
    below_five = np.mean(ages < 5.0)
    assert below_five == pytest.approx(FALLING.weight(4.0, 5.0) / FALLING.weight(4.0, 6.0), abs=0.015)


def test_a_rate_that_starts_from_nothing_draws_without_error():
    rise = AgeLaw((3.0, 10.0), (0.0, 1.0))
    rng = draw.Stream(1)
    ages = [rise.sample(3.0, 10.0, rng) for _ in range(2000)]
    assert 3.0 <= min(ages) and max(ages) <= 10.0
    assert np.mean(ages) > 6.5


def test_the_history_gives_a_law_for_each_disk_population_and_none_for_the_rest():
    assert set(LAWS) == set(tuning.DISK_AGE_LAW_POPULATIONS)
    for name, law in LAWS.items():
        assert law.window == tuple(tuning.STELLAR_POPULATION_AGE_RANGES_GY[name])
        assert max(law.rates) == 1.0
    assert disk_history.population_age_laws(SHAPE) is LAWS


def test_the_laws_probabilities_are_the_histories_shares_inside_the_window():
    """Within one window the law is the history's own rate: the share of a window's stars younger than a mid age is the
    history's `share_between`, to the accuracy of the law's points."""
    found = disk_history.history_of(SHAPE)
    for name, law in LAWS.items():
        low, high = law.window
        if name == "old":
            high = tuning.STELLAR_THICK_START_GY    # GEN.1025: from there the thick disc adds its stars to the thin disc's own
        middle = 0.5 * (low + high)
        expected = found.share_between(low, middle) / found.share_between(low, high)
        assert law.probability(low, middle) / law.probability(low, high) == pytest.approx(expected, abs=0.01)


def test_a_light_galaxy_forms_its_old_stars_later_than_a_heavy_one():
    """Downsizing inside the old window: the mean age of its stars falls from a heavy galaxy to a light one."""
    def mean_age(mass):
        law = disk_history.history_of(_shape(mass)).age_law(*tuning.THIN_DISC_OLD_WINDOW_GY)    # the thin disc's own (GEN.1025)
        ages = np.array(law.ages)
        rates = np.array(law.rates)
        return np.trapezoid(ages * rates, ages) / np.trapezoid(rates, ages)
    assert mean_age(1.0e10) < mean_age(5.0e10)


def test_a_living_star_of_a_front_loaded_law_is_younger_than_one_of_an_even_spread():
    law = AgeLaw((3.0, 10.0), (1.0, 0.0))
    rng_even, rng_law = draw.Stream(11), draw.Stream(11)
    even = [stellar_evolution.sample_living_star(rng=rng_even, population="old")[1] for _ in range(1500)]
    spread = [stellar_evolution.sample_living_star(rng=rng_law, population="old", age_law=law)[1] for _ in range(1500)]
    assert np.mean(spread) < np.mean(even) - 0.8


def test_a_fallback_to_the_whole_disk_drops_the_law():
    """A habitable host needs a star older than the young window: the draw falls back to the whole disk's ages, and
    the young window's law must not follow it there."""
    law = LAWS["young"]
    ages = [stellar_evolution.sample_living_star(rng=draw.Stream(i), population="young", habitable_host=True,
                                                 age_law=law)[1] for i in range(20)]
    assert min(ages) >= tuning.LIFE_MIN_STAR_AGE_GY


FRONT_LOADED = AgeLaw((0.1, 3.0), (1.0, 0.0))
"""AgeLaw: The intermediate window with its stars mostly young."""
BACK_LOADED = AgeLaw((0.1, 3.0), (0.0, 1.0))
"""AgeLaw: The intermediate window with its stars mostly old."""


def test_the_bright_fraction_with_a_law_is_the_share_the_whole_model_draw_gives():
    """`bright_star_fraction` weighs the bright phases with the law and the whole-model sampler draws the ages with it,
    so the two agree; a law other than the even one moves the fraction (most of the bright stars of this window are
    young main-sequence stars)."""
    threshold = star_population.MASS_PASS_MIN_LUMINOSITY_SOL
    laws = {"intermediate": FRONT_LOADED}
    even = star_population.bright_star_fraction(threshold, "intermediate")
    spread = star_population.bright_star_fraction(threshold, "intermediate", age_laws=laws)
    assert spread > 1.4 * even
    assert star_population.bright_star_fraction(
        threshold, "intermediate", age_laws={"intermediate": BACK_LOADED}) < 0.6 * even
    rng = draw.Stream(546)
    draws = 40_000
    bright = 0
    for _ in range(draws):
        mass, age, state = stellar_evolution.sample_living_star(rng=rng, population="intermediate",
                                                                age_law=FRONT_LOADED)
        bright += state["luminosity_sol"] >= threshold
    assert bright / draws == pytest.approx(spread, abs=4.0 * math.sqrt(spread / draws))


def test_the_bright_star_draw_follows_the_law():
    """The stars the scatter draws above a floor are the whole-model draw conditioned on the floor: a front-loaded law
    draws younger ones than a back-loaded one."""
    floor = star_population.MASS_PASS_MIN_LUMINOSITY_SOL
    young = star_population.sample_bright_stars(800, floor, "intermediate", draw.Stream(3),
                                                age_laws={"intermediate": FRONT_LOADED})
    old = star_population.sample_bright_stars(800, floor, "intermediate", draw.Stream(3),
                                              age_laws={"intermediate": BACK_LOADED})
    assert np.mean([star["age_gy"] for star in young]) < np.mean([star["age_gy"] for star in old]) - 0.1


def test_the_mass_band_draw_and_its_fraction_use_the_law():
    """A massive star is alive only early in the young window, so the law moves the share of the heavy band and the
    ages of the stars drawn in it."""
    front, back = AgeLaw((0.0, 0.1), (1.0, 0.0)), AgeLaw((0.0, 0.1), (0.0, 1.0))
    assert star_population.mass_band_fraction(8.0, 20.0, population="young", age_laws={"young": front}) > 1.5 * (
        star_population.mass_band_fraction(8.0, 20.0, population="young", age_laws={"young": back}))
    young = star_population.sample_mass_band_stars(800, 8.0, 20.0, "young", draw.Stream(2), age_laws={"young": front})
    old = star_population.sample_mass_band_stars(800, 8.0, 20.0, "young", draw.Stream(2), age_laws={"young": back})
    assert np.mean([star["age_gy"] for star in young]) < np.mean([star["age_gy"] for star in old])


TOY_SHAPE = toy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0, bulge_amplitude=2.0, arm_count=2,
    pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
EDGE_PC = 9.2
"""float: The toy galaxy's sector edge, parsecs (its cells are small, so a few shells hold bright stars)."""


@pytest.fixture
def recorded_laws(monkeypatch):
    """Every `age_laws` the scatter and the backfills hand the samplers, by sampler."""
    seen = {"bright": [], "mass": []}
    real_bright, real_mass = bright_stars.sample_bright_stars, bright_stars.sample_mass_band_stars

    def bright(*args, **kwargs):
        seen["bright"].append(kwargs.get("age_laws"))
        return real_bright(*args, **kwargs)

    def mass(*args, **kwargs):
        seen["mass"].append(kwargs.get("age_laws"))
        return real_mass(*args, **kwargs)

    monkeypatch.setattr(bright_stars, "sample_bright_stars", bright)
    monkeypatch.setattr(bright_stars, "sample_mass_band_stars", mass)
    return seen


def test_the_scatter_and_both_backfills_draw_with_the_galaxys_laws(recorded_laws):
    laws = bright_stars.age_laws_of(TOY_SHAPE)
    assert set(laws) == set(tuning.DISK_AGE_LAW_POPULATIONS)
    rows = list(bright_stars.scatter_shell(TOY_SHAPE, 3, 3, EDGE_PC, 2000.0, 100.0, 7))
    assert rows and recorded_laws["bright"] and all(found is laws for found in recorded_laws["bright"])
    recorded_laws["bright"].clear()
    addresses = [(3, 0, wedge) for wedge in range(12)]
    assert list(bright_stars.backfill_cells(TOY_SHAPE, addresses, EDGE_PC, 2000.0, 50.0, 100.0, 2))
    assert recorded_laws["bright"] and all(found is laws for found in recorded_laws["bright"])
    assert list(bright_stars.backfill_mass_cells(TOY_SHAPE, addresses, EDGE_PC, 2000.0, 5.0, None, None, 99,
                                                 bright_stars.mass_band_edges(None)))
    assert recorded_laws["mass"] and all(found is laws for found in recorded_laws["mass"])


def test_a_sector_fill_hands_its_dim_systems_the_law_of_their_population():
    from planetgen.galaxy.profile import profile_at
    from planetgen.galaxy.spherical import sector_position_pc
    fill = bright_stars.FillContext(profile_at(sector_position_pc(2, 0, 3, EDGE_PC), TOY_SHAPE, EDGE_PC, 2000.0),
                                    min_luminosity_sol=2500.0)
    laws = bright_stars.age_laws_of(TOY_SHAPE)
    seen = set()
    for seed in range(40):
        config = fill.apply(SystemConfig(), draw.Stream(seed))
        assert config.AGE_LAW is laws.get(config.POPULATION)
        seen.add(config.POPULATION)
    assert seen & set(laws)
