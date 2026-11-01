# tests/test_gen_system_mass.py -- GEN.545: the mean mass of a system is worked out from the generator's own draws, not
# recalled. Each piece is checked against a seeded draw made with the generator's own samplers.

import math
import statistics

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import density, disk_history, system_mass as galaxy_system_mass
from planetgen.generation import mass_ratio, multistar, system_mass
from planetgen.physics import imf, mathcheck, remnant_mass, remnant_progenitors, stellar_evolution
from planetgen.util import draw


def test_the_expected_star_count_is_the_mean_of_the_generators_own_draw():
    for mass in (0.1, 1.0, 5.0, 50.0):
        rng = draw.Stream(f"gen545-count-{mass}")
        counts = [multistar._draw_star_count(mass, multistar.MAX_MEMBERS, 1.0, rng) for _ in range(20000)]
        error = statistics.stdev(counts) / math.sqrt(len(counts))
        assert abs(statistics.fmean(counts) - multistar.expected_star_count(mass)) < 4.0 * error


def test_a_lower_cap_and_a_lower_prevalence_lower_the_expected_count():
    full = multistar.expected_star_count(50.0)
    assert multistar.expected_star_count(50.0, max_stars=7) < multistar.expected_star_count(50.0, max_stars=9) < full
    assert multistar.expected_star_count(50.0, prevalence=0.5) < full
    assert multistar.expected_star_count(50.0, max_stars=1) == pytest.approx(
        sum((n + 1) * w for n, w in enumerate(multistar.star_count_weights(50.0))) / sum(multistar.star_count_weights(50.0)))


@pytest.mark.parametrize("age_gy", [0.3, 3.0, 10.0])
def test_the_placed_stars_mean_mass_is_a_seeded_draw_of_the_initial_mass_function(age_gy):
    rng = draw.Stream(f"gen545-imf-{age_gy}")
    masses = imf.sample
    drawn = system_mass.present_mass_sol(np.array([masses(rng=rng) for _ in range(40000)]), age_gy)
    error = drawn.std() / math.sqrt(len(drawn))
    number, mean_mass, _others = system_mass._placed_bins()
    derived = float(np.sum(number * system_mass.present_mass_sol(mean_mass, age_gy)) / number.sum())
    assert abs(derived - drawn.mean()) < 4.0 * error


def test_the_placed_stars_mean_agrees_with_the_moments_the_perturber_mass_uses():
    # They differ only in where a star leaves the main sequence (the giant phases' end, not the main sequence's).
    for age_gy in (1.0, 5.0, 10.0):
        mass, _square = system_mass.present_day_moments(age_gy)
        old = mass / imf.integral(0.0, tuning.IMF_BREAKS_SOL[0], tuning.IMF_BREAKS_SOL[-1])
        number, mean_mass, _others = system_mass._placed_bins()
        new = float(np.sum(number * system_mass.present_mass_sol(mean_mass, age_gy)) / number.sum())
        # The primaries are heavier than the stars born by the light stars that are companions: 1.6 to 2.4 percent with a flat mass ratio,
        # 3.3 to 4.3 with the measured one (GEN.951: massive primaries' companions lean light, so fewer light stars are primaries).
        assert new == pytest.approx(old, rel=0.05)


@pytest.mark.parametrize("mass_sol, age_gy", [(0.3, 5.0), (1.0, 5.0), (1.0, 12.0), (3.0, 0.5), (3.0, 10.0), (8.0, 0.02)])
def test_the_other_stars_of_a_system_are_what_the_member_tree_draws(mass_sol, age_gy, monkeypatch):
    # The star-type draws of GEN.270 (a Be star's light stripped companion, a merger's missing star) are measured on
    # their own in test_close_binary_types.py.
    monkeypatch.setattr(multistar, "draw_be_star", lambda *args: (False, False))
    monkeypatch.setattr(multistar, "draw_merger_product", lambda *args: None)
    # The real trees (stability screens and dropped members included) against the expected count times the mean
    # present mass of a star born with 0.1 to 1 of the primary's mass.
    rng = draw.Stream(f"gen545-tree-{mass_sol}-{age_gy}")
    placed = float(system_mass.present_mass_sol(np.array([mass_sol]), age_gy)[0])
    others = [multistar.draw_member_tree(multistar.Placed("star", mass_sol), age_gy * 1.0e9, rng=rng).system_mass_solar - placed
              for _ in range(2500)]
    fractions = mass_ratio.Ratios(mass_sol).quantiles((np.arange(64) + 0.5) / 64)[0]
    derived = (multistar.expected_star_count(mass_sol) - 1.0) * float(
        system_mass.present_mass_sol(mass_sol * fractions, age_gy).mean())
    error = statistics.stdev(others) / math.sqrt(len(others))
    assert abs(statistics.fmean(others) - derived) < 4.0 * error + 0.01 * derived


def _derived_mean_from(age_gy, lowest_primary_sol):
    """The derived mean over the placed stars of `lowest_primary_sol` or more (the same sums as
    `system_mass.mean_system_mass_at_age_solar`, over the bins whose mean mass reaches it)."""
    number, mean_mass, others = system_mass._placed_bins()
    kept = mean_mass >= lowest_primary_sol
    fractions = system_mass.companion_fractions(mean_mass)
    further = system_mass.present_mass_sol(mean_mass[:, None] * fractions, age_gy).mean(axis=1)
    placed = system_mass.present_mass_sol(mean_mass, age_gy)
    return float(np.sum((number * (placed + others * further))[kept]) / np.sum(number[kept]))


@pytest.mark.parametrize("age_gy", [0.3, 3.0, 10.0])
def test_the_mean_system_mass_is_a_seeded_draw_of_systems(age_gy):
    # The generator's own samplers end to end: an initial mass, the member tree of a star that is alive or a white dwarf,
    # the remnant mass of a neutron star or a black hole (the scatter places those without companions).
    lowest = tuning.IMF_BREAKS_SOL[0]
    rng = draw.Stream(f"gen545-systems-{age_gy}")
    found = []
    for _ in range(4000):
        mass = imf.sample(rng=rng)
        if stellar_evolution.evolve_star(mass, age_gy, rng) is None:
            kind = "black-hole" if remnant_progenitors.fails(mass) else "neutron-star"
            found.append(remnant_mass.draw_mass_sol(kind, rng=rng))
            continue
        found.append(multistar.draw_member_tree(multistar.Placed("star", mass), age_gy * 1.0e9, rng=rng).system_mass_solar)
    error = statistics.stdev(found) / math.sqrt(len(found))
    assert abs(statistics.fmean(found) - _derived_mean_from(age_gy, lowest)) < 4.0 * error + 0.01 * statistics.fmean(found)


def test_the_mean_falls_with_age_and_sits_above_the_placed_star_alone():
    ages = (0.05, 0.5, 2.0, 6.0, 12.0)
    means = [system_mass.mean_system_mass_at_age_solar(age) for age in ages]
    assert means == sorted(means, reverse=True)
    for age, mean in zip(ages, means):
        number, mean_mass, _others = system_mass._placed_bins()
        alone = float(np.sum(number * system_mass.present_mass_sol(mean_mass, age)) / number.sum())
        assert mean > alone
    assert 0.45 < means[-1] < 0.55            # an old population: a placed star of 0.35 and 1.4 stars a system


def test_the_window_mean_is_the_mean_over_its_ages_and_a_law_weights_them():
    assert system_mass.window_mean_system_mass_solar(3.0, 10.0) == pytest.approx(
        system_mass.cohort_mean_system_mass_solar(3.0 + 7.0 * (np.arange(system_mass.COHORT_AGES) + 0.5) / system_mass.COHORT_AGES))
    young, old = system_mass.mean_system_mass_at_age_solar(0.5), system_mass.mean_system_mass_at_age_solar(8.0)
    assert system_mass.cohort_mean_system_mass_solar([0.5, 8.0], [3.0, 1.0]) == pytest.approx(0.75 * young + 0.25 * old)
    with pytest.raises(ValueError):
        system_mass.cohort_mean_system_mass_solar([1.0], [0.0])


def test_the_compact_companions_left_out_are_under_the_bound():
    # 1% of the stars of 2 to 14 solar masses get a white dwarf, 0.03% a neutron star and 0.02% a black hole.
    table = multistar.MEMBER_STATE_BY_PLACED["star-intermediate"]
    number = imf.integral(0.0, 2.0, 14.0) / imf.integral(
        0.0, tuning.IMF_BREAKS_SOL[0], tuning.IMF_BREAKS_SOL[-1])
    massive = multistar.MEMBER_STATE_BY_PLACED["star-massive"]
    heavy = imf.integral(0.0, 14.0, tuning.IMF_BREAKS_SOL[-1]) / imf.integral(
        0.0, tuning.IMF_BREAKS_SOL[0], tuning.IMF_BREAKS_SOL[-1])
    added = (number * (table["white-dwarf"][0] * 0.9 + table["neutron-star"][0] * 1.4 + table["black-hole"][0] * 9.0)
             + heavy * (massive["neutron-star"][0] * 1.4 + massive["black-hole"][0] * 9.0))
    assert added / system_mass.mean_system_mass_at_age_solar(5.0) < system_mass.COMPACT_COMPANION_SHARE_BOUND


def test_a_galaxys_mean_is_its_populations_averaged_by_their_systems():
    shape = mathcheck._reference_galaxy_shape()
    means = galaxy_system_mass.population_means_solar(shape)
    units = density.population_units(shape)
    expected = sum(units[name] * means[name] for name in units) / sum(units.values())
    assert density.mean_system_mass_solar(shape) == pytest.approx(expected)
    assert means["young"] > means["intermediate"] > means["old"] > means["bulge"] > means["halo"]
    assert 0.45 < density.mean_system_mass_solar(shape) < 0.56
    assert galaxy_system_mass.mean_over_solar({"old": 3.0, "bulge": 1.0}, shape) == pytest.approx(
        0.75 * means["old"] + 0.25 * means["bulge"])
    assert galaxy_system_mass.mean_over_solar({}, shape) == means["old"]


def test_the_galaxy_holds_the_systems_its_stellar_mass_makes_at_the_derived_mean():
    shape = mathcheck._reference_galaxy_shape()
    local = mathcheck.pc.LOCAL_STELLAR_DENSITY_LY3 / mathcheck.units.ly_to_pc(1.0) ** 3
    systems = shape.k_norm * local * density.total_mass(shape)
    assert systems == pytest.approx(shape.structure.stellar_mass_solar / density.mean_system_mass_solar(shape), rel=1e-9)


def test_the_history_does_not_depend_on_the_normalization():
    shape = mathcheck._reference_galaxy_shape()
    assert disk_history.history_of(shape) is disk_history.history_of(shape._replace(k_norm=shape.k_norm * 3.0))


@pytest.mark.parametrize("age_gy", [1.0e-4, 0.02, 0.3, 1.7, 4.0, 9.5, 13.8])
def test_the_cached_means_are_the_means_the_public_masses_give(age_gy):
    # PERF.146: the means a galaxy's thousands of clusters ask for share what no age changes (the remnants' masses and the
    # companions' birth masses); the value each age gets is the one `present_mass_sol` gives for it, to rounding (the table
    # of running sums adds in another order than the pass over every star did, which differs in the 14th digit).
    number, mean_mass, others = system_mass._placed_bins()
    fractions = system_mass.companion_fractions(mean_mass)
    further = system_mass.present_mass_sol(mean_mass[:, None] * fractions, age_gy).mean(axis=1)
    placed = float(np.sum(number * system_mass.present_mass_sol(mean_mass, age_gy)) / np.sum(number))
    assert system_mass.mean_system_mass_at_age_solar(age_gy) == pytest.approx(
        placed + float(np.sum(number * others * further) / np.sum(number)), rel=1e-12)
