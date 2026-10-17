"""
GEN.286 and GEN.287: the stellar evolution model has a horizontal branch and an asymptotic giant branch, and its red
giant branch spends its time at each luminosity as the real branch does (docs/design/realistic-galactic-core.md 5.3).

A star lighter than `HELIUM_FLASH_MASS_SOL` goes main sequence, subgiant, red giant, horizontal branch (core helium
burning; red, RR Lyrae, blue or extreme by [Fe/H] and age), asymptotic giant, white dwarf. The red giant branch rises
from its base to its tip with the time per unit of ln L falling as L^-0.78. The scatter's tables mirror those phases.
"""

import math
import random

import pytest

from planetgen import tuning
from planetgen.generation import multistar
from planetgen.generation import star_population as sp
from planetgen.generation.bright_stars import star_fe_h
from planetgen.physics import stellar_evolution as se

MASS = 0.95
"""A halo-age turn-off star: its main sequence lasts 11.4 Gy."""
T_MS = se.main_sequence_lifetime_gy(MASS)

FLUX_PER_YEAR_PER_LSUN = 1.7e-11
"""float: The specific evolutionary flux of an old population, stars leaving the main sequence each year for each
solar luminosity of the population's light: 1.7 +- 0.4e-11 (Buzzoni 1989, quoted by Renzini 1997, astro-ph/9711072)."""


def _state_at(fraction, **kwargs):
    return se.evolve_star(MASS, T_MS * fraction, random.Random(5), **kwargs)


def test_a_light_star_has_a_horizontal_branch_and_an_asymptotic_branch_between_its_giant_branch_and_its_white_dwarf():
    pc = tuning
    giant_end = pc.GIANT_PHASE_END_MS_FRACTION
    horizontal_end = giant_end + pc.HORIZONTAL_BRANCH_MS_FRACTION
    agb_end = horizontal_end + pc.AGB_MS_FRACTION
    sequence = [(giant_end - 0.01, "III", None), (giant_end + 0.1 * (horizontal_end - giant_end), "III", "horizontal-branch"),
                (horizontal_end + 0.5 * pc.AGB_MS_FRACTION, "III", "asymptotic-giant"), (agb_end + 0.001, "VII", None)]
    for fraction, yerkes, phase in sequence:
        state = _state_at(fraction)
        assert state["yerkes_class"] == yerkes and state.get("phase") == phase, fraction
    # Each phase ends where the next begins, and the star becomes a remnant at the end of the last one.
    assert _state_at(giant_end + 0.001)["phase_end_gy"] == pytest.approx(T_MS * horizontal_end)
    assert _state_at(horizontal_end + 0.001)["phase_end_gy"] == pytest.approx(T_MS * agb_end)
    assert _state_at(giant_end - 0.01)["lifespan_gy"] == pytest.approx(T_MS * agb_end)
    assert se.remnant_age_fraction(MASS) == pytest.approx(agb_end)
    # The branch lasts a hundred million years and the asymptotic giant branch twenty.
    assert 0.08 < T_MS * pc.HORIZONTAL_BRANCH_MS_FRACTION < 0.13
    assert 0.015 < T_MS * pc.AGB_MS_FRACTION < 0.03


def test_a_star_of_two_solar_masses_or_more_goes_straight_from_the_giant_branch_to_its_remnant():
    mass = tuning.HELIUM_FLASH_MASS_SOL
    t_ms = se.main_sequence_lifetime_gy(mass)
    rng = random.Random(1)
    assert not se.has_core_helium_phases(mass) and se.has_core_helium_phases(mass - 0.01)
    assert se.remnant_age_fraction(mass) == tuning.GIANT_PHASE_END_MS_FRACTION
    assert se.evolve_star(mass, t_ms * (tuning.GIANT_PHASE_END_MS_FRACTION - 0.001), rng)["yerkes_class"] == "III"
    assert se.evolve_star(mass, t_ms * (tuning.GIANT_PHASE_END_MS_FRACTION + 0.001), rng)["yerkes_class"] == "VII"


def test_the_colour_shares_follow_fe_h_and_sum_to_one():
    blue, rr_lyrae, red, extreme = zip(*(se.horizontal_branch_shares(fe_h, 12.0) for fe_h in (-2.2, -1.5, -0.1)))
    assert blue[0] > 0.75 > 0.3 > blue[2] and blue[0] > blue[1] > blue[2]
    assert max(rr_lyrae) == rr_lyrae[1] > 0.2 and rr_lyrae[2] == 0.0
    assert red[2] > 0.9 > 0.2 > red[0]
    assert all(0.0 < share < 0.15 for share in (extreme[0], extreme[1], extreme[2]))
    for fe_h in (-3.0, -2.2, -1.9, -1.0, -0.4, 0.5):
        assert sum(se.horizontal_branch_shares(fe_h, 12.0)) == pytest.approx(1.0)


def test_a_younger_branch_is_redder_at_the_same_fe_h():
    old, young = se.horizontal_branch_shares(-1.5, 13.0), se.horizontal_branch_shares(-1.5, 8.0)
    assert young[2] > old[2] and young[0] < old[0]


@pytest.mark.parametrize("fe_h", [-2.2, -1.5, -0.1])
def test_sampled_colours_match_the_shares(fe_h):
    rng = random.Random(7)
    draws = 4000
    counts = {colour: 0 for colour in se.HORIZONTAL_BRANCH_COLOURS}
    for _ in range(draws):
        counts[se.sample_horizontal_branch_colour(fe_h, 12.0, rng)] += 1
    for colour, share in zip(se.HORIZONTAL_BRANCH_COLOURS, se.horizontal_branch_shares(fe_h, 12.0)):
        sigma = math.sqrt(max(share * (1 - share), 1e-9) / draws)
        assert abs(counts[colour] / draws - share) < 4 * sigma + 1e-3, colour


def test_the_temperature_of_a_horizontal_branch_star_is_in_its_colours_range():
    rng = random.Random(9)
    for fe_h, colour in ((-2.2, "blue"), (-0.1, "red")):
        low, high = tuning.HORIZONTAL_BRANCH_COLOUR_TEMPERATURE_RANGE_K[colour]
        temperatures = [se.evolve_star(MASS, T_MS * 1.204, rng, fe_h=fe_h)["temperature_k"] for _ in range(400)]
        inside = sum(low <= value <= high for value in temperatures) / len(temperatures)
        assert inside > (0.7 if colour == "blue" else 0.9), (fe_h, inside)


def test_a_horizontal_branch_star_is_a_class_iii_giant_by_luminosity_and_says_so_in_its_type():
    state = _state_at(1.203, fe_h=-1.5)
    low, high = tuning.HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL
    assert state["yerkes_class"] == "III" and low <= state["luminosity_sol"] <= high
    params = se.star_params(MASS, T_MS * 1.203, state)
    assert "III " in params["type"] and "Horizontal Branch Giant Star" in params["type"]
    agb = se.star_params(MASS, T_MS * 1.2115, _state_at(1.2115))
    assert "Asymptotic Giant Star" in agb["type"] and agb["yerkes_class"] == "III"
    ordinary = se.star_params(MASS, T_MS * 1.15, _state_at(1.15))
    assert "III " in ordinary["type"] and "Giant Star" in ordinary["type"]
    assert "Horizontal" not in ordinary["type"] and "Asymptotic" not in ordinary["type"]


def test_a_horizontal_branch_star_without_an_fe_h_draws_one_and_keeps_it():
    drawn = _state_at(1.203, population="halo")
    assert "metallicity_fe_h" in drawn and -3.5 <= drawn["metallicity_fe_h"] <= 0.6
    given = _state_at(1.203, fe_h=-0.3)
    assert "metallicity_fe_h" not in given
    params = se.star_params(MASS, T_MS * 1.203, drawn)
    assert params["metallicity_fe_h"] == drawn["metallicity_fe_h"]
    # A scattered bright star stores the one its colour was drawn from, not a second draw.
    assert star_fe_h("halo", (0, 0, 0), 5, params) == params["metallicity_fe_h"]
    assert star_fe_h("cluster", (0, 0, 0), 5, params) is None


def test_no_horizontal_branch_star_is_brighter_than_the_top_of_its_range_or_a_floor_above_it():
    top = tuning.HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL[1]
    assert se.horizontal_branch_bright_chance(top + 1.0) == 0.0
    assert se.horizontal_branch_bright_chance(top - 5.0) > 0.0
    windows = sp._bright_windows(MASS, 100.0)
    horizontal_start = T_MS * tuning.GIANT_PHASE_END_MS_FRACTION
    assert all(not (start >= horizontal_start - 1e-9 and end <= horizontal_start + T_MS * tuning.HORIZONTAL_BRANCH_MS_FRACTION
                    + 1e-9) for start, end, _chance in windows)
    rng = random.Random(2)
    stars = sp.sample_bright_stars(300, 100.0, "halo", rng)
    assert not any("Horizontal" in star["type"] for star in stars)
    assert any("Asymptotic" in star["type"] for star in stars)


@pytest.mark.parametrize("population", ["halo", "old", "bulge"])
def test_a_horizontal_branch_star_sits_in_the_luminosity_range_below_every_scatter_floor(population):
    """The branch is 35 to 50 Lsun, so no scatter (floors of 50 Lsun and up) holds one; the star-by-star draws do."""
    low, high = tuning.HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL
    middle = tuning.GIANT_PHASE_END_MS_FRACTION + 0.5 * tuning.HORIZONTAL_BRANCH_MS_FRACTION
    rng = random.Random(4)
    for _ in range(50):
        state = se.evolve_star(MASS, T_MS * middle, rng, population=population)
        assert state["phase"] == "horizontal-branch" and low <= state["luminosity_sol"] <= high
    stars = sp.sample_bright_stars(600, 50.0, population, random.Random(4))
    assert not any("Horizontal" in star["type"] for star in stars)


def test_the_giant_branch_time_law_is_a_power_in_luminosity_between_its_base_and_its_tip():
    low, high, slope = se.giant_luminosity_law(MASS)
    assert low == pytest.approx(tuning.SUBGIANT_MAX_BRIGHTENING * se.main_sequence_luminosity_sol(MASS))
    assert high == tuning.GIANT_BRANCH_TIP_LUMINOSITY_SOL and slope == tuning.GIANT_BRANCH_TIME_SLOPE
    # The bright giants (4 solar masses and up) keep a log-uniform spread over their own range.
    assert se.giant_luminosity_law(5.0) == (*tuning.BRIGHT_GIANT_LUMINOSITY_RANGE_SOL, 0.0)
    # The share above a value falls monotonically from 1 at the base to 0 at the tip.
    values = [low * (high / low) ** (i / 10) for i in range(11)]
    shares = [se.giant_bright_chance(MASS, value) for value in values]
    assert shares[0] == 1.0 and shares[-1] == 0.0 and shares == sorted(shares, reverse=True)


def test_sampled_giant_luminosities_follow_the_law():
    rng = random.Random(3)
    draws = [se.sample_giant_luminosity_sol(MASS, rng) for _ in range(6000)]
    low, high, _slope = se.giant_luminosity_law(MASS)
    assert all(low <= value <= high for value in draws)
    for threshold in (10.0, 50.0, 200.0):
        expected = se.giant_bright_chance(MASS, threshold)
        found = sum(value >= threshold for value in draws) / len(draws)
        assert abs(found - expected) < 4 * math.sqrt(expected * (1 - expected) / len(draws)) + 1e-3
    # Most of the branch is spent near its base: the median giant is under 10 solar luminosities.
    assert sorted(draws)[len(draws) // 2] < 10.0
    above = [se.sample_giant_luminosity_sol(MASS, rng, 500.0) for _ in range(200)]
    assert min(above) >= 500.0 and max(above) <= high


def test_the_horizontal_branch_to_red_giants_above_it_ratio_is_the_observed_one():
    # Globular clusters have about 1.4 horizontal-branch stars for each red giant brighter than the branch
    # (Salaris et al. 2004) [R].
    level = sum(tuning.HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL) / 2.0
    above = (tuning.GIANT_PHASE_END_MS_FRACTION - tuning.SUBGIANT_PHASE_END_MS_FRACTION) * se.giant_bright_chance(MASS, level)
    ratio = tuning.HORIZONTAL_BRANCH_MS_FRACTION / above
    # The model sits a little under: its whole red giant branch (0.1 of the main-sequence lifetime, 1.1 Gy) is longer
    # than the observed branch, while the 114 Myr horizontal branch is the observed one.
    assert 0.8 < ratio < 2.0, ratio


@pytest.mark.parametrize("floor", [100.0, 500.0, 1000.0])
def test_the_halo_counts_match_the_evolutionary_flux_of_an_old_population(floor):
    """The number of stars at least `floor` bright for each solar luminosity of an old population's light is the
    evolutionary flux times the time a turn-off star spends above it, on the red giant and asymptotic giant branches
    (the horizontal branch is under 51 Lsun). The model's count must agree to within 35%."""
    rng = random.Random(11)
    total = sum(se.sample_living_star(rng=rng, population="halo")[2]["luminosity_sol"] for _ in range(30000))
    per_lsun = sp.bright_star_fraction(floor, "halo") / (total / 30000)
    giant_years = T_MS * 1.0e9 * (tuning.GIANT_PHASE_END_MS_FRACTION - tuning.SUBGIANT_PHASE_END_MS_FRACTION)
    years = giant_years * se.giant_bright_chance(MASS, floor) + (
        T_MS * 1.0e9 * tuning.AGB_MS_FRACTION * se.agb_bright_chance(floor))
    expected = FLUX_PER_YEAR_PER_LSUN * years
    assert expected * 0.65 < per_lsun < expected * 1.35, (per_lsun, expected)


def test_the_bright_fraction_still_matches_a_brute_force_draw_of_the_whole_population():
    """The scatter's table mirrors `evolve_star`: the share of living stars at least 50 Lsun bright, which includes
    the horizontal branch, equals the share a plain draw of the model gives."""
    rng = random.Random(21)
    draws = 150000
    found = sum(se.sample_living_star(rng=rng, population="old")[2]["luminosity_sol"] >= 50.0 for _ in range(draws))
    expected = sp.bright_star_fraction(50.0, "old")
    sigma = math.sqrt(expected * draws)
    assert abs(found - expected * draws) < 4 * sigma + 2


def test_the_youngest_progenitor_of_a_remnant_accounts_for_the_helium_burning_phases():
    for age_yr in (1.0e10, 1.3e10, 5.0e8):
        mass = multistar._youngest_progenitor_sol(age_yr)
        lifetime = se.main_sequence_lifetime_gy(mass) * se.remnant_age_fraction(mass) * 1.0e9
        assert lifetime == pytest.approx(age_yr, rel=1e-6)
