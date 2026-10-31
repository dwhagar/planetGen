"""
GEN.286 and GEN.287: the stellar evolution model has a horizontal branch and an asymptotic giant branch, and its red
giant branch spends its time at each luminosity as the real branch does (docs/design/galaxy/realistic-galactic-core.md 5.3).

A star lighter than `HELIUM_FLASH_MASS_SOL` goes main sequence, subgiant, red giant, horizontal branch (core helium
burning; red, RR Lyrae, blue or extreme by [Fe/H] and age), asymptotic giant, white dwarf. The red giant branch rises
from its base to its tip with the time per unit of ln L falling as L^-0.78 (the Hurley, Pols and Tout branch of
`giant_branch`, which the MIST track the star reads, GEN.935, agrees with). The scatter's tables mirror those phases.
"""

import math
import random

import pytest

from planetgen import tuning
from planetgen.generation import star_population as sp
from planetgen.generation.bright_stars import star_fe_h
from planetgen.physics import alpha_tracks, giant_branch, ms_lifetime, post_main_sequence
from planetgen.physics import stellar_evolution as se

MASS = 0.95
"""A halo-age turn-off star: its main sequence lasts 12.5 Gy."""
T_MS = se.main_sequence_lifetime_gy(MASS)
SUBGIANT_END = giant_branch.subgiant_end_fraction(MASS)
GIANT_END = giant_branch.giant_end_fraction(MASS)
HORIZONTAL_END = GIANT_END + tuning.HORIZONTAL_BRANCH_DURATION_GY / T_MS
"""float: The age each of these phases ends, as a multiple of the main-sequence lifetime (GEN.502: the giant branch ends
at the Hurley, Pols and Tout helium ignition, the horizontal branch lasts a fixed time)."""

SALARIS_R = 1.47
"""float: The weighted mean of the horizontal-branch to red-giant ratio R of 57 Galactic globular clusters, +-0.025 (its
standard error; the clusters scatter by 0.36 around it): Salaris et al. 2004, A&A 420, 911, Table 2
(astro-ph/0403600)."""

FLUX_PER_YEAR_PER_LSUN = 1.7e-11
"""float: The specific evolutionary flux of an old population, stars leaving the main sequence each year for each
solar luminosity of the population's light: 1.7 +- 0.4e-11 (Buzzoni 1989, quoted by Renzini 1997, astro-ph/9711072)."""


def _state_at(fraction, **kwargs):
    kwargs.setdefault("fe_h", 0.0)
    return se.evolve_star(MASS, T_MS * fraction, random.Random(5), **kwargs)


def _horizontal_branch_age(fe_h, share):
    """The age, Gy, `share` of the way through the horizontal branch of the test star of birth [Fe/H] `fe_h`
    (lifetimes follow the composition, GEN.651, on the track of the star's branch, GEN.1024: the piece of the age axis the age falls in)."""
    for start, end, track in alpha_tracks.track_segments(fe_h):
        t_ms = se.main_sequence_lifetime_gy(MASS, track)
        first = giant_branch.giant_end_fraction(MASS, track)
        age = t_ms * (first + share * tuning.HORIZONTAL_BRANCH_DURATION_GY / t_ms)
        if start <= age < end:
            return age
    raise AssertionError(fe_h)


def test_a_light_star_has_a_horizontal_branch_and_an_asymptotic_branch_between_its_giant_branch_and_its_white_dwarf():
    pc = tuning
    giant_end, horizontal_end = GIANT_END, HORIZONTAL_END
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
    assert 0.08 < pc.HORIZONTAL_BRANCH_DURATION_GY < 0.13
    assert 0.015 < T_MS * pc.AGB_MS_FRACTION < 0.03


def test_a_star_of_two_solar_masses_or_more_burns_its_helium_in_its_core_as_a_giant_and_not_on_a_horizontal_branch():
    mass = tuning.HELIUM_FLASH_MASS_SOL
    line = post_main_sequence.timeline(mass, 0.0)
    rng = random.Random(1)
    assert not se.has_core_helium_phases(mass) and se.has_core_helium_phases(mass - 0.01)
    assert [window.name for window in line.windows] == ["subgiant", "giant", "core-helium-burning", "asymptotic-giant"]
    assert se.remnant_age_fraction(mass) == pytest.approx(line.end_gy / line.main_sequence_gy)
    helium = line.window("core-helium-burning")
    state = se.evolve_star(mass, (helium.start_gy + helium.end_gy) / 2.0, rng, fe_h=0.0)
    assert state["yerkes_class"] == "III" and state["phase"] == "core-helium-burning"
    assert se.evolve_star(mass, line.end_gy * 1.001, rng, fe_h=0.0)["yerkes_class"] == "VII"


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
        age = _horizontal_branch_age(fe_h, 0.3)
        temperatures = [se.evolve_star(MASS, age, rng, fe_h=fe_h)["temperature_k"] for _ in range(400)]
        inside = sum(low <= value <= high for value in temperatures) / len(temperatures)
        assert inside > (0.6 if colour == "blue" else 0.9), (fe_h, inside)


def test_a_horizontal_branch_star_is_a_class_iii_giant_by_luminosity_and_says_so_in_its_type():
    early_age = _horizontal_branch_age(-1.5, 0.2)
    state = se.evolve_star(MASS, early_age, random.Random(5), fe_h=-1.5)
    low, high = tuning.HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL
    assert state["yerkes_class"] == "III" and low <= state["luminosity_sol"] <= high
    params = se.star_params(MASS, early_age, state)
    assert "III " in params["type"] and "Horizontal Branch Giant Star" in params["type"]
    asymptotic = HORIZONTAL_END + 0.5 * tuning.AGB_MS_FRACTION
    agb = se.star_params(MASS, T_MS * asymptotic, _state_at(asymptotic))
    assert "Asymptotic Giant Star" in agb["type"] and agb["yerkes_class"] == "III"
    ordinary_age = (SUBGIANT_END + GIANT_END) / 2.0
    ordinary = se.star_params(MASS, T_MS * ordinary_age, _state_at(ordinary_age))
    assert "III " in ordinary["type"] and "Giant Star" in ordinary["type"]
    assert "Horizontal" not in ordinary["type"] and "Asymptotic" not in ordinary["type"]


def test_a_star_without_an_fe_h_draws_one_first_and_keeps_it():
    drawn = se.evolve_star(MASS, 5.0, random.Random(5), population="halo")
    assert -3.5 <= drawn["metallicity_fe_h"] <= 0.6
    given = _state_at(0.2, fe_h=-0.3)
    assert given["metallicity_fe_h"] == -0.3
    params = se.star_params(MASS, 5.0, drawn)
    assert params["metallicity_fe_h"] == drawn["metallicity_fe_h"]
    # A scattered bright star stores the one its colour was drawn from, not a second draw.
    assert star_fe_h(params) == params["metallicity_fe_h"]


def test_no_horizontal_branch_star_is_brighter_than_the_top_of_its_range_or_a_floor_above_it():
    top = tuning.HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL[1]
    assert se.horizontal_branch_bright_chance(top + 1.0) == 0.0
    assert se.horizontal_branch_bright_chance(top - 5.0) > 0.0
    windows = sp._bright_windows(MASS, 100.0)
    horizontal_start = T_MS * GIANT_END
    assert all(not (start >= horizontal_start - 1e-9 and end <= horizontal_start + tuning.HORIZONTAL_BRANCH_DURATION_GY
                    + 1e-9) for start, end, _chance in windows)
    rng = random.Random(2)
    stars = sp.sample_bright_stars(300, 100.0, "halo", rng)
    assert not any("Horizontal" in star["type"] for star in stars)
    assert any("Asymptotic" in star["type"] for star in stars)


@pytest.mark.parametrize("population", ["halo", "old", "bulge"])
def test_a_horizontal_branch_star_sits_in_the_luminosity_range_below_every_scatter_floor(population):
    """The branch is 35 to 50 Lsun, so no scatter (floors of 50 Lsun and up) holds one; the star-by-star draws do."""
    low, high = tuning.HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL
    middle = (GIANT_END + HORIZONTAL_END) / 2.0
    rng = random.Random(4)
    for _ in range(50):
        state = se.evolve_star(MASS, T_MS * middle, rng, population=population, fe_h=0.0)
        assert state["phase"] == "horizontal-branch" and low <= state["luminosity_sol"] <= high
    stars = sp.sample_bright_stars(600, 50.0, population, random.Random(4))
    assert not any("Horizontal" in star["type"] for star in stars)


def test_the_giant_branch_is_the_hurley_pols_tout_branch_of_the_star():
    """Table 1 of Hurley, Pols and Tout 2000 gives the branch's duration, (t_HeI - t_BGB) / t_BGB, as 6.4e-2 at one solar
    mass, 2.0e-2 at two and 2.4e-3 at five (the two heavier ones only through the same formulae, which the module
    shares with the light star's)."""
    branch = giant_branch.branch(1.0)
    assert branch.duration_myr / branch.start_myr == pytest.approx(6.4e-2, rel=0.05)
    assert giant_branch.branch(2.0).duration_myr / giant_branch.branch(2.0).start_myr == pytest.approx(2.0e-2, rel=0.1)
    # The branch of the test star runs from its base to the tip of equation 49, and the star's windows are its times.
    branch = giant_branch.branch(MASS)
    assert branch.tip == pytest.approx(giant_branch.tip_luminosity_sol(MASS)) and 2000.0 < branch.tip < 3000.0
    assert branch.start_myr == pytest.approx(T_MS * 1.0e3 * SUBGIANT_END)
    assert branch.end_myr == pytest.approx(T_MS * 1.0e3 * GIANT_END)
    line = post_main_sequence.timeline(MASS, 0.0)
    assert line.window("subgiant").end_gy == pytest.approx(T_MS * SUBGIANT_END)
    assert line.window("giant").end_gy == pytest.approx(T_MS * GIANT_END)
    # The share above a value falls monotonically from 1 at the base to 0 at the tip.
    values = [branch.base * (branch.tip / branch.base) ** (i / 10) for i in range(11)]
    shares = [branch.share_above(value) for value in values]
    assert shares[0] == 1.0 and shares[-1] == 0.0 and shares == sorted(shares, reverse=True)
    # The star climbs quickly at the end: the time per unit of ln L falls from about L^-0.8 at the base to L^-0.6 at the tip.
    # Compare the time spent in a factor of 1.2 in luminosity near the base and near the tip.
    near_base = branch.share_above(branch.base * 2.0) - branch.share_above(branch.base * 2.4)
    near_tip = branch.share_above(branch.tip / 2.4) - branch.share_above(branch.tip / 2.0)
    assert near_base > 5.0 * near_tip > 0.0


def test_the_tip_is_lower_at_low_metallicity_as_the_models_of_a_metal_poor_star_have_it():
    """Serenelli et al. 2017 (arXiv:1706.07564) find the tip of an old, metal-poor star at 1,650 to 1,780 Lsun; the
    paper's fit lands within 20% of them at Z = 1e-4 ([Fe/H] = -2.19 against the reference Z of 0.0154)."""
    tip = giant_branch.tip_luminosity_sol(0.8, fe_h=math.log10(1.0e-4 / tuning.MS_LIFETIME_REFERENCE_Z))
    assert 1500.0 < tip < 2100.0
    assert tip < giant_branch.tip_luminosity_sol(0.8, fe_h=0.0)


def test_the_time_a_star_spends_above_a_luminosity_on_the_mist_branch_follows_the_hurley_branch():
    """The star reads its luminosity from the MIST track and its times from the Hurley, Pols and Tout branch; the share of the
    branch above a luminosity must agree to within a factor of about 1.6 (the MIST and the fitted branch differ in the base and
    the tip by up to 30%) at 10, 50 and 200 Lsun."""
    line = post_main_sequence.timeline(MASS, 0.0)
    branch = giant_branch.branch(MASS)
    giant = line.window("giant")
    for threshold in (10.0, 50.0, 200.0):
        time_above = sum(end - start for name, start, end in line.bright_pieces(threshold) if name == "giant")
        mist_share = time_above / (giant.end_gy - giant.start_gy)
        assert branch.share_above(threshold) / 1.6 < mist_share < branch.share_above(threshold) * 1.6, threshold


def test_the_horizontal_branch_to_red_giants_above_it_ratio_is_the_observed_one():
    """The horizontal branch lasts 100 Myr and the red giant branch spends 69 Myr above its luminosity, a ratio of
    1.45 for every star that ignites its helium in a flash and whose base of the giant branch lies below it; the 57
    Galactic globular clusters of Salaris et al. 2004 have a weighted mean of 1.47 +- 0.025."""
    level = sum(tuning.HORIZONTAL_BRANCH_LUMINOSITY_RANGE_SOL) / 2.0
    for mass in (0.8, MASS, 1.0):
        branch = giant_branch.branch(mass)
        above_myr = branch.share_above(level) * branch.duration_myr
        ratio = tuning.HORIZONTAL_BRANCH_DURATION_GY * 1.0e3 / above_myr
        assert abs(ratio - SALARIS_R) < 0.05, (mass, ratio)


HALO_FE_H = sp.reference_fe_h("halo")
"""float: The composition `bright_star_fraction` counts the halo at (-1.6 since GEN.651/GEN.662; the test once worked at the Sun's)."""


@pytest.fixture(scope="module")
def halo_mean_luminosity():
    """The mean luminosity of a halo star at the reference composition, from 120,000 draws: a few giants carry most of it,
    so a smaller sample moves it by a third (a mean of 40,000 draws has a spread of 25%, of 120,000 about 14%)."""
    total = 0.0
    for seed in (11, 12):
        rng = random.Random(seed)
        total += sum(se.sample_living_star(rng=rng, population="halo", fe_h=HALO_FE_H)[2]["luminosity_sol"] for _ in range(60000))
    return total / 120000


@pytest.mark.parametrize("floor", [100.0, 500.0, 1000.0])
def test_the_halo_counts_match_the_evolutionary_flux_of_an_old_population(floor, halo_mean_luminosity):
    """The number of stars at least `floor` bright for each solar luminosity of an old population's light is the
    evolutionary flux times the time a turn-off star spends above it, on the red giant and asymptotic giant branches
    (the horizontal branch is under 51 Lsun). The model's count must agree to within 35%."""
    # All at the reference composition, the one `bright_star_fraction` counts: the mean luminosity above, and the turn-off star, the mass
    # whose main sequence lasts as long as MASS's does at the Sun's composition (12.5 Gy; 0.82 solar masses at [Fe/H] -1.6). The
    # test took the Sun's composition for the mean and the track while the count moved to the halo's, which left the model's count
    # 1.2 to 1.5 times the flux's, and a mean of too few draws moved it as much again (TEST.481).
    per_lsun = sp.bright_star_fraction(floor, "halo") / halo_mean_luminosity
    track = alpha_tracks.track_fe_h(HALO_FE_H, T_MS)      # GEN.1024: the halo is alpha-enhanced and reads the track of its branch
    turnoff = ms_lifetime.mass_for_lifetime_gy(T_MS, track)
    line = post_main_sequence.timeline(turnoff, track)
    years = 1.0e9 * sum(end - start for _, start, end in line.bright_pieces(floor, skip=("horizontal-branch",)))
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
        mass = se.youngest_progenitor_sol(age_yr)
        lifetime = se.main_sequence_lifetime_gy(mass) * se.remnant_age_fraction(mass) * 1.0e9
        assert lifetime == pytest.approx(age_yr, rel=1e-6)


@pytest.mark.parametrize("fe_h", [-2.3, -0.5, 0.0, 0.3])
def test_the_metallicity_terms_worked_out_once_give_the_luminosities_the_plain_polynomials_do(fe_h):
    """PERF.146: `a27 .. a32` and `b9 .. b13` depend on [Fe/H] alone and are worked out once; the base luminosity is
    what the polynomials worked out each call give, bit for bit."""
    zeta = giant_branch.ms_lifetime.zeta(fe_h)
    a27, a28, a29p, a30, a31, a32 = (giant_branch._poly(giant_branch._BASE[n], zeta) for n in (27, 28, 29, 30, 31, 32))
    b9, b10 = (giant_branch._poly(giant_branch._TIP[n], zeta) for n in (9, 10))
    b11, b13 = giant_branch._poly(giant_branch._TIP[11], zeta) ** 2, giant_branch._poly(giant_branch._TIP[13], zeta) ** 2
    b12 = giant_branch._poly(giant_branch._TIP[12], zeta)
    assert giant_branch._base_terms(zeta) == (a27, a28, a29p, a30, a31, a32)
    assert giant_branch._tip_terms(zeta) == (b9, b10, b11, b12, b13)
    for m in (0.8, 1.5, 2.4, 5.0, 20.0):
        assert giant_branch.base_luminosity_sol(m, fe_h) == (
            (a27 * m ** a31 + a28 * m ** giant_branch._BASE_C2) / (a29p ** a32 + a30 * m ** giant_branch._BASE_C3 + m ** a32))
