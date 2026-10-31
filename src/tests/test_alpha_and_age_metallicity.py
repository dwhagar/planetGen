# tests/test_alpha_and_age_metallicity.py

"""GEN.744: the alpha-enhanced equivalent metallicity (`physics/alpha_tracks.py`) and the age-metallicity relation (`physics/age_metallicity.py`)."""

import random

import pytest

from planetgen import tuning
from planetgen.physics import age_metallicity, alpha_tracks, stellar_abundances

THIN, THICK = stellar_abundances.BRANCH_THIN, stellar_abundances.BRANCH_HIGH_ALPHA


def test_a_solar_scaled_star_reads_its_own_metallicity():
    assert alpha_tracks.metallicity_shift(0.0) == pytest.approx(0.0, abs=1e-12)


def test_the_salaris_shift_for_plus_four_tenths_is_about_three_tenths_of_a_dex():
    """Salaris et al. 1993: [alpha/Fe] = +0.4 is worth +0.29 dex in [M/H]."""
    assert alpha_tracks.metallicity_shift(0.4) == pytest.approx(0.293, abs=0.002)
    assert alpha_tracks.metallicity_shift(0.2) < alpha_tracks.metallicity_shift(0.4)


def test_the_high_alpha_branch_is_more_enhanced_than_the_thin_one_and_plateaus_when_metal_poor():
    assert alpha_tracks.alpha_fe(-0.5, THICK) > alpha_tracks.alpha_fe(-0.5, THIN)
    assert 0.3 < alpha_tracks.alpha_fe(-1.5, THICK) < 0.45
    assert alpha_tracks.alpha_fe(-1.5, THICK) == pytest.approx(alpha_tracks.alpha_fe(-1.2, THICK))
    assert alpha_tracks.alpha_fe(0.0, THICK) < 0.1


def test_a_thick_disc_star_reads_a_track_richer_in_iron_than_itself_and_a_thin_one_nearly_its_own():
    assert alpha_tracks.equivalent_fe_h(-0.8, THICK) == pytest.approx(-0.8 + 0.27, abs=0.04)
    assert alpha_tracks.equivalent_fe_h(0.0, THIN) == pytest.approx(0.0, abs=0.02)
    assert alpha_tracks.equivalent_fe_h(-0.8, THICK) > alpha_tracks.equivalent_fe_h(-0.8, THIN)


def test_the_thin_disc_gains_metals_at_the_slope_haywood_gives_and_the_thick_disc_six_times_faster():
    thin = age_metallicity.mean_fe_h(2.0, THIN) - age_metallicity.mean_fe_h(6.0, THIN)
    assert thin / 4.0 == pytest.approx(0.025)
    thick = age_metallicity.mean_fe_h(9.0, THICK) - age_metallicity.mean_fe_h(12.0, THICK)
    assert thick / 3.0 == pytest.approx(0.15)
    assert thick / 3.0 / (thin / 4.0) == pytest.approx(6.0)


def test_the_thick_disc_ends_at_its_anchor_and_is_held_outside_its_window():
    assert age_metallicity.mean_fe_h(tuning.STELLAR_THICK_START_GY, THICK) == pytest.approx(age_metallicity.THICK_ANCHOR_FE_H)
    assert age_metallicity.mean_fe_h(5.0, THICK) == age_metallicity.mean_fe_h(8.0, THICK)
    assert age_metallicity.mean_fe_h(20.0, THICK) == age_metallicity.mean_fe_h(13.0, THICK)
    assert -0.9 < age_metallicity.mean_fe_h(12.0, THICK) < -0.5


def test_the_relation_matches_the_populations_the_generator_has_today():
    """The young and intermediate disc average 0 in `tuning.METALLICITY_BY_POPULATION`, and the old (3 to 13 Gyr, the thick disc's 8 to 13 among them) its offset."""
    assert age_metallicity.mean_fe_h(1.5, THIN) == pytest.approx(0.0, abs=0.06)
    assert age_metallicity.mean_fe_h(6.5, THIN) == pytest.approx(-0.2, abs=0.1)
    assert tuning.METALLICITY_BY_POPULATION["old"]["offset"] == pytest.approx(-0.325)


def test_only_the_discs_populations_follow_the_relation():
    for population in ("young", "intermediate", "old"):
        assert age_metallicity.takes(population)
    for population in (None, "halo", "bulge", "cluster", "satellite", "open_cluster_03"):
        assert not age_metallicity.takes(population)


def test_the_thick_disc_exists_only_past_the_gate():
    assert age_metallicity.thick_share(7.9) == 0.0
    assert age_metallicity.thick_share(8.0) == age_metallicity.THICK_SHARE


def test_a_metal_poor_star_is_likelier_old_and_a_metal_rich_one_young():
    assert age_metallicity.log_density(-0.5, 9.0) > age_metallicity.log_density(-0.5, 1.0)
    assert age_metallicity.log_density(0.1, 1.0) > age_metallicity.log_density(0.1, 9.0)


def test_the_density_never_underflows_far_from_both_branches():
    value = age_metallicity.log_density(-3.0, 10.0)
    assert value < -100.0 and value > -1.0e6


def test_a_population_outside_the_disc_draws_one_age_and_no_more():
    rng = random.Random(1)
    calls = []
    age = age_metallicity.tilted_age(-1.5, "halo", lambda: calls.append(1) or 11.0, rng)
    assert age == 11.0 and len(calls) == 1


def test_the_tilt_moves_the_old_disc_ages_with_the_metallicity():
    rng = random.Random(7)

    def mean_age(fe_h):
        ages = [age_metallicity.tilted_age(fe_h, "old", lambda: rng.uniform(3.0, 10.0), rng) for _ in range(4000)]
        assert all(3.0 <= age <= 10.0 for age in ages)
        return sum(ages) / len(ages)

    poor, rich = mean_age(-0.5), mean_age(0.1)
    assert poor > 6.5 > rich
    assert poor - rich > 1.0


def test_the_tilt_leaves_a_flat_window_alone_where_the_relation_is_flat():
    """A 0.1 Gyr window changes the relation by 0.0025 dex: the mean stays at the window's middle."""
    rng = random.Random(3)
    ages = [age_metallicity.tilted_age(0.0, "young", lambda: rng.uniform(0.0, 0.1), rng) for _ in range(4000)]
    assert sum(ages) / len(ages) == pytest.approx(0.05, abs=0.005)


def test_an_unknown_branch_is_refused():
    with pytest.raises(ValueError):
        age_metallicity.mean_fe_h(4.0, "halo")


def test_a_drawn_old_disc_star_is_older_the_poorer_it_is():
    """The generator's own draw (`sample_living_star` with the birth [Fe/H] fixed) follows the tilt."""
    from planetgen.physics import stellar_evolution
    rng = random.Random(5)

    def mean_age(fe_h):
        ages = [stellar_evolution.sample_living_star(rng=rng, population="old", fe_h=fe_h)[1] for _ in range(150)]
        return sum(ages) / len(ages)

    assert mean_age(-0.5) - mean_age(0.1) > 1.5


def test_the_tilt_leaves_a_population_s_mean_age_where_it_was():
    """Each disc population's [Fe/H] is about the relation's mix over its window, so tilting the ages moves them with the metals and not as a whole."""
    from planetgen.physics import stellar_evolution
    rng = random.Random(11)
    for population in ("intermediate", "old"):
        low, high = stellar_evolution.population_age_range_gy(population)
        plain, tilted = [], []
        for _ in range(6000):
            fe_h = stellar_evolution.birth_fe_h(rng, population)
            plain.append(rng.uniform(low, high))
            tilted.append(age_metallicity.tilted_age(fe_h, population, lambda: rng.uniform(low, high), rng))
        assert sum(tilted) / len(tilted) == pytest.approx(sum(plain) / len(plain), abs=0.15)


def test_a_star_has_the_branch_its_birth_feh_and_age_give_and_keeps_it():
    """GEN.1024: the branch is a function of the two numbers a stored star keeps."""
    assert stellar_abundances.branch_of(-0.8, 10.0) == stellar_abundances.branch_of(-0.8, 10.0)
    assert all(stellar_abundances.branch_of(fe_h / 10.0, 7.9) == THIN for fe_h in range(-15, 5))


def test_the_branches_follow_the_logistic_in_feh_among_the_old():
    rng = random.Random(4)
    for fe_h, low, high in ((-1.2, 0.97, 1.0), (-0.35, 0.4, 0.6), (0.1, 0.0, 0.1)):
        share = sum(stellar_abundances.branch_of(fe_h + rng.uniform(-0.02, 0.02), 11.0) == THICK for _ in range(2000)) / 2000
        assert low <= share <= high, (fe_h, share)


def test_a_thick_disc_star_reads_a_richer_track_than_a_thin_one_of_the_same_feh():
    old = [fe_h for fe_h in (-0.9 - 0.001 * n for n in range(200)) if stellar_abundances.branch_of(fe_h, 11.0) == THICK]
    assert old
    for fe_h in old[:5]:
        assert alpha_tracks.track_fe_h(fe_h, 11.0) > alpha_tracks.track_fe_h(fe_h, 3.0) + 0.05    # a younger star is thin-disc at any metallicity
        assert alpha_tracks.track_fe_h(fe_h, 3.0) == alpha_tracks.equivalent_fe_h(fe_h, THIN)


def test_a_star_is_read_off_the_track_of_its_branch():
    """The lifetime and the state `evolve_star` gives are those of the equivalent [Fe/H], not of the star's own."""
    from planetgen.physics import post_main_sequence, stellar_evolution
    fe_h = next(f for f in (-1.0 - 0.001 * n for n in range(300)) if stellar_abundances.branch_of(f, 11.0) == THICK)
    state = stellar_evolution.evolve_star(0.9, 11.0, random.Random(1), fe_h=fe_h)
    assert state["metallicity_fe_h"] == fe_h
    assert state["lifespan_gy"] == pytest.approx(post_main_sequence.timeline(0.9, alpha_tracks.track_fe_h(fe_h, 11.0)).end_gy)
    assert state["lifespan_gy"] != pytest.approx(post_main_sequence.timeline(0.9, fe_h).end_gy, rel=1e-4)


def test_a_thin_disc_star_of_solar_feh_sits_on_the_reference_track():
    """The grid is solar-scaled: the fits' own offset at the Sun (-0.017 dex) is the zero of the shift, so the Sun's [Fe/H] is read as it is."""
    assert alpha_tracks.track_fe_h(0.0, 4.5) == pytest.approx(0.0, abs=1e-12)
    assert alpha_tracks.track_fe_h(0.0, 4.5) == alpha_tracks.equivalent_fe_h(0.0, THIN)


def test_the_track_pieces_of_a_star_are_the_tracks_it_reads_at_every_age():
    """GEN.1024: `track_segments` is what the scatter's windows are drawn on, so it must agree with `track_fe_h` at any age."""
    thick = next(f for f in (-1.0 - 0.001 * n for n in range(300)) if stellar_abundances.branch_of(f, 11.0) == THICK)
    for fe_h in (-1.2, -0.5, 0.0, 0.25, thick):
        pieces = alpha_tracks.track_segments(fe_h)
        assert pieces[0][0] == 0.0 and pieces[-1][1] == float("inf")
        assert all(a[1] == b[0] for a, b in zip(pieces, pieces[1:]))
        for age in (0.01, 1.0, 4.5, 7.9, 8.1, 11.0, 13.5):
            track = next(track for start, end, track in pieces if start <= age < end)
            assert track == pytest.approx(alpha_tracks.track_fe_h(fe_h, age), abs=1e-12)
