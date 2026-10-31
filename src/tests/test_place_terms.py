# tests/test_place_terms.py -- GEN.103: the close-companion factor that follows a star's [Fe/H], the widest axis a wide
# pair keeps at its place, and how `StarSystem._companion` uses them.

import math
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.galaxy import system_mass as galaxy_system_mass
from planetgen.generation import close_binary_trend, companion_classes, mass_ratio, place_terms, system_mass
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import imf, remnant_progenitors
from planetgen.util import draw
from tests.toy_shape import toy_shape

POINTS = sorted(tuning.CLOSE_BINARY_FRACTION_BY_FE_H.items())


def test_the_factor_is_the_pooled_line_over_its_solar_value():
    f0 = 0.20  # the line at [Fe/H] 0: between the -0.2 and +0.5 points
    assert close_binary_trend.close_companion_factor(0.0) == pytest.approx(1.0)
    for fe_h, fraction in POINTS:
        assert close_binary_trend.close_companion_factor(fe_h) == pytest.approx(fraction / f0)
    # Flat past the ends, and 1 for a star with no [Fe/H].
    assert close_binary_trend.close_companion_factor(-9.0) == close_binary_trend.close_companion_factor(-3.0)
    assert close_binary_trend.close_companion_factor(3.0) == close_binary_trend.close_companion_factor(0.5)
    assert close_binary_trend.close_companion_factor(None) == 1.0
    # Straight between two points.
    assert close_binary_trend.close_companion_factor(-0.6) == pytest.approx((0.40 + 0.24) / 2 / f0)


def test_the_factor_moves_the_odds_and_never_passes_one():
    # Boss, 7:19 PM EDT: the odds p / (1 - p), not p.
    assert close_binary_trend.scaled_close_chance(0.22, 0.0) == pytest.approx(0.22)
    for chance in (0.05, 0.22, 0.45, 0.9, 0.99):
        for fe_h in (-5.0, -3.0, -1.0, -0.2, 0.0, 0.5, 2.0):
            moved = close_binary_trend.scaled_close_chance(chance, fe_h)
            assert 0.0 < moved < 1.0
            odds = chance / (1 - chance) * close_binary_trend.close_companion_factor(fe_h)
            assert moved / (1 - moved) == pytest.approx(odds)
    assert close_binary_trend.scaled_close_chance(0.0, -3.0) == 0.0 and close_binary_trend.scaled_close_chance(1.0, 0.5) == 1.0
    # A star already likely to have a companion rises less in absolute terms.
    assert (close_binary_trend.scaled_close_chance(0.45, -1.0) - 0.45) < (1 - close_binary_trend.scaled_close_chance(0.1, -1.0)) * 10
    assert close_binary_trend.scaled_close_chance(0.45, -0.7) / 0.45 == pytest.approx(1.29, abs=0.02)   # inside "within about 30%"


def test_the_factor_falls_as_metallicity_rises():
    values = [close_binary_trend.close_companion_factor(x / 10) for x in range(-40, 10)]
    assert all(a >= b for a, b in zip(values, values[1:]))


def test_the_encounter_limit_reproduces_jiang_and_tremaines_worked_value():
    # 3.1e4 AU for M1 + M2 = 2, sigma = 50 km/s, rho_2 = 0.03 and t = 10 Gyr (design note 3.3). Their sigma is the
    # dispersion of the *relative* velocity, root 2 times the field's one-dimensional dispersion (their section 5:
    # 29 km/s becoming 40 km/s), and rho_2 the mass density times the mass-weighted perturber mass (GEN.447).
    one_dimensional = 50.0 / math.sqrt(2.0)
    terms = place_terms.PlaceTerms(mass_density_solar_pc3=0.075, dispersion_kms=one_dimensional,
                                   perturber_mass_solar=0.4)
    assert place_terms.encounter_limit_au(terms, 2.0, 10.0) == pytest.approx(3.1e4, rel=0.01)
    assert place_terms.RELATIVE_VELOCITY_FACTOR == pytest.approx(math.sqrt(2.0))
    # Jiang and Tremaine's own solar neighbourhood: 29 km/s one-dimensional, rho_2 0.029 (their equation 47).
    local = place_terms.PlaceTerms(mass_density_solar_pc3=0.029 / 0.4, dispersion_kms=29.0, perturber_mass_solar=0.4)
    assert place_terms.encounter_limit_au(local, 2.0, 10.0) == pytest.approx(
        3.1e4 * (29.0 * math.sqrt(2.0) / 50.0) * (0.03 / 0.029), rel=0.01)


def test_the_mass_weighted_perturber_mass_of_the_stars_matches_jiang_and_tremaines_neighbourhood():
    # Their equation 47: rho_2 / rho_1 = 0.029 / 0.045 = 0.64 solar masses, from the observed stars alone. The
    # generator's initial mass function, aged over the disk's star formation history with no remnants, gives 0.59.
    disk = [(1.0, tuning.STAR_FORMATION_AGE_RANGE_GY)]
    stars = place_terms.mass_weighted_perturber_mass_solar(disk, remnants=False)
    assert stars == pytest.approx(0.029 / 0.045, rel=0.15)


def test_the_remnants_raise_the_perturber_mass_and_the_black_holes_do_most_of_it():
    disk = [(1.0, tuning.STAR_FORMATION_AGE_RANGE_GY)]
    stars = place_terms.mass_weighted_perturber_mass_solar(disk, remnants=False)
    everything = place_terms.mass_weighted_perturber_mass_solar(disk)
    assert 1.3 < everything / stars < 2.5
    assert place_terms.PlaceTerms(1.0, 30.0).perturber_mass_solar == pytest.approx(everything)


def test_the_black_hole_inputs_lie_inside_the_published_ranges():
    """GEN.541 (docs/design/generation/placement-by-population.md): the share of core-collapse stars that leave a hole
    (the failed share, GEN.554: inside the verified 0 to 45 percent, Sukhbold et al. 2016, Table 4), the mean hole mass (Sukhbold et al.: 9 to 14 solar masses;
    Olejak et al. 2020: 11.7 bulge, 14.2 disk) and the number in a galaxy of 5e10 solar masses of stars (Olejak et al.:
    1.3e8; the older estimates they list: 1e8 to 1e9)."""
    low, high = tuning.CORE_COLLAPSE_MASS_RANGE_SOL
    collapse = imf.integral(0.0, low, high)
    holes = remnant_progenitors.number("black-hole")
    assert holes / collapse == pytest.approx(remnant_progenitors.reference_failed_share(), rel=1e-4)
    assert 0.0 <= holes / collapse <= 0.45
    assert 9.0 <= sum(tuning.BLACK_HOLE_MASS_RANGE_SOLAR) / 2.0 <= 14.2
    per_solar_mass_formed = holes / imf.integral(1.0, tuning.IMF_BREAKS_SOL[0], high)
    widest = remnant_progenitors.number("black-hole", 0.0, tuning.BLACK_HOLE_FAILED_SHARE_SOLAR_RANGE[1])
    assert 1.0e8 <= 5.0e10 * widest / imf.integral(1.0, tuning.IMF_BREAKS_SOL[0], high) <= 1.0e9
    assert 1.0e8 <= 5.0e10 * per_solar_mass_formed < 1.0e9


def test_the_black_holes_carry_the_share_of_the_second_moment_the_design_note_records():
    disk = [(1.0, tuning.STAR_FORMATION_AGE_RANGE_GY)]
    stars = place_terms.mass_weighted_perturber_mass_solar(disk, remnants=False)
    everything = place_terms.mass_weighted_perturber_mass_solar(disk)
    assert stars == pytest.approx(0.57, abs=0.06)
    assert everything == pytest.approx(1.19, abs=0.06)


def test_a_young_population_has_the_heavier_perturbers():
    young = place_terms.mass_weighted_perturber_mass_solar([(1.0, tuning.STELLAR_POPULATION_AGE_RANGES_GY["young"])])
    old = place_terms.mass_weighted_perturber_mass_solar([(1.0, tuning.STELLAR_POPULATION_AGE_RANGES_GY["old"])])
    assert young > old > 0.5
    mixed = place_terms.mass_weighted_perturber_mass_solar([(0.5, (0.0, 0.1)), (0.5, (3.0, 10.0))])
    assert old < mixed < young


def test_a_cell_with_no_stars_falls_back_to_the_mean_system_mass_of_the_disk():
    assert place_terms.mass_weighted_perturber_mass_solar([]) == system_mass.window_mean_system_mass_solar(
        *tuning.STAR_FORMATION_AGE_RANGE_GY)


def test_the_limits_scale_the_way_the_formulas_do():
    terms = place_terms.PlaceTerms(0.05, 30.0)
    base = place_terms.encounter_limit_au(terms, 2.0, 5.0)
    assert place_terms.encounter_limit_au(terms, 2.0, 10.0) == pytest.approx(base / 2)           # 1 / t
    assert place_terms.encounter_limit_au(terms, 4.0, 5.0) == pytest.approx(base * 2)            # M
    assert place_terms.encounter_limit_au(place_terms.PlaceTerms(0.1, 30.0), 2.0, 5.0) == pytest.approx(base / 2)
    tide = place_terms.tidal_limit_au(terms, 2.0)
    assert place_terms.tidal_limit_au(place_terms.PlaceTerms(0.4, 30.0), 2.0) == pytest.approx(tide / 2)  # rho^(-1/3)
    assert place_terms.wide_axis_limit_au(terms, 2.0, 5.0) == min(base, tide)


def test_the_design_notes_table_follows():
    # docs/design/generation/placement-by-population.md 3.3, 2 solar masses, 10 Gyr, the default perturber mass (<M^2> / <M> of the
    # whole disk with remnants).
    for density, sigma, encounter, tide in ((0.048, 30, 13_900, 496_000), (6.9, 120, 387, 94_600),
                                            (77, 120, 35, 42_400), (660, 120, 4.0, 20_700)):
        terms = place_terms.PlaceTerms(density, float(sigma))
        assert place_terms.encounter_limit_au(terms, 2.0, 10.0) == pytest.approx(encounter, rel=0.1)
        assert place_terms.tidal_limit_au(terms, 2.0) == pytest.approx(tide, rel=0.1)


def test_no_density_or_no_age_leaves_nothing_to_cut_the_axis():
    zero = place_terms.PlaceTerms(0.0, 30.0)
    assert math.isinf(place_terms.wide_axis_limit_au(zero, 2.0, 10.0))
    assert math.isinf(place_terms.encounter_limit_au(place_terms.PlaceTerms(0.1, 30.0), 2.0, 0.0))


def test_terms_from_a_profile_read_its_density_and_the_weighted_dispersion():
    def pop(share, sigma, age_gy=(3.0, 10.0)):
        return SimpleNamespace(share=share, dispersion_kms=sigma, age_gy=age_gy)

    shape = toy_shape(disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
                      bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    profile = SimpleNamespace(expected_systems=50.0, volume_pc3=100.0, shape=shape,
                              populations={"old": pop(0.75, 20.0), "bulge": pop(0.25, 60.0)})
    terms = place_terms.from_profile(profile)
    mean_mass = galaxy_system_mass.mean_over_solar({"old": 0.75, "bulge": 0.25}, shape)
    assert terms.mass_density_solar_pc3 == pytest.approx(0.5 * mean_mass)
    # The profile's dispersions are three-axis totals; the terms hold the one-axis value (GEN.284).
    assert terms.dispersion_kms == pytest.approx(math.sqrt(0.75 * 400 + 0.25 * 3600) / math.sqrt(3.0))
    assert terms.perturber_mass_solar == pytest.approx(
        place_terms.mass_weighted_perturber_mass_solar([(1.0, (3.0, 10.0))]))
    empty = SimpleNamespace(expected_systems=0.0, volume_pc3=100.0, shape=shape, populations={"old": pop(0.0, 20.0)})
    assert place_terms.from_profile(empty).dispersion_kms == tuning.SOLAR_NEIGHBOURHOOD_DISPERSION_KMS


# --- The tree of stars a StarSystem draws --------------------------------------------------------------------------------------------

def _systems(terms=None, count=150, seed=11, **flags):
    """`count` real systems (no planets) of a place's `terms` and the config `flags`."""
    systems = []
    with draw.bound(seed):
        for _ in range(count):
            config = SystemConfig()
            config.PLANETS = False
            config.PLACE_TERMS = terms
            for name, value in flags.items():
                setattr(config, name, value)
            systems.append(StarSystem(system_config=config))
    return systems


@pytest.mark.parametrize("fe_h", [-3.0, -1.0, -0.2, 0.5])
def test_a_g_dwarfs_chance_inside_10_au_follows_the_pooled_line_on_the_odds(fe_h):
    # The factor F(fe_h) / F(0) multiplies the odds of a companion inside 10 AU (GEN.103's rule, now on the close and
    # inner ranges together, GEN.448): odds(p(fe_h)) / odds(p(0)) = F(fe_h) / F(0).
    def inside(x):
        return companion_classes.odds(1.0, 0.44, x, 11_400.0).within_10_au

    def odds_of(p):
        return p / (1.0 - p)
    assert odds_of(inside(fe_h)) / odds_of(inside(0.0)) == pytest.approx(close_binary_trend.close_companion_factor(fe_h),
                                                                         rel=1e-6)


def test_a_g_dwarfs_chance_inside_10_au_falls_with_metallicity_and_stays_below_moe_kratter_badenes_ranges():
    # Moe, Kratter & Badenes 2019: 53, 40, 24 and 10% at [Fe/H] -3.0, -1.0, -0.2 and +0.5 (+-12, 6, 4, 3). Moe & Di
    # Stefano's Table 13 puts 0.15 +-0.03 companions inside logP 3.7 for solar-type stars, so the level is about 0.6 of
    # the pooled line's (both are in the design note); the shape follows it.
    found = [companion_classes.odds(1.0, 0.44, x, 11_400.0).within_10_au for x in (-3.0, -1.0, -0.2, 0.5)]
    assert found == sorted(found, reverse=True)
    for value, paper in zip(found, (0.53, 0.40, 0.24, 0.10)):
        assert 0.4 * paper < value < paper


def test_the_wide_share_outside_10_au_does_not_follow_metallicity():
    low, high = (companion_classes.odds(1.0, 0.44, x, 11_400.0) for x in (-1.0, 0.5))
    outer = lambda o: o.wide - (o.within_10_au - o.close)    # noqa: E731
    assert outer(low) == pytest.approx(outer(high), rel=1e-6)


def test_the_chances_never_pass_one_and_a_massive_star_changes_by_the_odds_not_the_chance():
    saturated = companion_classes.odds(30.0, 0.9, -3.0, 11_400.0)
    assert saturated.close + saturated.wide <= 1.0 + 1e-12
    solar = companion_classes.odds(30.0, 0.9, 0.0, 11_400.0).within_10_au
    poor = companion_classes.odds(30.0, 0.9, -0.7, 11_400.0).within_10_au
    assert poor / solar < 1.38        # the O and early-B surveys find no change within about 30%
    m_low, m_high = (companion_classes.odds(0.3, 0.26, x, 11_400.0).within_10_au for x in (-1.0, 0.5))
    assert m_low > 2 * m_high


def _widest_limit(terms, system):
    """The widest axis the place keeps for any star of `system` as the placed one (the heavier companion may take the primary's
    place, the limit was drawn for the star placed)."""
    def pair_mass(star):
        birth = star.initial_mass_sol or star.mass / 1.98841e30
        return birth * (1.0 + mass_ratio.mean_ratio(birth))

    return max(place_terms.wide_axis_limit_au(terms, pair_mass(star), system.primary_star.age) for star in system.stars)


def test_a_place_that_keeps_no_wide_pair_loses_it_from_the_chances():
    dense = place_terms.PlaceTerms(660.0, 120.0)       # the limit is under 5 AU: the 0.25 AU edge is all that is left
    here = companion_classes.odds(1.0, 0.44, 0.0, 0.1)
    assert here.wide == 0.0 and here.close < 0.44      # nothing above 0.25 AU, and the lost chances are lost
    assert here.close == pytest.approx(companion_classes.odds(1.0, 0.44, 0.0, 11_400.0).close)
    # The tree of stars keeps no axis the place cannot (GEN.631): a dense place has a pair or none, and no wide one beyond.
    pairs = [system for system in _systems(terms=dense, count=400) if system.binary_type == "wide"]
    assert all(system.wide_binary.separation_au <= _widest_limit(dense, system) for system in pairs)
    # A forced companion in a place that keeps little still keeps the place's limit.
    crowded = place_terms.PlaceTerms(60_000.0, 150.0)
    forced = [system for system in _systems(terms=crowded, BINARY_SYSTEM=True, count=60) if system.binary_type == "wide"]
    assert all(system.wide_binary.separation_au <= _widest_limit(crowded, system) for system in forced)


def test_the_ranges_of_separation_share_the_companions_by_moe_and_di_stefanos_lines():
    # A solar-type primary with a 1.55 solar mass pair: the lines give about 5% of its companions inside 0.25 AU, 95%
    # beyond; a massive primary has a larger close share (the lines are higher at log P 1).
    solar = companion_classes.odds(1.0, 0.44, 0.0, 1.0e6)
    massive = companion_classes.odds(20.0, 0.9, 0.0, 1.0e6)
    assert solar.close / (solar.close + solar.wide) == pytest.approx(0.055, abs=0.02)
    assert massive.close / (massive.close + massive.wide) > solar.close / (solar.close + solar.wide)


def test_the_binary_flags_still_force_the_outcome():
    assert {len(system.stars) for system in _systems(BINARY_SYSTEM=False, count=30)} == {1}
    forced = _systems(BINARY_SYSTEM=True, WIDE_BINARY=False, count=30)
    assert {system.binary_type for system in forced} <= {"close", None}    # none: a star with no lighter star to pair with
    assert any(system.binary_type == "close" for system in forced)
