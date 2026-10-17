# tests/test_place_terms.py -- GEN.103: the close-companion factor that follows a star's [Fe/H], the widest axis a wide
# pair keeps at its place, and how `StarSystem._companion` uses them.

import math
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation import place_terms
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

POINTS = sorted(tuning.CLOSE_BINARY_FRACTION_BY_FE_H.items())


def test_the_factor_is_the_pooled_line_over_its_solar_value():
    f0 = 0.20  # the line at [Fe/H] 0: between the -0.2 and +0.5 points
    assert place_terms.close_companion_factor(0.0) == pytest.approx(1.0)
    for fe_h, fraction in POINTS:
        assert place_terms.close_companion_factor(fe_h) == pytest.approx(fraction / f0)
    # Flat past the ends, and 1 for a star with no [Fe/H].
    assert place_terms.close_companion_factor(-9.0) == place_terms.close_companion_factor(-3.0)
    assert place_terms.close_companion_factor(3.0) == place_terms.close_companion_factor(0.5)
    assert place_terms.close_companion_factor(None) == 1.0
    # Straight between two points.
    assert place_terms.close_companion_factor(-0.6) == pytest.approx((0.40 + 0.24) / 2 / f0)


def test_the_factor_moves_the_odds_and_never_passes_one():
    # Boss, 7:19 PM EDT: the odds p / (1 - p), not p.
    assert place_terms.scaled_close_chance(0.22, 0.0) == pytest.approx(0.22)
    for chance in (0.05, 0.22, 0.45, 0.9, 0.99):
        for fe_h in (-5.0, -3.0, -1.0, -0.2, 0.0, 0.5, 2.0):
            moved = place_terms.scaled_close_chance(chance, fe_h)
            assert 0.0 < moved < 1.0
            odds = chance / (1 - chance) * place_terms.close_companion_factor(fe_h)
            assert moved / (1 - moved) == pytest.approx(odds)
    assert place_terms.scaled_close_chance(0.0, -3.0) == 0.0 and place_terms.scaled_close_chance(1.0, 0.5) == 1.0
    # A star already likely to have a companion rises less in absolute terms.
    assert (place_terms.scaled_close_chance(0.45, -1.0) - 0.45) < (1 - place_terms.scaled_close_chance(0.1, -1.0)) * 10
    assert place_terms.scaled_close_chance(0.45, -0.7) / 0.45 == pytest.approx(1.29, abs=0.02)   # inside "within about 30%"


def test_the_factor_falls_as_metallicity_rises():
    values = [place_terms.close_companion_factor(x / 10) for x in range(-40, 10)]
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
    assert 2.0 < everything / stars < 4.0
    assert place_terms.PlaceTerms(1.0, 30.0).perturber_mass_solar == pytest.approx(everything)


def test_a_young_population_has_the_heavier_perturbers():
    young = place_terms.mass_weighted_perturber_mass_solar([(1.0, tuning.STELLAR_POPULATION_AGE_RANGES_GY["young"])])
    old = place_terms.mass_weighted_perturber_mass_solar([(1.0, tuning.STELLAR_POPULATION_AGE_RANGES_GY["old"])])
    assert young > old > 1.0
    mixed = place_terms.mass_weighted_perturber_mass_solar([(0.5, (0.0, 0.1)), (0.5, (3.0, 10.0))])
    assert old < mixed < young


def test_a_cell_with_no_stars_falls_back_to_the_mean_system_mass():
    assert place_terms.mass_weighted_perturber_mass_solar([]) == tuning.MEAN_SYSTEM_MASS_SOLAR


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
    # docs/design/placement-by-population.md 3.3, 2 solar masses, 10 Gyr, the default perturber mass (<M^2> / <M> of the
    # whole disk with remnants).
    for density, sigma, encounter, tide in ((0.048, 30, 11_400, 500_000), (6.9, 120, 317, 95_000),
                                            (77, 120, 28, 42_000), (660, 120, 3.3, 21_000)):
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

    profile = SimpleNamespace(expected_systems=50.0, volume_pc3=100.0,
                              populations={"a": pop(0.75, 20.0), "b": pop(0.25, 60.0)})
    terms = place_terms.from_profile(profile)
    assert terms.mass_density_solar_pc3 == pytest.approx(0.5 * tuning.MEAN_SYSTEM_MASS_SOLAR)
    assert terms.dispersion_kms == pytest.approx(math.sqrt(0.75 * 400 + 0.25 * 3600))
    assert terms.perturber_mass_solar == pytest.approx(
        place_terms.mass_weighted_perturber_mass_solar([(1.0, (3.0, 10.0))]))
    empty = SimpleNamespace(expected_systems=0.0, volume_pc3=100.0, populations={"a": pop(0.0, 20.0)})
    assert place_terms.from_profile(empty).dispersion_kms == tuning.SOLAR_NEIGHBOURHOOD_DISPERSION_KMS


# --- StarSystem._companion --------------------------------------------------------------------------------------------

def _system(letter="G", fe_h=0.0, terms=None, age_gy=5.0, mass_sol=1.0, **flags):
    system = object.__new__(StarSystem)
    config = SystemConfig()
    config.PLACE_TERMS = terms
    for name, value in flags.items():
        setattr(config, name, value)
    system.system_config = config
    system.primary_star = SimpleNamespace(type=letter + "2V", metallicity_fe_h=fe_h, age=age_gy,
                                          initial_mass_sol=mass_sol, mass=mass_sol * 1.98841e30)
    return system


def _kinds(system, n=3000, seed=11):
    counts = {"close": 0, "wide": 0, None: 0}
    with draw.bound(seed):
        for _ in range(n):
            counts[system._companion()[0]] += 1
    return {kind: count / n for kind, count in counts.items()}


@pytest.mark.parametrize("fe_h, expected", [(-3.0, 0.43), (-1.0, 0.36), (-0.2, 0.25), (0.5, 0.12)])
def test_a_g_dwarfs_close_share_follows_the_pooled_line(fe_h, expected):
    # Today's table (G 0.44, split half and half) with the factor on the odds: 0.22 moved by F / F(0), against Moe,
    # Kratter and Badenes' 53, 40, 24 and 10% (+-12, 6, 4, 3).
    shares = _kinds(_system(fe_h=fe_h, terms=place_terms.PlaceTerms(0.05, 30.0)))
    assert shares["close"] == pytest.approx(expected, abs=0.035)


def test_the_wide_share_does_not_follow_metallicity():
    terms = place_terms.PlaceTerms(0.05, 30.0)
    low, high = (_kinds(_system(fe_h=x, terms=terms))["wide"] for x in (-1.0, 0.5))
    assert low == pytest.approx(0.22, abs=0.035) and high == pytest.approx(0.22, abs=0.035)


def test_every_mass_follows_the_line_and_the_chances_never_pass_one():
    terms = place_terms.PlaceTerms(0.05, 30.0)
    o_star = _kinds(_system("O", fe_h=-3.0, terms=terms, mass_sol=30.0, age_gy=0.005))
    assert o_star[None] == pytest.approx(0.0, abs=0.01)         # the companion chance saturates: one every time
    assert o_star["close"] == pytest.approx(0.68, abs=0.035)    # the odds form: up from 0.45, not to 1
    assert o_star["wide"] == pytest.approx(0.32, abs=0.035)     # the wide share gives way where the two pass 1
    # O stars at [Fe/H] -0.7 change by about 30% at most (the O and early-B surveys find no trend within that).
    solar = _kinds(_system("O", fe_h=0.0, terms=terms, mass_sol=30.0, age_gy=0.005))["close"]
    poor = _kinds(_system("O", fe_h=-0.7, terms=terms, mass_sol=30.0, age_gy=0.005))["close"]
    assert poor / solar < 1.38
    m_dwarf_low, m_dwarf_high = (_kinds(_system("M", fe_h=x, terms=terms, mass_sol=0.3))["close"] for x in (-1.0, 0.5))
    assert m_dwarf_low > 3 * m_dwarf_high


def test_a_place_that_keeps_no_wide_pair_dissolves_it():
    dense = place_terms.PlaceTerms(660.0, 120.0)       # the limit is about 9 AU: under the 50 AU lower edge
    shares = _kinds(_system(terms=dense))
    assert shares["wide"] == 0.0
    # A forced companion there is a close one; a pair forced wide is wide all the same.
    with draw.bound(3):
        assert {_system(terms=dense, BINARY_SYSTEM=True)._companion()[0] for _ in range(200)} == {"close"}
        assert {_system(terms=dense, WIDE_BINARY=True, BINARY_SYSTEM=True)._companion()[0] for _ in range(20)} == {"wide"}


def test_the_binary_flags_still_force_the_outcome():
    with draw.bound(5):
        assert _system(BINARY_SYSTEM=False)._companion() == (None, None, None)
        for _ in range(50):
            kind, ratio, _reason = _system(BINARY_SYSTEM=True, WIDE_BINARY=False)._companion()
            assert kind == "close" and tuning.BINARY_MASS_RATIO_RANGE[0] <= ratio <= tuning.BINARY_MASS_RATIO_RANGE[1]
