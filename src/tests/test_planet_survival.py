# tests/test_planet_survival.py -- GEN.325: the orbit beyond which passing systems have stripped a planet
# (docs/design/galaxy/realistic-galactic-core.md section 7, docs/design/galaxy/globular-clusters.md section 7).

import math
from types import SimpleNamespace

import pytest

from planetgen.generation import place_terms
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.generation.validation import reach_au
from planetgen.physics import metallicity, planet_survival
from planetgen.util import draw

# The table of section 7: place, systems per pc^3, dispersion (km/s), largest surviving orbit (AU, as the note rounds it).
NUCLEUS_TABLE = [
    ("nuclear cluster, 0.5 pc", 2.4e6, 212.0, 0.03),
    ("1 pc", 7.1e5, 173.0, 0.06),
    ("4 pc", 2.3e4, 128.0, 0.4),
    ("10 pc", 1.2e3, 96.0, 2.0),
    ("nuclear disc, centre", 540.0, 90.0, 3.0),
    ("100 pc", 180.0, 80.0, 5.6),
    ("200 pc", 58.0, 70.0, 11.0),
]

# The cluster note's four rows: systems per pc^3, dispersion (km/s), the hard-binary boundary a_h (AU) and the expected
# ionising passes of an orbit at a_h that the note reads (188, 11, 2.6 and 2.1).
CLUSTER_TABLE = [(4.8e5, 11.5, 3.4, 188.0), (1.0e4, 8.0, 6.9, 11.0), (1.0e3, 6.0, 12.0, 2.6), (1.0e2, 3.0, 49.0, 2.1)]


@pytest.mark.parametrize("place, density, sigma, note", NUCLEUS_TABLE, ids=[row[0] for row in NUCLEUS_TABLE])
def test_the_largest_surviving_orbit_is_the_note_s_table(place, density, sigma, note):
    limit = planet_survival.planet_survival_limit_au(density, sigma)
    assert limit == pytest.approx(note, rel=0.08)      # the note rounds to one or two figures
    assert planet_survival.expected_passes(density, sigma, limit) == pytest.approx(1.0, rel=1e-9)


@pytest.mark.parametrize("density, sigma, a_h, note", CLUSTER_TABLE)
def test_the_cluster_table_reproduces_the_expected_passes_at_the_hard_binary_boundary(density, sigma, a_h, note):
    passes = planet_survival.expected_passes(density, sigma, a_h)
    assert passes == pytest.approx(note, rel=0.03)
    assert math.exp(-passes) < 0.13 or passes < 3.0      # beyond a_h essentially nothing survives a dense core


def test_the_limit_falls_as_the_density_rises_and_takes_any_density_and_dispersion():
    densities = (0.01, 0.2, 5.0, 300.0, 1.0e4, 1.0e6)
    for sigma in (1.0, 10.0, 60.0, 200.0):
        limits = [planet_survival.planet_survival_limit_au(n, sigma) for n in densities]
        assert all(a > b > 0.0 for a, b in zip(limits, limits[1:]))
    # Wide of 2 G M / v^2 the cross section decides and a hotter field strips more; inside it the focusing decides
    # (it goes as 1 / v) and a hotter field strips less.
    wide = [planet_survival.planet_survival_limit_au(0.2, s) for s in (60.0, 100.0, 200.0)]
    assert wide[0] > wide[1] > wide[2]
    focused = [planet_survival.planet_survival_limit_au(1.0e4, s) for s in (0.5, 1.0, 2.0)]
    assert focused[0] < focused[1] < focused[2]


def test_the_limit_is_the_root_of_the_passes_and_older_planets_keep_less():
    for age in (0.5, 4.0, 12.0, 13.8):
        limit = planet_survival.planet_survival_limit_au(500.0, 70.0, age_gy=age)
        assert planet_survival.expected_passes(500.0, 70.0, limit, age_gy=age) == pytest.approx(1.0)
    assert (planet_survival.planet_survival_limit_au(500.0, 70.0, age_gy=1.0)
            > planet_survival.planet_survival_limit_au(500.0, 70.0, age_gy=10.0))
    # A heavier host focuses the passers-by less than it holds its planets: more mass, more passes at one orbit.
    assert (planet_survival.expected_passes(500.0, 70.0, 5.0, host_mass_solar=1.0)
            > planet_survival.expected_passes(500.0, 70.0, 5.0, host_mass_solar=0.2))


def test_a_place_with_nothing_passing_by_cuts_nothing():
    assert math.isinf(planet_survival.planet_survival_limit_au(0.0, 30.0))
    assert math.isinf(planet_survival.planet_survival_limit_au(0.2, 0.0))
    assert math.isinf(planet_survival.planet_survival_limit_au(0.2, 30.0, age_gy=0.0))


def test_the_suns_neighbourhood_keeps_planets_out_to_a_few_hundred_au():
    terms = place_terms.SOLAR_NEIGHBOURHOOD
    assert terms.systems_per_pc3 == pytest.approx(terms.mass_density_solar_pc3 / terms.mean_system_mass_solar)
    limit = place_terms.planet_limit_au(terms)
    assert 150.0 < limit < 600.0
    assert place_terms.planet_limit_au(terms, 100.0) == limit                       # older than 12 Gy is exposed 12 Gy
    assert place_terms.planet_limit_au(terms, 2.0) > limit                          # a young star's planets are safer
    assert place_terms.planet_limit_au(terms, 2.0) == planet_survival.planet_survival_limit_au(
        terms.systems_per_pc3, terms.dispersion_kms, age_gy=2.0)


def test_the_profile_s_systems_per_volume_reach_the_limit():
    class _Profile:
        expected_systems, volume_pc3 = 3.0e4, 1.0
        populations = {"bulge": SimpleNamespace(share=1.0, dispersion_kms=128.0 * math.sqrt(3.0), age_gy=(8.0, 12.0))}   # three-axis total
        shape = None

    from planetgen.galaxy import system_mass as galaxy_system_mass
    original = galaxy_system_mass.mean_over_solar
    galaxy_system_mass.mean_over_solar = lambda weights, shape: 0.5
    try:
        terms = place_terms.from_profile(_Profile())
    finally:
        galaxy_system_mass.mean_over_solar = original
    assert terms.systems_per_pc3 == pytest.approx(3.0e4)
    assert place_terms.planet_limit_au(terms) == pytest.approx(planet_survival.planet_survival_limit_au(3.0e4, 128.0))


# --- The cull in a system ------------------------------------------------------------------------------------------------


def _g_dwarf(terms, seed):
    with draw.bound(f"gen325:{seed}"):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        cfg.PLACE_TERMS = terms
        return StarSystem(system_config=cfg)


def test_the_ceiling_is_the_tighter_of_the_stable_orbit_and_the_survival_limit():
    dense = place_terms.PlaceTerms(1.2e3 * 0.5, 96.0, mean_system_mass_solar=0.5)
    system = _g_dwarf(dense, 0)
    star = system.primary_star
    expected = place_terms.planet_limit_au(dense, star.age)
    assert expected < 5.0
    assert system._orbit_ceiling_au(star) == pytest.approx(expected)
    sparse = _g_dwarf(None, 1)
    stable = sparse._orbit_ceiling_au(sparse.primary_star)
    assert stable == pytest.approx(min(place_terms.planet_limit_au(place_terms.SOLAR_NEIGHBOURHOOD, sparse.primary_star.age),
                                       stable))


def test_no_planet_of_a_dense_core_lies_past_the_survival_limit_and_the_sun_s_neighbourhood_keeps_far_ones():
    dense = place_terms.PlaceTerms(2.3e4 * 0.5, 128.0, mean_system_mass_solar=0.5)
    far = 0
    for seed in range(12):
        system = _g_dwarf(dense, seed)
        limit = place_terms.planet_limit_au(dense, system.primary_star.age)
        assert all(reach_au(planet) <= limit * (1.0 + 1e-9) for planet in system.planets), seed
        open_system = _g_dwarf(None, seed)
        far += sum(1 for planet in open_system.planets if reach_au(planet) > limit)
    assert far > 0       # the same stars with the Sun's neighbourhood keep planets past the dense core's limit


def test_a_metal_poor_cluster_thins_its_giants_by_the_metallicity_factor_as_well():
    # GEN.236 (the factor), the metallicity exponent of the gas outcome (GEN.828) and GEN.325 (the cut) act together: the giants of a
    # metal-poor star are fewer, by about a hundredth at [Fe/H] -1 (giants go as 10^(2 [Fe/H]), [R1377]). The same streams
    # draw the same cores, so a star that keeps a giant at lower [Fe/H] keeps it at the Sun's too.
    from planetgen.physics import planet_mass

    def giants(fe_h):
        return sum(planet_mass.planet_mass(3.0, 1.989e30, 2.7, fe_h, 1.0, draw.Stream(seed)).runaway for seed in range(4000))

    solar = giants(0.0)
    assert solar > 150
    for fe_h in (-1.0, -1.3):
        assert giants(fe_h) / solar <= 3.0 * metallicity.giant_occurrence_factor(fe_h)
        assert metallicity.giant_occurrence_factor(fe_h) < metallicity.planet_occurrence_factor(fe_h) < 1.0
    assert giants(-2.2) <= 2
