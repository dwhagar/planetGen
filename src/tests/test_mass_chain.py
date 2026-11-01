# tests/test_mass_chain.py

"""GEN.400: the galaxy's mass chain (stellar mass function, stellar-to-halo relation, virial radius, concentration),
drawn from the galaxy seed (the Hill spheres read the mass enclosed inside their radius, `test_enclosed_mass.py`)."""

import math

import pytest

from planetgen.galaxy import halo_function, keepout, mass_chain, seed as galaxySeed
from planetgen.physics import constants
from planetgen.physics.orbits import calculate_hill_sphere


def _seed(number):
    return number.to_bytes(16, "big")


LOW, HIGH = mass_chain.AUTO_LOG_MASS_RANGE


def test_the_same_seed_gives_the_same_chain_and_two_seeds_differ():
    assert mass_chain.draw_mass_chain(_seed(7)) == mass_chain.draw_mass_chain(_seed(7))
    assert mass_chain.draw_mass_chain(_seed(7)) != mass_chain.draw_mass_chain(_seed(8))
    assert mass_chain.draw_mass_chain(None) == mass_chain.draw_mass_chain(mass_chain.FALLBACK_SEED)


def test_every_draw_stays_in_the_auto_window():
    for number in range(1, 2001):
        chain = mass_chain.draw_mass_chain(_seed(number))
        assert LOW - 1e-9 <= math.log10(chain.stellar_mass_solar) <= HIGH + 1e-9
        assert chain.halo_mass_solar > chain.stellar_mass_solar
        assert 1.0 < chain.concentration < 40.0
        assert chain.virial_radius_kpc > 30.0


def test_the_stellar_masses_of_10000_seeds_follow_the_double_schechter_function():
    """A histogram of the drawn masses against the function's own integral, within sampling error (4 sigma a bin)."""
    edges = [LOW + (HIGH - LOW) * i / 8 for i in range(9)]
    counts = [0] * 8
    draws = 10000
    for number in range(draws):
        value = mass_chain.log_stellar_mass_at(mass_chain._stream(_seed(number + 1), "stellar-mass").random())
        counts[min(int((value - LOW) / (HIGH - LOW) * 8), 7)] += 1
    # The expected share of each bin: the function integrated finely over it.
    steps = 400
    weights = []
    for i in range(8):
        step = (edges[i + 1] - edges[i]) / steps
        weights.append(sum(mass_chain.stellar_mass_density_per_dex(edges[i] + step * (k + 0.5)) for k in range(steps)) * step)
    total = sum(weights)
    for count, weight in zip(counts, weights):
        expected = draws * weight / total
        assert abs(count - expected) <= 4.0 * math.sqrt(expected * (1.0 - weight / total)), (counts, expected)


def test_the_function_is_the_papers_per_dex_density():
    """[R0238: Eq. 6] per dex: at M = M* the density is ln(10) e^-1 (phi1 + phi2)."""
    at_knee = mass_chain.stellar_mass_density_per_dex(mass_chain.LOG_M_CHARACTERISTIC)
    assert at_knee == pytest.approx(math.log(10.0) * math.exp(-1.0) * (mass_chain.PHI_1 + mass_chain.PHI_2))


@pytest.mark.parametrize("stellar, halo, ratio, c200, r200", [
    (3e9, 1.62e11, 0.0185, 10.1, 116), (1e10, 3.11e11, 0.0322, 9.4, 144), (3e10, 7.92e11, 0.0379, 8.6, 196),
    (5e10, 1.57e12, 0.0319, 8.0, 247), (1e11, 5.69e12, 0.0176, 7.0, 379),
])
def test_the_worked_table_of_the_design_note_is_reproduced(stellar, halo, ratio, c200, r200):
    median = mass_chain.median_halo_mass(stellar)
    assert median == pytest.approx(halo, rel=0.01)
    assert stellar / median == pytest.approx(ratio, rel=0.02)
    assert mass_chain.median_concentration(median) == pytest.approx(c200, abs=0.05)
    assert mass_chain.virial_radius_kpc(median) == pytest.approx(r200, abs=1.0)


def test_the_conditional_halo_is_the_median_relation_inverted_where_the_scatter_is_small():
    """Up to 3e10 solar masses of stars the first build's median inverse was within 0.03 dex of the conditional
    (GEN.399); above it the conditional is lighter by the Eddington bias of the steep halo mass function."""
    for stellar in (3e9, 1e10, 2e10):
        inverse = _inverse_of_the_median_relation(stellar)
        assert abs(math.log10(mass_chain.median_halo_mass(stellar) / inverse)) < 0.03
    gaps = [math.log10(_inverse_of_the_median_relation(s) / mass_chain.median_halo_mass(s)) for s in (6e10, 1e11, 2e11)]
    assert gaps[0] > 0.15 and gaps == sorted(gaps) and 0.3 < gaps[2] < 0.4


def _inverse_of_the_median_relation(stellar):
    low, high = 8.0, 17.0
    for _ in range(60):
        middle = 0.5 * (low + high)
        if math.log10(mass_chain.median_stellar_mass(10.0 ** middle)) < math.log10(stellar):
            low = middle
        else:
            high = middle
    return 10.0 ** (0.5 * (low + high))


def test_the_halo_quantiles_rise_with_the_share_and_the_stellar_mass():
    for stellar in (3e9, 5e10, 2e11):
        halos = [mass_chain.conditional_halo_mass(stellar, share) for share in (0.01, 0.16, 0.5, 0.84, 0.99)]
        assert halos == sorted(halos) and halos[0] < halos[-1]
    medians = [mass_chain.median_halo_mass(10.0 ** log_stellar) for log_stellar in [9.5 + 0.2 * i for i in range(10)]]
    assert medians == sorted(medians)
    with pytest.raises(ValueError):
        mass_chain.median_halo_mass(1.0e20)


def test_the_conditional_is_the_halo_function_times_the_relations_log_normal():
    """The share at a halo mass is the cumulative integral of n(Mh) x LogNormal on a finer grid than the table's."""
    stellar, log_stellar = 5e10, math.log10(5e10)
    delta = mass_chain.HALO_DELTA
    step = 0.002
    grid = [10.0 + step * i for i in range(int(round(6.0 / step)) + 1)]
    weights = [halo_function.mass_function_per_dex(g, delta)
               * math.exp(-0.5 * ((log_stellar - math.log10(mass_chain.median_stellar_mass(10.0 ** g)))
                                  / mass_chain.STELLAR_SCATTER_DEX) ** 2) for g in grid]
    total = math.fsum(weights)
    running, expected = 0.0, {}
    for g, w in zip(grid, weights):
        running += w
        for share in (0.16, 0.5, 0.84):
            if share not in expected and running / total >= share:
                expected[share] = g
    for share, log_halo in expected.items():
        assert math.log10(mass_chain.conditional_halo_mass(stellar, share)) == pytest.approx(log_halo, abs=0.015)


def test_the_halo_function_is_tinkers_and_the_variance_is_normalised_to_sigma_8():
    assert halo_function.sigma(4.0 / 3.0 * math.pi * (8.0 / halo_function.HUBBLE_H) ** 3
                               * halo_function._mean_matter_density_msun_mpc3()) == pytest.approx(
        halo_function.SIGMA_8, rel=1e-9)
    assert halo_function.tinker_parameters(200.0) == (0.186, 1.47, 2.57, 1.19)
    assert halo_function.tinker_parameters(3200.0) == (0.260, 2.66, 1.41, 2.44)
    a_norm, _, _, _ = halo_function.tinker_parameters(630.0)
    assert 0.218 < a_norm < 0.248
    with pytest.raises(ValueError):
        halo_function.tinker_parameters(100.0)
    # the number of halos per Mpc^3 per dex: about 4e-3 at 1e12 and 4e-5 at 1e14 solar masses ([R4367: Fig. 5])
    assert halo_function.mass_function_per_dex(12.0, 630.0) == pytest.approx(3.6e-3, rel=0.1)
    assert halo_function.mass_function_per_dex(14.0, 630.0) == pytest.approx(3.8e-5, rel=0.15)
    assert halo_function.sigma(1.0e12) == pytest.approx(2.29, abs=0.03)


def test_the_halo_and_concentration_scatter_follow_the_stated_widths():
    stellar = 5e10
    quantiles = {share: math.log10(mass_chain.conditional_halo_mass(stellar, share)) for share in (0.16, 0.5, 0.84)}
    halos, concentrations = [], []
    for number in range(1, 2001):
        chain = mass_chain.chain_from_stellar_mass(stellar, _seed(number))
        halos.append(math.log10(chain.halo_mass_solar))
        concentrations.append(math.log10(chain.concentration) - math.log10(mass_chain.median_concentration(chain.halo_mass_solar)))
    halos.sort()
    for share, expected in quantiles.items():
        assert halos[int(share * len(halos))] == pytest.approx(expected, abs=0.04)
    mean = sum(concentrations) / len(concentrations)
    spread = math.sqrt(sum((v - mean) ** 2 for v in concentrations) / (len(concentrations) - 1))
    assert abs(mean) < 4.0 * mass_chain.CONCENTRATION_SCATTER_DEX / math.sqrt(len(concentrations))
    assert spread == pytest.approx(mass_chain.CONCENTRATION_SCATTER_DEX, rel=0.08)


def test_the_virial_radius_is_the_radius_of_two_hundred_critical_densities():
    mass_kg = 1.0e12 * constants.SOLAR_MASS_TO_KG
    radius_m = mass_chain.virial_radius_kpc(1.0e12) * 1000.0 * constants.PARSEC_M
    mean = mass_kg / (4.0 / 3.0 * math.pi * radius_m ** 3)
    assert mean == pytest.approx(200.0 * mass_chain.CRITICAL_DENSITY_KG_M3, rel=1e-9)
    assert mass_chain.CRITICAL_DENSITY_KG_M3 == pytest.approx(1.878e-26 * mass_chain.HUBBLE_H ** 2, rel=1e-3)


def test_the_unit_seed_is_restored_after_a_unit():
    assert galaxySeed.current_galaxy_seed() is None
    with galaxySeed.seeded(_seed(3), "sector", "1/1/1"):
        assert galaxySeed.current_galaxy_seed() == _seed(3)
        with galaxySeed.seeded(_seed(4), "nebula-cell", "0/0/0"):
            assert galaxySeed.current_galaxy_seed() == _seed(4)
        assert galaxySeed.current_galaxy_seed() == _seed(3)
    assert galaxySeed.current_galaxy_seed() is None


def test_the_galactic_hill_radius_matches_the_hill_sphere():
    assert keepout.galactic_hill_radius_km(1e30, 2000.0, 2e42) == pytest.approx(
        calculate_hill_sphere(2000.0 * constants.PARSEC_M, 1e30, 2e42) / 1000.0)
