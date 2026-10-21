# tests/test_mass_chain.py

"""GEN.400: the galaxy's mass chain (stellar mass function, stellar-to-halo relation, virial radius, concentration),
drawn from the galaxy seed (the Hill spheres read the mass enclosed inside their radius, `test_enclosed_mass.py`)."""

import math

import pytest

from planetgen.galaxy import keepout, mass_chain, seed as galaxySeed
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
    """Baldry et al. 2012 Eq. 6 per dex: at M = M* the density is ln(10) e^-1 (phi1 + phi2)."""
    at_knee = mass_chain.stellar_mass_density_per_dex(mass_chain.LOG_M_CHARACTERISTIC)
    assert at_knee == pytest.approx(math.log(10.0) * math.exp(-1.0) * (mass_chain.PHI_1 + mass_chain.PHI_2))


@pytest.mark.parametrize("stellar, halo, ratio, c200, r200", [
    (3e9, 1.65e11, 0.018, 10.0, 116), (1e10, 3.15e11, 0.032, 9.4, 144), (3e10, 8.27e11, 0.036, 8.5, 199),
    (5e10, 2.00e12, 0.025, 7.8, 268), (1e11, 1.07e13, 0.0093, 6.6, 468),
])
def test_the_worked_table_of_the_design_note_is_reproduced(stellar, halo, ratio, c200, r200):
    median = mass_chain.median_halo_mass(stellar)
    assert median == pytest.approx(halo, rel=0.01)
    assert stellar / median == pytest.approx(ratio, rel=0.02)
    assert mass_chain.median_concentration(median) == pytest.approx(c200, abs=0.05)
    assert mass_chain.virial_radius_kpc(median) == pytest.approx(r200, abs=1.0)


def test_the_inverse_is_exact_and_monotonic():
    last = 0.0
    for log_stellar in [9.5 + 0.1 * i for i in range(19)]:
        halo = mass_chain.median_halo_mass(10.0 ** log_stellar)
        assert mass_chain.median_stellar_mass(halo) == pytest.approx(10.0 ** log_stellar, rel=1e-9)
        assert halo > last
        last = halo
    with pytest.raises(ValueError):
        mass_chain.median_halo_mass(1.0e20)


def test_the_effective_scatter_is_the_relations_scatter_over_its_slope():
    low, high = mass_chain.stellar_halo_slope(1e11), mass_chain.stellar_halo_slope(1e14)
    assert low == pytest.approx(1.0 + mass_chain.BETA, abs=0.2)  # steep at low mass
    assert high == pytest.approx(1.0 - mass_chain.GAMMA, abs=0.05)  # shallow at high mass
    assert mass_chain.effective_halo_scatter_dex(1e14) == pytest.approx(0.15 / high)


def test_the_halo_and_concentration_scatter_follow_the_stated_widths():
    stellar = 5e10
    median = mass_chain.median_halo_mass(stellar)
    sigma = mass_chain.effective_halo_scatter_dex(median)
    halos, concentrations = [], []
    for number in range(1, 2001):
        chain = mass_chain.chain_from_stellar_mass(stellar, _seed(number))
        halos.append(math.log10(chain.halo_mass_solar) - math.log10(median))
        concentrations.append(math.log10(chain.concentration) - math.log10(mass_chain.median_concentration(chain.halo_mass_solar)))
    for values, expected in ((halos, sigma), (concentrations, mass_chain.CONCENTRATION_SCATTER_DEX)):
        mean = sum(values) / len(values)
        spread = math.sqrt(sum((v - mean) ** 2 for v in values) / (len(values) - 1))
        assert abs(mean) < 4.0 * expected / math.sqrt(len(values))
        assert spread == pytest.approx(expected, rel=0.08)


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


def test_the_old_milky_way_constant_is_gone():
    assert not hasattr(constants, "MILKY_WAY_MASS")
    assert keepout.galactic_hill_radius_km(1e30, 2000.0, 2e42) == pytest.approx(
        calculate_hill_sphere(2000.0 * constants.PARSEC_M, 1e30, 2e42) / 1000.0)
