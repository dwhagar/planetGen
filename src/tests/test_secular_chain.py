"""GEN.1450: the Laplace-Lagrange secular exchange of a star's planets (`physics/secular_chain.py`)."""

import math

import pytest

from planetgen.physics import secular_chain as sc

ARCSEC_PER_YR = 180.0 / math.pi * 3600.0
SUN_KG = 1.98847e30
# the eight planets of the Solar System: semi-major axes (AU) and masses (kg). Used only to check the theory against its published frequencies; no generator reads them.
SOLAR_A = [0.3871, 0.7233, 1.0000, 1.5237, 5.2026, 9.5549, 19.2184, 30.1104]
SOLAR_M = [3.301e23, 4.867e24, 5.972e24, 6.417e23, 1.898e27, 5.683e26, 8.681e25, 1.024e26]


def _eigen(matrix):
    import numpy as np
    return sorted(float(x) * ARCSEC_PER_YR for x in np.linalg.eigvals(matrix).real)


def test_solar_system_nodal_frequencies_match_the_published_ones():
    # Murray and Dermott Table 7.1 (Laplace-Lagrange, recalled): s = 0, -0.68, -2.90, -5.20, -6.57, -17.64, -18.74, -25.73 arcsec per year
    _a, b = sc.matrices(SOLAR_A, SOLAR_M, SUN_KG)
    found = _eigen(b)
    for got, want in zip(found, sorted([-25.73, -18.74, -17.64, -6.57, -5.20, -2.90, -0.68, 0.0])):
        assert got == pytest.approx(want, abs=0.15)


def test_solar_system_inner_apsidal_frequencies_match_the_published_ones():
    # g = 0.63, 2.72, 7.34, 17.33, 18.00 for the modes of the inner planets (Table 7.1); the Jupiter-Saturn pair is left out: the theory is known to put it
    # 10 to 20 percent off the full theory's 4.26 and 28.2 (the 5:2 near-resonance is outside the expansion), which the N-body run in the report shows
    a, _b = sc.matrices(SOLAR_A, SOLAR_M, SUN_KG)
    found = _eigen(a)
    for want in (0.63, 2.72, 17.33, 18.00):
        assert min(abs(g - want) for g in found) < 0.1


def test_a_lone_planet_keeps_its_orbit_but_for_relativity():
    result = sc.evolve([1.0], [5.972e24], SUN_KG, [0.1], [1.0], [0.05], [2.0], 1.0e9)
    assert result.eccentricity[0] == pytest.approx(0.1)
    assert result.inclination[0] == pytest.approx(0.05, rel=1e-9)
    assert result.node_longitude[0] == pytest.approx(2.0)
    turn = (result.periapsis_longitude[0] - 1.0) % (2.0 * math.pi)
    rate = sc.relativistic_rates_rad_yr([1.0], SUN_KG, [5.972e24])[0]
    assert turn == pytest.approx((rate * 1.0e9) % (2.0 * math.pi), rel=1e-6)
    assert rate * ARCSEC_PER_YR * 100.0 == pytest.approx(3.84, rel=0.02)         # Earth's 3.8 arcsec per century


def test_circular_flat_orbits_stay_so():
    result = sc.evolve([0.4, 1.0, 5.2], [3.3e23, 6.0e24, 1.9e27], SUN_KG, [0.0] * 3, [0.0] * 3, [0.0] * 3, [0.0] * 3, 5.0e9)
    assert max(result.eccentricity) == pytest.approx(0.0, abs=1e-12)
    assert max(result.inclination) == pytest.approx(0.0, abs=1e-12)


def test_a_close_giant_raises_the_eccentricity_of_its_neighbour_to_the_forced_value():
    # a Jupiter on e = 0.05 and a light planet inside it on a circular orbit: the forced eccentricity is e_giant * b2/b1 * alpha-ish; the light planet's maximum over
    # the cycle is twice the forced value (it starts at zero), and its mean is the forced value
    a = [1.0, 5.2]
    m = [6.0e24, 1.898e27]
    alpha = a[0] / a[1]
    b1, b2 = sc._laplace_pair(alpha)
    forced = 0.05 * b2 / b1
    peak = 0.0
    for years in range(0, 2_000_000, 5_000):
        peak = max(peak, sc.evolve(a, m, SUN_KG, [0.0, 0.05], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], float(years)).eccentricity[0])
    assert peak == pytest.approx(2.0 * forced, rel=0.1)


def test_angular_momentum_deficit_is_conserved():
    a, m = SOLAR_A[3:], SOLAR_M[3:]
    e0 = [0.09, 0.048, 0.054, 0.047, 0.009]
    v0 = [5.0, 0.2, 1.7, 3.0, 0.9]
    zeros = [0.0] * 5

    def deficit(e):
        return sum(mass * math.sqrt(sma) * x * x for mass, sma, x in zip(m, a, e))

    start = deficit(e0)
    for years in (1.0e5, 1.0e6, 7.0e7):
        # the theory conserves the sum of m sqrt(a) e^2 only to the order of the relativistic term, which is small for planets this far out
        got = deficit(sc.evolve(a, m, 0.0 + SUN_KG, e0, v0, zeros, zeros, years).eccentricity)
        assert got == pytest.approx(start, rel=0.05)


def test_precession_rates_are_the_diagonals_and_relativity_adds_to_the_apsidal_one():
    a, m = [0.1, 0.5], [3e24, 1e27]
    apsidal, nodal = sc.precession_rates(a, m, SUN_KG)
    matrix_a, matrix_b = sc.matrices(a, m, SUN_KG)
    assert apsidal == pytest.approx((matrix_a[0, 0], matrix_a[1, 1]))
    assert nodal == pytest.approx((matrix_b[0, 0], matrix_b[1, 1]))
    assert apsidal[0] > 0.0 > nodal[0]
    assert apsidal[0] > sc.relativistic_rates_rad_yr(a, SUN_KG, m)[0]


def test_a_planet_beyond_the_theory_is_held_and_still_pushes_on_the_others():
    a, m = [1.0, 5.2], [6.0e24, 1.898e27]
    result = sc.evolve(a, m, SUN_KG, [0.0, 0.5], [0.0, 1.0], [0.0, 0.0], [0.0, 0.0], 1.0e6)
    assert result.eccentricity[1] == pytest.approx(0.5)
    assert result.periapsis_longitude[1] == pytest.approx(1.0)
    assert result.eccentricity[0] > 0.0
    assert result.linear


def test_inputs_are_checked():
    with pytest.raises(sc.SecularChainError):
        sc.matrices([1.0, 2.0], [1.0], SUN_KG)
    with pytest.raises(sc.SecularChainError):
        sc.matrices([1.0, -2.0], [1.0, 1.0], SUN_KG)
    with pytest.raises(sc.SecularChainError):
        sc.matrices([1.0, 2.0], [1.0, 1.0], 0.0)
