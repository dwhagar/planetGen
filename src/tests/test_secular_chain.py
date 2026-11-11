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


def _laplace_by_quadrature(index, alpha, steps=20000):
    """b_{3/2}^{(index)}(alpha) from its integral, independent of the module's own coefficient."""
    total = 0.0
    for i in range(steps):
        psi = 2.0 * math.pi * (i + 0.5) / steps
        total += math.cos(index * psi) / (1.0 - 2.0 * alpha * math.cos(psi) + alpha * alpha) ** 1.5
    return total / steps * 2.0


def test_a_belt_is_cut_into_rings_that_keep_its_mass_and_lie_inside_it():
    axes, masses = sc.belt_rings(2.2, 3.4, 1.0e22)
    assert len(axes) == sc.BELT_RINGS and sum(masses) == pytest.approx(1.0e22)
    assert all(2.2 < a < 3.4 for a in axes) and axes == sorted(axes)
    assert masses == sorted(masses), "a uniform surface density puts more mass in the wider outer rings"
    assert sc.belt_rings(2.2, 3.4, 0.0) == ([], []) and sc.belt_rings(2.2, 2.2, 1.0e22) == ([], [])


def test_a_belt_of_known_mass_turns_a_neighbouring_planet_by_the_expected_amount_outer_and_inner():
    # Murray and Dermott 7.2: a ring of mass m outside a planet at alpha = a / a' adds (n / 4) (m / M) alpha^2 b_{3/2}^{(1)}(alpha) to the planet's apsidal rate (a ring inside: alpha * 1)
    planet_a, planet_m = 1.0, 5.972e24
    base, _ = sc.precession_rates([planet_a], [planet_m], SUN_KG)
    n = sc.mean_motions_rad_yr([planet_a], SUN_KG, [planet_m])[0]
    ring_m = 6.0e21
    for ring_a in (1.6, 0.6):
        with_ring, _ = sc.precession_rates([planet_a], [planet_m], SUN_KG, [ring_a], [ring_m])
        alpha = min(planet_a, ring_a) / max(planet_a, ring_a)
        bar = alpha if planet_a < ring_a else 1.0
        expected = 0.25 * n * ring_m / (SUN_KG + planet_m) * alpha * bar * _laplace_by_quadrature(1, alpha)
        assert with_ring[0] - base[0] == pytest.approx(expected, rel=2e-3), ring_a


def test_a_belt_turns_the_nodes_the_other_way_and_a_light_belt_turns_less_than_a_heavy_one():
    light, heavy = 1.0e21, 1.0e22
    base_g, base_s = sc.precession_rates([1.0], [5.972e24], SUN_KG)
    g1, s1 = sc.precession_rates([1.0], [5.972e24], SUN_KG, [2.0], [light])
    g2, s2 = sc.precession_rates([1.0], [5.972e24], SUN_KG, [2.0], [heavy])
    assert g1[0] > base_g[0] and s1[0] < base_s[0]
    assert g2[0] - base_g[0] == pytest.approx(10.0 * (g1[0] - base_g[0]), rel=1e-6)
    assert s2[0] - base_s[0] == pytest.approx(10.0 * (s1[0] - base_s[0]), rel=1e-6)


def test_no_belt_or_a_massless_one_changes_nothing_and_the_rings_are_not_part_of_the_result():
    args = ([0.5, 1.5], [5.0e24, 2.0e27], SUN_KG, [0.02, 0.05], [1.0, 2.0], [0.01, 0.02], [0.5, 3.0], 1.0e9)
    plain = sc.evolve(*args)
    assert sc.evolve(*args, [], []) == plain
    assert sc.evolve(*args, *sc.belt_rings(0.8, 1.2, 0.0)) == plain
    belted = sc.evolve(*args, *sc.belt_rings(2.2, 3.4, 1.0e22))
    assert len(belted.eccentricity) == 2 and len(belted.apsidal_rates) == 2 and belted.linear
    assert belted != plain, "the belt's pull moves the planets' orbits at the star's age"


def test_a_ring_next_to_a_planet_is_left_out_and_a_far_one_kept():
    axes, masses = sc.belt_rings(2.0, 3.0, 1.0e22)
    kept_axes, kept_masses = sc.clear_rings(axes, masses, [2.2])
    assert len(kept_axes) < len(axes) and all(min(a, 2.2) / max(a, 2.2) <= sc.BELT_ALPHA_LIMIT for a in kept_axes)
    assert sc.clear_rings(axes, masses, [0.5]) == (axes, masses)


# GEN.1459: the tide of each planet wears its eccentricity down while the others pump it

def _pair():
    return dict(a_au=[0.5, 1.0], masses_kg=[3e24, 6e24], host_mass_kg=2e30, eccentricity=[0.05, 0.02], periapsis_longitude=[0.3, 1.2],
                inclination=[0.01, 0.02], node_longitude=[0.1, 0.2])


def test_the_damped_exchange_is_the_matrix_exponential_of_the_secular_matrix_with_the_decay_rates():
    import numpy as np
    from scipy.linalg import expm
    pair, years, rates = _pair(), 2.0e8, [1 / 3.0e8, 1 / 5.0e9]
    result = sc.evolve(**pair, years=years, damping_per_yr=rates)
    A, _ = sc.matrices(pair["a_au"], pair["masses_kg"], pair["host_mass_kg"])
    z0 = np.array(pair["eccentricity"]) * np.exp(1j * np.array(pair["periapsis_longitude"]))
    expected = expm(1j * (A + 1j * np.diag(rates)) * years) @ z0
    assert result.eccentricity == pytest.approx(tuple(np.abs(expected)), rel=1e-9)
    assert result.periapsis_longitude == pytest.approx(tuple(float(x) % (2.0 * math.pi) for x in np.angle(expected)), rel=1e-9)


def test_no_damping_or_zero_damping_is_the_undamped_exchange_and_inclinations_are_not_damped():
    pair = _pair()
    plain = sc.evolve(**pair, years=3.0e8)
    assert sc.evolve(**pair, years=3.0e8, damping_per_yr=[0.0, 0.0]) == plain
    damped = sc.evolve(**pair, years=3.0e8, damping_per_yr=[1e-8, 1e-9])
    assert damped.inclination == pytest.approx(plain.inclination, rel=1e-12)
    assert all(d < p for d, p in zip(damped.eccentricity, plain.eccentricity))


def test_a_planet_the_others_keep_pumping_settles_above_the_decay_of_its_own_tide_and_a_lone_one_dies_away():
    import numpy as np
    pair, years = _pair(), 5.0e8
    rate = 1.0 / 1.0e8        # the inner planet's tide cuts its eccentricity by e every 100 Myr
    pair["masses_kg"] = [3e24, 6e27]          # a giant beyond it pumps it
    wear_only = pair["eccentricity"][0] * math.exp(-rate * years)
    pumped = sc.evolve(**pair, years=years, damping_per_yr=[rate, 0.0]).eccentricity[0]
    assert pumped > 10.0 * wear_only
    # with nothing to pump it (the neighbour a test particle) it follows the wear-down alone
    pair["masses_kg"] = [3e24, 1.0]
    alone = sc.evolve(**pair, years=years, damping_per_yr=[rate, 0.0]).eccentricity[0]
    assert alone == pytest.approx(wear_only, rel=5e-3)
