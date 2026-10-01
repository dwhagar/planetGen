"""
Kepler solver extremes (TODO TEST.34).

Direct tests of `keplerMotion.solve_eccentric_anomaly` past e = 0.99, with
mean anomalies below zero and above 2*pi, with too few iterations for
Newton-Raphson to converge (it falls back to bisection rather than
returning the last unconverged iterate), and of `_real_cube_root` at 0 and
for negative values. Reference values belong to the physics reference
gate; these check behavior only.
"""
import math

import pytest

from stellarObjects import keplerMotion as km

TWO_PI = 2 * math.pi


def kepler_residual(eccentric_anomaly, mean_anomaly, eccentricity):
    """How far `E - e*sin(E)` is from the wrapped mean anomaly."""
    return eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly) - (mean_anomaly % TWO_PI)


@pytest.mark.parametrize("eccentricity", [0.991, 0.995, 0.999, 0.9999, 0.999999, 1 - 1e-12])
@pytest.mark.parametrize("mean_anomaly", [1e-9, 1e-4, 0.01, 0.5, math.pi, 5.0, TWO_PI - 1e-9])
def test_eccentricity_above_0_99_satisfies_keplers_equation(eccentricity, mean_anomaly):
    e_anomaly = km.solve_eccentric_anomaly(mean_anomaly, eccentricity)
    assert math.isfinite(e_anomaly)
    assert 0 <= e_anomaly <= TWO_PI
    assert abs(kepler_residual(e_anomaly, mean_anomaly, eccentricity)) < 1e-9


@pytest.mark.parametrize("eccentricity", [0.0, 0.5, 0.95, 0.999])
@pytest.mark.parametrize("mean_anomaly", [0.3, 2.0, 4.5])
@pytest.mark.parametrize("turns", [-3, -1, 1, 7])
def test_mean_anomaly_outside_0_to_2pi_wraps(eccentricity, mean_anomaly, turns):
    """A negative mean anomaly, or one past 2*pi, gives the same E as
    its wrapped value."""
    base = km.solve_eccentric_anomaly(mean_anomaly, eccentricity)
    shifted = km.solve_eccentric_anomaly(mean_anomaly + turns * TWO_PI, eccentricity)
    assert shifted == pytest.approx(base, abs=1e-9)


@pytest.mark.parametrize("eccentricity", [0.2, 0.9, 0.999])
def test_small_negative_mean_anomaly_lands_just_below_2pi(eccentricity):
    e_anomaly = km.solve_eccentric_anomaly(-1e-3, eccentricity)
    assert math.pi < e_anomaly < TWO_PI
    assert abs(kepler_residual(e_anomaly, -1e-3, eccentricity)) < 1e-9


@pytest.mark.parametrize("mean_anomaly", [-1e6, -100.0, 100.0, 1e6])
def test_huge_mean_anomalies_still_converge(mean_anomaly):
    e_anomaly = km.solve_eccentric_anomaly(mean_anomaly, 0.995)
    assert abs(kepler_residual(e_anomaly, mean_anomaly, 0.995)) < 1e-8


@pytest.mark.parametrize("eccentricity", [0.1, 0.7, 0.99, 0.9999])
@pytest.mark.parametrize("max_iterations", [0, 1, 2])
def test_non_convergence_is_not_silently_returned(eccentricity, max_iterations):
    """Too few Newton steps used to hand back the last, unconverged
    iterate; the solver now falls back to bisection, so the answer still
    satisfies Kepler's equation."""
    mean_anomaly = 0.1
    e_anomaly = km.solve_eccentric_anomaly(mean_anomaly, eccentricity, max_iterations=max_iterations)
    assert abs(kepler_residual(e_anomaly, mean_anomaly, eccentricity)) < 1e-9


def test_zero_tolerance_never_breaks_early_but_still_converges():
    e_anomaly = km.solve_eccentric_anomaly(2.0, 0.9, tolerance=0.0)
    assert abs(kepler_residual(e_anomaly, 2.0, 0.9)) < 1e-12


def test_bisection_fallback_brackets_and_converges():
    for eccentricity in (0.0, 0.5, 0.999999):
        for mean_anomaly in (0.0, 1e-12, 1.0, math.pi, TWO_PI - 1e-12):
            e_anomaly = km._bisect_eccentric_anomaly(mean_anomaly, eccentricity, 1e-12)
            assert 0.0 <= e_anomaly <= TWO_PI
            assert abs(kepler_residual(e_anomaly, mean_anomaly, eccentricity)) < 1e-9


@pytest.mark.parametrize("eccentricity", [1.0, 1.5, -0.01, -1.0])
def test_solver_rejects_eccentricity_outside_0_to_1(eccentricity):
    with pytest.raises(ValueError, match="eccentricity"):
        km.solve_eccentric_anomaly(1.0, eccentricity)


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
def test_solver_rejects_non_finite_mean_anomaly(value):
    with pytest.raises(ValueError, match="must be a finite number"):
        km.solve_eccentric_anomaly(value, 0.5)


@pytest.mark.parametrize("eccentricity", [0.991, 0.999, 0.999999])
def test_elliptical_distance_stays_between_periapsis_and_apoapsis(eccentricity):
    for step in range(-12, 25):
        mean_anomaly = step * 0.55
        nu, r = km.true_anomaly_and_distance_elliptical(mean_anomaly, eccentricity, 10.0)
        assert -1e-12 <= nu <= TWO_PI + 1e-12
        assert 10.0 * (1 - eccentricity) * (1 - 1e-9) <= r <= 10.0 * (1 + eccentricity) * (1 + 1e-9)


def test_elliptical_negative_mean_anomaly_mirrors_positive():
    nu_pos, r_pos = km.true_anomaly_and_distance_elliptical(0.4, 0.995, 3.0)
    nu_neg, r_neg = km.true_anomaly_and_distance_elliptical(-0.4, 0.995, 3.0)
    assert nu_neg == pytest.approx(TWO_PI - nu_pos, abs=1e-9)
    assert r_neg == pytest.approx(r_pos, rel=1e-9)


# --- _real_cube_root ------------------------------------------------------

def test_real_cube_root_of_zero_is_zero():
    assert km._real_cube_root(0.0) == 0.0
    assert km._real_cube_root(0) == 0.0


def test_real_cube_root_keeps_the_sign_of_negative_zero():
    assert math.copysign(1.0, km._real_cube_root(-0.0)) == -1.0


@pytest.mark.parametrize("value", [-1.0, -8.0, -27.0, -1e-30, -1e30, -0.125])
def test_real_cube_root_of_negative_is_real_and_negative(value):
    root = km._real_cube_root(value)
    assert isinstance(root, float)
    assert root < 0
    assert root ** 3 == pytest.approx(value, rel=1e-12)
    assert root == pytest.approx(-km._real_cube_root(-value), rel=1e-15)


@pytest.mark.parametrize("value", [1e-30, 0.125, 1.0, 8.0, 1e30])
def test_real_cube_root_of_positive_matches_power(value):
    assert km._real_cube_root(value) == pytest.approx(value ** (1 / 3), rel=1e-15)


def test_real_cube_root_is_monotonic_through_zero():
    values = [-1e6, -8.0, -1.0, -1e-9, 0.0, 1e-9, 1.0, 8.0, 1e6]
    roots = [km._real_cube_root(v) for v in values]
    assert roots == sorted(roots)


def test_barker_solution_is_finite_at_extreme_negative_anomalies():
    for mp in (-1e200, -1e12, -1.0, -1e-12):
        d = km.solve_barker_equation(mp)
        assert math.isfinite(d) and d <= 0
