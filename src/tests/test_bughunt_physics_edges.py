# tests/test_bughunt_physics_edges.py

"""
Tier 1 bug-hunt coverage: degenerate/out-of-contract inputs to the pure
physics functions in `planetPhysics.py`/`keplerMotion.py` -- zero/negative/
near-zero mass, radius, distance, and out-of-range eccentricity. These
functions are never called with such inputs by any real generation path
(every caller's own value is drawn from an always-positive/always-in-range
program_constants table), so a crash here is not a live bug for a normal
run -- but two real, if latent, defects were found and fixed while writing
this file (see the two "regression" tests below): a `ZeroDivisionError`/
`math domain error` where the module's own established convention
(`calculate_surface_gravity`) is a clean `ValueError`, and, worse, an
out-of-range eccentricity silently returning a plausible-looking but
physically meaningless number instead of failing loudly at all. Kept as a
permanent regression suite so a future refactor that starts calling these
functions with less-trusted input (e.g. a `--system-file`-driven orbital
override) inherits the same clean failure mode instead of quietly
corrupting output.
"""

import math

import pytest

from stellarObjects import keplerMotion as km
from stellarObjects import planetPhysics as pp


# --- calculate_orbital_period_years -----------------------------------

def test_orbital_period_zero_mass_raises_clean_value_error():
    """Regression: previously a bare ZeroDivisionError; see the guard
    clause added to calculate_orbital_period_years."""
    with pytest.raises(ValueError):
        pp.calculate_orbital_period_years(1.0, 0.0)


def test_orbital_period_negative_mass_raises_clean_value_error():
    with pytest.raises(ValueError):
        pp.calculate_orbital_period_years(1.0, -1.989e30)


def test_orbital_period_negative_distance_raises_clean_value_error():
    """Regression: previously a bare `ValueError: math domain error` from
    `math.sqrt` of a negative `distance_au ** 3` -- still a ValueError, but
    now with a message that actually names the problem."""
    with pytest.raises(ValueError):
        pp.calculate_orbital_period_years(-1.0, 1.989e30)


def test_orbital_period_zero_distance_raises_clean_value_error():
    with pytest.raises(ValueError):
        pp.calculate_orbital_period_years(0.0, 1.989e30)


def test_orbital_period_valid_inputs_still_finite_and_positive():
    period = pp.calculate_orbital_period_years(1.0, 1.989e30)
    assert period > 0 and math.isfinite(period)


# --- calculate_surface_gravity ------------------------------------------

def _make_planet_stub(**overrides):
    """A minimal object exposing only the attributes
    calculate_surface_gravity/calculate_atmospheric_conditions actually
    read/write -- avoids constructing a full Planet (which requires a
    whole SystemConfig/Star/generation pipeline) just to fuzz these two
    pure-ish calculation functions in isolation."""
    from types import SimpleNamespace
    defaults = dict(radius=6371.0, mass=5.972e24, gravity=None)
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


def test_surface_gravity_zero_radius_raises():
    planet = _make_planet_stub(radius=0.0)
    with pytest.raises((ValueError, ZeroDivisionError)):
        pp.calculate_surface_gravity(planet)


def test_surface_gravity_zero_mass_raises_value_error():
    """Already guarded (`raise ValueError('Invalid value for gravity.')`)
    -- confirms it stays that way."""
    planet = _make_planet_stub(mass=0.0)
    with pytest.raises(ValueError):
        pp.calculate_surface_gravity(planet)


def test_surface_gravity_negative_mass_raises_value_error():
    planet = _make_planet_stub(mass=-5.972e24)
    with pytest.raises(ValueError):
        pp.calculate_surface_gravity(planet)


def test_surface_gravity_earth_like_is_close_to_one_g():
    planet = _make_planet_stub()
    pp.calculate_surface_gravity(planet)
    assert 0.9 < planet.gravity < 1.1


# --- keplerMotion: mean_motion_per_year ---------------------------------

def test_mean_motion_zero_mass_raises():
    with pytest.raises((ValueError, ZeroDivisionError)):
        km.mean_motion_per_year(1.0, 0.0)


def test_mean_motion_negative_mass_raises():
    with pytest.raises((ValueError, ZeroDivisionError)):
        km.mean_motion_per_year(1.0, -1.0)


def test_mean_motion_valid_input_is_finite_positive():
    n = km.mean_motion_per_year(1.0, 1.0)
    assert n > 0 and math.isfinite(n)


# --- keplerMotion: solve_eccentric_anomaly / true_anomaly_and_distance_elliptical --

@pytest.mark.parametrize("bad_eccentricity", [1.0, 1.2, -0.1, -1.0])
def test_true_anomaly_elliptical_rejects_out_of_range_eccentricity(bad_eccentricity):
    """Regression: e=1.0 (exactly parabolic) and e=-0.1 (unphysical)
    previously returned a silently-wrong-but-plausible-looking
    (true_anomaly, distance) pair instead of failing -- the worse kind of
    bug, since nothing downstream would ever notice. See the guard clause
    added to true_anomaly_and_distance_elliptical."""
    with pytest.raises(ValueError):
        km.true_anomaly_and_distance_elliptical(0.5, bad_eccentricity, 1.0)


@pytest.mark.parametrize("eccentricity", [0.0, 0.3, 0.7, 0.9, 0.999, 0.999999])
def test_true_anomaly_elliptical_valid_eccentricities_are_finite(eccentricity):
    true_anomaly, distance = km.true_anomaly_and_distance_elliptical(0.5, eccentricity, 1.0)
    assert math.isfinite(true_anomaly)
    assert math.isfinite(distance) and distance > 0


def test_solve_eccentric_anomaly_converges_across_mean_anomaly_range():
    """Sweep the full [0, 2*pi) mean-anomaly range at a high (but valid)
    eccentricity -- Kepler's equation solver convergence degrades as
    e -> 1 (see the function's own docstring), so this is the region most
    likely to silently fail to converge within max_iterations."""
    for mean_anomaly in [i * 0.1 for i in range(0, 63)]:
        e_anomaly = km.solve_eccentric_anomaly(mean_anomaly, 0.999)
        assert math.isfinite(e_anomaly)
        # Kepler's equation itself, satisfied to a loose tolerance --
        # confirms convergence, not just finiteness.
        residual = e_anomaly - 0.999 * math.sin(e_anomaly) - (mean_anomaly % km.TWO_PI)
        assert abs(residual) < 1e-6, f"mean_anomaly={mean_anomaly}: residual {residual} too large, solver did not converge"


# --- keplerMotion: parabolic_mean_anomaly / vis_viva_speed_kms ----------

def test_parabolic_mean_anomaly_zero_perihelion_raises():
    with pytest.raises((ValueError, ZeroDivisionError)):
        km.parabolic_mean_anomaly(1.0, 0.0, 1.0)


def test_parabolic_mean_anomaly_negative_perihelion_raises():
    with pytest.raises((ValueError, ZeroDivisionError)):
        km.parabolic_mean_anomaly(1.0, -1.0, 1.0)


def test_vis_viva_zero_distance_raises():
    with pytest.raises((ValueError, ZeroDivisionError)):
        km.vis_viva_speed_kms(0.0, 1.0, 1.0)


def test_vis_viva_distance_beyond_apoapsis_raises_not_silently_wrong():
    """r > 2a is unreachable for any genuinely bound (e < 1) orbit, but if
    it's ever passed anyway (a caller bug elsewhere), the vis-viva formula
    should fail loudly (negative sqrt argument -> ValueError), not return
    a silently-wrong complex-discarded value."""
    with pytest.raises(ValueError):
        km.vis_viva_speed_kms(100.0, 1.0, 1.0)


def test_vis_viva_parabolic_semi_major_axis_infinite_is_finite_speed():
    speed = km.vis_viva_speed_kms(1.0, math.inf, 1.0)
    assert speed > 0 and math.isfinite(speed)


# --- comet_orbital_state: dispatch validation ----------------------------

def test_comet_orbital_state_unknown_orbit_type_raises_value_error():
    with pytest.raises(ValueError):
        km.comet_orbital_state(
            "hyperbolic", 1.0, 0.5, 0.0, 0.0, 0.0, 1.0, mean_anomaly_rad=0.0,
        )


def test_comet_orbital_state_elliptical_missing_mean_anomaly_raises():
    with pytest.raises(ValueError):
        km.comet_orbital_state(
            "elliptical", 1.0, 0.5, 0.0, 0.0, 0.0, 1.0,
        )


def test_comet_orbital_state_parabolic_missing_anomaly_raises():
    with pytest.raises(ValueError):
        km.comet_orbital_state(
            "parabolic", 1.0, 0.999, 0.0, 0.0, 0.0, 1.0,
        )
