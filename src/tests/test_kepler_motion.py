"""
Kepler/Barker orbital-mechanics regression tests
=================================================

Covers `stellarObjects/keplerMotion.py`'s two-body propagation for
eccentric (elliptical) and parabolic orbits -- the physics
`cometData.Comet` relies on for realistic, non-uniform-angular-speed
orbital motion (see docs/design/comet-orbital-realism.md).

Run with: pytest src/tests/test_kepler_motion.py
"""

import math

import pytest

from stellarObjects import keplerMotion


# ---------------------------------------------------------------------------
# Kepler's equation (elliptical).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("eccentricity", [0.0, 0.3, 0.7, 0.9, 0.99])
@pytest.mark.parametrize("mean_anomaly_deg", [0.0, 45.0, 90.0, 180.0, 270.0, 359.0])
def test_solve_eccentric_anomaly_satisfies_keplers_equation(eccentricity, mean_anomaly_deg):
    mean_anomaly_rad = math.radians(mean_anomaly_deg)
    eccentric_anomaly = keplerMotion.solve_eccentric_anomaly(mean_anomaly_rad, eccentricity)
    # M = E - e*sin(E), by definition -- the solved E must satisfy this to
    # within the solver's own tolerance.
    reconstructed_m = eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly)
    assert math.isclose(reconstructed_m % (2 * math.pi), mean_anomaly_rad % (2 * math.pi), abs_tol=1e-8)


def test_elliptical_at_perihelion_gives_zero_true_anomaly_and_min_distance():
    semi_major_axis_au = 10.0
    eccentricity = 0.8
    true_anomaly, distance_au = keplerMotion.true_anomaly_and_distance_elliptical(0.0, eccentricity, semi_major_axis_au)
    assert math.isclose(true_anomaly, 0.0, abs_tol=1e-9)
    assert math.isclose(distance_au, semi_major_axis_au * (1 - eccentricity), rel_tol=1e-9)


def test_elliptical_at_aphelion_gives_max_distance():
    semi_major_axis_au = 10.0
    eccentricity = 0.8
    true_anomaly, distance_au = keplerMotion.true_anomaly_and_distance_elliptical(math.pi, eccentricity, semi_major_axis_au)
    assert math.isclose(true_anomaly, math.pi, abs_tol=1e-9)
    assert math.isclose(distance_au, semi_major_axis_au * (1 + eccentricity), rel_tol=1e-9)


def test_circular_orbit_moves_at_uniform_angular_speed():
    # At e=0, Kepler's equation degenerates to true_anomaly == mean_anomaly
    # exactly (no equation-of-center correction), and distance is constant.
    semi_major_axis_au = 5.0
    for mean_anomaly_deg in (0, 60, 130, 250, 340):
        mean_anomaly_rad = math.radians(mean_anomaly_deg)
        true_anomaly, distance_au = keplerMotion.true_anomaly_and_distance_elliptical(
            mean_anomaly_rad, 0.0, semi_major_axis_au
        )
        assert math.isclose(true_anomaly % (2 * math.pi), mean_anomaly_rad % (2 * math.pi), abs_tol=1e-9)
        assert math.isclose(distance_au, semi_major_axis_au, rel_tol=1e-9)


def test_eccentric_orbit_moves_faster_near_perihelion_than_aphelion():
    # Kepler's second law: equal areas in equal time, so a small step in
    # MEAN anomaly (linear in time) produces a much larger step in TRUE
    # anomaly (actual angular position) near perihelion than near aphelion.
    semi_major_axis_au = 10.0
    eccentricity = 0.9
    step_rad = math.radians(1.0)

    true_near_perihelion_start, _ = keplerMotion.true_anomaly_and_distance_elliptical(0.0, eccentricity, semi_major_axis_au)
    true_near_perihelion_end, _ = keplerMotion.true_anomaly_and_distance_elliptical(step_rad, eccentricity, semi_major_axis_au)
    delta_near_perihelion = true_near_perihelion_end - true_near_perihelion_start

    true_near_aphelion_start, _ = keplerMotion.true_anomaly_and_distance_elliptical(math.pi, eccentricity, semi_major_axis_au)
    true_near_aphelion_end, _ = keplerMotion.true_anomaly_and_distance_elliptical(math.pi + step_rad, eccentricity, semi_major_axis_au)
    delta_near_aphelion = true_near_aphelion_end - true_near_aphelion_start

    assert delta_near_perihelion > delta_near_aphelion > 0


# ---------------------------------------------------------------------------
# Barker's equation (parabolic).
# ---------------------------------------------------------------------------

def test_barker_equation_at_perihelion_gives_zero():
    assert math.isclose(keplerMotion.solve_barker_equation(0.0), 0.0, abs_tol=1e-12)


@pytest.mark.parametrize("parabolic_mean_anomaly", [-5.0, -1.0, -0.1, 0.1, 1.0, 5.0])
def test_barker_equation_satisfies_its_own_cubic(parabolic_mean_anomaly):
    d = keplerMotion.solve_barker_equation(parabolic_mean_anomaly)
    reconstructed = d ** 3 + 3 * d
    assert math.isclose(reconstructed, 3 * parabolic_mean_anomaly, rel_tol=1e-9, abs_tol=1e-9)


def test_barker_equation_is_antisymmetric():
    # D(-Mp) == -D(Mp): a parabolic orbit is symmetric in time around
    # perihelion passage.
    for mp in (0.5, 2.0, 10.0):
        assert math.isclose(keplerMotion.solve_barker_equation(-mp), -keplerMotion.solve_barker_equation(mp), rel_tol=1e-9)


def test_parabolic_at_perihelion_gives_zero_true_anomaly_and_min_distance():
    perihelion_distance_au = 0.5
    true_anomaly, distance_au = keplerMotion.true_anomaly_and_distance_parabolic(0.0, perihelion_distance_au)
    assert math.isclose(true_anomaly, 0.0, abs_tol=1e-9)
    assert math.isclose(distance_au, perihelion_distance_au, rel_tol=1e-9)


def test_parabolic_distance_increases_away_from_perihelion():
    perihelion_distance_au = 0.5
    _, distance_at_perihelion = keplerMotion.true_anomaly_and_distance_parabolic(0.0, perihelion_distance_au)
    _, distance_before = keplerMotion.true_anomaly_and_distance_parabolic(-2.0, perihelion_distance_au)
    _, distance_after = keplerMotion.true_anomaly_and_distance_parabolic(2.0, perihelion_distance_au)
    assert distance_before > distance_at_perihelion
    assert distance_after > distance_at_perihelion
    # Symmetric around perihelion passage.
    assert math.isclose(distance_before, distance_after, rel_tol=1e-9)


def test_parabolic_mean_anomaly_matches_definition():
    primary_mass_solar = 1.0
    perihelion_distance_au = 1.0
    time_years = 10.0
    mp = keplerMotion.parabolic_mean_anomaly(time_years, perihelion_distance_au, primary_mass_solar)
    mu = keplerMotion.gravitational_parameter_au3_yr2(primary_mass_solar)
    expected = math.sqrt(mu / (2 * perihelion_distance_au ** 3)) * time_years
    assert math.isclose(mp, expected, rel_tol=1e-12)


# ---------------------------------------------------------------------------
# Vis-viva speed.
# ---------------------------------------------------------------------------

def test_vis_viva_speed_matches_circular_orbital_speed():
    from stellarObjects.utils import circular_orbital_speed_kms

    semi_major_axis_au = 3.0
    primary_mass_solar = 1.0
    period_years = math.sqrt(semi_major_axis_au ** 3 / primary_mass_solar)
    expected_speed_kms = circular_orbital_speed_kms(semi_major_axis_au, period_years)

    speed_kms = keplerMotion.vis_viva_speed_kms(semi_major_axis_au, semi_major_axis_au, primary_mass_solar)
    assert math.isclose(speed_kms, expected_speed_kms, rel_tol=1e-6)


def test_vis_viva_speed_is_faster_at_perihelion_than_aphelion():
    semi_major_axis_au = 10.0
    eccentricity = 0.8
    primary_mass_solar = 1.0
    perihelion_au = semi_major_axis_au * (1 - eccentricity)
    aphelion_au = semi_major_axis_au * (1 + eccentricity)

    speed_at_perihelion = keplerMotion.vis_viva_speed_kms(perihelion_au, semi_major_axis_au, primary_mass_solar)
    speed_at_aphelion = keplerMotion.vis_viva_speed_kms(aphelion_au, semi_major_axis_au, primary_mass_solar)
    assert speed_at_perihelion > speed_at_aphelion > 0


def test_parabolic_semi_major_axis_of_infinity_gives_escape_speed():
    distance_au = 1.0
    primary_mass_solar = 1.0
    speed_kms = keplerMotion.vis_viva_speed_kms(distance_au, math.inf, primary_mass_solar)
    mu = keplerMotion.gravitational_parameter_au3_yr2(primary_mass_solar)
    expected_au_per_year = math.sqrt(2 * mu / distance_au)
    from stellarObjects import physical_constants
    expected_kms = expected_au_per_year * physical_constants.AU_TO_KM / physical_constants.SECONDS_PER_YEAR
    assert math.isclose(speed_kms, expected_kms, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# comet_orbital_state -- full position/speed dispatch.
# ---------------------------------------------------------------------------

def test_comet_orbital_state_elliptical_never_closer_than_perihelion():
    perihelion_distance_au = 1.0
    eccentricity = 0.6
    for mean_anomaly_deg in range(0, 360, 15):
        state = keplerMotion.comet_orbital_state(
            "elliptical", perihelion_distance_au, eccentricity, inclination_deg=0.0,
            arg_periapsis_deg=0.0, ascending_node_deg=0.0, primary_mass_solar=1.0,
            mean_anomaly_rad=math.radians(mean_anomaly_deg),
        )
        assert state["distance_au"] >= perihelion_distance_au - 1e-9


def test_comet_orbital_state_parabolic_never_closer_than_perihelion():
    perihelion_distance_au = 0.3
    for mp in (-5, -1, 0, 1, 5):
        state = keplerMotion.comet_orbital_state(
            "parabolic", perihelion_distance_au, eccentricity=1.0, inclination_deg=20.0,
            arg_periapsis_deg=40.0, ascending_node_deg=80.0, primary_mass_solar=1.0,
            parabolic_mean_anomaly_value=mp,
        )
        assert state["distance_au"] >= perihelion_distance_au - 1e-9


def test_comet_orbital_state_position_magnitude_matches_distance():
    state = keplerMotion.comet_orbital_state(
        "elliptical", perihelion_distance_au=2.0, eccentricity=0.5, inclination_deg=30.0,
        arg_periapsis_deg=50.0, ascending_node_deg=70.0, primary_mass_solar=1.2,
        mean_anomaly_rad=math.radians(123.0),
    )
    magnitude = math.sqrt(state["position_x_au"] ** 2 + state["position_y_au"] ** 2 + state["position_z_au"] ** 2)
    assert math.isclose(magnitude, state["distance_au"], rel_tol=1e-9)


def test_comet_orbital_state_requires_matching_anomaly_argument():
    with pytest.raises(ValueError):
        keplerMotion.comet_orbital_state(
            "elliptical", perihelion_distance_au=1.0, eccentricity=0.5, inclination_deg=0.0,
            arg_periapsis_deg=0.0, ascending_node_deg=0.0, primary_mass_solar=1.0,
        )
    with pytest.raises(ValueError):
        keplerMotion.comet_orbital_state(
            "parabolic", perihelion_distance_au=1.0, eccentricity=1.0, inclination_deg=0.0,
            arg_periapsis_deg=0.0, ascending_node_deg=0.0, primary_mass_solar=1.0,
        )


def test_comet_orbital_state_rejects_unknown_orbit_type():
    with pytest.raises(ValueError):
        keplerMotion.comet_orbital_state(
            "hyperbolic", perihelion_distance_au=1.0, eccentricity=1.5, inclination_deg=0.0,
            arg_periapsis_deg=0.0, ascending_node_deg=0.0, primary_mass_solar=1.0,
        )
