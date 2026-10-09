# tests/test_orbital_limits.py

"""
GEN.108: each place the orbital maths breaks down has a guard, and each
guard is hit here (docs/design/orbital-updates.md section 7).
"""

import math

import pytest

from planetgen.db import store
from planetgen.galaxy.galactic_orbit import calculate_galactic_orbit
from planetgen.physics import constants as pc
from planetgen.physics import kepler
from planetgen.physics import state_vectors as sv
from planetgen.physics.orbits import classify_encounter, roche_limit_m
from planetgen.physics.units import ly_to_pc

MU_SUN = 4 * math.pi ** 2  # AU^3 / yr^2


def _state(q, e, nu, mu=MU_SUN):
    return sv.state_from_elements(mu, q, e, 0.4, 1.1, 2.3, nu)


def _norm(v):
    return math.sqrt(sum(c * c for c in v))


# --- Near-parabolic orbits ---------------------------------------------------------------

def test_a_near_parabolic_radius_keeps_its_digits_at_perihelion():
    q, e = 0.5, 1 - 1e-12
    a = q / (1 - e)
    eccentric = 1e-7
    mean = eccentric - e * math.sin(eccentric)
    _nu, r = kepler.true_anomaly_and_distance_elliptical(mean, e, a, eccentric_anomaly_rad=eccentric)
    reference = q + a * e * (eccentric ** 2 / 2 - eccentric ** 4 / 24)
    assert r == pytest.approx(reference, rel=1e-12)
    naive = a * (1 - e * math.cos(eccentric))
    assert abs(naive - reference) / reference > 1e-6  # what the old form lost


@pytest.mark.parametrize("q, e, nu, dt", [
    (1.0, 0.0167, 0.1, 0.37),        # Earth-like
    (1.0, 0.5, 2.0, -3.3),           # backwards
    (0.5, 0.999, 0.01, 0.2),         # near-parabolic ellipse
])
def test_the_universal_step_matches_keplers_equation_on_ellipses(q, e, nu, dt):
    position, velocity = _state(q, e, nu)
    moved, _ = kepler.universal_step(position, velocity, MU_SUN, dt)
    a = q / (1 - e)
    mean = sv.mean_anomaly_from_true(nu, e) + math.sqrt(MU_SUN / a ** 3) * dt
    nu_after, _ = kepler.true_anomaly_and_distance_elliptical(mean, e, a)
    expected, _ = _state(q, e, nu_after)
    assert max(abs(x - y) for x, y in zip(moved, expected)) < 1e-9 * max(1.0, _norm(expected))


def test_past_the_kepler_limit_a_comet_moves_by_the_universal_variable():
    q, e = 0.5, 1 - 1e-9
    assert e > kepler.KEPLER_MAX_ECCENTRICITY
    a = q / (1 - e)
    for since_perihelion in (-0.3, 1e-4, 0.05, 2.0):
        mean = since_perihelion * math.sqrt(MU_SUN / a ** 3)
        state = kepler.comet_orbital_state("elliptical", q, e, 0.0, 0.0, 0.0, 1.0, mean_anomaly_rad=mean)
        # So close to a parabola, Barker's equation is right to about 1e-9.
        _nu, parabolic_r = kepler.true_anomaly_and_distance_parabolic(
            kepler.parabolic_mean_anomaly(since_perihelion, q, 1.0), q)
        assert state["distance_au"] == pytest.approx(parabolic_r, rel=1e-6)
        position, _ = sv.state_from_elements(MU_SUN, q, e, 0.0, 0.0, 0.0, 0.0)
        moved, _ = kepler.universal_step(position, (0.0, math.sqrt(MU_SUN * (1 + e) / q), 0.0), MU_SUN,
                                         since_perihelion)
        assert state["position_x_au"] == pytest.approx(moved[0], rel=1e-9, abs=1e-12)


@pytest.mark.parametrize("nu, dt", [(0.3, 2.0), (-1.0, -2.0), (0.0, 1e-6)])
def test_the_universal_step_matches_barkers_equation_on_a_parabola(nu, dt):
    q = 0.5
    position, velocity = _state(q, 1.0, nu)
    moved, _ = kepler.universal_step(position, velocity, MU_SUN, dt)
    d = math.tan(nu / 2)
    since_perihelion = (d + d ** 3 / 3) / math.sqrt(MU_SUN / (2 * q ** 3))
    nu_after, _ = kepler.true_anomaly_and_distance_parabolic(
        kepler.parabolic_mean_anomaly(since_perihelion + dt, q, 1.0), q)
    expected, _ = _state(q, 1.0, nu_after)
    assert max(abs(x - y) for x, y in zip(moved, expected)) < 1e-10 * max(1.0, _norm(expected))


@pytest.mark.parametrize("e, nu, dt", [(1.5, 0.3, 5.0), (3.0, -0.5, -1.0), (1 + 1e-9, 0.2, 0.3)])
def test_the_universal_step_keeps_a_hyperbola_on_its_orbit(e, nu, dt):
    position, velocity = _state(0.5, e, nu)
    moved, moved_velocity = kepler.universal_step(position, velocity, MU_SUN, dt)
    after = sv.elements_from_state(moved, moved_velocity, MU_SUN)
    assert after["eccentricity"] == pytest.approx(e, rel=1e-9)
    assert after["periapsis_distance"] == pytest.approx(0.5, rel=1e-9)
    back, _ = kepler.universal_step(moved, moved_velocity, MU_SUN, -dt)
    assert max(abs(x - y) for x, y in zip(back, position)) < 1e-9


def test_the_universal_step_is_a_proper_two_body_map():
    """f g' - f' g = 1 (the Lagrange coefficients' Wronskian), which the
    research document's uncorrected equation broke by up to 15%."""
    position, velocity = _state(1.0, 0.6, 0.7)
    dt = 0.9
    moved, moved_velocity = kepler.universal_step(position, velocity, MU_SUN, dt)
    h0 = sv._cross(position, velocity)
    h1 = sv._cross(moved, moved_velocity)
    assert h1 == pytest.approx(h0, rel=1e-12)  # angular momentum kept, so the map is area-preserving
    energy = [0.5 * sv._dot(v, v) - MU_SUN / _norm(r) for r, v in ((position, velocity), (moved, moved_velocity))]
    assert energy[1] == pytest.approx(energy[0], rel=1e-12)


def test_the_universal_step_refuses_what_it_cannot_move():
    with pytest.raises(ValueError):
        kepler.universal_step((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), MU_SUN, 1.0)
    with pytest.raises(ValueError):
        kepler.universal_step((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), 0.0, 1.0)
    with pytest.raises(ValueError):
        kepler.universal_step((1.0, 0.0, 0.0), (0.0, math.nan, 0.0), MU_SUN, 1.0)


def test_the_stumpff_series_joins_the_closed_forms():
    for z in (1e-3, -1e-3):
        series = kepler.stumpff_c2_c3(z * 0.999999)
        closed = kepler.stumpff_c2_c3(z * 1.000001)
        assert series == pytest.approx(closed, rel=1e-8)


# --- Circular and equatorial orbits -----------------------------------------------------

@pytest.mark.parametrize("e, i", [(0.0, 0.0), (1e-12, 1e-12), (0.0, 1.2), (0.3, 0.0), (0.3, 0.5), (1.4, 0.1)])
def test_equinoctial_elements_stay_exact_where_classical_ones_are_undefined(e, i):
    position, velocity = sv.state_from_elements(1.0, 1.0, e, i, 0.8, 1.9, 0.6)
    elements = sv.equinoctial_from_state(position, velocity, 1.0)
    again = sv.state_from_equinoctial(elements["p"], elements["f"], elements["g"], elements["h"], elements["k"],
                                      elements["L"], 1.0)
    assert max(abs(x - y) for x, y in zip(position + velocity, again[0] + again[1])) < 1e-14
    assert math.hypot(elements["f"], elements["g"]) == pytest.approx(e, abs=1e-14)


def test_equinoctial_elements_refuse_a_retrograde_equatorial_orbit():
    position, velocity = sv.state_from_elements(1.0, 1.0, 0.1, math.pi, 0.0, 0.0, 0.3)
    with pytest.raises(ValueError):
        sv.equinoctial_from_state(position, velocity, 1.0)


# --- Steps longer than an orbit -------------------------------------------------------------

def test_a_galactic_step_of_many_orbits_wraps_exactly():
    period_gy = 0.25
    rest = 1.234e6
    many = 4e6 * period_gy * 1e9 + rest
    assert store._galactic_turn(many, period_gy) == pytest.approx(store._galactic_turn(rest, period_gy), rel=1e-6)
    assert store._galactic_turn(many, 0.0) == 0.0


def test_the_sql_phase_turn_drops_whole_orbits(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        period = 0.001  # a close moon's year
        orbits = 1e9  # a million years
        end = orbits * period * pc.SECONDS_PER_YEAR + 0.25 * period * pc.SECONDS_PER_YEAR
        row = conn.execute(f"SELECT {store._phase_turn_sql('period_years')} AS turn"
                           " FROM (SELECT 0 AS epoch_unix, ? AS period_years) t", (end, period)).fetchone()
        assert float(row["turn"]) == pytest.approx(90.0, abs=1e-4)
    finally:
        conn.close()


def test_a_universal_step_of_many_orbits_lands_where_the_remainder_does():
    position, velocity = _state(1.0, 0.3, 0.4)
    period = 2 * math.pi * math.sqrt((1.0 / 0.7) ** 3 / MU_SUN)
    far, _ = kepler.universal_step(position, velocity, MU_SUN, 1000 * period + 0.3)
    near, _ = kepler.universal_step(position, velocity, MU_SUN, 0.3)
    assert max(abs(x - y) for x, y in zip(far, near)) < 1e-9


# --- Close passes ------------------------------------------------------------------------

def test_a_close_pass_is_a_collision_a_disruption_a_capture_or_a_flyby():
    earth_radius, moon_radius = 6.371e6, 1.7374e6
    roche = roche_limit_m(earth_radius, 5514, 3344)
    assert roche == pytest.approx(1.84e7, rel=0.01)  # 18,400 km for the Moon
    hill = 1.5e9
    assert classify_encounter(7.0e6, earth_radius + moon_radius, roche, hill) == "collision"
    assert classify_encounter(1.0e7, earth_radius + moon_radius, roche, hill) == "disruption"
    assert classify_encounter(3.84e8, earth_radius + moon_radius, roche, hill) == "captured"
    assert classify_encounter(5.0e9, earth_radius + moon_radius, roche, hill) == "flyby"
    with pytest.raises(ValueError):
        classify_encounter(-1.0, 1.0, 1.0, 1.0)


# --- The galactic centre -------------------------------------------------------------------

def test_the_galactic_centre_gives_no_infinite_turn():
    core_pc = pc.GALACTIC_ROTATION_CORE_RADIUS_PC
    speed = pc.GALACTIC_ROTATION_FLAT_VELOCITY_KMS
    core_period_gy = 2 * math.pi * core_pc * pc.PARSEC_M / (speed * 1000) / pc.SECONDS_PER_YEAR / 1e9
    for distance_ly in (1e-6, 1e-3, 1.0):
        speed_kms, period_gy = calculate_galactic_orbit(distance_ly)
        assert speed_kms < speed * ly_to_pc(distance_ly) / core_pc * 1.0001
        assert period_gy == pytest.approx(core_period_gy, rel=1e-6)
    assert calculate_galactic_orbit(0.0) == (0.0, 0.0)
    assert store._galactic_turn(1e6, 0.0) == 0.0  # the exact centre doesn't turn
