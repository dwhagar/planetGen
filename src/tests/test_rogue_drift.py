# tests/test_rogue_drift.py

"""
GEN.1114: the closed-form drift of an ejected planet (`planetgen.physics.rogue_drift`): the Hill solution against a direct nonlinear orbit integration in a flat-rotation-curve
potential, the pieces of the solution (bounded radial and vertical kicks, the slide along the orbit), the velocity as the derivative of the position, the ring form for large
drifts, the design note's median displacements, and the inputs refused.
"""

import math

import pytest

from planetgen.physics import rogue_drift as rd

KM_S_PC_MYR = rd._KM_S_TO_PC_MYR
R0_PC = 8000.0
V0_KM_S = 230.0


def _flat_galaxy():
    omega = V0_KM_S / R0_PC * 1000.0
    return rd.Galaxy.from_oort(omega / 2.0, -omega / 2.0, 76.0)


def _integrate(kick, time_myr, step=0.05):
    """Star and rogue in a flat-rotation-curve potential (`a = -v0^2 r / r^2`), RK4, in pc and Myr; the rogue's offset in the star's rotating frame, pc."""
    v0 = V0_KM_S * KM_S_PC_MYR

    def derivative(s):
        out = []
        for i in (0, 4):
            x, y, vx, vy = s[i:i + 4]
            f = -v0 * v0 / (x * x + y * y)
            out += [vx, vy, f * x, f * y]
        return out

    u, v = kick
    state = [R0_PC, 0.0, 0.0, v0, R0_PC, 0.0, u * KM_S_PC_MYR, (V0_KM_S + v) * KM_S_PC_MYR]
    for _ in range(round(time_myr / step)):
        k1 = derivative(state)
        k2 = derivative([s + 0.5 * step * k for s, k in zip(state, k1)])
        k3 = derivative([s + 0.5 * step * k for s, k in zip(state, k2)])
        k4 = derivative([s + step * k for s, k in zip(state, k3)])
        state = [s + step / 6.0 * (a + 2 * b + 2 * c + d) for s, a, b, c, d in zip(state, k1, k2, k3, k4)]
    angle = math.atan2(state[1], state[0])
    rel = (state[4] - state[0], state[5] - state[1])
    return (rel[0] * math.cos(angle) + rel[1] * math.sin(angle), -rel[0] * math.sin(angle) + rel[1] * math.cos(angle))


@pytest.mark.parametrize("kick", [(0.0, 3.0), (3.0, 0.0), (2.0, -2.0), (1.5, 2.5), (-1.0, 1.0)])
@pytest.mark.parametrize("time_myr", [10.0, 50.0, 100.0, 200.0])
def test_the_closed_form_agrees_with_a_direct_orbit_integration_inside_a_kiloparsec(kick, time_myr):
    x, y = _integrate(kick, time_myr)
    here = rd.drift((kick[0], kick[1], 0.0), time_myr, _flat_galaxy())
    assert here.separation_pc < 1000.0
    assert here.x_pc == pytest.approx(x, abs=0.1 * max(here.separation_pc, 10.0))
    assert here.y_pc == pytest.approx(y, abs=0.1 * max(here.separation_pc, 10.0))
    assert math.hypot(here.x_pc - x, here.y_pc - y) < 0.08 * max(here.separation_pc, 10.0)


def test_nothing_has_moved_at_the_moment_of_the_kick_and_the_velocity_is_the_kick():
    here = rd.drift((1.0, 2.0, 3.0), 0.0)
    assert (here.x_pc, here.y_pc, here.z_pc) == (0.0, 0.0, 0.0)
    assert (here.vx_km_s, here.vy_km_s, here.vz_km_s) == (1.0, 2.0, 3.0)


def test_the_velocity_is_the_time_derivative_of_the_position():
    kick, t, h = (1.5, 2.5, 1.0), 60.0, 1e-3
    here = rd.drift(kick, t)
    ahead, behind = rd.drift(kick, t + h), rd.drift(kick, t - h)
    for pos, vel in (("x_pc", "vx_km_s"), ("y_pc", "vy_km_s"), ("z_pc", "vz_km_s")):
        rate = (getattr(ahead, pos) - getattr(behind, pos)) / (2.0 * h) / KM_S_PC_MYR
        assert rate == pytest.approx(getattr(here, vel), abs=1e-4)


def test_a_radial_or_vertical_kick_stays_bounded_at_its_epicycle_amplitude():
    galaxy = rd.SOLAR_NEIGHBOURHOOD
    radial = [rd.drift((3.0, 0.0, 0.0), t, galaxy) for t in range(0, 5000, 25)]
    assert max(abs(d.x_pc) for d in radial) <= 1000.0 * 3.0 / galaxy.kappa * 1.001
    assert max(abs(d.x_pc) for d in radial) > 0.99 * 1000.0 * 3.0 / galaxy.kappa * 0.9
    vertical = [rd.drift((0.0, 0.0, 3.0), t, galaxy) for t in range(0, 5000, 25)]
    assert max(abs(d.z_pc) for d in vertical) <= 1000.0 * 3.0 / galaxy.nu * 1.001
    assert all(d.x_pc == 0.0 or abs(d.x_pc) < 1e-9 for d in vertical)


def test_an_along_orbit_kick_slides_at_one_point_two_eight_five_times_the_kick():
    galaxy = rd.SOLAR_NEIGHBOURHOOD
    assert galaxy.slide_ratio == pytest.approx(1.2857, abs=1e-3)
    assert rd.along_orbit_speed_pc_myr(3.0) == pytest.approx(-1.2857 * 3.0 * 1.0227, rel=1e-3)             # 3.9 pc/Myr for 3 km/s
    late = rd.drift((0.0, 3.0, 0.0), 2000.0, galaxy)
    assert late.y_pc == pytest.approx(rd.along_orbit_speed_pc_myr(3.0) * 2000.0, rel=0.05)      # the 4 Omega^2 v / kappa^3 oscillation is 2 percent of it
    # a retrograde kick slides the other way
    assert rd.drift((0.0, -3.0, 0.0), 2000.0, galaxy).y_pc > 0.0


def test_the_guiding_radius_shift_is_two_omega_over_kappa_squared_times_the_kick():
    galaxy = rd.SOLAR_NEIGHBOURHOOD
    ring = rd.ring_drift((0.0, 3.0, 0.0), 10000.0, 8000.0, galaxy)
    assert ring.guiding_radius_shift_pc == pytest.approx(1000.0 * 3.0 * 2.0 * galaxy.omega / galaxy.kappa ** 2)
    # the radial oscillation is about that guiding radius: its mean is the shift
    mean = sum(rd.drift((0.0, 3.0, 0.0), t * 0.5, galaxy).x_pc for t in range(0, 4000)) / 4000.0
    assert mean == pytest.approx(ring.guiding_radius_shift_pc, rel=0.05)


def test_beyond_the_linear_limit_the_drift_is_flagged_and_the_ring_form_wraps_the_azimuth():
    near = rd.drift((0.0, 3.0, 0.0), 100.0)
    far = rd.drift((0.0, 3.0, 0.0), 1000.0)
    assert near.linear and not far.linear
    ring = rd.ring_drift((0.0, 3.0, 0.0), 10000.0, 8000.0)
    assert -math.pi < ring.azimuth_rad <= math.pi
    unwrapped = rd.drift((0.0, 3.0, 0.0), 10000.0).y_pc / 8000.0
    assert ring.azimuth_rad + 2.0 * math.pi * ring.turns == pytest.approx(unwrapped)
    assert ring.turns != 0
    assert rd.ring_drift((0.0, 3.0, 0.0), 1.0, 8000.0).turns == 0


def test_the_donor_lies_the_opposite_way_from_the_rogue():
    kick = (1.0, 2.0, -1.0)
    here = rd.drift(kick, 80.0)
    assert rd.donor_offset_pc(kick, 80.0) == (-here.x_pc, -here.y_pc, -here.z_pc)


def test_the_design_notes_median_displacements_for_a_three_kilometre_a_second_kick():
    # rogue-planet-drift.md 3.2: isotropic kick, median speed 3 km/s: about 30 pc after 10 Myr, 270 pc after 100 Myr (median over directions)
    import random
    rng = random.Random(7)
    for t, low, high in ((10.0, 20.0, 40.0), (100.0, 200.0, 350.0)):
        sep = []
        for _ in range(4000):
            d = [rng.gauss(0, 1) for _ in range(3)]
            sep.append(rd.drift(rd.kick_from_direction(3.0, d), t).separation_pc)
        sep.sort()
        assert low < sep[len(sep) // 2] < high


def test_the_frequencies_of_the_potential_can_stand_in_for_the_oort_constants():
    galaxy = rd.Galaxy.from_frequencies(28.21, 38.89, 74.29)               # galactic-potential.md section 7.1 at 8.128 kpc
    assert galaxy.oort_a == pytest.approx(28.21 - 38.89 ** 2 / (4 * 28.21))
    assert galaxy.slide_ratio > 1.0
    same = rd.Galaxy.from_oort(rd.OORT_A_KM_S_KPC, rd.OORT_B_KM_S_KPC, rd.VERTICAL_FREQUENCY_KM_S_KPC)
    assert same == rd.SOLAR_NEIGHBOURHOOD
    assert rd.SOLAR_NEIGHBOURHOOD.kappa == pytest.approx(36.0, abs=0.05)


def test_a_kick_is_a_speed_in_a_direction():
    assert rd.kick_from_direction(5.0, (0, 3, 4)) == pytest.approx((0.0, 3.0, 4.0))
    assert rd.kick_from_direction(0.0, (1, 0, 0)) == (0.0, 0.0, 0.0)


def test_inputs_the_model_cannot_use_are_refused():
    with pytest.raises(rd.RogueDriftError):
        rd.drift((1.0, 2.0, 3.0), -1.0)
    with pytest.raises(rd.RogueDriftError):
        rd.drift((1.0, 2.0), 1.0)
    with pytest.raises(rd.RogueDriftError):
        rd.drift((1.0, math.nan, 0.0), 1.0)
    with pytest.raises(rd.RogueDriftError):
        rd.ring_drift((1.0, 2.0, 3.0), 1.0, 0.0)
    with pytest.raises(rd.RogueDriftError):
        rd.Galaxy.from_oort(15.3, 11.9)
    with pytest.raises(rd.RogueDriftError):
        rd.Galaxy.from_frequencies(0.0, 30.0, 70.0)
    with pytest.raises(rd.RogueDriftError):
        rd.kick_from_direction(-1.0, (1, 0, 0))
    with pytest.raises(rd.RogueDriftError):
        rd.kick_from_direction(1.0, (0, 0, 0))
