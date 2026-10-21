# tests/test_closest_approach.py

"""GEN.334: the pair closest-approach library (`physics/closest_approach.py`): the sweep, the pair window, the window cap
and the chunking, against hand cases, exact rational arithmetic and dense time sampling."""

import decimal
import math
import random
from decimal import Decimal
from fractions import Fraction

import pytest

from planetgen.physics import closest_approach as ca


def test_a_head_on_pair_hits_at_the_time_the_gap_closes():
    hit, t = ca.sweep((1000.0, 0.0, 0.0), (-10.0, 0.0, 0.0), 100.0, 200.0)
    assert hit and t == pytest.approx(90.0)


def test_a_pair_that_arrives_after_the_window_misses_and_returns_the_window_length():
    assert ca.sweep((1000.0, 0.0, 0.0), (-10.0, 0.0, 0.0), 100.0, 50.0) == (False, 50.0)


def test_a_hit_at_the_window_end_counts():
    hit, t = ca.sweep((1000.0, 0.0, 0.0), (-10.0, 0.0, 0.0), 100.0, 90.0)
    assert hit and t == pytest.approx(90.0)
    assert ca.sweep((1000.0, 0.0, 0.0), (-10.0, 0.0, 0.0), 100.0, 89.9)[0] is False


def test_a_separating_pair_misses():
    assert ca.sweep((1000.0, 0.0, 0.0), (10.0, 0.0, 0.0), 100.0, 1e9) == (False, 1e9)


def test_a_pair_already_overlapping_hits_at_zero_even_when_separating():
    assert ca.sweep((50.0, 0.0, 0.0), (10.0, 0.0, 0.0), 100.0, 10.0) == (True, 0.0)
    assert ca.sweep((100.0, 0.0, 0.0), (0.0, 0.0, 0.0), 100.0, 10.0) == (True, 0.0)


def test_a_tangent_graze_is_a_hit_and_a_miss_by_a_hair_is_not():
    hit, t = ca.sweep((1000.0, 100.0, 0.0), (-10.0, 0.0, 0.0), 100.0, 200.0)
    assert hit and t == pytest.approx(100.0)
    assert ca.sweep((1000.0, 100.0 * (1.0 + 1e-7), 0.0), (-10.0, 0.0, 0.0), 100.0, 200.0)[0] is False


def test_no_relative_speed_and_no_window_are_misses_unless_overlapping():
    assert ca.sweep((1000.0, 0.0, 0.0), (0.0, 0.0, 0.0), 100.0, 1e6) == (False, 1e6)
    assert ca.sweep((1000.0, 0.0, 0.0), (-10.0, 0.0, 0.0), 100.0, 0.0) == (False, 0.0)
    assert ca.sweep((10.0, 0.0, 0.0), (-10.0, 0.0, 0.0), 100.0, 0.0) == (True, 0.0)


def test_bad_input_is_refused():
    for args in (((math.nan, 0, 0), (0, 0, 0), 1.0, 1.0), ((1, 0, 0), (math.inf, 0, 0), 1.0, 1.0),
                 ((1, 0, 0), (0, 0, 0), -1.0, 1.0), ((1, 0, 0), (0, 0, 0), 1.0, -1.0),
                 ((1, 0, 0), (0, 0, 0), math.nan, 1.0)):
        with pytest.raises(ValueError):
            ca.sweep(*args)


def _exact_hit(dp, dv, R, dt):
    """The sweep in exact arithmetic: the first time `|dp + dv t| = R` in `[0, dt]` (80 digits for the one root), or
    None."""
    decimal.getcontext().prec = 80
    dp = [Fraction(v) for v in dp]
    dv = [Fraction(v) for v in dv]
    R, dt = Fraction(R), Fraction(dt)
    d2 = sum(v * v for v in dp)
    if d2 <= R * R:
        return Decimal(0)
    a = sum(v * v for v in dv)
    if a == 0:
        return None
    b = sum(p * v for p, v in zip(dp, dv))
    t_min = -b / a
    if t_min < 0:
        return None
    m2 = d2 - b * b / a
    if m2 > R * R:
        return None
    gap = (R * R - m2) / a
    t_hit = Decimal(t_min.numerator) / Decimal(t_min.denominator) - (Decimal(gap.numerator) / Decimal(gap.denominator)).sqrt()
    return t_hit if t_hit <= Decimal(dt.numerator) / Decimal(dt.denominator) else None


def test_far_starts_where_the_swept_sphere_form_loses_the_miss_distance_are_right():
    """Section 2.1: two Earth-size bodies at 30 km/s over 10,000 years start 0.3 pc away; the closest-approach form
    must agree with exact arithmetic on hit or miss and on the time."""
    rng = random.Random(334)
    r_sum, speed, dt = 1.27e7, 3.0e4, 1.0e4 * 365.25 * 86400.0
    wrong = 0
    for _ in range(300):
        direction = [rng.gauss(0, 1) for _ in range(3)]
        norm = math.sqrt(sum(c * c for c in direction))
        velocity = [-speed * c / norm for c in direction]
        # aim at the target with an offset of 0 to 2 R_sum across the path, then start far enough back
        offset = [rng.uniform(-1, 1) * r_sum for _ in range(3)]
        start = [-v / speed * (speed * dt * rng.uniform(0.2, 0.95)) + o for v, o in zip(velocity, offset)]
        hit, t = ca.sweep(start, velocity, r_sum, dt)
        exact = _exact_hit(start, velocity, r_sum, dt)
        if (exact is not None) != hit or (hit and abs(t - float(exact)) > 1e-9 * dt):
            wrong += 1
    assert wrong == 0


def test_twenty_thousand_random_cases_agree_with_dense_time_sampling():
    rng = random.Random(20000)
    samples = 4000
    for _ in range(20000):
        dp = [rng.uniform(-50.0, 50.0) for _ in range(3)]
        dv = [rng.uniform(-5.0, 5.0) for _ in range(3)]
        R, dt = rng.uniform(0.5, 6.0), rng.uniform(0.5, 20.0)
        hit, t_hit = ca.sweep(dp, dv, R, dt)
        speed = math.sqrt(sum(v * v for v in dv))
        step = dt / samples
        resolution = speed * step
        first = None
        for i in range(samples + 1):
            t = i * step
            if math.dist((dp[0] + dv[0] * t, dp[1] + dv[1] * t, dp[2] + dv[2] * t), (0.0, 0.0, 0.0)) <= R:
                first = t
                break
        if (first is not None) != hit:
            # only a graze finer than the sampling can differ: the closest approach inside the window is within one step
            a = speed * speed
            t_min = 0.0 if a == 0.0 else min(max(-sum(p * v for p, v in zip(dp, dv)) / a, 0.0), dt)
            closest = math.dist([p + v * t_min for p, v in zip(dp, dv)], (0.0, 0.0, 0.0))
            assert abs(closest - R) <= resolution + 1e-9, (dp, dv, R, dt)
        elif hit:
            assert abs(t_hit - first) <= step + 1e-9, (dp, dv, R, dt, t_hit, first)


def test_the_window_is_from_the_later_epoch_to_the_updates_end():
    y = ca.Segment((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 100.0)
    assert ca.pair_window(y, 50.0, 200.0) == (100.0, 200.0)
    assert ca.pair_window(ca.Segment((0, 0, 0), (0, 0, 0), 30.0), 50.0, 200.0) == (50.0, 200.0)
    assert ca.pair_window(ca.Segment((0, 0, 0), (0, 0, 0), 200.0), 50.0, 200.0) == (200.0, 200.0)
    assert ca.pair_window(ca.Segment((0, 0, 0), (0, 0, 0), 201.0), 50.0, 200.0) is None
    with pytest.raises(ValueError):
        ca.pair_window(y, 200.0, 50.0)


def test_two_stored_segments_meet_inside_the_window_from_their_own_epochs():
    # X is at (1000, 0, 0) at t = 100 moving at -10; Y sat at the origin at t = 40 (so its segment is extended to t = 100)
    x = ca.Segment((1000.0, 0.0, 0.0), (-10.0, 0.0, 0.0), 100.0)
    y = ca.Segment((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 40.0)
    start, end = ca.pair_window(y, 100.0, 400.0)
    found = ca.first_hit(ca.linear_state_at(x, y), start, end, 100.0)
    assert found.hit and found.t_hit_s == pytest.approx(190.0)
    short = ca.first_hit(ca.linear_state_at(x, y), start, 150.0, 100.0)
    assert not short.hit and short.t_hit_s == 150.0


def test_positions_are_differenced_from_sector_relative_values():
    """Two bodies 2 km apart across a sector face at 1.5e20 m from the centre, where absolute positions are only good to
    33 km: differencing the origins first (exactly) and then the sector-relative positions keeps the 2 km."""
    edge, far = 2.0 ** 57, 2.0 ** 67
    x = ca.Segment((-edge / 2 + 1000.0, 0.0, 0.0), (-1.0, 0.0, 0.0), 0.0, origin_m=(far + edge, 0.0, 0.0))
    y = ca.Segment((edge / 2 - 1000.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0, origin_m=(far, 0.0, 0.0))
    dp, dv = ca.linear_state_at(x, y)(0.0)
    assert abs(dp[0] - 2000.0) <= 32.0 and dp[1:] == (0.0, 0.0) and dv == (-1.0, 0.0, 0.0)
    naive = ((far + edge) + x.position_m[0]) - (far + y.position_m[0])
    assert abs(naive - 2000.0) > 500.0
    found = ca.first_hit(ca.linear_state_at(x, y), 0.0, 5000.0, 10.0)
    assert found.hit and found.t_hit_s == pytest.approx(1990.0, abs=33.0)


def test_the_window_cap_is_the_tide_or_one_sector():
    year = 365.25 * 86400.0
    omega = 220e3 / (8128.0 * 3.0857e16)            # a flat 220 km/s curve at 8.1 kpc, rad/s
    earth = ca.window_cap(1.27e7, 3.0e4, omega, 1.23e17)
    assert earth == pytest.approx(1800.0 * year, rel=0.35)         # section 5.3: about 1,800 years for two Earths
    small = ca.window_cap(2000.0, 3.0e4, omega, 1.23e17)
    assert 50.0 * year < small < 200.0 * year                      # about 100 years for two 1 km bodies
    assert ca.window_cap(1.27e7, 1.0e9, 0.0, 1.23e17) == pytest.approx(1.23e17 / 1.0e9)      # no tide: one sector
    assert ca.window_cap(1.27e7, 5.6e4, omega * 1e-6, 1.23e17) == pytest.approx(1.23e17 / 5.6e4)   # the sector binds
    assert ca.window_cap(1.27e7, 0.0, omega, 1.23e17) == math.inf
    for bad in ((0.0, 1.0, 1.0, 1.0), (1.0, -1.0, 1.0, 1.0), (1.0, 1.0, -1.0, 1.0), (1.0, 1.0, 1.0, 0.0),
                (math.nan, 1.0, 1.0, 1.0)):
        with pytest.raises(ValueError):
            ca.window_cap(*bad)


def test_a_long_step_is_cut_at_the_cap_and_the_curvature_chunk():
    assert ca.chunks(0.0, 10.0, 100.0, 100.0) == [(0.0, 10.0)]
    assert ca.chunks(5.0, 5.0, 100.0) == []
    pieces = ca.chunks(0.0, 10.0, 3.0, 100.0)
    assert len(pieces) == 4 and pieces[0][0] == 0.0 and pieces[-1][1] == 10.0
    assert all(b - a <= 3.0 + 1e-12 for a, b in pieces) and all(p[1] == q[0] for p, q in zip(pieces, pieces[1:]))
    year = 365.25 * 86400.0
    long = ca.chunks(0.0, 0.35e6 * year, math.inf)
    assert len(long) == 4 and all(b - a <= ca.CURVATURE_CHUNK_S * (1 + 1e-12) for a, b in long)
    for bad in ((10.0, 0.0, 1.0), (0.0, 1.0, 0.0), (0.0, 1.0, -1.0), (math.nan, 1.0, 1.0)):
        with pytest.raises(ValueError):
            ca.chunks(*bad)


def test_a_curved_path_is_followed_chunk_by_chunk():
    """A caller that re-propagates (the state at each chunk's start) sees a hit the straight path misses."""
    straight = ca.first_hit(lambda t: ((1000.0, 500.0, 0.0), (-10.0, 0.0, 0.0)), 0.0, 200.0, 100.0)
    assert not straight.hit

    def bent(t):
        # the relative path turns toward the origin after t = 50
        return ((1000.0 - 10.0 * t, 500.0, 0.0), (-10.0, 0.0, 0.0)) if t < 50.0 else (
            (500.0 - 10.0 * (t - 50.0), 500.0 - 10.0 * (t - 50.0), 0.0), (-10.0, -10.0, 0.0))

    found = ca.first_hit(bent, 0.0, 200.0, 100.0, cap_s=50.0)
    assert found.hit and found.t_hit_s == pytest.approx(50.0 + (500.0 * math.sqrt(2.0) - 100.0) / (10.0 * math.sqrt(2.0)))


def test_a_window_of_no_length_is_a_hit_only_when_already_overlapping():
    assert ca.first_hit(lambda t: ((10.0, 0.0, 0.0), (0.0, 0.0, 0.0)), 7.0, 7.0, 100.0) == ca.Approach(True, 7.0)
    assert ca.first_hit(lambda t: ((1000.0, 0.0, 0.0), (0.0, 0.0, 0.0)), 7.0, 7.0, 100.0) == ca.Approach(False, 7.0)
    with pytest.raises(ValueError):
        ca.first_hit(lambda t: ((0, 0, 0), (0, 0, 0)), 8.0, 7.0, 1.0)
