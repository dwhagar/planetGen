# planetgen/physics/closest_approach.py

"""
The pair closest-approach library (GEN.334; docs/design/orbits/collisions-and-mergers.md section 5)
============================================================================================

One question, asked of any two paths: do they come within `R` of each other in this window? For `R = R1 + R2` it is a
collision test, for a system's radius an encounter test. Nothing here decides what a collision does (GEN.345 and GEN.110)
or finds the candidates (the broad phase of section 5.4); it answers the question for one pair.

- `sweep(dp, dv, R, dt)` is the closest-approach form of section 5.2: the time of closest approach and the miss distance
  as a vector first, then the entry time, so nothing large is taken from something nearly equal (the swept-sphere form
  in the source document loses the miss distance to rounding once the start is more than a few thousand times `R`
  away, section 2.1). A pair already overlapping at the start is a hit at 0: the caller leaves out pairs made in the
  same run.
- `Segment` and `linear_state_at` are the pair-window rule of section 5.3: each body stores a position, a velocity and an
  epoch valid until its next update, and an update of X tests it against Y's stored segment over `[max(e_X_prev, e_Y),
  e_X_now]`. Positions are differenced from sector-relative values: a segment's `origin` is its sector's offset and its
  `position` is relative to it, so two bodies in one sector never subtract two large numbers.
- `window_cap` is the longest window one straight relative path may cover: the tidal curvature of the galaxy bends the
  relative path by about `0.03 Omega^2 v_rel t^3`, held under 1 percent of `R` (`CURVATURE_FRACTION`), and a pair must
  not move more than one sector in a window. `chunks` cuts a longer span at the cap (and at most
  `CURVATURE_CHUNK_S`, 0.1 Myr, the chunk length the galactic curvature of a long step is carried at).
- `first_hit` walks the chunks of a window, asking `state_at(t)` for the relative state at each chunk's start, so a caller
  that re-propagates a curved path gives its own states and a straight pair uses `linear_state_at`.

All lengths are metres, times seconds, speeds m/s. Non-finite input, a negative radius or a negative window raises
`ValueError`.
"""

import dataclasses
import math

COEFFICIENT = 0.03
"""float: The curvature fit of section 5.3: the straight relative path is off by about `0.03 Omega^2 v_rel t^3`."""
CURVATURE_FRACTION = 0.01
"""float: The share of `R` the straight path may be off by at the end of a window."""
CURVATURE_CHUNK_S = 0.1e6 * 365.25 * 86400.0
"""float: The longest chunk a long step is cut into, 0.1 million years (section 5.3's chunk for galactic curvature)."""


def _finite(*values):
    for value in values:
        if not math.isfinite(value):
            raise ValueError(f"closest approach needs finite numbers, got {value!r}")


def sweep(dp, dv, R, dt):
    """
    Whether a pair comes within `R` in the window `[0, dt]`, and when.

    Args:
        dp (sequence): The relative position `(x, y, z)` at the window start, metres, differenced from sector-relative
            values (never from absolute galactic coordinates).
        dv (sequence): The relative velocity, m/s.
        R (float): The distance that counts as a hit (`R1 + R2` for a collision), metres, `>= 0`.
        dt (float): The window length, seconds, `>= 0`.

    Returns:
        tuple: `(hit, t_hit)`: `t_hit` seconds after the start (0 for a pair already within `R`), or `dt` when there is
            no hit. A pair that only touches (a miss distance exactly `R`) is a hit.

    Raises:
        ValueError: For a non-finite input, a negative `R` or a negative `dt`.
    """
    _finite(*dp, *dv, R, dt)
    if R < 0.0 or dt < 0.0:
        raise ValueError(f"R and dt must be >= 0, got {R} and {dt}")
    d2 = dp[0] ** 2 + dp[1] ** 2 + dp[2] ** 2
    r2 = R * R
    if d2 <= r2:
        return True, 0.0
    a = dv[0] ** 2 + dv[1] ** 2 + dv[2] ** 2
    if a == 0.0:
        return False, dt
    t_min = -(dp[0] * dv[0] + dp[1] * dv[1] + dp[2] * dv[2]) / a
    if t_min < 0.0:                         # separating
        return False, dt
    mx, my, mz = (dp[0] + dv[0] * t_min, dp[1] + dv[1] * t_min, dp[2] + dv[2] * t_min)
    m2 = mx * mx + my * my + mz * mz        # the miss distance squared, from the vector
    if m2 > r2:
        return False, dt
    t_hit = t_min - math.sqrt((r2 - m2) / a)
    if t_hit <= dt:
        return True, max(t_hit, 0.0)
    return False, dt


@dataclasses.dataclass(frozen=True)
class Segment:
    """
    A body's stored path: at `epoch_s` it is at `origin_m + position_m` (`origin_m` its sector's offset, `position_m`
    relative to it, metres) moving at `velocity_ms` until its next update.
    """

    position_m: tuple
    velocity_ms: tuple
    epoch_s: float
    origin_m: tuple = (0.0, 0.0, 0.0)

    def local_at(self, t_s):
        """The position at `t_s` relative to this segment's own origin, the straight path extended."""
        elapsed = t_s - self.epoch_s
        return tuple(p + v * elapsed for p, v in zip(self.position_m, self.velocity_ms))


def linear_state_at(x, y):
    """
    `state_at(t)` for `first_hit`: the relative position `x - y` and velocity at `t` of two segments' straight paths,
    the two origins differenced first (exactly, when both bodies share a sector) and then the sector-relative
    positions.
    """
    shift = tuple(a - b for a, b in zip(x.origin_m, y.origin_m))
    dv = tuple(a - b for a, b in zip(x.velocity_ms, y.velocity_ms))

    def state_at(t_s):
        px, py = x.local_at(t_s), y.local_at(t_s)
        return tuple(s + (a - b) for s, a, b in zip(shift, px, py)), dv

    return state_at


def window_cap(R, v_rel, omega_s, edge_m, fraction=CURVATURE_FRACTION):
    """
    The longest window, seconds, one straight relative path may cover: the smaller of the time for the galactic tide to
    bend it by `fraction` of `R`, `(fraction R / (0.03 Omega^2 v_rel))^(1/3)`, and the time to cross one sector,
    `edge_m / v_rel`. Infinite for a pair at rest relative to each other (nothing moves, nothing bends).

    Args:
        R (float): The hit distance, metres, `> 0`.
        v_rel (float): The pair's relative speed, m/s, `>= 0`.
        omega_s (float): The galaxy's angular speed at the pair, rad/s, `>= 0` (0 for no tide).
        edge_m (float): The sector edge, metres, `> 0`.
    """
    _finite(R, v_rel, omega_s, edge_m)
    if R <= 0.0 or edge_m <= 0.0 or v_rel < 0.0 or omega_s < 0.0:
        raise ValueError(f"window_cap needs R > 0, edge > 0, v_rel >= 0 and omega >= 0, got {R}, {edge_m}, {v_rel}, {omega_s}")
    if v_rel == 0.0:
        return math.inf
    reach = edge_m / v_rel
    if omega_s == 0.0:
        return reach
    return min((fraction * R / (COEFFICIENT * omega_s ** 2 * v_rel)) ** (1.0 / 3.0), reach)


def chunks(t0, t1, cap_s, chunk_s=CURVATURE_CHUNK_S):
    """
    `[t0, t1]` cut into equal windows no longer than `min(cap_s, chunk_s)`: a list of `(start, end)`, the last ending at
    `t1` exactly. Empty for a window of no length (`t1 == t0`).

    Raises:
        ValueError: For `t1 < t0`, a non-finite bound or a cap that is not positive.
    """
    _finite(t0, t1)
    if t1 < t0:
        raise ValueError(f"the window ends before it starts: {t0} to {t1}")
    limit = min(cap_s, chunk_s)
    if not limit > 0.0:
        raise ValueError(f"the chunk length must be > 0, got {limit}")
    if t1 == t0:
        return []
    count = max(1, math.ceil((t1 - t0) / limit))
    edges = [t0 + (t1 - t0) * i / count for i in range(count)] + [t1]
    return list(zip(edges[:-1], edges[1:]))


@dataclasses.dataclass(frozen=True)
class Approach:
    """The answer for one pair and window: `hit`, and `t_hit_s` (on the caller's clock; the window's end when no hit)."""

    hit: bool
    t_hit_s: float


def first_hit(state_at, t0, t1, R, cap_s=math.inf, chunk_s=CURVATURE_CHUNK_S):
    """
    Whether a pair comes within `R` in `[t0, t1]`, and the first time it does.

    Args:
        state_at (callable): `state_at(t)` gives `(dp, dv)`, the pair's relative position and velocity at time `t`; the
            window's chunks start at `t0` and at each cut. `linear_state_at` gives straight paths; a caller that
            re-propagates the paths over a long span gives its own.
        t0, t1 (float): The window, `t0 <= t1`: `max(e_X_prev, e_Y)` to `e_X_now` for X's update (section 5.3).
        R (float): The hit distance, metres.
        cap_s (float): The longest straight window (`window_cap`); the window is cut at it and at `chunk_s`.

    Returns:
        Approach: `hit` and the first hit's time; for a window of no length, a hit only if the pair already overlaps.
    """
    _finite(t0, t1, R)
    if t1 < t0:
        raise ValueError(f"the window ends before it starts: {t0} to {t1}")
    if t1 == t0:
        dp, dv = state_at(t0)
        hit, _ = sweep(dp, dv, R, 0.0)
        return Approach(hit, t0)
    for start, end in chunks(t0, t1, cap_s, chunk_s):
        dp, dv = state_at(start)
        hit, t_hit = sweep(dp, dv, R, end - start)
        if hit:
            return Approach(True, start + t_hit)
    return Approach(False, t1)


def pair_window(y, previous_epoch_x_s, now_x_s):
    """
    The window `[max(e_X_prev, e_Y), e_X_now]` of section 5.3, for body X updating from `previous_epoch_x_s` to
    `now_x_s` against Y's stored segment `y`: `(start, end)`, or `None` when Y's segment begins after the update's end
    (Y's path then has nothing in X's window; the pair is tested when Y updates).
    """
    _finite(previous_epoch_x_s, now_x_s, y.epoch_s)
    if now_x_s < previous_epoch_x_s:
        raise ValueError(f"the update ends before it starts: {previous_epoch_x_s} to {now_x_s}")
    if y.epoch_s > now_x_s:
        return None
    return max(previous_epoch_x_s, y.epoch_s), now_x_s
