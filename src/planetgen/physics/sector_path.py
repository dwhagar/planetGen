# planetgen/physics/sector_path.py

"""
The path of a body through a sector (GEN.123)
=============================================

A body with no closed orbit (a star, a rogue planet, a hyperbolic or
parabolic comet, an interstellar object) crosses a sector on a path that
is not a straight line when it passes a heavy mass. `integrate_path`
follows a test particle from where it enters the sector, with the velocity
it has there, against the sector's point masses (and the neighbours' that
matter), until it leaves the sector's cell. The masses stay where they are
for the crossing and the particle does not pull back.

The path is kept as **cubic Hermite spline knots** (time, position,
velocity). Two knots when it is nearly straight; more where it bends,
placed where the spline is furthest from the integration, up to a hard
cap. The exit knot is the next sector's entry, so paths chain across
sectors: pass `path.end` of one as the start of the next.

Everything here is SI (metres, seconds, kilograms) and in the galactic
axes. Positions are tuples, not `SpatialPosition3D`s, so the module needs
no database or sector.

Point masses are softened like the galaxy's own (`Plummer`, 1 pc by
default, see `docs/design/orbital-updates.md` section 10.1), so the pull
near a black hole stays finite. A single dominant flyby could use the
hyperbolic deflection formula instead; the integration covers it without a
second path.
"""

import bisect
import heapq
import math
from dataclasses import dataclass, field

import numpy as np

from planetgen.galaxy import geometry
from planetgen.physics import constants

DEFAULT_SOFTENING_M = constants.PARSEC_M
"""float: Plummer softening length of a point mass, 1 pc."""

DEFAULT_TOLERANCE_FRACTION = 2.0e-3
"""float: How far, as a fraction of the sector edge, the spline may stray
from the integrated path before it gets another knot."""

DEFAULT_MAX_KNOTS = 48
"""int: The most knots one sector's path keeps."""

DEFAULT_MAX_STEPS = 20000
"""int: The most integration steps before the path is cut off."""

DEFAULT_RELEVANT_MASSES = 32
"""int: The most point masses one path is integrated against."""

STEP_FRACTION = 0.02
"""float: A step is this fraction of the shortest dynamical or crossing time
of any mass that matters."""

CROSSINGS_ALLOWED = 100.0
"""float: A path is cut off after this many straight-line crossings of the
sector's edge (a body held by a mass never leaves)."""


@dataclass(frozen=True)
class PointMass:
    """A mass the path bends round: its position (m, galactic axes) and mass (kg)."""

    position: tuple
    mass_kg: float
    key: object = field(default=None, compare=False)
    """Names the mass (say `("star_systems", 12)`), so the body it is does
    not feel its own pull: pass the body's key as `exclude_key`."""

    @property
    def mu(self):
        """The gravitational parameter G M, m^3/s^2."""
        return constants.G * self.mass_kg


class MassTable:
    """
    The point masses of a sector and its neighbours, kept as arrays so that
    picking the few that matter for each of thousands of bodies stays cheap
    (`relevant_masses`). Build it once and reuse it for every body.
    """

    def __init__(self, masses=()):
        self.masses = list(masses)
        self.positions = np.array([m.position for m in self.masses], dtype=float).reshape(-1, 3)
        self.mus = np.array([m.mu for m in self.masses], dtype=float)
        self.keys = [m.key for m in self.masses]

    def __len__(self):
        return len(self.masses)


@dataclass(frozen=True)
class PathKnot:
    """A spline knot: seconds after the path starts, position (m) and velocity (m/s)."""

    t_s: float
    position: tuple
    velocity: tuple


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a):
    return math.sqrt(_dot(a, a))


def hermite(a, b, t_s):
    """
    The cubic Hermite spline between knots `a` and `b` at time `t_s`
    (clamped to the knots): `(position, velocity)`.
    """
    h = b.t_s - a.t_s
    if h <= 0.0:
        return a.position, a.velocity
    s = min(1.0, max(0.0, (t_s - a.t_s) / h))
    s2, s3 = s * s, s * s * s
    h00, h10, h01, h11 = 2 * s3 - 3 * s2 + 1, s3 - 2 * s2 + s, -2 * s3 + 3 * s2, s3 - s2
    d00, d10, d01, d11 = 6 * s2 - 6 * s, 3 * s2 - 4 * s + 1, -6 * s2 + 6 * s, 3 * s2 - 2 * s
    position = tuple(h00 * pa + h10 * h * va + h01 * pb + h11 * h * vb
                     for pa, va, pb, vb in zip(a.position, a.velocity, b.position, b.velocity))
    velocity = tuple((d00 * pa + d10 * h * va + d01 * pb + d11 * h * vb) / h
                     for pa, va, pb, vb in zip(a.position, a.velocity, b.position, b.velocity))
    return position, velocity


class SectorPath:
    """
    A body's path through one sector: `knots` (a Hermite spline in time),
    and whether it `exited` the cell (rather than being cut off).
    """

    def __init__(self, knots, exited):
        if not knots:
            raise ValueError("a path needs at least one knot")
        self.knots = list(knots)
        self.exited = bool(exited)
        self._times = [knot.t_s for knot in self.knots]

    @property
    def start(self):
        """The `PathKnot` where the path enters."""
        return self.knots[0]

    @property
    def end(self):
        """The `PathKnot` where it leaves (or is cut off): the next
        sector's start, once its time is set to 0 by `restart`."""
        return self.knots[-1]

    @property
    def duration_s(self):
        return self.knots[-1].t_s - self.knots[0].t_s

    def state_at(self, t_s):
        """`(position, velocity)` at `t_s` seconds after the path starts (clamped)."""
        if len(self.knots) == 1:
            return self.knots[0].position, self.knots[0].velocity
        i = min(max(bisect.bisect_right(self._times, t_s) - 1, 0), len(self.knots) - 2)
        return hermite(self.knots[i], self.knots[i + 1], t_s)

    def position_at(self, t_s):
        return self.state_at(t_s)[0]

    def velocity_at(self, t_s):
        return self.state_at(t_s)[1]

    def points(self, count=64):
        """`count` positions evenly spaced in time along the spline, to draw."""
        if count < 2:
            raise ValueError(f"a path needs at least 2 points, got {count!r}")
        t0, span = self.knots[0].t_s, self.duration_s
        return [self.position_at(t0 + span * k / (count - 1)) for k in range(count)]

    def restart(self):
        """The end of this path as the start of the next sector's: the same
        position and velocity at time 0."""
        return PathKnot(0.0, self.end.position, self.end.velocity)


def sector_inside(address, edge_pc):
    """
    A predicate `inside(position_m)` that is true while the point (metres,
    galactic axes) is in the sector cell `address` of the grid with edge
    `edge_pc` (`galaxy.geometry`).
    """
    to_pc = 1.0 / constants.PARSEC_M

    def inside(position_m):
        return geometry.sector_address_at(tuple(c * to_pc for c in position_m), edge_pc) == tuple(address)

    return inside


def relevant_masses(position, velocity, masses, edge_m, softening_m=DEFAULT_SOFTENING_M, tolerance_m=None,
                    limit=DEFAULT_RELEVANT_MASSES, exclude_key=None):
    """
    The point masses that bend a path from `position` along `velocity`
    through the next two sector edges by enough to matter, strongest first
    and at most `limit` of them. A mass bends a passing body by about
    `2 mu / (b v)` in velocity, which over the rest of the crossing moves it
    by `2 mu t / (b v)`; masses that move it by far less than the spline's
    tolerance are dropped, and so is the one named `exclude_key` (the
    body's own).

    Args:
        masses (MassTable or iterable of PointMass): A `MassTable` is
            cheaper when many bodies share the masses.
    """
    table = masses if isinstance(masses, MassTable) else MassTable(masses)
    speed = _norm(velocity)
    if speed == 0.0 or not len(table):
        return []
    if tolerance_m is None:
        tolerance_m = DEFAULT_TOLERANCE_FRACTION * edge_m
    length = 2.0 * edge_m
    direction = np.array(velocity, dtype=float) / speed
    offsets = table.positions - np.array(position, dtype=float)
    along = np.clip(offsets @ direction, 0.0, length)
    closest = np.linalg.norm(offsets - np.outer(along, direction), axis=1)
    impulse = np.minimum(speed, 2.0 * table.mus / (speed * np.maximum(closest, softening_m)))
    displacement = impulse * (length - along) / speed
    wanted = np.nonzero(displacement >= 0.1 * tolerance_m)[0]
    wanted = wanted[np.argsort(-displacement[wanted], kind="stable")]
    chosen = []
    for index in wanted:
        if exclude_key is not None and table.keys[index] == exclude_key:
            continue
        chosen.append(table.masses[index])
        if len(chosen) >= limit:
            break
    return chosen


def _acceleration(position, masses, softening_sq):
    ax = ay = az = 0.0
    for mass in masses:
        dx, dy, dz = mass.position[0] - position[0], mass.position[1] - position[1], mass.position[2] - position[2]
        r2 = dx * dx + dy * dy + dz * dz + softening_sq
        k = mass.mu / (r2 * math.sqrt(r2))
        ax += k * dx
        ay += k * dy
        az += k * dz
    return (ax, ay, az)


def _step_size(position, speed, masses, softening_sq, max_step):
    step = max_step
    for mass in masses:
        dx, dy, dz = mass.position[0] - position[0], mass.position[1] - position[1], mass.position[2] - position[2]
        r2 = dx * dx + dy * dy + dz * dz + softening_sq
        r = math.sqrt(r2)
        step = min(step, STEP_FRACTION * math.sqrt(r * r2 / mass.mu), STEP_FRACTION * r / speed)
    return max(step, max_step * 1.0e-6)


def _exit_state(inside, t0, p0, v0, t1, p1, v1):
    """Where the step from `(p0, v0)` to `(p1, v1)` leaves the cell, found by
    bisecting the step's own Hermite curve."""
    a, b = PathKnot(t0, p0, v0), PathKnot(t1, p1, v1)
    low, high = t0, t1
    for _ in range(40):
        middle = 0.5 * (low + high)
        if inside(hermite(a, b, middle)[0]):
            low = middle
        else:
            high = middle
    position, velocity = hermite(a, b, high)
    return PathKnot(high, position, velocity)


def _simplify(samples, tolerance_m, max_knots):
    """Picks the knots: the ends, then whichever integrated sample the
    spline is furthest from, until it is within `tolerance_m` or `max_knots`."""
    chosen = {0, len(samples) - 1}

    def worst(i, j):
        a, b = samples[i], samples[j]
        best_error, best_k = 0.0, None
        for k in range(i + 1, j):
            position, _velocity = hermite(a, b, samples[k].t_s)
            error = _norm(_sub(position, samples[k].position))
            if error > best_error:
                best_error, best_k = error, k
        return best_error, best_k

    heap = []
    if len(samples) > 2:
        error, k = worst(0, len(samples) - 1)
        heap.append((-error, 0, len(samples) - 1, k))
    while heap and len(chosen) < max_knots:
        negative_error, i, j, k = heapq.heappop(heap)
        if -negative_error <= tolerance_m or k is None:
            break
        chosen.add(k)
        for lo, hi in ((i, k), (k, j)):
            if hi - lo > 1:
                error, mid = worst(lo, hi)
                heapq.heappush(heap, (-error, lo, hi, mid))
    return [samples[index] for index in sorted(chosen)]


def integrate_path(position, velocity, masses, inside, edge_m, softening_m=DEFAULT_SOFTENING_M, tolerance_m=None,
                   max_knots=DEFAULT_MAX_KNOTS, max_steps=DEFAULT_MAX_STEPS, exclude_key=None):
    """
    The path of a body that enters a sector at `position` (m) with `velocity`
    (m/s), until it leaves it.

    Args:
        position, velocity (tuple): The entry state, galactic axes.
        masses (MassTable or iterable of PointMass): The sector's point
            masses and any neighbours'; those that cannot bend the path by
            the tolerance are skipped (`relevant_masses`).
        inside (callable): `inside(position_m)` is true while the body is in
            the sector (`sector_inside`).
        edge_m (float): The sector's edge, m, which sets the step limit and
            the spline's tolerance.
        softening_m (float): Plummer softening of every mass.
        tolerance_m (float or None): How far the spline may stray from the
            integration; `None` is `DEFAULT_TOLERANCE_FRACTION` of the edge.
        max_knots (int): The most knots kept (at least 2).
        max_steps (int): The most steps before the path is cut off.
        exclude_key: The key of the body's own `PointMass`, if it is among
            `masses`.

    Returns:
        SectorPath: Two knots when the path is nearly straight. `exited` is
        false when the body never left (it is at rest, or held by a mass, or
        the step limit ran out).

    Raises:
        ValueError: For fewer than 2 knots, or a start outside the sector.
    """
    if max_knots < 2:
        raise ValueError(f"a path needs at least 2 knots, got {max_knots!r}")
    position, velocity = tuple(map(float, position)), tuple(map(float, velocity))
    speed = _norm(velocity)
    start = PathKnot(0.0, position, velocity)
    if not inside(position):
        raise ValueError("the path starts outside the sector")
    if speed == 0.0:
        return SectorPath([start], exited=False)
    if tolerance_m is None:
        tolerance_m = DEFAULT_TOLERANCE_FRACTION * edge_m
    softening_sq = softening_m * softening_m
    active = relevant_masses(position, velocity, masses, edge_m, softening_m, tolerance_m, exclude_key=exclude_key)
    max_step = edge_m / (speed * 20.0)
    max_time = CROSSINGS_ALLOWED * edge_m / speed

    samples = [start]
    t, p, v = 0.0, position, velocity
    a = _acceleration(p, active, softening_sq)
    exited = False
    for _ in range(max_steps):
        h = _step_size(p, _norm(v), active, softening_sq, max_step)
        v_half = tuple(vi + 0.5 * h * ai for vi, ai in zip(v, a))
        p_next = tuple(pi + h * vi for pi, vi in zip(p, v_half))
        a_next = _acceleration(p_next, active, softening_sq)
        v_next = tuple(vi + 0.5 * h * ai for vi, ai in zip(v_half, a_next))
        if not inside(p_next):
            samples.append(_exit_state(inside, t, p, v, t + h, p_next, v_next))
            exited = True
            break
        t, p, v, a = t + h, p_next, v_next, a_next
        samples.append(PathKnot(t, p, v))
        if t > max_time:
            break
    return SectorPath(_simplify(samples, tolerance_m, max_knots), exited)
