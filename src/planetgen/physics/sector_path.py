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

The integration itself runs in `physics/field_kernel.py`, compiled by numba (PERF.88); this module holds what
the kernel needs described in plain terms (`Cell`, `MassTable`, `PointMass`) and builds the `SectorPath`s from what
it returns. It imports the kernel inside the functions that call it, so the web app, which reads paths but never
integrates them, never imports numba.

Point masses are softened like the galaxy's own (`Plummer`, 1 pc by
default, see `docs/design/orbital-updates.md` section 10.1), so the pull
near a black hole stays finite. A single dominant flyby could use the
hyperbolic deflection formula instead; the integration covers it without a
second path.
"""

import bisect
import math
from dataclasses import dataclass, field

import numpy as np

from planetgen.galaxy import spherical
from planetgen.physics import constants

DEFAULT_SOFTENING_M = constants.PARSEC_M
"""float: Plummer softening length of a point mass, 1 pc."""

DEFAULT_TOLERANCE_FRACTION = 2.0e-3
"""float: How far, as a fraction of the sector edge, the spline may stray
from the integrated path before it gets another knot."""

DEFAULT_MAX_KNOTS = 48
"""int: The most knots one sector's path keeps."""

DEFAULT_MAX_STEPS = 2000
"""int: The most integration steps before the path is cut off. A step is
a fraction of the shortest dynamical or crossing time, so a path that has
not left by now is held by a mass (see `MAX_TURN_RAD`) or is plunging
through a very dense stretch; it is the hard ceiling on a body's cost."""

MAX_TURN_RAD = 2.0 * math.pi
"""float: A path is cut off once its direction of travel has turned through
this much in all, one full revolution: a body held by a mass has come all
the way round, and the spline (`DEFAULT_MAX_KNOTS` knots) cannot follow it
further anyway. A flyby turns by less than half of it. Without this cap a
body in the well of a heavy black hole ran every one of the thousands of
steps for each of the thousands of bodies in a core sector (GEN.232)."""

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
        self._by_key = {}
        for index, key in enumerate(self.keys):
            if key is not None:
                self._by_key.setdefault(key, []).append(index)

    def __len__(self):
        return len(self.masses)

    def indices_of(self, key):
        """The indices of the masses named `key` (none for `None`), as the kernel's int64 array."""
        return np.array(self._by_key.get(key, ()) if key is not None else (), dtype=np.int64)


@dataclass(frozen=True)
class PathKnot:
    """A spline knot: seconds after the path starts, position (m) and velocity (m/s)."""

    t_s: float
    position: tuple
    velocity: tuple


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

    def __init__(self, knots, exited, steps=0):
        if not knots:
            raise ValueError("a path needs at least one knot")
        self.knots = list(knots)
        self.exited = bool(exited)
        self.steps = int(steps)
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


def _kernel():
    """The compiled kernel (`field_kernel`), imported when first needed so the web app never loads numba."""
    from planetgen.physics import field_kernel

    return field_kernel


@dataclass(frozen=True)
class Cell:
    """
    The region a path is followed through, described by numbers so the compiled integrator can test it
    (`field_kernel.inside`): build one with `box_cell`, `ball_cell` or `sector_inside`. Calling it with a position
    (metres, galactic axes) says whether the position is in the region.
    """

    kind: int
    params: tuple

    def contains(self, position):
        return bool(_kernel().inside(self.kind, self.padded(), float(position[0]), float(position[1]),
                                     float(position[2])))

    __call__ = contains

    def padded(self):
        """`params` as the float array the kernel takes."""
        kernel = _kernel()
        params = np.zeros(kernel.PARAM_COUNT)
        params[:len(self.params)] = self.params
        return params


def box_cell(low, high):
    """A box from `low` to `high` (each `(x, y, z)`, metres), the faces included."""
    return Cell(0, tuple(map(float, low)) + tuple(map(float, high)))


def ball_cell(center, radius):
    """A ball round `center` (metres) of `radius` metres, the surface left out."""
    return Cell(1, tuple(map(float, center)) + (float(radius),))


def sector_inside(address, edge_pc):
    """
    The `Cell` of the sector `address` of the grid with edge `edge_pc` (`galaxy.spherical`): true exactly where
    `spherical.sector_address_at` names that address, for a position in metres.
    """
    address = tuple(address)
    shell, band, wedge = address
    to_pc = 1.0 / constants.PARSEC_M
    spherical.sector_cell_test(address, edge_pc)   # raises for an address that is not a cell
    if shell == 0:
        return Cell(2, (edge_pc / 2, to_pc))
    n = spherical.wedge_count(shell, band)
    inner, outer = spherical.shell_bounds_pc(shell, edge_pc)
    north, south = spherical.band_colatitude_bounds(shell, band)
    return Cell(3, (inner, outer, north, south, float(n), float(wedge), spherical.TWO_PI / n, to_pc))


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
    if tolerance_m is None:
        tolerance_m = DEFAULT_TOLERANCE_FRACTION * edge_m
    chosen = _kernel().select(float(position[0]), float(position[1]), float(position[2]),
                              float(velocity[0]), float(velocity[1]), float(velocity[2]),
                              table.positions, table.mus, table.indices_of(exclude_key), float(edge_m),
                              float(softening_m), float(tolerance_m), int(limit))
    return [table.masses[index] for index in chosen]


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
        inside (Cell): The region the body is followed through
            (`sector_inside`, `box_cell`, `ball_cell`).
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
        the step limit ran out); `steps` is how many integration steps it took.

    Raises:
        ValueError: For fewer than 2 knots, or a start outside the sector.
    """
    if max_knots < 2:
        raise ValueError(f"a path needs at least 2 knots, got {max_knots!r}")
    position, velocity = tuple(map(float, position)), tuple(map(float, velocity))
    speed = math.sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2])
    start = PathKnot(0.0, position, velocity)
    if not inside.contains(position):
        raise ValueError("the path starts outside the sector")
    if speed == 0.0:
        return SectorPath([start], exited=False)
    table = masses if isinstance(masses, MassTable) else MassTable(masses)
    if tolerance_m is None:
        tolerance_m = DEFAULT_TOLERANCE_FRACTION * edge_m
    times, positions, velocities, exited, steps = _kernel().integrate(
        np.array(position), np.array(velocity), table.positions, table.mus, table.indices_of(exclude_key),
        inside.kind, inside.padded(), float(edge_m), float(softening_m), float(tolerance_m), int(max_knots),
        int(max_steps), DEFAULT_RELEVANT_MASSES, STEP_FRACTION, CROSSINGS_ALLOWED, MAX_TURN_RAD)
    knots = [PathKnot(t, tuple(p), tuple(v)) for t, p, v in zip(times.tolist(), positions.tolist(), velocities.tolist())]
    return SectorPath(knots, bool(exited), steps=int(steps))
