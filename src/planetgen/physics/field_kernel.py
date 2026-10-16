# planetgen/physics/field_kernel.py

"""
The compiled gravity kernel (PERF.88)
=====================================

Everything here is compiled by numba (`docs/design/math-and-library-efficiency.md`), so a sector path that took a
millisecond of Python takes about forty microseconds. There is no fallback without numba (OPS.42).

The module is **kept out of the web app's import path**: numba's default cache needs a writable code tree or home
folder, which Apache's web user has neither of. `sector_path` imports this module inside the functions that call it, and
`test_the_web_app_never_imports_numba` holds that.

Two layers:

- **The field**: `accelerations(points, masses, softening)` and `potentials(points, masses, softening)` over the
  `MassTable` of `physics/sector_path.py`, with Plummer softening as an argument. Positions are metres, the masses'
  `mus` are G M in m^3/s^2, the results m/s^2 and m^2/s^2 (the potential per unit mass, negative). `smooth`, when
  given, is a callable `smooth(points) -> (n, 3) array` (acceleration) or `(n,)` array (potential) added to the sum:
  the hook for the galaxy's smooth term (GEN.115). The sums are the same arithmetic in the same order as a plain loop,
  and equal a NumPy sum to better than 1e-12 relative.
- **The sector-path integrator**: `integrate` follows one body through a cell (`inside`), picking the masses that
  matter (`select`) and keeping the spline knots (`simplify`) in the one compiled call. `sector_path.integrate_path`
  calls it; the algorithm is the one the module documents (kick-drift-kick steps, a step a fraction of the shortest
  dynamical or crossing time, an exit found by bisecting the step's own Hermite curve), and numba, without
  `fastmath`, uses the same arithmetic in the same order as the Python it replaced, so paths agree to 4e-16.

A cell is described by numbers, not a callable, because a compiled loop cannot call back into Python: `inside(kind,
params, x, y, z)` with `kind` one of `BOX`, `BALL`, `CORE` or `SHELL` (see `sector_path.Cell`).
"""

import math

import numpy as np
from numba import njit

BOX = 0
"""int: A box, `params = [x0, y0, z0, x1, y1, z1]`, ends included."""

BALL = 1
"""int: A ball, `params = [x, y, z, radius]`, surface excluded."""

CORE = 2
"""int: The core ball of the sector grid, `params = [half_edge_pc, metres_to_pc]` (shell 0)."""

SHELL = 3
"""int: One cell of a sector-grid shell, `params = [r_in, r_out, north, south, wedges, wedge, wedge_step, metres_to_pc]`
(radii in pc, colatitudes in radians; `galaxy.spherical.sector_cell_test`'s arithmetic)."""

PARAM_COUNT = 8
"""int: The length every cell's `params` is padded to."""

TWO_PI = 2.0 * math.pi


@njit(cache=True)
def inside(kind, params, x, y, z):
    """Whether the point (metres) is in the cell `kind` and `params` describe."""
    if kind == BOX:
        return (params[0] <= x <= params[3]) and (params[1] <= y <= params[4]) and (params[2] <= z <= params[5])
    if kind == BALL:
        dx, dy, dz = x - params[0], y - params[1], z - params[2]
        return math.sqrt(dx * dx + dy * dy + dz * dz) < params[3]
    to_pc = params[1] if kind == CORE else params[7]
    x, y, z = x * to_pc, y * to_pc, z * to_pc
    if kind == CORE:
        return math.sqrt(x * x + y * y + z * z) < params[0]
    planar = math.sqrt(x * x + y * y)
    radius = math.sqrt(planar * planar + z * z)
    if radius < params[0] or radius >= params[1]:
        return False
    colatitude = math.atan2(planar, z)
    if colatitude < params[2]:
        return False
    if colatitude >= params[3] and params[3] < math.pi:
        return False
    n = int(params[4])
    if n == 1:
        return True
    wedge = int(params[5])
    step = params[6]
    theta = math.atan2(y, x) % TWO_PI
    found = min(n - 1, int(theta * n / TWO_PI))
    if found > 0 and theta < found * step:
        found -= 1
    elif found < n - 1 and theta >= (found + 1) * step:
        found += 1
    return found == wedge


# --- The field ---------------------------------------------------------------------------------------------------


@njit(cache=True)
def _accelerations(points, positions, mus, softening_sq):
    out = np.zeros((points.shape[0], 3))
    for i in range(points.shape[0]):
        ax = ay = az = 0.0
        px, py, pz = points[i, 0], points[i, 1], points[i, 2]
        for j in range(mus.shape[0]):
            dx, dy, dz = positions[j, 0] - px, positions[j, 1] - py, positions[j, 2] - pz
            r2 = dx * dx + dy * dy + dz * dz + softening_sq
            k = mus[j] / (r2 * math.sqrt(r2))
            ax += k * dx
            ay += k * dy
            az += k * dz
        out[i, 0], out[i, 1], out[i, 2] = ax, ay, az
    return out


@njit(cache=True)
def _potentials(points, positions, mus, softening_sq):
    out = np.zeros(points.shape[0])
    for i in range(points.shape[0]):
        total = 0.0
        px, py, pz = points[i, 0], points[i, 1], points[i, 2]
        for j in range(mus.shape[0]):
            dx, dy, dz = positions[j, 0] - px, positions[j, 1] - py, positions[j, 2] - pz
            total -= mus[j] / math.sqrt(dx * dx + dy * dy + dz * dz + softening_sq)
        out[i] = total
    return out


def _points(points):
    array = np.ascontiguousarray(points, dtype=np.float64)
    if array.ndim == 1:
        array = array.reshape(1, -1)
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError(f"points must be (n, 3) metres, got shape {array.shape}")
    return array


def accelerations(points, masses, softening, smooth=None):
    """
    The pull on each point from the point masses, Plummer-softened.

    Args:
        points: `(n, 3)` positions, m (a single `(3,)` point is one row).
        masses (MassTable): The masses (`sector_path.MassTable`).
        softening (float): Plummer softening length, m.
        smooth (callable, optional): `smooth(points) -> (n, 3)`, the galaxy's smooth term, added.

    Returns:
        numpy.ndarray: `(n, 3)` accelerations, m/s^2.
    """
    array = _points(points)
    out = _accelerations(array, masses.positions, masses.mus, float(softening) ** 2)
    if smooth is not None:
        out = out + np.asarray(smooth(array), dtype=np.float64).reshape(-1, 3)
    return out


def potentials(points, masses, softening, smooth=None):
    """
    The potential per unit mass at each point from the point masses, Plummer-softened (negative, m^2/s^2).

    Args:
        points: `(n, 3)` positions, m.
        masses (MassTable): The masses.
        softening (float): Plummer softening length, m.
        smooth (callable, optional): `smooth(points) -> (n,)`, the galaxy's smooth term, added.

    Returns:
        numpy.ndarray: `(n,)` potentials.
    """
    array = _points(points)
    out = _potentials(array, masses.positions, masses.mus, float(softening) ** 2)
    if smooth is not None:
        out = out + np.asarray(smooth(array), dtype=np.float64).reshape(-1)
    return out


# --- The sector-path integrator ------------------------------------------------------------------------------------


@njit(cache=True)
def select(px, py, pz, vx, vy, vz, positions, mus, exclude, edge_m, softening_m, tolerance_m, limit):
    """
    The indices of the masses that bend a path from `(px, py, pz)` along `(vx, vy, vz)` through the next two sector
    edges by enough to matter, strongest first, at most `limit`, leaving out the indices in `exclude`
    (`sector_path.relevant_masses` has the reasoning).
    """
    n = mus.shape[0]
    speed = math.sqrt(vx * vx + vy * vy + vz * vz)
    if speed == 0.0 or n == 0:
        return np.empty(0, dtype=np.int64)
    length = 2.0 * edge_m
    dx, dy, dz = vx / speed, vy / speed, vz / speed
    displacement = np.empty(n)
    for i in range(n):
        ox, oy, oz = positions[i, 0] - px, positions[i, 1] - py, positions[i, 2] - pz
        along = min(max(ox * dx + oy * dy + oz * dz, 0.0), length)
        rx, ry, rz = ox - along * dx, oy - along * dy, oz - along * dz
        closest = math.sqrt(rx * rx + ry * ry + rz * rz)
        impulse = min(speed, 2.0 * mus[i] / (speed * max(closest, softening_m)))
        displacement[i] = impulse * (length - along) / speed
    wanted = np.nonzero(displacement >= 0.1 * tolerance_m)[0]
    order = np.argsort(-displacement[wanted], kind="mergesort")
    chosen = np.empty(min(limit, wanted.shape[0]), dtype=np.int64)
    count = 0
    for o in range(order.shape[0]):
        index = wanted[order[o]]
        skip = False
        for e in range(exclude.shape[0]):
            if exclude[e] == index:
                skip = True
        if skip:
            continue
        if count >= chosen.shape[0]:
            break
        chosen[count] = index
        count += 1
    return chosen[:count]


@njit(cache=True)
def _pull(px, py, pz, positions, mus, softening_sq, step_fraction):
    """The pull at a point and the two limits it sets on the next step: `(ax, ay, az, dynamical_step, nearest_m)`."""
    ax = ay = az = 0.0
    dynamical = math.inf
    nearest_sq = math.inf
    for j in range(mus.shape[0]):
        dx, dy, dz = positions[j, 0] - px, positions[j, 1] - py, positions[j, 2] - pz
        r2 = dx * dx + dy * dy + dz * dz + softening_sq
        r = math.sqrt(r2)
        k = mus[j] / (r2 * r)
        ax += k * dx
        ay += k * dy
        az += k * dz
        dynamical = min(dynamical, step_fraction * math.sqrt(r * r2 / mus[j]))
        nearest_sq = min(nearest_sq, r2)
    return ax, ay, az, dynamical, math.sqrt(nearest_sq)


@njit(cache=True)
def _hermite(t0, p0, v0, t1, p1, v1, t, out_p, out_v):
    """`sector_path.hermite`'s arithmetic for one spline segment, into `out_p` and `out_v`."""
    h = t1 - t0
    if h <= 0.0:
        for c in range(3):
            out_p[c] = p0[c]
            out_v[c] = v0[c]
        return
    s = min(1.0, max(0.0, (t - t0) / h))
    s2, s3 = s * s, s * s * s
    h00, h10, h01, h11 = 2 * s3 - 3 * s2 + 1, s3 - 2 * s2 + s, -2 * s3 + 3 * s2, s3 - s2
    d00, d10, d01, d11 = 6 * s2 - 6 * s, 3 * s2 - 4 * s + 1, -6 * s2 + 6 * s, 3 * s2 - 2 * s
    for c in range(3):
        out_p[c] = h00 * p0[c] + h10 * h * v0[c] + h01 * p1[c] + h11 * h * v1[c]
        out_v[c] = (d00 * p0[c] + d10 * h * v0[c] + d01 * p1[c] + d11 * h * v1[c]) / h


@njit(cache=True)
def _worst(st, sp, sv, i, j, scratch_p, scratch_v):
    """The sample strictly between `i` and `j` the spline `i`-`j` is furthest from: `(error, index or -1)`."""
    best_error = 0.0
    best_k = -1
    for k in range(i + 1, j):
        _hermite(st[i], sp[i], sv[i], st[j], sp[j], sv[j], st[k], scratch_p, scratch_v)
        dx, dy, dz = scratch_p[0] - sp[k, 0], scratch_p[1] - sp[k, 1], scratch_p[2] - sp[k, 2]
        error = math.sqrt(dx * dx + dy * dy + dz * dz)
        if error > best_error:
            best_error = error
            best_k = k
    return best_error, best_k


@njit(cache=True)
def simplify(st, sp, sv, count, tolerance_m, max_knots):
    """
    The sample indices to keep as knots: the ends, then whichever sample the spline is furthest from, until it is
    within `tolerance_m` or there are `max_knots` (ties go to the segment that starts first, as a heap of
    `(-error, start, end)` would).
    """
    chosen = np.zeros(count, dtype=np.bool_)
    chosen[0] = True
    chosen[count - 1] = True
    kept = 2
    scratch_p = np.empty(3)
    scratch_v = np.empty(3)
    limit = 2 * max_knots + 2
    seg_error = np.empty(limit)
    seg_i = np.empty(limit, dtype=np.int64)
    seg_j = np.empty(limit, dtype=np.int64)
    seg_k = np.empty(limit, dtype=np.int64)
    segments = 0
    if count > 2:
        error, k = _worst(st, sp, sv, 0, count - 1, scratch_p, scratch_v)
        seg_error[0], seg_i[0], seg_j[0], seg_k[0] = error, 0, count - 1, k
        segments = 1
    while segments > 0 and kept < max_knots:
        best = 0
        for s in range(1, segments):
            if seg_error[s] > seg_error[best] or (seg_error[s] == seg_error[best] and (
                    seg_i[s] < seg_i[best] or (seg_i[s] == seg_i[best] and seg_j[s] < seg_j[best]))):
                best = s
        error, i, j, k = seg_error[best], seg_i[best], seg_j[best], seg_k[best]
        seg_error[best], seg_i[best], seg_j[best], seg_k[best] = (
            seg_error[segments - 1], seg_i[segments - 1], seg_j[segments - 1], seg_k[segments - 1])
        segments -= 1
        if error <= tolerance_m or k < 0:
            break
        chosen[k] = True
        kept += 1
        for lo, hi in ((i, k), (k, j)):
            if hi - lo > 1:
                e, mid = _worst(st, sp, sv, lo, hi, scratch_p, scratch_v)
                seg_error[segments], seg_i[segments], seg_j[segments], seg_k[segments] = e, lo, hi, mid
                segments += 1
    return np.nonzero(chosen)[0]


@njit(cache=True)
def integrate(position, velocity, positions, mus, exclude, kind, params, edge_m, softening_m, tolerance_m, max_knots,
              max_steps, limit, step_fraction, crossings_allowed, max_turn):
    """
    One body's path through a cell, `sector_path.integrate_path`'s algorithm in one compiled call.

    Returns:
        tuple: `(times (k,), positions (k, 3), velocities (k, 3), exited, steps)`, the kept knots.
    """
    px, py, pz = position[0], position[1], position[2]
    vx, vy, vz = velocity[0], velocity[1], velocity[2]
    softening_sq = softening_m * softening_m
    chosen = select(px, py, pz, vx, vy, vz, positions, mus, exclude, edge_m, softening_m, tolerance_m, limit)
    active_positions = np.empty((chosen.shape[0], 3))
    active_mus = np.empty(chosen.shape[0])
    for c in range(chosen.shape[0]):
        active_positions[c] = positions[chosen[c]]
        active_mus[c] = mus[chosen[c]]
    speed0 = math.sqrt(vx * vx + vy * vy + vz * vz)
    max_step = edge_m / (speed0 * 20.0)
    max_time = crossings_allowed * edge_m / speed0

    st = np.empty(max_steps + 2)
    sp = np.empty((max_steps + 2, 3))
    sv = np.empty((max_steps + 2, 3))
    count = 1
    st[0] = 0.0
    sp[0, 0], sp[0, 1], sp[0, 2] = px, py, pz
    sv[0, 0], sv[0, 1], sv[0, 2] = vx, vy, vz
    t = 0.0
    ax, ay, az, dynamical, nearest = _pull(px, py, pz, active_positions, active_mus, softening_sq, step_fraction)
    exited = False
    turned = 0.0
    steps = 0
    for _ in range(max_steps):
        steps += 1
        speed_prev = math.sqrt(vx * vx + vy * vy + vz * vz)
        h = min(max_step, dynamical, step_fraction * nearest / speed_prev)
        h = max(h, max_step * 1.0e-6)
        hx, hy, hz = vx + 0.5 * h * ax, vy + 0.5 * h * ay, vz + 0.5 * h * az
        nx, ny, nz = px + h * hx, py + h * hy, pz + h * hz
        bx, by, bz, dynamical_next, nearest_next = _pull(nx, ny, nz, active_positions, active_mus, softening_sq,
                                                         step_fraction)
        wx, wy, wz = hx + 0.5 * h * bx, hy + 0.5 * h * by, hz + 0.5 * h * bz
        if not inside(kind, params, nx, ny, nz):
            # Where the step leaves the cell: bisect the step's own Hermite curve.
            a_p = np.array([px, py, pz])
            a_v = np.array([vx, vy, vz])
            b_p = np.array([nx, ny, nz])
            b_v = np.array([wx, wy, wz])
            out_p = np.empty(3)
            out_v = np.empty(3)
            low, high = t, t + h
            for _bisect in range(40):
                middle = 0.5 * (low + high)
                _hermite(t, a_p, a_v, t + h, b_p, b_v, middle, out_p, out_v)
                if inside(kind, params, out_p[0], out_p[1], out_p[2]):
                    low = middle
                else:
                    high = middle
            _hermite(t, a_p, a_v, t + h, b_p, b_v, high, out_p, out_v)
            st[count] = high
            for c in range(3):
                sp[count, c] = out_p[c]
                sv[count, c] = out_v[c]
            count += 1
            exited = True
            break
        speed_next = math.sqrt(wx * wx + wy * wy + wz * wz)
        if speed_next > 0.0:
            cosine = (vx * wx + vy * wy + vz * wz) / (speed_prev * speed_next)
            turned += math.acos(max(-1.0, min(1.0, cosine)))
        t = t + h
        px, py, pz = nx, ny, nz
        vx, vy, vz = wx, wy, wz
        ax, ay, az, dynamical, nearest = bx, by, bz, dynamical_next, nearest_next
        st[count] = t
        sp[count, 0], sp[count, 1], sp[count, 2] = px, py, pz
        sv[count, 0], sv[count, 1], sv[count, 2] = vx, vy, vz
        count += 1
        if t > max_time or turned > max_turn:
            break
    kept = simplify(st, sp, sv, count, tolerance_m, max_knots)
    times = np.empty(kept.shape[0])
    out_positions = np.empty((kept.shape[0], 3))
    out_velocities = np.empty((kept.shape[0], 3))
    for c in range(kept.shape[0]):
        times[c] = st[kept[c]]
        out_positions[c] = sp[kept[c]]
        out_velocities[c] = sv[kept[c]]
    return times, out_positions, out_velocities, exited, steps
