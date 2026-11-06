# planetgen/galaxy/gravity_field.py

"""
The gravity map's field as the map shows it: total or local, one number per zone, and its saddle points (API.24; docs/design/orbits/gravity-map.md)
=====================================================================================================================================================

`db.gravity_grid` stores the **local** field of a sector (the generated stars less the smooth galaxy's own stars, MAP.383) as `(zones^3, 10)` numbers in
`gravity_kernel`'s units. The map shows the **total** field (Boss, 2026-11-05: the generated stars plus the smooth galaxy added), so this module adds the smooth
galaxy's pull, potential and tidal tensor at the zone centers (`galaxy.potential`, closed forms and a table, no compiled code, so the web app can call it),
and turns the ten numbers into the three scalars of the map's modes:

- **pull**: the size of the acceleration, m/s^2;
- **depth**: how far below the shallowest zone of the grid a zone's potential lies, m^2/s^2 (the potential has no zero of its own; the map shows wells, so the
  zone with the highest potential is the zero);
- **tidal**: the largest eigenvalue of the tidal tensor, 1/s^2: the stretch (positive) or, where every eigenvalue is negative, the squeeze of the least
  compressed direction.

The smooth galaxy's tidal tensor is the derivative of its pull, taken by a central difference over `TENSOR_STEP_PC` (the smooth field changes over kiloparsecs,
so the difference is exact to the table's own accuracy). Why a `centered` reading (Boss, 2026-11-06; `reports/gravity-range-2026-11-06/`): the smooth galaxy's slope across a sector is 700 to 900 times the stars' wells in the depth, so
colours over the total depth are one ramp; with the smooth galaxy held at the sector's center the colours show the wells, and the pull and tidal strength differ from the total's
by the galaxy's own change across 4 parsecs only. The **saddle points** are the zones where the pull is smallest against its neighbours' and the tidal
tensor has both signs: the places a body would sit between wells (the Lagrange points of the map).

The reading of the local field alone (`which="local"`) is kept for the cases where the smooth galaxy swamps the stars: it is 99.9 percent of the pull and all of
the depth inside a sector, and about 1/750 of the tidal strength within 0.25 pc of a star (`reports/map383-background-rule-2026-11-05.md`).
"""

import numpy as np

from planetgen.galaxy import potential as galactic
from planetgen.physics import constants

PULL_X, PULL_Y, PULL_Z, POTENTIAL, T_XX, T_YY, T_ZZ, T_XY, T_XZ, T_YZ = range(10)
"""int: The columns of a grid row (`galaxy.gravity_kernel`'s; repeated here because that module imports numba)."""
MODES = ("pull", "depth", "tidal")
"""tuple: The scalars the map draws."""
UNITS = {"pull": "m/s^2", "depth": "m^2/s^2", "tidal": "1/s^2"}
"""dict: The unit of each mode's numbers."""
TENSOR_STEP_PC = 0.5
"""float: The step the smooth galaxy's pull is differenced over for its tidal tensor, parsecs."""
MAX_SADDLES = 64
"""int: The most saddle points an answer lists (the ones nearest zero pull first)."""
POTENTIAL_TO_SI = 1.0e6
"""float: One (km/s)^2 in m^2/s^2."""
TENSOR_TO_SI = 1.0e6 / constants.PARSEC_M ** 2
"""float: One (km/s)^2 per pc^2 in 1/s^2."""


def smooth_fields(points, pot):
    """The smooth galaxy's ten numbers at `points` `(n, 3)` in `gravity_kernel`'s units: the pull, the potential and the tidal tensor `d a_i / d x_j`."""
    points = np.asarray(points, dtype=np.float64)
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    ax, ay, az = galactic.acceleration(pot, x, y, z)
    out = np.empty((points.shape[0], 10))
    out[:, PULL_X], out[:, PULL_Y], out[:, PULL_Z] = ax, ay, az
    out[:, POTENTIAL] = galactic.potential(pot, x, y, z)
    step = TENSOR_STEP_PC
    columns = {}
    for axis, offset in enumerate(np.eye(3) * step):
        plus = np.stack(galactic.acceleration(pot, x + offset[0], y + offset[1], z + offset[2]), axis=1)
        minus = np.stack(galactic.acceleration(pot, x - offset[0], y - offset[1], z - offset[2]), axis=1)
        columns[axis] = (plus - minus) / (2.0 * step)           # d a / d x_axis: column `axis` of the tensor, shape (n, 3)
    # the derivative of a pull is symmetric; the two off-diagonal estimates are averaged so a rounding difference does not make it otherwise
    out[:, T_XX], out[:, T_YY], out[:, T_ZZ] = columns[0][:, 0], columns[1][:, 1], columns[2][:, 2]
    out[:, T_XY] = 0.5 * (columns[1][:, 0] + columns[0][:, 1])
    out[:, T_XZ] = 0.5 * (columns[2][:, 0] + columns[0][:, 2])
    out[:, T_YZ] = 0.5 * (columns[2][:, 1] + columns[1][:, 2])
    return out


FIELDS = ("total", "centered", "local")
"""tuple: The readings of the field: the `total` (the smooth galaxy added at each zone), the `centered` (the smooth galaxy added as it is at the sector's center, the same for every zone, so the colours show what the stars add: the reading the map colours by) and the `local` (the stars alone)."""


def field(local, points, pot, which="total"):
    """The field the map shows, `(n, 10)` float64: `local` (a grid's numbers) plus the smooth galaxy at each of `points` (`which` "total"), at their mean point (`"centered"`)
    or not at all (`"local"`)."""
    local = np.asarray(local, dtype=np.float64)
    if which == "local":
        return local
    if which == "centered":
        points = np.asarray(points, dtype=np.float64)
        return local + smooth_fields(points.mean(axis=0, keepdims=True), pot)
    if which != "total":
        raise ValueError(f"unknown field {which!r}; one of {', '.join(FIELDS)}")
    return local + smooth_fields(points, pot)


def tensors(result):
    """The tidal tensors of a grid's rows as an `(n, 3, 3)` array (symmetric)."""
    t = np.empty((result.shape[0], 3, 3))
    t[:, 0, 0], t[:, 1, 1], t[:, 2, 2] = result[:, T_XX], result[:, T_YY], result[:, T_ZZ]
    t[:, 0, 1] = t[:, 1, 0] = result[:, T_XY]
    t[:, 0, 2] = t[:, 2, 0] = result[:, T_XZ]
    t[:, 1, 2] = t[:, 2, 1] = result[:, T_YZ]
    return t


def pull_size(result):
    """The size of the pull, (km/s)^2 per pc, at each row of a grid (`gravity_kernel`'s units)."""
    return np.sqrt(result[:, PULL_X] ** 2 + result[:, PULL_Y] ** 2 + result[:, PULL_Z] ** 2)


def tidal_strength(result):
    """The largest eigenvalue of the tidal tensor, (km/s)^2 per pc^2, at each row of a grid: the stretch the map's Tidal strength mode draws."""
    return np.linalg.eigvalsh(tensors(result))[:, 2] if result.shape[0] else np.empty(0)


def scalar(result, mode):
    """One number per row of `result` for `mode` (see the module doc), in the mode's SI unit (`UNITS`)."""
    if mode == "pull":
        return pull_size(result) * galactic.ACCELERATION_TO_SI
    if mode == "depth":
        potential = result[:, POTENTIAL]
        return (potential.max() - potential) * POTENTIAL_TO_SI
    if mode == "tidal":
        return tidal_strength(result) * TENSOR_TO_SI
    raise ValueError(f"unknown mode {mode!r}; one of {', '.join(MODES)}")


def saddles(result, zones):
    """
    The saddle points of the grid `result` (`zones^3` rows, zone `(i, j, k)` at row `(i * zones + j) * zones + k`): the zones whose pull is smaller than each of
    their 26 neighbours' and whose tidal tensor has an eigenvalue of each sign, the ones with the smallest pull against their tidal strength first, at most
    `MAX_SADDLES`. Returns a list of `(i, j, k)`.
    """
    pull = pull_size(result).reshape(zones, zones, zones)
    padded = np.pad(pull, 1, mode="constant", constant_values=np.inf)
    minimum = np.ones(pull.shape, dtype=bool)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if (dx, dy, dz) != (0, 0, 0):
                    neighbor = padded[1 + dx:1 + dx + zones, 1 + dy:1 + dy + zones, 1 + dz:1 + dz + zones]
                    # of two zones with exactly the same pull the one earlier in the order is the minimum, so a tie names one zone, not both
                    minimum &= (pull <= neighbor) if (dx, dy, dz) > (0, 0, 0) else (pull < neighbor)
    eigenvalues = np.linalg.eigvalsh(tensors(result)).reshape(zones, zones, zones, 3)
    indefinite = (eigenvalues[..., 0] < 0.0) & (eigenvalues[..., 2] > 0.0)
    found = np.argwhere(minimum & indefinite)
    if not len(found):
        return []
    strength = np.abs(eigenvalues).max(axis=3)
    ranking = sorted(range(len(found)), key=lambda n: pull[tuple(found[n])] / max(strength[tuple(found[n])], 1.0e-300))
    return [tuple(int(v) for v in found[n]) for n in ranking[:MAX_SADDLES]]
