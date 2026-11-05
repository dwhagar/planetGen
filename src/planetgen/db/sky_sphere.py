# planetgen/db/sky_sphere.py

"""
The stars of a sphere above a luminosity floor, as arrays (VIEW.22, System sky S4)
==================================================================================

The sky a system sees (`reports/system-sky-from-the-system-2026-11-04.md`) is made of the stars around it that are bright enough to be seen.
`stars_in_sphere(conn, center, radius_pc, min_luminosity_sol)` is the one read of the stored stars that the sky needs: every star and phenomenon
inside a sphere whose luminosity is at least the floor, from the four places the project keeps them, with the cell levels the sphere reached.

The four sources, as `SOURCES` numbers them:

* `CHARTED` (0): the stars of charted systems (`stars`), at the point of their `star_systems` row and with its velocity and epoch. A system is
  read only when one of its stars is above the floor. A star of a cluster carries its `cluster_id`.
* `BRIGHT` (1): the pre-placed stars no system has been built for (`bright_stars` with no `star_system_id`). One that has been built is a
  charted star and is read there, so no star is counted twice. A cluster member carries its `cluster_id`. The velocity columns are NULL until
  the star's system is built; they come back as `nan`.
* `SCATTER` (2): the scattered phenomena not built yet (`phenomenon_scatter` with no `built_at`) that have a map luminosity: neutron stars, black holes,
  quasars and the rest. `kind` names the phenomenon.

Each source is one statement on `db/region.py` (at most 256 runs of cell keys, a covering index): the sphere is turned into key runs once and
the same runs serve the three tables. The floor is part of the statement, so a row below it is never read. The SQL test of the distance is a little
loose, so the exact distance is measured here with `detmath` and a row outside the radius is dropped.

**The cell levels.** A sector's stars are complete only down to the level its backfill reached (`sector_stats.bright_level_sol`: -1 untouched, so the
galaxy scatter's own level; a positive number the dimmest luminosity drawn; 0 filled) and, for the mass rings, `bright_mass_sol`.
`CellLevels` lists each sector that can hold part of the sphere (`spherical.cells_touching_sphere`) with those two numbers and the luminosity above which
its stars are complete (`complete_above_sol`, `inf` when nothing was drawn). `complete_above_sol` of the whole sphere is the largest of them: a star
brighter than it is stored wherever it lies in the sphere, a fainter one may not be (the sky tiers, VIEW.24 and VIEW.25, fill the rest).

The result is plain data, sorted by distance and then by source and ID, so the same database always gives the same arrays. The sky catalogue
(`galaxy/sky.py`, VIEW.23) joins it with the photometry (VIEW.20) and the light-travel positions (VIEW.21). Reads only; no row is changed.
"""

import dataclasses
import math

import numpy as np

from planetgen.db import region, store
from planetgen.galaxy import cellkey, spherical
from planetgen.physics import unit_registry as units
from planetgen.util import detmath

CHARTED, BRIGHT, SCATTER = 0, 1, 2
SOURCES = {CHARTED: "charted", BRIGHT: "bright_star", SCATTER: "scatter"}
"""dict: The source of each row, as the `source` array numbers them."""

CHUNK = 500
"""int: How many system IDs one statement reads the stars of."""

_to_watts = units.converter("Lsun", "W")
_to_solar = units.converter("W", "Lsun")

_FIELDS = ("source", "ident", "system_id", "x_pc", "y_pc", "z_pc", "distance_pc", "vx_kms", "vy_kms", "vz_kms", "epoch_unix", "luminosity_sol",
           "temperature_k", "age_gy", "lifespan_gy", "cluster_id", "kind")


class SkySphereError(ValueError):
    """A sphere or a floor that cannot be searched."""


@dataclasses.dataclass(frozen=True)
class CellLevel:
    """One sector that can hold part of the sphere: its address, what its row of `sector_stats` says, and the luminosity above which its stars are complete."""

    address: tuple
    level_sol: float
    mass_sol: object
    complete_above_sol: float


@dataclasses.dataclass(frozen=True)
class CellLevels:
    """The sectors a sphere reached and how far each was drawn. `complete_above_sol` is the largest of the sectors' (`inf` when there are none to
    judge by, or when a sector was never drawn and the galaxy has no scatter level)."""

    cells: tuple
    scatter_level_sol: object
    complete_above_sol: float

    @property
    def untouched(self):
        """How many of the sectors no backfill has drawn (their stars are the galaxy scatter's alone)."""
        return sum(1 for cell in self.cells if cell.level_sol == store.UNTOUCHED_LEVEL)

    @property
    def filled(self):
        """How many of the sectors are filled (every star down to the smallest was drawn)."""
        return sum(1 for cell in self.cells if cell.level_sol == store.FILLED_LEVEL)


@dataclasses.dataclass(frozen=True)
class SphereStars:
    """
    The stars of a sphere as arrays of one length, in the order of `distance_pc` (then source, then ID). Positions are the galaxy frame, parsecs;
    velocities km/s and `epoch_unix` seconds, `nan` where nothing is stored; `luminosity_sol` in solar luminosities; `age_gy` and `lifespan_gy` in Gy
    (`nan` where none; a star that never ends has `lifespan_gy` of `inf`); `cluster_id` 0 for a star of no cluster; `kind` the phenomenon of a scatter
    row and `""` for a star. `levels` is the `CellLevels` of the sphere.
    """

    source: np.ndarray
    ident: np.ndarray
    system_id: np.ndarray
    x_pc: np.ndarray
    y_pc: np.ndarray
    z_pc: np.ndarray
    distance_pc: np.ndarray
    vx_kms: np.ndarray
    vy_kms: np.ndarray
    vz_kms: np.ndarray
    epoch_unix: np.ndarray
    luminosity_sol: np.ndarray
    temperature_k: np.ndarray
    age_gy: np.ndarray
    lifespan_gy: np.ndarray
    cluster_id: np.ndarray
    kind: np.ndarray
    levels: CellLevels
    center: tuple
    radius_pc: float
    min_luminosity_sol: float

    def __len__(self):
        return len(self.source)

    def count(self, source):
        """How many rows come from `source` (`CHARTED`, `BRIGHT` or `SCATTER`)."""
        return int(np.count_nonzero(self.source == source))

    def positions_pc(self):
        """The positions as one `n` by 3 array."""
        return np.column_stack((self.x_pc, self.y_pc, self.z_pc))

    def velocities_kms(self):
        """The velocities as one `n` by 3 array (`nan` where none is stored)."""
        return np.column_stack((self.vx_kms, self.vy_kms, self.vz_kms))


def _check(center, radius_pc, min_luminosity_sol):
    try:
        point = tuple(float(value) for value in center)
        radius = float(radius_pc)
        floor = float(min_luminosity_sol)
    except (TypeError, ValueError) as exc:
        raise SkySphereError("the center, the radius and the floor must be numbers") from exc
    if len(point) != 3 or not all(math.isfinite(value) for value in point):
        raise SkySphereError("the center is three finite numbers of parsecs")
    if not math.isfinite(radius) or radius < 0.0:
        raise SkySphereError("the radius is a finite number of parsecs, zero or more")
    if not math.isfinite(floor) or floor < 0.0:
        raise SkySphereError("the luminosity floor is a finite number of solar luminosities, zero or more")
    return point, radius, floor


def _number(value, empty=math.nan):
    return empty if value is None else float(value)


def _chunks(items, size=CHUNK):
    items = list(items)
    for start in range(0, len(items), size):
        yield items[start:start + size]


def _charted(conn, center, radius, floor_w, runs):
    """`(rows, runs)`: one row per star of a system in the sphere that is at least `floor_w` W."""
    columns = "p.id, p.galaxy_x_pc, p.galaxy_y_pc, p.galaxy_z_pc, p.velocity_x_kms, p.velocity_y_kms, p.velocity_z_kms, p.epoch_unix"
    position = ("p.galaxy_x_pc", "p.galaxy_y_pc", "p.galaxy_z_pc")
    where = "EXISTS (SELECT 1 FROM stars s WHERE s.star_system_id = p.id AND s.luminosity_w >= ?)"
    systems = {}
    for row in region.region_rows(conn, "star_systems", columns, position, center, radius, where, (floor_w,), runs):
        point = (row["galaxy_x_pc"], row["galaxy_y_pc"], row["galaxy_z_pc"])
        away = detmath.dist(center, point)
        if away <= radius:
            systems[row["id"]] = (point, away, row)
    rows = []
    for chunk in _chunks(sorted(systems)):
        marks = ", ".join("?" for _ in chunk)
        for star in conn.execute(
                "SELECT id, star_system_id, luminosity_w, temperature_k, age_gy, lifespan_gy, cluster_id FROM stars"
                f" WHERE star_system_id IN ({marks}) AND luminosity_w >= ?", (*chunk, floor_w)).fetchall():
            point, away, system = systems[star["star_system_id"]]
            rows.append((CHARTED, star["id"], star["star_system_id"], *point, away, _number(system["velocity_x_kms"]), _number(system["velocity_y_kms"]),
                         _number(system["velocity_z_kms"]), _number(system["epoch_unix"]), _to_solar(star["luminosity_w"]),
                         _number(star["temperature_k"]), _number(star["age_gy"]), _number(star["lifespan_gy"], math.inf), star["cluster_id"] or 0, ""))
    return rows


def _mpc(row, axis):
    return row[axis] / 1000.0


def _bright(conn, center, radius, floor_w, runs):
    """The unbuilt pre-placed stars in the sphere at least `floor_w` W."""
    position = ("p.position_x_mpc / 1000.0", "p.position_y_mpc / 1000.0", "p.position_z_mpc / 1000.0")
    columns = ("p.id, p.position_x_mpc, p.position_y_mpc, p.position_z_mpc, p.velocity_x_kms, p.velocity_y_kms, p.velocity_z_kms, p.epoch_unix, "
               "p.luminosity_w, p.temperature_k, p.age_gy, p.lifespan_gy, p.cluster_id")
    rows = []
    for row in region.region_rows(conn, "bright_stars", columns, position, center, radius,
                                  "p.star_system_id IS NULL AND p.luminosity_w >= ?", (floor_w,), runs):
        point = (_mpc(row, "position_x_mpc"), _mpc(row, "position_y_mpc"), _mpc(row, "position_z_mpc"))
        away = detmath.dist(center, point)
        if away <= radius:
            rows.append((BRIGHT, row["id"], 0, *point, away, _number(row["velocity_x_kms"]), _number(row["velocity_y_kms"]),
                         _number(row["velocity_z_kms"]), _number(row["epoch_unix"]), _to_solar(row["luminosity_w"]), _number(row["temperature_k"]),
                         _number(row["age_gy"]), _number(row["lifespan_gy"], math.inf), row["cluster_id"] or 0, ""))
    return rows


def _scatter(conn, center, radius, floor_sol, runs):
    """The unbuilt scattered phenomena in the sphere with a map luminosity of at least `floor_sol`."""
    position = ("p.position_x_mpc / 1000.0", "p.position_y_mpc / 1000.0", "p.position_z_mpc / 1000.0")
    columns = ("p.id, p.kind, p.position_x_mpc, p.position_y_mpc, p.position_z_mpc, p.velocity_x_kms, p.velocity_y_kms, p.velocity_z_kms, "
               "p.epoch_unix, p.map_luminosity_sol")
    rows = []
    for row in region.region_rows(conn, "phenomenon_scatter", columns, position, center, radius,
                                  "p.built_at IS NULL AND p.map_luminosity_sol >= ?", (floor_sol,), runs):
        point = (_mpc(row, "position_x_mpc"), _mpc(row, "position_y_mpc"), _mpc(row, "position_z_mpc"))
        away = detmath.dist(center, point)
        if away <= radius:
            rows.append((SCATTER, row["id"], 0, *point, away, _number(row["velocity_x_kms"]), _number(row["velocity_y_kms"]),
                         _number(row["velocity_z_kms"]), _number(row["epoch_unix"]), float(row["map_luminosity_sol"]), math.nan, math.nan,
                         math.nan, 0, row["kind"]))
    return rows


def _edge_pc(conn):
    shape = store.get_galaxy_shape(conn, light=True)
    return float(shape.edge_pc) if shape is not None else store._edge_pc(conn)


def cell_levels(conn, center, radius_pc):
    """
    The `CellLevels` of the sphere: each sector that can hold part of it, with its `sector_stats` level and mass and the luminosity above which its
    stars are complete. A galaxy with no sectors gives no cells and `inf`.
    """
    edge = _edge_pc(conn)
    scatter = store.bright_star_scatter_settings(conn)
    scatter_level = None if scatter is None else float(scatter[0])
    if edge <= 0.0:
        return CellLevels((), scatter_level, math.inf)
    addresses = sorted((shell, band, wedge) for shell, band, wedge, *_rest in spherical.cells_touching_sphere(center, radius_pc, edge))
    levels = store.sector_bright_levels(conn, addresses)
    masses = store.sector_bright_masses(conn, addresses)
    cells = []
    for address in addresses:
        level = levels.get(address, store.UNTOUCHED_LEVEL)
        if level == store.UNTOUCHED_LEVEL:
            complete = math.inf if scatter_level is None else scatter_level
        else:
            complete = float(level)
        cells.append(CellLevel(address, float(level), masses.get(address), complete))
    return CellLevels(tuple(cells), scatter_level, max((cell.complete_above_sol for cell in cells), default=math.inf))


def stars_in_sphere(conn, center, radius_pc, min_luminosity_sol):
    """
    Every star and phenomenon within `radius_pc` of `center` whose luminosity is at least `min_luminosity_sol`, as a `SphereStars`.

    Args:
        conn (Connection): An open connection (a read-only one will do).
        center (sequence): `(x, y, z)`, parsecs, galaxy frame.
        radius_pc (float): Zero or more.
        min_luminosity_sol (float): The floor, solar luminosities, zero or more. A star at exactly the floor is in.

    Returns:
        SphereStars: the rows of the four sources, nearest first, and the cell levels the sphere reached.

    Raises:
        SkySphereError: A center, radius, or floor that is not a finite number.
    """
    center, radius, floor = _check(center, radius_pc, min_luminosity_sol)
    floor_w = _to_watts(floor)
    runs = cellkey.sphere_runs(center, radius)
    rows = (_charted(conn, center, radius, floor_w, runs) + _bright(conn, center, radius, floor_w, runs)
            + _scatter(conn, center, radius, floor, runs))
    rows.sort(key=lambda row: (row[6], row[0], row[1]))
    columns = list(zip(*rows)) if rows else [()] * len(_FIELDS)
    dtypes = {"source": np.int8, "ident": np.int64, "system_id": np.int64, "cluster_id": np.int64, "kind": object}
    arrays = {name: np.array(values, dtype=dtypes.get(name, np.float64)) for name, values in zip(_FIELDS, columns)}
    return SphereStars(levels=cell_levels(conn, center, radius), center=center, radius_pc=radius, min_luminosity_sol=floor, **arrays)
