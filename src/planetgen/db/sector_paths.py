# planetgen/db/sector_paths.py

"""
Saved sector paths (GEN.123)
============================

`compute_sector_paths` works out, for every star system, rogue planet and
interstellar comet in a sector, the path it takes from where it is now to
where it leaves the sector (`physics/sector_path.py`), against the point
masses of the sector and the sectors around it, and saves the spline knots
in `sector_paths` / `sector_path_knots`. `load_sector_paths` reads them
back as `SectorPath`s.

A path starts at the body's stored position and velocity, so it is as old
as the last time paths were computed for its sector: the orbit update
(`planetgen.cli.orbits`) recomputes them for every sector holding a body,
and a `galaxy` run settles the sectors it created and the sectors around
them as its last step (`settle_sectors`, GEN.126). They are never computed
sector by sector during the run: a path depends on the neighbouring sectors
that exist at the time, which varies with the order workers fill sectors
in, and generated data must not.

Point masses are the stars of each star system (summed per system), black
holes and neutron stars. The stored velocity of a star system is its
galactic velocity (`star_systems.velocity_*_kms`, v61); a rogue planet or
comet moves at its stored vector.
"""

import concurrent.futures
import multiprocessing

from planetgen.db import kind_tables, point_columns
from planetgen.galaxy import spherical as geometry
from planetgen.physics import sector_path as sp
from planetgen.physics import unit_registry as units

_msun_to_kg = units.converter("Msun", "kg")

_pc_to_m = units.converter("pc", "m")
_yr_to_s = units.converter("yr", "s")

PATH_TABLES = kind_tables.tables("path")
"""tuple: The tables whose rows get a path (the `sector_paths.object_table` values)."""

NEIGHBOUR_REACH_EDGES = 2.5
"""float: Point masses are gathered from sectors whose center is within this
many edges of the sector's, along each axis."""

_PC = _pc_to_m(1.0)
_YEAR = _yr_to_s(1.0)


def _sector(conn, sector_id):
    row = conn.execute(
        "SELECT id, edge_mpc, shell_index, band_index, wedge_index, center_x_pc, center_y_pc, center_z_pc"
        " FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    if row is None:
        raise ValueError(f"no sectors row with id {sector_id}")
    return row


def _center(row):
    return (row["center_x_pc"], row["center_y_pc"], row["center_z_pc"])


def _neighbour_centers(conn, sector, edge_pc):
    """`{sector id: center_pc}` for the sector and every placed sector near it."""
    reach = NEIGHBOUR_REACH_EDGES * edge_pc
    cx, cy, cz = _center(sector)
    rows = conn.execute(
        "SELECT id, center_x_pc, center_y_pc, center_z_pc FROM sectors"
        " WHERE center_x_pc BETWEEN ? AND ? AND center_y_pc BETWEEN ? AND ? AND center_z_pc BETWEEN ? AND ?",
        (cx - reach, cx + reach, cy - reach, cy + reach, cz - reach, cz + reach)).fetchall()
    return {row["id"]: _center(row) for row in rows}


def _point_masses(conn, centers):
    """Every star system, black hole and neutron star in the sectors `centers`
    names, as `PointMass`es (metres, kilograms) keyed `(table, id)`."""
    ids = sorted(centers)
    marks = ", ".join("?" for _ in ids)
    masses = []
    for row in conn.execute(
        f"SELECT ss.id, ss.sector_id, ss.position_x_mpc AS x, ss.position_y_mpc AS y, ss.position_z_mpc AS z,"
        f" SUM(s.mass_kg) AS mass_kg FROM star_systems ss JOIN stars s ON s.star_system_id = ss.id"
        f" WHERE ss.sector_id IN ({marks}) AND ss.position_x_mpc IS NOT NULL GROUP BY ss.id, ss.sector_id, ss.position_x_mpc, ss.position_y_mpc, ss.position_z_mpc", ids
    ).fetchall():
        point = geometry.local_to_galaxy_pc(centers[row["sector_id"]], (row["x"] / 1000.0, row["y"] / 1000.0,
                                                                         row["z"] / 1000.0))
        masses.append(sp.PointMass(tuple(c * _PC for c in point), row["mass_kg"], ("star_systems", row["id"])))
    for table in ("black_holes", "neutron_stars"):
        for row in conn.execute(
            f"SELECT id, center_x_pc, center_y_pc, center_z_pc, mass_solar FROM {table}"
            f" WHERE sector_id IN ({marks}) AND center_x_pc IS NOT NULL", ids
        ).fetchall():
            masses.append(sp.PointMass(tuple(row[f"center_{axis}_pc"] * _PC for axis in "xyz"),
                                       _msun_to_kg(row["mass_solar"]), (table, row["id"])))
    return sp.MassTable(masses)


def _bodies(conn, sector_id, center):
    """`[(table, id, position_m, velocity_m_s)]` for everything in the sector that gets a path."""
    bodies = []
    for row in conn.execute(
        f"SELECT id, {point_columns.LAYOUTS['star_systems'].select('position', 'velocity')}"
        " FROM star_systems WHERE sector_id = ? AND position_x_mpc IS NOT NULL", (sector_id,)
    ).fetchall():
        point = geometry.local_to_galaxy_pc(center, (row["position_x_mpc"] / 1000.0, row["position_y_mpc"] / 1000.0,
                                                      row["position_z_mpc"] / 1000.0))
        bodies.append(("star_systems", row["id"], tuple(c * _PC for c in point),
                       tuple(row[f"velocity_{axis}_kms"] * 1000.0 for axis in "xyz")))
    for table in (name for name in kind_tables.tables("path") if name != "star_systems"):
        for row in conn.execute(
            f"SELECT id, center_x_pc, center_y_pc, center_z_pc, velocity_x_kms, velocity_y_kms, velocity_z_kms"
            f" FROM {table}"
            f" WHERE sector_id = ? AND center_x_pc IS NOT NULL", (sector_id,)
        ).fetchall():
            point = tuple(row[f"center_{axis}_pc"] for axis in "xyz")
            bodies.append((table, row["id"], tuple(c * _PC for c in point),
                           tuple(row[f"velocity_{axis}_kms"] * 1000.0 for axis in "xyz")))
    return bodies


PARALLEL_CHUNK = 250
"""int: Bodies one worker process integrates per task (PERF.80); a sector with
fewer than twice this many is integrated in the calling process."""


def _integrate_bodies(bodies, masses, address, edge_pc):
    """`[(table, object_id, SectorPath)]` for `bodies` (`_bodies`' rows), in their
    order: each against `masses`, skipping any that is filed here but past the
    cell's face (a rotated-frame rounding), whose own sector's pass handles it."""
    inside = sp.sector_inside(address, edge_pc)
    edge_m = edge_pc * _PC
    results = []
    for table, object_id, position, velocity in bodies:
        if not inside(position):
            continue
        results.append((table, object_id, sp.integrate_path(position, velocity, masses, inside, edge_m,
                                                           exclude_key=(table, object_id))))
    return results


def _integrate_bodies_in_pool(pool, bodies, masses, address, edge_pc):
    """`_integrate_bodies`, split into `PARALLEL_CHUNK`-body tasks for the
    process pool `pool` when there are enough of them; the same results in the
    same order either way (nothing here is random or reads the clock)."""
    if pool is None or len(bodies) < 2 * PARALLEL_CHUNK:
        return _integrate_bodies(bodies, masses, address, edge_pc)
    futures = [pool.submit(_integrate_bodies, bodies[start:start + PARALLEL_CHUNK], masses, address, edge_pc)
               for start in range(0, len(bodies), PARALLEL_CHUNK)]
    return [result for future in futures for result in future.result()]


def compute_sector_paths(conn, sector_id, pool=None):
    """
    Rewrites the saved path of every star system, rogue planet and
    interstellar comet in sector `sector_id`: the sector's earlier paths
    are deleted first, so a body that left or was deleted loses its path.

    Args:
        conn (Connection): An open, schema-initialized, read-write connection.
        sector_id (int): The `sectors.id`.
        pool (concurrent.futures.Executor, optional): A process pool to
            integrate a crowded sector's bodies in (`settle_sectors` makes
            one for `workers` above 1); the paths saved are the same.

    Returns:
        int: How many paths were saved (0 for a sector that is not a placed
        cell of the galaxy's grid).

    Raises:
        ValueError: If there is no such sector.
    """
    sector = _sector(conn, sector_id)
    _delete_paths(conn, "SELECT id FROM sector_paths WHERE sector_id = ?", (sector_id,))
    if sector["center_x_pc"] is None or sector["shell_index"] is None or not sector["edge_mpc"]:
        return 0
    edge_pc = sector["edge_mpc"] / 1000.0
    centers = _neighbour_centers(conn, sector, edge_pc)
    masses = _point_masses(conn, centers)
    address = (sector["shell_index"], sector["band_index"], sector["wedge_index"])
    results = _integrate_bodies_in_pool(pool, _bodies(conn, sector_id, _center(sector)), masses, address, edge_pc)
    _save_paths(conn, sector_id, results)
    return len(results)


def _delete_paths(conn, select_sql, params):
    """Deletes the paths `select_sql` finds, by primary key: a DELETE with a
    non-key WHERE on a missing row takes gap locks, and workers filling
    neighbouring sectors at once deadlock on them."""
    ids = [row["id"] for row in conn.execute(select_sql, params).fetchall()]
    for path_id in ids:
        conn.execute("DELETE FROM sector_paths WHERE id = ?", (path_id,))


_ROWS_PER_LOOKUP = 500
"""int: Bodies named in one lookup of their old paths."""


def _save_paths(conn, sector_id, results):
    """
    Saves `results`, `[(table, object_id, SectorPath)]`, as sector
    `sector_id`'s paths in a few statements for the whole sector (PERF.80: one
    lookup, delete and pair of inserts per body made the saving a third of
    the settle step): any path the same bodies had in another sector (a body
    that moved) goes first, by primary key, then the `sector_paths` rows go in
    together, their ids are read back, and the knots go in together.
    """
    if not results:
        return
    for table in PATH_TABLES:
        object_ids = [object_id for t, object_id, _path in results if t == table]
        for start in range(0, len(object_ids), _ROWS_PER_LOOKUP):
            chunk = object_ids[start:start + _ROWS_PER_LOOKUP]
            marks = ", ".join("?" for _ in chunk)
            stale = [row["id"] for row in conn.execute(
                f"SELECT id FROM sector_paths WHERE object_table = ? AND object_id IN ({marks})",
                (table, *chunk)).fetchall()]
            for path_id in stale:
                conn.execute("DELETE FROM sector_paths WHERE id = ?", (path_id,))
    conn.executemany(
        "INSERT INTO sector_paths (sector_id, object_table, object_id, exited, duration_years)"
        " VALUES (?, ?, ?, ?, ?)",
        [(sector_id, table, object_id, int(path.exited), path.duration_s / _YEAR) for table, object_id, path in results])
    path_ids = {(row["object_table"], row["object_id"]): row["id"] for row in conn.execute(
        "SELECT id, object_table, object_id FROM sector_paths WHERE sector_id = ?", (sector_id,)).fetchall()}
    conn.executemany(
        "INSERT INTO sector_path_knots (path_id, position, t_years, x_pc, y_pc, z_pc, vx_kms, vy_kms, vz_kms)"
        " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [(path_ids[(table, object_id)], index, knot.t_s / _YEAR, *(c / _PC for c in knot.position),
          *(v / 1000.0 for v in knot.velocity))
         for table, object_id, path in results for index, knot in enumerate(path.knots)])


def sectors_to_settle(conn, sector_ids):
    """
    The sectors whose saved paths a change to `sector_ids` can alter: the
    sectors themselves and every sector within the neighbour reach of one
    (a body's path sees the point masses of the sectors around its own).

    Returns:
        list: Sector ids, ascending. Ids that are not placed galaxy sectors
        are kept (their paths clear to nothing).
    """
    wanted = set(sector_ids)
    for sector_id in sector_ids:
        sector = _sector(conn, sector_id)
        if sector["center_x_pc"] is None or not sector["edge_mpc"]:
            continue
        wanted.update(_neighbour_centers(conn, sector, sector["edge_mpc"] / 1000.0))
    return sorted(wanted)


def settle_sectors(conn, sector_ids, on_progress=None, expand=True, workers=1):
    """
    Recomputes the saved paths of `sector_ids` and the sectors around them
    (`sectors_to_settle`), committing every sector as it finishes: the last
    step of a generation run or of an edit that changes the masses in a
    sector (GEN.126), when the neighbour set is final.
    Nothing here is random or depends on the clock, so it leaves a
    reproducible galaxy reproducible.

    Args:
        conn (Connection): An open, schema-initialized, read-write connection.
        sector_ids (iterable): The sectors a run created or changed.
        on_progress (callable, optional): `on_progress(done, total)` after each sector.
        expand (bool, optional): Also settle the sectors around `sector_ids`
            (the default). With `False`, exactly `sector_ids` (those that still exist),
            for a caller that worked the neighbours out before deleting a sector.
        workers (int, optional): Above 1, a crowded sector's bodies are
            integrated in that many worker processes (PERF.80); the paths
            saved are the same, the saving and the order stay in this process.

    Returns:
        int: How many paths were saved.
    """
    if expand:
        sectors = sectors_to_settle(conn, list(sector_ids))
    else:
        sectors = sorted(set(sector_ids))
        existing = {row["id"] for row in conn.execute(
            f"SELECT id FROM sectors WHERE id IN ({', '.join('?' for _ in sectors)})", sectors).fetchall()} if sectors else set()
        sectors = [sector_id for sector_id in sectors if sector_id in existing]
    saved = 0
    pool = None
    if workers > 1 and len(sectors) > 0:
        pool = concurrent.futures.ProcessPoolExecutor(max_workers=workers,
                                                      mp_context=multiprocessing.get_context("spawn"))
    try:
        for done, sector_id in enumerate(sectors, start=1):
            saved += compute_sector_paths(conn, sector_id, pool)
            conn.commit()
            if on_progress is not None:
                on_progress(done, len(sectors))
    finally:
        if pool is not None:
            pool.shutdown()
    return saved


def sector_ids_since(conn, since):
    """The ids of every galaxy sector created at or after `since` (`database_now`), ascending."""
    return [row["id"] for row in conn.execute(
        "SELECT id FROM sectors WHERE created_at >= ? AND shell_index IS NOT NULL ORDER BY id", (since,)).fetchall()]


def sector_addresses_since(conn, since):
    """The `(shell, band, wedge)` of every galaxy sector created at or after `since` (`database_now`), in id order."""
    return [(row["shell_index"], row["band_index"], row["wedge_index"]) for row in conn.execute(
        "SELECT shell_index, band_index, wedge_index FROM sectors WHERE created_at >= ? AND shell_index IS NOT NULL"
        " AND band_index IS NOT NULL AND wedge_index IS NOT NULL ORDER BY id", (since,)).fetchall()]


def load_sector_paths(conn, sector_id):
    """
    The saved paths of sector `sector_id`.

    Returns:
        dict: `(object table, object id)` -> `SectorPath` (metres, seconds,
        galactic axes), for every body that has one.
    """
    paths = {}
    rows = conn.execute(
        "SELECT p.object_table, p.object_id, p.exited, k.position, k.t_years, k.x_pc, k.y_pc, k.z_pc,"
        " k.vx_kms, k.vy_kms, k.vz_kms FROM sector_paths p JOIN sector_path_knots k ON k.path_id = p.id"
        " WHERE p.sector_id = ? ORDER BY p.id, k.position", (sector_id,)).fetchall()
    knots, exited = {}, {}
    for row in rows:
        key = (row["object_table"], row["object_id"])
        exited[key] = bool(row["exited"])
        knots.setdefault(key, []).append(sp.PathKnot(
            row["t_years"] * _YEAR, (row["x_pc"] * _PC, row["y_pc"] * _PC, row["z_pc"] * _PC),
            (row["vx_kms"] * 1000.0, row["vy_kms"] * 1000.0, row["vz_kms"] * 1000.0)))
    for key, key_knots in knots.items():
        paths[key] = sp.SectorPath(key_knots, exited[key])
    return paths
