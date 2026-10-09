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
(`planetgen.cli.orbits`) recomputes them for every sector holding a body.
Generation does not compute them: a path depends on the
neighbouring sectors that exist at the time, which varies with the order
workers fill sectors in, and generated data must not.

Point masses are the stars of each star system (summed per system), black
holes and neutron stars. The stored velocity of a star system is its
galactic velocity (`star_systems.velocity_*_kms`, v61); a rogue planet or
comet follows the galaxy's rotation at its own `galactic_orbital_speed_kms`.
"""

from planetgen.galaxy import geometry
from planetgen.galaxy.system_position import galactic_velocity_ms
from planetgen.physics import constants
from planetgen.physics import sector_path as sp

PATH_TABLES = ("star_systems", "rogue_planets", "interstellar_comets")
"""tuple: The tables whose rows get a path (the `sector_paths.object_table` values)."""

NEIGHBOUR_REACH_EDGES = 2.5
"""float: Point masses are gathered from sectors whose center is within this
many edges of the sector's, along each axis."""

_PC = constants.PARSEC_M
_YEAR = constants.SECONDS_PER_YEAR


def _sector(conn, sector_id):
    row = conn.execute(
        "SELECT id, edge_mpc, ring_index, layer_index, ring_slot_index, center_x_pc, center_y_pc, center_z_pc"
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
                                       row["mass_solar"] * constants.SOLAR_MASS_TO_KG, (table, row["id"])))
    return sp.MassTable(masses)


def _bodies(conn, sector_id, center):
    """`[(table, id, position_m, velocity_m_s)]` for everything in the sector that gets a path."""
    bodies = []
    for row in conn.execute(
        "SELECT id, position_x_mpc, position_y_mpc, position_z_mpc, velocity_x_kms, velocity_y_kms, velocity_z_kms"
        " FROM star_systems WHERE sector_id = ? AND position_x_mpc IS NOT NULL", (sector_id,)
    ).fetchall():
        point = geometry.local_to_galaxy_pc(center, (row["position_x_mpc"] / 1000.0, row["position_y_mpc"] / 1000.0,
                                                      row["position_z_mpc"] / 1000.0))
        bodies.append(("star_systems", row["id"], tuple(c * _PC for c in point),
                       tuple(row[f"velocity_{axis}_kms"] * 1000.0 for axis in "xyz")))
    for table in ("rogue_planets", "interstellar_comets"):
        for row in conn.execute(
            f"SELECT id, center_x_pc, center_y_pc, center_z_pc, galactic_orbital_speed_kms FROM {table}"
            f" WHERE sector_id = ? AND center_x_pc IS NOT NULL", (sector_id,)
        ).fetchall():
            point = tuple(row[f"center_{axis}_pc"] for axis in "xyz")
            bodies.append((table, row["id"], tuple(c * _PC for c in point),
                           galactic_velocity_ms(point, row["galactic_orbital_speed_kms"])))
    return bodies


def compute_sector_paths(conn, sector_id):
    """
    Rewrites the saved path of every star system, rogue planet and
    interstellar comet in sector `sector_id`: the sector's earlier paths
    are deleted first, so a body that left or was deleted loses its path.

    Args:
        conn (Connection): An open, schema-initialized, read-write connection.
        sector_id (int): The `sectors.id`.

    Returns:
        int: How many paths were saved (0 for a sector that is not a placed
        cell of the galaxy's grid).

    Raises:
        ValueError: If there is no such sector.
    """
    sector = _sector(conn, sector_id)
    _delete_paths(conn, "SELECT id FROM sector_paths WHERE sector_id = ?", (sector_id,))
    if sector["center_x_pc"] is None or sector["ring_index"] is None or not sector["edge_mpc"]:
        return 0
    edge_pc = sector["edge_mpc"] / 1000.0
    edge_m = edge_pc * _PC
    centers = _neighbour_centers(conn, sector, edge_pc)
    masses = _point_masses(conn, centers)
    inside = sp.sector_inside((sector["ring_index"], sector["layer_index"], sector["ring_slot_index"]), edge_pc)
    saved = 0
    for table, object_id, position, velocity in _bodies(conn, sector_id, _center(sector)):
        if not inside(position):
            continue  # filed here but past the cell's face (a rotated-frame rounding): its own sector's pass handles it
        path = sp.integrate_path(position, velocity, masses, inside, edge_m, exclude_key=(table, object_id))
        _save_path(conn, sector_id, table, object_id, path)
        saved += 1
    return saved


def _delete_paths(conn, select_sql, params):
    """Deletes the paths `select_sql` finds, by primary key: a DELETE with a
    non-key WHERE on a missing row takes gap locks, and workers filling
    neighbouring sectors at once deadlock on them."""
    ids = [row["id"] for row in conn.execute(select_sql, params).fetchall()]
    for path_id in ids:
        conn.execute("DELETE FROM sector_paths WHERE id = ?", (path_id,))


def _save_path(conn, sector_id, table, object_id, path):
    _delete_paths(conn, "SELECT id FROM sector_paths WHERE object_table = ? AND object_id = ?", (table, object_id))
    path_id = conn.execute(
        "INSERT INTO sector_paths (sector_id, object_table, object_id, exited, duration_years)"
        " VALUES (?, ?, ?, ?, ?)",
        (sector_id, table, object_id, int(path.exited), path.duration_s / _YEAR)).lastrowid
    conn.executemany(
        "INSERT INTO sector_path_knots (path_id, position, t_years, x_pc, y_pc, z_pc, vx_kms, vy_kms, vz_kms)"
        " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [(path_id, index, knot.t_s / _YEAR, *(c / _PC for c in knot.position), *(v / 1000.0 for v in knot.velocity))
         for index, knot in enumerate(path.knots)])


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
