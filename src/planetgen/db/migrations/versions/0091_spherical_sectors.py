"""Spherical shell sectors (DB.32, docs/design/spherical-sector-coordinates.md).

Schema v91. A sector's address `(ring_index, layer_index, ring_slot_index)` becomes `(shell_index, band_index,
wedge_index)` in `sectors`, `bright_stars`, `sector_stats` and `phenomenon_scatter`; `galaxy_shape.outer_ring_index`
becomes `outer_shell_index`; `galaxy_layer` and `galaxy_column` are replaced by `galaxy_shell`; `galaxy_block_stats` goes.

The database is CONVERTED, nothing is deleted. Every address means something else on the new grid, so the data moves:

* Every star system, phenomenon and facility keeps its place in the galaxy, its velocity and its epoch. A system's
  stored offset is worked out again from its place and its new sector's centre.
* The sector rows are rebuilt on shells: one for every cell that holds an object, and one for every cell whose centre
  lay inside a sector that was generated (so a generated, empty sector stays generated). The objects, and the planets,
  moons and belts under their stars, are linked to the new sector that holds them. A new sector takes the name of the
  old sector holding its centre (the second and later ones in the same old sector get " 2", " 3" ...); a system's
  location text takes the new sector's name.
* Object IDs are renumbered for the new birth sector (generated serials by system then phenomenon id; run-time births,
  and field clouds, which are now run-time births, from the new sector's counter). Body numbers inside a system are
  kept. The sector IDs are the new designations.
* The pre-placed bright stars and phenomenon scatter get the address of the cell holding their position. A star built
  into a system takes that system's cell. `sector_stats` is rebuilt for the generated sectors from their stars; the
  untouched and partly backfilled rows are dropped (a cell with no row is untouched), as a new scatter does.
* `sector_paths` of the old sectors go with them and come back with the next orbit update; `nearest_systems` follow
  their objects; `sector_system_counts` is counted again; the timing history of a deleted sector loses its sector.
* The outline (`galaxy_shell`) is built again from the stored shape.

The statements are in two groups: the table changes, which are safe to run again, and the conversion, which runs as one
transaction (it is skipped once `galaxy_layer` is gone, which is its last step) so a failure leaves the database as it was.
The old cylindrical grid is repeated here (`_old_*`) so the revision keeps meaning what it meant when it was written.
"""

import math

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner

revision = "0091"
down_revision = "0090"
branch_labels = None
depends_on = None

_PHENOMENA = ("black_holes", "neutron_stars", "nebulae", "supernova_remnants", "quasars", "rogue_planets",
              "interstellar_comets", "asteroid_fields")
_BODIES = ("stars", "planets", "moons", "asteroid_belts", "comets")
_CHUNK = 20000
_READDRESS_CHUNK = 100000
_SOLAR_LUMINOSITY_W = 3.828e26

# The object ID layout of galaxy/object_uid.py: 12/12/16 bits of sector, 28 of serial (2 of kind), 12 of body.
_BAND_BIAS = 1 << 11
_SERIAL_COUNT_BITS = 26
_GENERATED, _RUNTIME = 0, 1
_NO_SECTOR = (0, 0, 0)


# ---------------------------------------------------------------------------------------------------------------------
# The old cylindrical grid (galaxy/geometry.py as of v90), only what the conversion asks of it.

def _old_master_count(ring):
    c = 2 * math.pi * (ring + 0.5)
    master = 3
    while c >= 2 * master * 8:
        master *= 2
    return master


def _old_slot_count(ring):
    master = _old_master_count(ring)
    return max(master, master * round(2 * math.pi * (ring + 0.5) / master))


def _old_settle(value, guess, bounds):
    low, high = bounds(guess)
    if value < low:
        return guess - 1
    if value >= high:
        return guess + 1
    return guess


def _old_address_at(point, edge):
    x, y, z = point
    radius = math.sqrt(x * x + y * y)
    ring = _old_settle(radius, int(math.floor(radius / edge)), lambda i: (i * edge, (i + 1) * edge))
    layer = _old_settle(z, int(math.floor(z / edge + 0.5)), lambda j: ((j - 0.5) * edge, (j + 0.5) * edge))
    n = _old_slot_count(ring)
    step = 2 * math.pi / n
    theta = math.atan2(y, x) % (2 * math.pi)
    slot = min(n - 1, int(theta * n / (2 * math.pi)))
    if slot > 0 and theta < slot * step:
        slot -= 1
    elif slot < n - 1 and theta >= (slot + 1) * step:
        slot += 1
    return ring, layer, slot


# The local axes of a sector (galaxy/frame.py): +X radially outward, +Y toward increasing longitude, +Z galactic north.

def _axes(center):
    r = math.hypot(center[0], center[1])
    if r < 1e-12:
        return (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
    c, s = center[0] / r, center[1] / r
    return (c, s, 0.0), (-s, c, 0.0), (0.0, 0.0, 1.0)


def _to_galaxy(center, offset):
    ax, ay, az = _axes(center)
    return tuple(center[i] + offset[0] * ax[i] + offset[1] * ay[i] + offset[2] * az[i] for i in range(3))


def _to_local(center, point):
    delta = tuple(point[i] - center[i] for i in range(3))
    return tuple(sum(delta[i] * axis[i] for i in range(3)) for axis in _axes(center))


# ---------------------------------------------------------------------------------------------------------------------
# Small helpers

def _has_column(connection, table, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :t AND COLUMN_NAME = :c"), {"t": table, "c": column}).scalar()


def _has_table(connection, table):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = :t"),
        {"t": table}).scalar()


def _touch(connection, table):
    """`, modified_at = modified_at` for a table that has the column: an UPDATE then leaves the time alone."""
    return ", modified_at = modified_at" if _has_column(connection, table, "modified_at") else ""


def _has_index(connection, table, index):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.STATISTICS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :t AND INDEX_NAME = :i"), {"t": table, "i": index}).scalar()


_ADDRESS_INDEXED = ("bright_stars", "phenomenon_scatter")


def _drop_address_indexes(connection):
    """The two big tables' address indexes are dropped for the bulk update and built again after (far cheaper than
    keeping them up to date row by row); a DDL statement ends the conversion's transaction, so both are outside it."""
    for table in _ADDRESS_INDEXED:
        if _has_index(connection, table, f"idx_{table}_address"):
            connection.execute(sa.text(f"ALTER TABLE {table} DROP INDEX idx_{table}_address"))


def _add_address_indexes(connection):
    for table in _ADDRESS_INDEXED:
        if not _has_index(connection, table, f"idx_{table}_address"):
            connection.execute(sa.text(
                f"ALTER TABLE {table} ADD KEY idx_{table}_address (shell_index, band_index, wedge_index)"))


def _id_ranges(connection, table, where=""):
    """`(first, last)` id ranges of `_CHUNK` rows of `table`, for updates in pieces."""
    low, high = connection.execute(sa.text(f"SELECT MIN(id), MAX(id) FROM {table} {where}")).fetchone()
    if low is None:
        return
    for first in range(low, high + 1, _CHUNK):
        yield first, min(first + _CHUNK - 1, high)


# ---------------------------------------------------------------------------------------------------------------------
# Table changes

def _change_tables(connection):
    for table, kind, band_type in (("sectors", "NULL", "INT"), ("bright_stars", "NOT NULL", "SMALLINT"),
                                   ("sector_stats", "NOT NULL", "SMALLINT"),
                                   ("phenomenon_scatter", "NOT NULL", "SMALLINT")):
        if _has_column(connection, table, "ring_index"):
            connection.execute(sa.text(
                f"ALTER TABLE {table} CHANGE ring_index shell_index INT {kind},"
                f" CHANGE layer_index band_index {band_type} {kind},"
                f" CHANGE ring_slot_index wedge_index INT {kind}"))
    if _has_column(connection, "galaxy_shape", "outer_ring_index"):
        connection.execute(sa.text("ALTER TABLE galaxy_shape CHANGE outer_ring_index outer_shell_index INT NOT NULL"))
    connection.execute(sa.text(
        "CREATE TABLE IF NOT EXISTS galaxy_shell ("
        " shell_index INT NOT NULL PRIMARY KEY, high_band_index INT NOT NULL, cell_count BIGINT UNSIGNED NOT NULL)"
        " ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"))


# ---------------------------------------------------------------------------------------------------------------------
# The conversion

class _Conversion:
    """One run of the conversion on `connection`; see the module docstring."""

    def __init__(self, connection, edge_pc):
        from planetgen.galaxy import spherical

        self.c = connection
        self.edge = edge_pc
        self.spherical = spherical
        self.done = 0

    def sql(self, text, params=None):
        return self.c.execute(sa.text(text), params or {})

    def step(self, label):
        alembic_runner.report_progress(f"Converting to shells: {label}", self.done, 12)
        self.done += 1

    # -- the old sectors --------------------------------------------------------------------------------------------

    def read_old_sectors(self):
        self.sql("CREATE TEMPORARY TABLE old_sector ("
                 " id BIGINT UNSIGNED NOT NULL PRIMARY KEY, name VARCHAR(255) NOT NULL,"
                 " x DOUBLE NOT NULL, y DOUBLE NOT NULL, z DOUBLE NOT NULL, ring INT NOT NULL, layer INT NOT NULL,"
                 " slot INT NOT NULL, wiki_url VARCHAR(2048), version_key CHAR(22), planetgen_version VARCHAR(32),"
                 " python_version VARCHAR(32), platform VARCHAR(64)) ENGINE=InnoDB")
        self.sql("INSERT INTO old_sector SELECT id, name, center_x_pc, center_y_pc, center_z_pc, shell_index,"
                 " band_index, wedge_index, wiki_url, version_key, planetgen_version, python_version, platform"
                 " FROM sectors WHERE center_x_pc IS NOT NULL AND shell_index IS NOT NULL")
        self.first_new_id = (self.sql("SELECT COALESCE(MAX(id), 0) + 1 FROM sectors").scalar())
        # The old rows give up their addresses and IDs, which the new rows are about to take.
        self.sql("UPDATE sectors s JOIN old_sector o ON o.id = s.id"
                 " SET s.shell_index = NULL, s.band_index = NULL, s.wedge_index = NULL, s.uid = NULL")
        self.sql("CREATE TEMPORARY TABLE new_cell ("
                 " shell INT NOT NULL, band INT NOT NULL, wedge INT NOT NULL, container BIGINT UNSIGNED NULL,"
                 " PRIMARY KEY (shell, band, wedge)) ENGINE=InnoDB")
        self.sql("CREATE TEMPORARY TABLE target ("
                 " tbl VARCHAR(24) NOT NULL, id BIGINT UNSIGNED NOT NULL, shell INT NOT NULL, band INT NOT NULL,"
                 " wedge INT NOT NULL, px DOUBLE NULL, py DOUBLE NULL, pz DOUBLE NULL, PRIMARY KEY (tbl, id))"
                 " ENGINE=InnoDB")

    # -- where every object goes ------------------------------------------------------------------------------------

    def _place(self, table, rows):
        """`rows` are `(id, old sector id, point or None, offset or None, old centre)`; writes the targets and
        the cells they need."""
        targets, cells = [], {}
        for row_id, sector_id, point, offset, center in rows:
            address = self.spherical.sector_address_at(point if point is not None else center, self.edge)
            px = py = pz = None
            if offset is not None:
                new_center = self.spherical.sector_position_pc(*address, self.edge)
                px, py, pz = (round(v * 1000.0) for v in _to_local(new_center, point))
            targets.append({"t": table, "i": row_id, "a": address[0], "b": address[1], "w": address[2],
                            "x": px, "y": py, "z": pz})
            cells.setdefault(address, sector_id)
        if targets:
            self.sql("INSERT INTO target (tbl, id, shell, band, wedge, px, py, pz)"
                     " VALUES (:t, :i, :a, :b, :w, :x, :y, :z)", targets)
            self.sql("INSERT IGNORE INTO new_cell (shell, band, wedge, container) VALUES (:a, :b, :w, :c)",
                     [{"a": a[0], "b": a[1], "w": a[2], "c": s} for a, s in cells.items()])

    def place_systems(self):
        last = 0
        while True:
            rows = self.sql(
                "SELECT ss.id, ss.sector_id, ss.position_x_mpc, ss.position_y_mpc, ss.position_z_mpc, o.x, o.y, o.z"
                " FROM star_systems ss JOIN old_sector o ON o.id = ss.sector_id WHERE ss.id > :last"
                " ORDER BY ss.id LIMIT :n", {"last": last, "n": _CHUNK}).fetchall()
            if not rows:
                return
            batch = []
            for row_id, sector_id, px, py, pz, x, y, z in rows:
                center = (x, y, z)
                if px is None or py is None or pz is None:
                    batch.append((row_id, sector_id, None, None, center))
                else:
                    batch.append((row_id, sector_id, _to_galaxy(center, (px / 1000.0, py / 1000.0, pz / 1000.0)),
                                  True, center))
            self._place("star_systems", batch)
            last = rows[-1][0]

    def place_centered(self, table, where=""):
        """Phenomena, and facilities that are not in a system: placed by the centre they store."""
        last = 0
        while True:
            rows = self.sql(
                f"SELECT t.id, t.sector_id, t.center_x_pc, t.center_y_pc, t.center_z_pc, o.x, o.y, o.z"
                f" FROM {table} t JOIN old_sector o ON o.id = t.sector_id WHERE t.id > :last {where}"
                f" ORDER BY t.id LIMIT :n", {"last": last, "n": _CHUNK}).fetchall()
            if not rows:
                return
            batch = []
            for row_id, sector_id, cx, cy, cz, x, y, z in rows:
                point = None if cx is None or cy is None or cz is None else (cx, cy, cz)
                batch.append((row_id, sector_id, point, None, (x, y, z)))
            self._place(table, batch)
            last = rows[-1][0]

    def cover_old_sectors(self):
        """Every cell whose centre lies in an old sector is a sector too, with that old sector as its container."""
        last = 0
        while True:
            rows = self.sql("SELECT id, x, y, z, ring, layer, slot FROM old_sector WHERE id > :last ORDER BY id"
                            " LIMIT :n", {"last": last, "n": 2000}).fetchall()
            if not rows:
                return
            cells = []
            for old_id, x, y, z, ring, layer, slot in rows:
                for shell, band, wedge, cx, cy, cz, _distance in self.spherical.enumerate_sectors_within_radius(
                        (x, y, z), self.edge * 1.05, self.edge):
                    if _old_address_at((cx, cy, cz), self.edge) == (ring, layer, slot):
                        cells.append({"a": shell, "b": band, "w": wedge, "c": old_id})
            if cells:
                self.sql("INSERT INTO new_cell (shell, band, wedge, container) VALUES (:a, :b, :w, :c)"
                         " ON DUPLICATE KEY UPDATE container = VALUES(container)", cells)
            last = rows[-1][0]

    # -- the new sectors --------------------------------------------------------------------------------------------

    def insert_sectors(self):
        spherical = self.spherical
        edge_mpc = round(self.edge * 1000.0)
        used = {}
        first = (-(2 ** 31), 0, 0)
        while True:
            rows = self.sql(
                "SELECT n.shell, n.band, n.wedge, n.container, o.name, o.wiki_url, o.version_key,"
                " o.planetgen_version, o.python_version, o.platform FROM new_cell n"
                " LEFT JOIN old_sector o ON o.id = n.container"
                " WHERE (n.shell, n.band, n.wedge) > (:a, :b, :w) ORDER BY n.shell, n.band, n.wedge LIMIT :n",
                {"a": first[0], "b": first[1], "w": first[2], "n": _CHUNK}).fetchall()
            if not rows:
                return
            batch = []
            for shell, band, wedge, container, name, wiki, key, version, python, platform in rows:
                x, y, z = spherical.sector_position_pc(shell, band, wedge, self.edge)
                designation = spherical.provisional_sector_designation(shell, band, wedge)
                count = used.get(container, 0)
                used[container] = count + 1
                if name is None:
                    name = f"Sector {designation}"
                elif count:
                    name = f"{name} {count + 1}"
                batch.append({"name": name, "edge": edge_mpc, "x": x, "y": y, "z": z, "r": math.sqrt(x * x + y * y + z * z),
                              "a": shell, "b": band, "w": wedge, "wiki": wiki if not count else None,
                              "uid": int(designation, 16), "key": key, "ver": version, "py": python, "pl": platform})
            self.sql("INSERT INTO sectors (name, edge_mpc, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc,"
                     " shell_index, band_index, wedge_index, wiki_url, uid, version_key, planetgen_version,"
                     " python_version, platform) VALUES (:name, :edge, :x, :y, :z, :r, :a, :b, :w, :wiki, :uid,"
                     " :key, :ver, :py, :pl)", batch)
            first = tuple(rows[-1][:3])

    # -- linking ----------------------------------------------------------------------------------------------------

    def relink(self):
        join = ("JOIN target t ON t.tbl = '{table}' AND t.id = x.id"
                " JOIN sectors s ON s.shell_index = t.shell AND s.band_index = t.band AND s.wedge_index = t.wedge")
        for first, last in _id_ranges(self.c, "star_systems"):
            self.sql(
                "UPDATE star_systems x " + join.format(table="star_systems") +
                " SET x.sector_id = s.id, x.position_x_mpc = t.px, x.position_y_mpc = t.py, x.position_z_mpc = t.pz,"
                " x.location = IF(x.location IS NULL, NULL, IF(LOCATE(' -- nearest: ', x.location) > 0,"
                " CONCAT(s.name, SUBSTRING(x.location, LOCATE(' -- nearest: ', x.location))), s.name)),"
                " x.modified_at = x.modified_at WHERE x.id BETWEEN :a AND :b", {"a": first, "b": last})
        for table in _PHENOMENA + ("facilities",):
            touch = _touch(self.c, table).replace("modified_at = modified_at", "x.modified_at = x.modified_at")
            for first, last in _id_ranges(self.c, table):
                self.sql(f"UPDATE {table} x " + join.format(table=table) +
                         f" SET x.sector_id = s.id{touch} WHERE x.id BETWEEN :a AND :b", {"a": first, "b": last})
        # A facility in a system is in the system's sector.
        self.sql("UPDATE facilities f JOIN star_systems ss ON ss.id = f.star_system_id"
                 " SET f.sector_id = ss.sector_id, f.modified_at = f.modified_at"
                 " WHERE f.sector_id IS NOT NULL AND ss.sector_id IS NOT NULL")
        # Rows that follow their object.
        self.sql("UPDATE nearest_systems n JOIN star_systems ss ON n.object_table = 'star_systems'"
                 " AND n.object_id = ss.id SET n.sector_id = ss.sector_id WHERE ss.sector_id IS NOT NULL")
        for table in _PHENOMENA:
            self.sql(f"UPDATE nearest_systems n JOIN {table} x ON n.object_table = '{table}' AND n.object_id = x.id"
                     " SET n.sector_id = x.sector_id WHERE x.sector_id IS NOT NULL")
        self.sql("DELETE FROM sector_system_counts")
        self.sql("INSERT INTO sector_system_counts (sector_id, system_count)"
                 " SELECT sector_id, COUNT(*) FROM star_systems WHERE sector_id IS NOT NULL GROUP BY sector_id")
        self.sql("UPDATE object_timing_stats t JOIN old_sector o ON o.id = t.sector_id SET t.sector_id = NULL")
        self.sql("UPDATE object_timing_samples t JOIN old_sector o ON o.id = t.sector_id SET t.sector_id = NULL")

    def forget_old_sectors(self):
        """Deletes the old rows, after checking nothing that is not rebuilt data still points at one."""
        for table in ("star_systems", "facilities") + _PHENOMENA:
            left = self.sql(f"SELECT COUNT(*) FROM {table} x JOIN old_sector o ON o.id = x.sector_id").scalar()
            if left:
                raise RuntimeError(f"{left} rows of {table} are still in an old sector; nothing was changed")
        self.sql("DELETE sp FROM sector_paths sp JOIN old_sector o ON o.id = sp.sector_id")
        self.sql("DELETE n FROM nearest_systems n JOIN old_sector o ON o.id = n.sector_id")
        self.sql("DELETE s FROM sectors s JOIN old_sector o ON o.id = s.id")
        # Every name is its own base name now; an old sector's registration went with its row.
        self.sql("INSERT IGNORE INTO sector_name_registry (base_name, occurrence_count, first_sector_id)"
                 " SELECT name, 1, id FROM sectors WHERE id >= :first", {"first": self.first_new_id})

    # -- scatter and stats ------------------------------------------------------------------------------------------

    def readdress(self, table):
        self.sql("CREATE TEMPORARY TABLE IF NOT EXISTS moved (id BIGINT UNSIGNED NOT NULL PRIMARY KEY,"
                 " shell INT NOT NULL, band INT NOT NULL, wedge INT NOT NULL) ENGINE=InnoDB")
        last = 0
        while True:
            rows = self.sql(f"SELECT id, position_x_mpc, position_y_mpc, position_z_mpc FROM {table}"
                            " WHERE id > :last ORDER BY id LIMIT :n", {"last": last, "n": _READDRESS_CHUNK}).fetchall()
            if not rows:
                break
            moved = []
            for row_id, px, py, pz in rows:
                address = self.spherical.sector_address_at((px / 1000.0, py / 1000.0, pz / 1000.0), self.edge)
                moved.append({"i": row_id, "a": address[0], "b": address[1], "w": address[2]})
            self.sql("TRUNCATE TABLE moved")
            self.sql("INSERT INTO moved (id, shell, band, wedge) VALUES (:i, :a, :b, :w)", moved)
            self.sql(f"UPDATE {table} t JOIN moved m ON m.id = t.id"
                     " SET t.shell_index = m.shell, t.band_index = m.band, t.wedge_index = m.wedge")
            last = rows[-1][0]
        self.sql("DROP TEMPORARY TABLE moved")

    def link_bright_stars(self):
        self.sql("UPDATE bright_stars b JOIN star_systems ss ON ss.id = b.star_system_id"
                 " JOIN sectors s ON s.id = ss.sector_id"
                 " SET b.shell_index = s.shell_index, b.band_index = s.band_index, b.wedge_index = s.wedge_index")

    def rebuild_stats(self):
        self.sql("DELETE FROM sector_stats")
        self.sql(
            "INSERT INTO sector_stats (shell_index, band_index, wedge_index, bright_level_sol, level_before_fill_sol,"
            " actual_systems, actual_stars, mean_temperature_k, mean_luminosity_sol, mean_age_gy,"
            " total_luminosity_sol, filled_at)"
            " SELECT s.shell_index, s.band_index, s.wedge_index, 0, -1, COUNT(DISTINCT ss.id), COUNT(st.id),"
            " AVG(st.temperature_k), AVG(st.luminosity_w) / :lum, AVG(st.age_gy), SUM(st.luminosity_w) / :lum,"
            " CURRENT_TIMESTAMP(3) FROM sectors s LEFT JOIN star_systems ss ON ss.sector_id = s.id"
            " LEFT JOIN stars st ON st.star_system_id = ss.id"
            " WHERE s.id >= :first AND s.shell_index IS NOT NULL GROUP BY s.id",
            {"lum": _SOLAR_LUMINOSITY_W, "first": self.first_new_id})

    # -- object IDs -------------------------------------------------------------------------------------------------

    def renumber(self):
        """New IDs for every object (see the module docstring); returns nothing."""
        self.sql("CREATE TEMPORARY TABLE uid_map (tbl VARCHAR(24) NOT NULL, id BIGINT UNSIGNED NOT NULL,"
                 " old_uid BINARY(10) NULL, uid BINARY(10) NOT NULL, PRIMARY KEY (tbl, id)) ENGINE=InnoDB")
        generated, runtime = {}, {}
        for table in ("star_systems",) + _PHENOMENA + ("facilities",):
            last = 0
            while True:
                rows = self.sql(
                    f"SELECT x.id, x.uid, s.shell_index, s.band_index, s.wedge_index FROM {table} x"
                    " LEFT JOIN sectors s ON s.id = x.sector_id WHERE x.id > :last ORDER BY x.id LIMIT :n",
                    {"last": last, "n": _CHUNK}).fetchall()
                if not rows:
                    break
                batch = []
                for row_id, old, shell, band, wedge in rows:
                    address = _NO_SECTOR if shell is None else (shell, band, wedge)
                    kind = _RUNTIME
                    if table != "facilities" and old is not None and shell is not None:
                        kind = (int.from_bytes(bytes(old), "big") >> 12 >> _SERIAL_COUNT_BITS) & 3
                    if kind == _GENERATED:
                        serial = generated.get(address, 0)
                        generated[address] = serial + 1
                    else:
                        kind, serial = _RUNTIME, runtime.get(address, 0)
                        runtime[address] = serial + 1
                    batch.append({"t": table, "i": row_id, "o": old, "u": _pack(address, kind, serial)})
                self.sql("INSERT INTO uid_map (tbl, id, old_uid, uid) VALUES (:t, :i, :o, :u)", batch)
                last = rows[-1][0]
        for table in _BODIES:
            last = 0
            while True:
                rows = self.sql(
                    f"SELECT b.id, b.uid, m.uid FROM {table} b JOIN uid_map m ON m.tbl = 'star_systems'"
                    " AND m.id = b.star_system_id WHERE b.id > :last AND b.uid IS NOT NULL ORDER BY b.id LIMIT :n",
                    {"last": last, "n": _CHUNK}).fetchall()
                if not rows:
                    break
                batch = []
                for row_id, old, parent in rows:
                    body = int.from_bytes(bytes(old), "big") & 0xFFF
                    batch.append({"t": table, "i": row_id, "o": old,
                                  "u": ((int.from_bytes(bytes(parent), "big") & ~0xFFF) | body).to_bytes(10, "big")})
                self.sql("INSERT INTO uid_map (tbl, id, old_uid, uid) VALUES (:t, :i, :o, :u)", batch)
                last = rows[-1][0]
        for table in ("star_systems",) + _BODIES + _PHENOMENA + ("facilities",):
            touch = _touch(self.c, table)
            self.sql(f"UPDATE {table} SET uid = NULL{touch}")
            self.sql(f"UPDATE {table} x JOIN uid_map m ON m.tbl = '{table}' AND m.id = x.id"
                     f" SET x.uid = m.uid{touch.replace('modified_at = modified_at', 'x.modified_at = x.modified_at')}")
        # A stored field cloud keeps pointing at its nebula.
        self.sql("UPDATE nebula_field f JOIN uid_map m ON m.tbl = 'nebulae' AND m.old_uid = f.uid SET f.uid = m.uid"
                 " WHERE f.uid IS NOT NULL")
        # Counters: a system's follows it to its new ID, a sector's start again past the run-time IDs just given.
        self.sql("CREATE TEMPORARY TABLE kept_counter (scope VARBINARY(16) NOT NULL PRIMARY KEY,"
                 " next_value BIGINT UNSIGNED NOT NULL) ENGINE=InnoDB")
        self.sql("INSERT INTO kept_counter SELECT m.uid, c.next_value FROM id_counters c JOIN uid_map m"
                 " ON m.tbl = 'star_systems' AND m.old_uid = c.scope WHERE c.kind = 'system'")
        self.sql("DELETE FROM id_counters WHERE kind IN ('system', 'sector')")
        self.sql("INSERT INTO id_counters (kind, scope, next_value) SELECT 'system', scope, next_value FROM kept_counter")
        self.sql("INSERT INTO id_counters (kind, scope, next_value) VALUES ('sdel', '', 1)"
                 " ON DUPLICATE KEY UPDATE next_value = next_value + 1")
        counters = [{"s": _sector_scope(address), "n": count} for address, count in runtime.items()]
        for first in range(0, len(counters), _CHUNK):
            self.sql("INSERT INTO id_counters (kind, scope, next_value) VALUES ('sector', :s, :n)",
                     counters[first:first + _CHUNK])

    # -- the outline ------------------------------------------------------------------------------------------------

    def rebuild_outline(self):
        from planetgen.galaxy import skeleton
        from planetgen.galaxy.density import GalaxyShape

        row = self.sql("SELECT * FROM galaxy_shape WHERE id = 1").mappings().fetchone()
        if row is None:
            return
        shape = GalaxyShape(
            disk_scale_length_pc=row["disk_scale_length_pc"], disk_scale_height_pc=row["disk_scale_height_pc"],
            bulge_scale_radius_pc=row["bulge_scale_radius_pc"], bulge_amplitude=row["bulge_amplitude"],
            arm_count=row["arm_count"], pitch_angle_rad=row["pitch_angle_rad"], arm_amplitude=row["arm_amplitude"],
            arm_level=row["arm_level"], core_amplitude=row["core_amplitude"],
            spiral_reference_radius_pc=row["spiral_reference_radius_pc"],
            spiral_reference_angle_rad=row["spiral_reference_angle_rad"], k_norm=row["k_norm"])
        extents, outer, _confirmed = skeleton.build_shell_extents(
            shape, row["edge_pc"], 1.0 / row["expected_system_count_at_density_1"])
        self.sql("DELETE FROM galaxy_shell")
        rows = [{"s": shell, "h": high, "n": count} for shell, high, count in skeleton.shell_candidate_counts(extents)]
        for first in range(0, len(rows), _CHUNK):
            self.sql("INSERT INTO galaxy_shell (shell_index, high_band_index, cell_count) VALUES (:s, :h, :n)",
                     rows[first:first + _CHUNK])
        self.sql("UPDATE galaxy_shape SET outer_shell_index = :o WHERE id = 1", {"o": outer})

    # -- all of it --------------------------------------------------------------------------------------------------

    def run(self):
        self.step("reading the old sectors")
        self.read_old_sectors()
        self.step("placing systems")
        self.place_systems()
        self.step("placing phenomena and facilities")
        for table in _PHENOMENA:
            self.place_centered(table)
        self.place_centered("facilities", "AND t.star_system_id IS NULL")
        self.step("covering the generated sectors")
        self.cover_old_sectors()
        self.step("making the new sectors")
        self.insert_sectors()
        self.step("linking objects")
        self.relink()
        self.step("removing the old sectors")
        self.forget_old_sectors()
        self.step("bright stars")
        self.readdress("bright_stars")
        self.link_bright_stars()
        self.step("phenomenon scatter")
        self.readdress("phenomenon_scatter")
        self.step("sector stats")
        self.rebuild_stats()
        self.step("object IDs")
        self.renumber()
        self.step("outline")
        self.rebuild_outline()
        self.step("done")


def _sector_scope(address):
    shell, band, wedge = address
    return ((((shell << 12) | (band + _BAND_BIAS)) << 16) | wedge).to_bytes(5, "big")


def _pack(address, kind, serial):
    shell, band, wedge = address
    sector = (((shell << 12) | (band + _BAND_BIAS)) << 16) | wedge
    return ((((sector << 28) | (kind << _SERIAL_COUNT_BITS) | serial) << 12)).to_bytes(10, "big")


def upgrade():
    connection = op.get_bind()
    _change_tables(connection)
    if _has_table(connection, "galaxy_layer"):
        _drop_address_indexes(connection)
        edge = connection.execute(sa.text("SELECT edge_pc FROM galaxy_shape WHERE id = 1")).scalar()
        if edge is None:
            edge = connection.execute(sa.text(
                "SELECT edge_mpc / 1000 FROM sectors WHERE center_x_pc IS NOT NULL LIMIT 1")).scalar()
        if edge is not None:
            _Conversion(connection, float(edge)).run()
    for table in ("galaxy_layer", "galaxy_column", "galaxy_block_stats"):
        connection.execute(sa.text(f"DROP TABLE IF EXISTS {table}"))
    _add_address_indexes(connection)
