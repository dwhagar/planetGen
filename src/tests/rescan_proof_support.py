# tests/rescan_proof_support.py

"""
The harness for the proof tests of the generation-rescan chain (GEN.1476 running block totals, GEN.1477 the shared neighbourhood read, GEN.1479 derived stamps; Boss,
8:50 PM EST on 2026-11-08, "all lanes work the rescan chain"; GEN.1482 the link stage). Each item replaces a pass that read the whole galaxy by something that keeps or shares what the pass would
have found, so each is proved the same way: the new way gives, to the last bit, what the old way gave. Three pieces, none of them built for one item:

* `fingerprint` / `assert_same_rows`: a table's rows as bytes (a float as its eight bytes, never as text), in a stated order, and the first difference named. Exact, so a
  sum that differs in the last place fails.
* `block_light_rows` / `assert_totals_equal_a_rebuild`: the stored block lights as they stand, then `block_lights.rebuild` (the check the chain keeps), then both compared.
  GEN.1476's running totals are proved by calling it after each way they change: the scatter, a sector delete, a regenerate.
* `TouchLog`: row triggers on chosen tables that record every row an INSERT, UPDATE or DELETE statement processes (a statement that writes the value a column already holds
  fires them too), so "touches zero rows" is a count of rows the job processed, not a guess from the values.
* `link_rows`: the nearest-system lists and the columns the link stage fills (GEN.1482), as bytes in an order that does not depend on row ids.
* `settled_rows`: the paths `settle_sectors` saves for some sectors as bytes, in an order that does not depend on row ids, for the settle chain (GEN.1477): the paths of two
  ways of reading the neighbourhood compare with `assert_same_rows`.

`small_galaxy` builds what they run on: three neighbouring sectors, each with systems and a black hole, with pre-placed bright stars and scattered phenomena in the blocks
that hold them. It reuses the builders the existing tests use (`test_sector_paths_db._save`, `test_block_lights._store_objects`), and seeds every system draw, so a run is the
same every time.
"""

import contextlib
import struct

from planetgen.db import store
from planetgen.db import block_lights
from tests.fuzz_support import deterministic_entropy

BLOCK_LIGHT_ORDER = "m, shell_index, band_index, wedge_index"

SMALL_GALAXY_ADDRESSES = ((1500, 2, 700), (1500, 2, 701), (1500, 3, 700))
"""tuple: The shell, band and wedge of the three sectors of the small galaxy: neighbours, so settling one covers the others and their blocks share parents."""

SMALL_GALAXY_SEED = 20261109


def bits(value):
    """A value as something that compares exactly: a float as its eight bytes (so -0.0 differs from 0.0 and a NaN from a number), the rest as it is."""
    if isinstance(value, float):
        return struct.pack("<d", value)
    if isinstance(value, bytearray):
        return bytes(value)
    return value


def fingerprint(conn, sql, params=()):
    """The rows of `sql`, each as a tuple of `bits` in the column order of the query; `sql` states its own ORDER BY, so two fingerprints of the same query compare row by row."""
    return [tuple(bits(value) for value in row.values()) for row in conn.execute(sql, params).fetchall()]


def first_difference(left, right, columns=None):
    """A sentence naming the first row that differs between two fingerprints (or the first one missing from the shorter), `None` when they are the same."""
    for index, (a, b) in enumerate(zip(left, right)):
        if a != b:
            where = [columns[i] if columns else i for i, (x, y) in enumerate(zip(a, b)) if x != y]
            return f"row {index} differs in {where}: {a!r} against {b!r}"
    if len(left) != len(right):
        longer, side = (left, "left") if len(left) > len(right) else (right, "right")
        return f"{len(left)} rows against {len(right)}: the first one only the {side} has is {longer[min(len(left), len(right))]!r}"
    return None


def assert_same_rows(left, right, what, columns=None):
    """Fails naming `what` and the first row that differs, when the two fingerprints are not the same bytes."""
    difference = first_difference(left, right, columns)
    assert difference is None, f"{what}: {difference}"


def block_light_rows(conn):
    """The stored `galaxy_block_lights` as they stand, every column, in the table's key order."""
    return fingerprint(conn, f"SELECT * FROM galaxy_block_lights ORDER BY {BLOCK_LIGHT_ORDER}")


def assert_totals_equal_a_rebuild(config, what):
    """
    The block lights as they stand against the block lights `block_lights.rebuild` makes from the objects (the check the running totals keep, GEN.1476): the same rows with the
    same bytes. The rebuild replaces the table, so the database ends holding the rebuilt rows; call it last in a step, or again after the next change. Returns the number of rows.
    """
    conn = store.get_connection(config)
    try:
        running = block_light_rows(conn)
        block_lights.rebuild(conn)
        conn.commit()
        rebuilt = block_light_rows(conn)
    finally:
        conn.close()
    assert_same_rows(running, rebuilt, f"{what}: the stored block lights against a rebuild")
    return len(rebuilt)


class TouchLog:
    """
    Records every row that an INSERT, UPDATE or DELETE processes on `tables` while it is open: `touched()` is a list of `(table, operation, id)`. It is made of triggers in the
    scratch database and table `touch_log`, dropped on exit. A statement that sets a column to the value it holds is processed too, which is the point: a job that "touches nothing"
    must not write at all. The tables need an `id` column.
    """

    def __init__(self, config, tables):
        self.config = config
        self.tables = tuple(tables)

    def __enter__(self):
        conn = store.get_connection(self.config)
        try:
            conn.execute("DROP TABLE IF EXISTS touch_log")
            conn.execute("CREATE TABLE touch_log (n BIGINT AUTO_INCREMENT PRIMARY KEY, table_name VARCHAR(64) NOT NULL, operation CHAR(1) NOT NULL, row_id BIGINT NOT NULL)")
            for table in self.tables:
                for operation, event, row in (("I", "INSERT", "NEW"), ("U", "UPDATE", "NEW"), ("D", "DELETE", "OLD")):
                    conn.execute(f"DROP TRIGGER IF EXISTS touch_{table}_{operation}")
                    conn.execute(f"CREATE TRIGGER touch_{table}_{operation} AFTER {event} ON {table} FOR EACH ROW "
                                 f"INSERT INTO touch_log (table_name, operation, row_id) VALUES ('{table}', '{operation}', {row}.id)")
            conn.commit()
        finally:
            conn.close()
        return self

    def __exit__(self, *exc):
        conn = store.get_connection(self.config)
        try:
            for table in self.tables:
                for operation in "IUD":
                    conn.execute(f"DROP TRIGGER IF EXISTS touch_{table}_{operation}")
            conn.execute("DROP TABLE IF EXISTS touch_log")
            conn.commit()
        finally:
            conn.close()
        return False

    def touched(self):
        conn = store.get_connection(self.config)
        try:
            return [(row["table_name"], row["operation"], row["row_id"]) for row in conn.execute("SELECT table_name, operation, row_id FROM touch_log ORDER BY n").fetchall()]
        finally:
            conn.close()

    def clear(self):
        conn = store.get_connection(self.config)
        try:
            conn.execute("DELETE FROM touch_log")
            conn.commit()
        finally:
            conn.close()


def settled_rows(conn, sector_ids):
    """
    The saved paths of the systems in `sector_ids`, every knot, as bytes, ordered by what the path belongs to and not by row ids (a path settled again gets a new id). Each row is
    `(sector_id, object_table, object_id, exited, duration_years, position, t_years, x, y, z, vx, vy, vz)`.
    """
    marks = ", ".join("?" for _ in sector_ids)
    return fingerprint(
        conn,
        "SELECT p.sector_id, p.object_table, p.object_id, p.exited, p.duration_years, k.position, k.t_years, k.x_pc, k.y_pc, k.z_pc, k.vx_kms, k.vy_kms, k.vz_kms "
        f"FROM sector_paths p JOIN sector_path_knots k ON k.path_id = p.id WHERE p.sector_id IN ({marks}) "
        "ORDER BY p.sector_id, p.object_table, p.object_id, k.position", tuple(sector_ids))


def link_rows(conn):
    """
    What the link stage (GEN.1482) writes, as bytes and in an order that does not depend on row ids: the nearest-system lists (`nearest_systems`: sector, object, rank, neighbour and
    distance) and the system columns it fills (`quadrant`, `inside_nebula_id`, `inside_remnant_id`, by the system's unique id). The phenomenon tables' `inside_*` and `quadrant`
    columns are Lane 14's to add here with the item, for they know which tables the stage fills. Compare the same sectors linked in batches of 2, 20 and 100, and a galaxy grown by a
    second run, with `assert_same_rows` (the pattern of `test_db_persistence.test_linking_in_small_batches_gives_the_same_lists_as_one_batch`).
    """
    return {
        "nearest_systems": fingerprint(
            conn, "SELECT sector_id, object_table, object_id, neighbor_rank, neighbor_system_id, distance_pc FROM nearest_systems "
                  "ORDER BY sector_id, object_table, object_id, neighbor_rank"),
        "star_systems": fingerprint(conn, "SELECT uid, quadrant, inside_nebula_id, inside_remnant_id FROM star_systems ORDER BY uid"),
    }


def small_galaxy(config, seed=SMALL_GALAXY_SEED):
    """
    Builds the small galaxy in `config`'s empty database and returns `{"sector_ids": [...], "addresses": [...]}`: three neighbouring sectors of two seeded G2V systems and a black
    hole each (`SMALL_GALAXY_ADDRESSES`), with a pre-placed bright star and a scattered phenomenon in the block of each, and the block lights made by a rebuild. The same draws every time.
    """
    from tests.test_block_lights import _store_objects
    from tests.test_sector_paths_db import _save

    sector_ids = []
    for index, address in enumerate(SMALL_GALAXY_ADDRESSES):
        with deterministic_entropy(seed + index):
            sector_id, _center = _save(config, systems=2, address=address)
        sector_ids.append(sector_id)
    conn = store.get_connection(config)
    try:
        stars = [(address, 2.0e6 + 1.0e5 * i, -1.0e6, 3.0e5, 1.0e3 * (i + 1), 8000.0 + 500.0 * i) for i, address in enumerate(SMALL_GALAXY_ADDRESSES)]
        points = [(address, -2.0e6, 1.0e6 + 1.0e5 * i, -3.0e5, float(10 * (i + 1)), False) for i, address in enumerate(SMALL_GALAXY_ADDRESSES)]
        _store_objects(conn, stars, points)
        block_lights.rebuild(conn)
        conn.commit()
    finally:
        conn.close()
    return {"sector_ids": sector_ids, "addresses": list(SMALL_GALAXY_ADDRESSES)}


@contextlib.contextmanager
def connection(config):
    """A connection that is closed on exit."""
    conn = store.get_connection(config)
    try:
        yield conn
    finally:
        conn.close()
