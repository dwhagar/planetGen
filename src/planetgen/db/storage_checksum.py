# planetgen/db/storage_checksum.py

"""
The storage checksum of a sector and the flag that says it is out of date (DB.45; docs/design/operations/db-check-and-parity-repair.md
sections 5.2 and 5.3).

**The checksum** is the SHA-256 of a lossless, id-preserving dump of everything a sector holds as the database stores it: its `sectors`
row, every object filed in it (star systems, their stars, planets, moons, belts and comets, the phenomena filed by `sector_id`,
the facilities) and those objects' child rows (gases, interior layers, reflection spectra, compositions, nebula balls, orbit nodes). Every
column of every row is in it, `id`, foreign keys and `uid` included, with a value written exactly: an integer as its digits, a float as
its 64 bits (so a flipped low bit of a DOUBLE changes the checksum, which GEN.58's nine-digit rounding would hide), bytes in hex, a
string as its UTF-8 bytes. The rows come in table order and then in primary key order, so the order they were saved in does not count.

Left out (`VOLATILE_COLUMNS`): what the daily orbit update rewrites (`store.advance_orbital_phases`, `advance_galactic_positions`,
`advance_facility_orbits`, `advance_comet_orbits`: phases, positions, velocities, reflex offsets, binary positions, comet anomalies, the
`location` text and each row's clocks), so a tick of the simulation never makes a sector look damaged; and the three columns that cannot
describe themselves (`storage_checksum`, `parity_dirty`, `leaf_digest`), and `modified_at`, which any write bumps. Not in it either: the
address-keyed rows (`bright_stars`, `phenomenon_scatter`), which exist without a sector, and the rebuilt tables (stats, paths, nearest
systems, the object index).

**The flag** is `sectors.parity_dirty` (DB.44): 1 means the stored checksum, and the parity record made from the same export (DB.9), may
no longer describe the sector. Every writer of a column the checksum covers sets it in its own transaction (`mark_sector`,
`mark_system`, `mark_systems`, `mark_star`, `mark_rows`), without touching `modified_at`. A sector that exists before the flag, and a
new one, is 1. `seal` computes the checksum from the rows as they read back, stores it and clears the flag: the parity writer calls it
once the parity record is current. `check_sector` tells the three states apart, which is what DB.8's `storage` check reports: `ok`, `unprotected`
(the flag is set, or no checksum is stored yet: not damage) and `damaged` (the flag is clear and the rows no longer hash to the checksum).

`tests/test_storage_checksum.py` pins the orbit update's columns against the real update and lists every function that writes a covered
table, so a new writer cannot be added without deciding how it marks the sector.
"""

import hashlib
import math
import struct
from datetime import date, datetime, time, timedelta
from decimal import Decimal

from planetgen.db import fingerprint

SECTOR_TABLES = fingerprint.SECTOR_TABLES
SYSTEM_TABLES = fingerprint.SYSTEM_TABLES
STAR_TABLES = fingerprint.STAR_TABLES
CHILD_TABLES = fingerprint.CHILD_TABLES

CONTENT_TABLES = ("sectors",) + tuple(dict.fromkeys(SECTOR_TABLES + SYSTEM_TABLES + STAR_TABLES)) + tuple(
    table for table, _column, _parent in CHILD_TABLES)
"""tuple: Every table a sector's checksum reads (the fingerprint's inventory without the address-keyed and plan tables)."""

_POSITION = ("position_x_km", "position_y_km", "position_z_km", "velocity_x_kms", "velocity_y_kms", "velocity_z_kms")
_GALACTIC = ("center_x_pc", "center_y_pc", "center_z_pc", "galactic_radius_pc", "quadrant", "velocity_x_kms", "velocity_y_kms",
             "velocity_z_kms", "galactic_orbital_phase_deg")
_BINARY = ("binary_galactic_orbital_phase_deg", "binary_mutual_orbital_phase_deg", "binary_mutual_orbital_speed_kms",
           "binary_mutual_position_x_km", "binary_mutual_position_y_km", "binary_mutual_position_z_km",
           "binary_primary_position_x_km", "binary_primary_position_y_km", "binary_primary_position_z_km",
           "binary_secondary_position_x_km", "binary_secondary_position_y_km", "binary_secondary_position_z_km",
           "binary_planetary_wobble_x_km", "binary_planetary_wobble_y_km", "binary_planetary_wobble_z_km")
_CLOCKS = ("epoch_unix", "next_update_due", "binary_epoch_unix", "binary_next_update_due")

VOLATILE_COLUMNS = {
    "sectors": ("storage_checksum", "parity_dirty", "leaf_digest", "modified_at"),
    "star_systems": ("sector_id", "position_x_mpc", "position_y_mpc", "position_z_mpc", "galaxy_x_pc", "galaxy_y_pc", "galaxy_z_pc", "quadrant", "location", "velocity_x_kms",
                     "velocity_y_kms", "velocity_z_kms", "modified_at") + _BINARY + _CLOCKS,
    "stars": ("galactic_orbital_phase_deg", "reflex_offset_x_km", "reflex_offset_y_km", "reflex_offset_z_km") + _CLOCKS,
    "planets": _POSITION + ("distance_km", "orbital_speed_kms", "reflex_offset_x_km", "reflex_offset_y_km", "reflex_offset_z_km") + _CLOCKS,
    "moons": _POSITION + ("distance_km", "orbital_speed_kms") + _CLOCKS,
    "asteroid_belts": _POSITION + _CLOCKS,
    "comets": _POSITION + ("mean_anomaly_deg", "parabolic_mean_anomaly", "distance_km", "orbital_speed_kms") + _CLOCKS,
    "facilities": ("sector_id", "orbit_phase_deg", "orbit_distance_km", "orbital_speed_kms", "modified_at") + _GALACTIC + _CLOCKS,
    **{table: ("sector_id", "modified_at") + _GALACTIC + _CLOCKS for table in (
        "black_holes", "neutron_stars", "nebulae", "supernova_remnants", "quasars", "rogue_planets", "interstellar_comets",
        "asteroid_fields")},
}
"""dict: Per table, the columns the checksum leaves out: what the orbit update rewrites (the design's section 5.2 table, found by
`tests/test_storage_checksum.py` against the real update), a refiled object's `sector_id`, and the columns that cannot describe themselves.
A column named here that a table does not have is ignored."""

_IN_BATCH = 500


class StorageChecksumError(ValueError):
    """A sector that has no row."""


def _encode(value):
    """One stored value as bytes that tell it from every other value of every type (a tag, then its exact form)."""
    if value is None:
        return b"N"
    if isinstance(value, bool):
        value = int(value)
    if isinstance(value, int):
        return b"I" + str(value).encode("ascii")
    if isinstance(value, float):
        return b"F" + struct.pack("<d", value).hex().encode("ascii")
    if isinstance(value, Decimal):
        return b"D" + str(value).encode("ascii")
    if isinstance(value, (bytes, bytearray, memoryview)):
        return b"B" + bytes(value).hex().encode("ascii")
    if isinstance(value, str):
        return b"S" + value.encode("utf-8")
    if isinstance(value, datetime):
        return b"T" + value.isoformat(timespec="microseconds").encode("ascii")
    if isinstance(value, date):
        return b"A" + value.isoformat().encode("ascii")
    if isinstance(value, time):
        return b"M" + value.isoformat(timespec="microseconds").encode("ascii")
    if isinstance(value, timedelta):
        return b"L" + str(value.total_seconds()).encode("ascii")
    raise TypeError(f"the storage checksum has no encoding for a {type(value).__name__}")


def _feed(digest, parts):
    """Adds each part of a record, with its length, so no two records can run together into the same bytes."""
    for part in parts:
        digest.update(struct.pack("<I", len(part)))
        digest.update(part)


class _Reader:
    """Reads one sector's rows, remembering each table's covered columns and primary key."""

    def __init__(self, conn):
        self.conn = conn
        self._columns = {}

    def columns(self, table):
        """`(covered columns, order-by columns)` of `table`: every column but the volatile ones, in table order; the primary key
        (else every covered column) to sort by."""
        if table not in self._columns:
            rows = self.conn.execute(
                "SELECT column_name AS name, column_key AS k, extra AS x FROM information_schema.columns"
                " WHERE table_schema = DATABASE() AND table_name = ? ORDER BY ordinal_position", (table,)).fetchall()
            left_out = set(VOLATILE_COLUMNS.get(table, ()))
            # A generated column (`cell_key`) is worked out from the others, so it adds nothing a flip could hide in.
            covered = [row["name"] for row in rows if row["name"] not in left_out and "GENERATED" not in (row["x"] or "").upper()]
            key = [row["name"] for row in rows if row["k"] == "PRI"]
            self._columns[table] = (covered, [c for c in key if c in covered] or covered)
        return self._columns[table]

    def select(self, table, where, params=()):
        covered, order = self.columns(table)
        names = ", ".join(f"`{c}`" for c in covered)
        sort = ", ".join(f"`{c}`" for c in order)
        return self.conn.execute(f"SELECT {names} FROM `{table}` WHERE {where} ORDER BY {sort}", tuple(params)).fetchall()

    def select_in(self, table, column, ids):
        ids = sorted(set(ids))
        rows = []
        for first in range(0, len(ids), _IN_BATCH):
            chunk = ids[first:first + _IN_BATCH]
            rows.extend(self.select(table, f"`{column}` IN ({', '.join('?' * len(chunk))})", chunk))
        return rows


def _ids(rows):
    return [row["id"] for row in rows if "id" in row.keys()]


def sector_rows(conn, sector_id):
    """
    The rows a sector's checksum covers: `[(table, covered columns, [row, ...])]` in `CONTENT_TABLES` order, each table's rows in primary key order.

    Raises:
        StorageChecksumError: The sector has no row.
    """
    reader = _Reader(conn)
    found = {"sectors": reader.select("sectors", "id = ?", (sector_id,))}
    if not found["sectors"]:
        raise StorageChecksumError(f"sector {sector_id} has no row")
    seen = {}

    def keep(table, rows):
        have = seen.setdefault(table, {})
        for row in rows:
            have[tuple(row[c] for c in reader.columns(table)[1])] = row
        found[table] = [have[key] for key in sorted(have, key=repr)]

    for table in SECTOR_TABLES:
        keep(table, reader.select(table, "sector_id = ?", (sector_id,)))
    systems = _ids(found.get("star_systems", []))
    for table in SYSTEM_TABLES:
        keep(table, reader.select_in(table, "star_system_id", systems))
    stars = _ids(found.get("stars", []))
    for table in STAR_TABLES:
        keep(table, reader.select_in(table, "star_id", stars))
    # DB.73: a moon of a rogue planet has no system or planet; it is filed with the rogue planet it orbits.
    keep("moons", reader.select_in("moons", "rogue_planet_id", _ids(found.get("rogue_planets", []))))
    for table, column, parent in CHILD_TABLES:
        keep(table, reader.select_in(table, column, _ids(found.get(parent, []))))
    return [(table, reader.columns(table)[0], found.get(table, [])) for table in CONTENT_TABLES]


def checksum(conn, sector_id):
    """The storage checksum of a sector, 32 bytes, from the rows as they are stored now."""
    digest = hashlib.sha256()
    for table, covered, rows in sector_rows(conn, sector_id):
        _feed(digest, [b"table", table.encode("ascii"), str(len(rows)).encode("ascii")])
        for row in rows:
            _feed(digest, [_encode(row[column]) for column in covered])
    return digest.digest()


def seal(conn, sector_id):
    """
    Stores the checksum of a sector, worked out from the rows as they read back, and clears its `parity_dirty` flag (the parity
    record has been made from the same rows). Returns the checksum. Call it in the transaction that wrote the rows, after the write.
    """
    value = checksum(conn, sector_id)
    conn.execute("UPDATE sectors SET storage_checksum = ?, parity_dirty = 0, modified_at = modified_at WHERE id = ?",
                 (value, sector_id))
    return value


OK = "ok"
UNPROTECTED = "unprotected"
DAMAGED = "damaged"


def check_sector(conn, sector_id):
    """
    `ok`, `unprotected` or `damaged` for one sector: a sector whose flag is set, or that has no checksum yet, is `unprotected` (a writer
    changed it since the last seal, or it was never sealed: not damage); one whose flag is clear is `ok` when the rows still hash to the
    stored checksum and `damaged` when they do not.
    """
    row = conn.execute("SELECT storage_checksum AS c, parity_dirty AS d FROM sectors WHERE id = ?", (sector_id,)).fetchone()
    if row is None:
        raise StorageChecksumError(f"sector {sector_id} has no row")
    if row["d"] or row["c"] is None:
        return UNPROTECTED
    return OK if bytes(row["c"]) == checksum(conn, sector_id) else DAMAGED


MARK_SECTOR_SQL = "UPDATE sectors SET parity_dirty = 1, modified_at = modified_at WHERE parity_dirty = 0 AND id = ?"
"""str: Marks one sector (`?` its id), leaving `modified_at` as it was."""

MARK_SYSTEM_SQL = ("UPDATE sectors SET parity_dirty = 1, modified_at = modified_at WHERE parity_dirty = 0"
                   " AND id = (SELECT sector_id FROM star_systems WHERE id = ?)")
"""str: Marks the sector of one star system (`?` its id)."""


def _mark(conn, where, params):
    conn.execute(f"UPDATE sectors SET parity_dirty = 1, modified_at = modified_at WHERE parity_dirty = 0 AND ({where})", tuple(params))


def mark_sector(conn, sector_id):
    """Marks a sector's storage checksum out of date. A `None` sector (an object not filed in one) is nothing to mark."""
    if sector_id is not None:
        _mark(conn, "id = ?", (sector_id,))


def mark_sectors(conn, sector_ids):
    """`mark_sector` for each of `sector_ids` (a refile marks the sector an object left and the one it entered)."""
    ids = sorted({int(i) for i in sector_ids if i is not None})
    for first in range(0, len(ids), _IN_BATCH):
        chunk = ids[first:first + _IN_BATCH]
        _mark(conn, f"id IN ({', '.join('?' * len(chunk))})", chunk)


def mark_system(conn, star_system_id):
    """Marks the sector that holds a star system, for a change to the system or to anything under it (its stars, planets, moons,
    belts, comets, facilities and their child rows)."""
    if star_system_id is not None:
        _mark(conn, "id = (SELECT sector_id FROM star_systems WHERE id = ?)", (star_system_id,))


def mark_systems(conn, star_system_ids):
    """`mark_system` for each of `star_system_ids`."""
    ids = sorted({int(i) for i in star_system_ids if i is not None})
    for first in range(0, len(ids), _IN_BATCH):
        chunk = ids[first:first + _IN_BATCH]
        _mark(conn, f"id IN (SELECT sector_id FROM star_systems WHERE id IN ({', '.join('?' * len(chunk))}))", chunk)


def mark_star(conn, star_id):
    """Marks the sector that holds a star (through its system, else through its own `sector_id` for a compact remnant filed there)."""
    if star_id is not None:
        _mark(conn, "id = (SELECT s.sector_id FROM stars t JOIN star_systems s ON s.id = t.star_system_id WHERE t.id = ?)", (star_id,))


def mark_rows(conn, table, row_ids):
    """Marks the sectors of rows of `table` that carry a `sector_id` (a phenomenon, a facility) or a `star_system_id`."""
    ids = sorted({int(i) for i in row_ids if i is not None})
    if table not in CONTENT_TABLES:
        raise ValueError(f"{table} is not a table a sector's checksum covers")
    for first in range(0, len(ids), _IN_BATCH):
        chunk = ids[first:first + _IN_BATCH]
        marks = ", ".join("?" * len(chunk))
        if table in SECTOR_TABLES:
            _mark(conn, f"id IN (SELECT sector_id FROM `{table}` WHERE id IN ({marks}))", chunk)
        if table in SYSTEM_TABLES:
            _mark(conn, f"id IN (SELECT sector_id FROM star_systems WHERE id IN"
                        f" (SELECT star_system_id FROM `{table}` WHERE id IN ({marks})))", chunk)
