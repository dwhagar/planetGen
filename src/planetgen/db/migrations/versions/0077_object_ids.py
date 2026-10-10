"""Object IDs in the schema (DB.20).

Schema v77. `uid` becomes `BINARY(10)` and UNIQUE on its own on every object
table (`star_systems`, `stars`, `planets`, `moons`, `asteroid_belts`,
`comets`, the eight phenomenon tables) and is added to `facilities`; a new
`id_counters` table holds the counters of run-time births. Existing rows are
numbered by row order: a sector's systems, then its phenomena, take
generated serials 0, 1, 2 ... in id order; a system's stars, planets, moons,
belts and comets take body numbers 1, 2, 3 ... in that order. Rows with no
sector address are run-time births at the no-sector address. The layout is
`galaxy/object_uid.py`'s 80-bit default, repeated here so the revision keeps
meaning what it meant when it was written.
"""

import sqlalchemy as sa
from alembic import op

revision = "0077"
down_revision = "0076"
branch_labels = None
depends_on = None

_TOP_LEVEL = ("star_systems",)
_PHENOMENA = ("black_holes", "neutron_stars", "nebulae", "supernova_remnants", "quasars", "rogue_planets",
              "interstellar_comets", "asteroid_fields")
_BODIES = ("stars", "planets", "moons", "asteroid_belts", "comets")

# 80 bits: ring 12, biased layer 12, slot 16 | kind 2 + serial 26 | body 12
_LAYER_BIAS = 1 << 11
_SERIAL_COUNT_BITS = 26
_GENERATED, _RUNTIME = 0, 1
_NO_SECTOR = (0, 0, 0)
_CHUNK = 1000


def _has_column(connection, table, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :t AND COLUMN_NAME = :c"), {"t": table, "c": column}).scalar()


def _has_index(connection, table, index):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.STATISTICS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :t AND INDEX_NAME = :i"), {"t": table, "i": index}).scalar()


def _pack(address, kind, serial, body=0):
    ring, layer, slot = address
    sector = (((ring << 12) | (layer + _LAYER_BIAS)) << 16) | slot
    return (((sector << 28) | (kind << _SERIAL_COUNT_BITS) | serial) << 12) | body


def _bytes(value):
    return value.to_bytes(10, "big")


def _sector_scope(address):
    ring, layer, slot = address
    return ((((ring << 12) | (layer + _LAYER_BIAS)) << 16) | slot).to_bytes(5, "big")


def _write(connection, table, pairs):
    touch = "" if not _has_column(connection, table, "modified_at") else ", modified_at = modified_at"
    statement = sa.text(f"UPDATE {table} SET uid = :u{touch} WHERE id = :i")
    for first in range(0, len(pairs), _CHUNK):
        connection.execute(statement, [{"u": uid, "i": row_id} for row_id, uid in pairs[first:first + _CHUNK]])


def _reshape(connection, table):
    if _has_column(connection, table, "uid"):
        touch = "" if not _has_column(connection, table, "modified_at") else ", modified_at = modified_at"
        connection.execute(sa.text(f"UPDATE {table} SET uid = NULL{touch}"))
        if _has_index(connection, table, f"uq_{table}_uid"):
            connection.execute(sa.text(f"ALTER TABLE {table} DROP INDEX uq_{table}_uid"))
        connection.execute(sa.text(f"ALTER TABLE {table} MODIFY COLUMN uid BINARY(10) NULL"))
    else:
        connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN uid BINARY(10) NULL"))
    connection.execute(sa.text(f"ALTER TABLE {table} ADD UNIQUE KEY uq_{table}_uid (uid)"))


def _number_top_level(connection):
    """Generated serials by sector, systems first; returns `{system id: ID}`."""
    addresses = {}
    for row in connection.execute(sa.text("SELECT id, ring_index, layer_index, ring_slot_index FROM sectors")):
        if None not in (row[1], row[2], row[3]):
            addresses[row[0]] = (row[1], row[2], row[3])
    next_serial = {}
    runtime = 0
    system_uid = {}
    for table in _TOP_LEVEL + _PHENOMENA:
        pairs = []
        for row_id, sector_id in connection.execute(sa.text(f"SELECT id, sector_id FROM {table} ORDER BY id")):
            address = addresses.get(sector_id)
            if address is None:
                value = _pack(_NO_SECTOR, _RUNTIME, runtime)
                runtime += 1
            else:
                serial = next_serial.get(sector_id, 0)
                next_serial[sector_id] = serial + 1
                value = _pack(address, _GENERATED, serial)
            pairs.append((row_id, _bytes(value)))
            if table == "star_systems":
                system_uid[row_id] = value
        _write(connection, table, pairs)
    if runtime:
        connection.execute(sa.text("INSERT INTO id_counters (kind, scope, next_value) VALUES ('sector', :s, :n)"
                                   " ON DUPLICATE KEY UPDATE next_value = GREATEST(next_value, :n)"),
                           {"s": _sector_scope(_NO_SECTOR), "n": runtime})
    return system_uid


def _number_bodies(connection, system_uid):
    count = {}
    for table in _BODIES:
        pairs = []
        for row_id, system_id in connection.execute(sa.text(
                f"SELECT id, star_system_id FROM {table} ORDER BY star_system_id, id")):
            parent = system_uid.get(system_id)
            number = count.get(system_id, 0) + 1
            if parent is None or number > 4095:
                continue
            count[system_id] = number
            pairs.append((row_id, _bytes(parent | number)))
        _write(connection, table, pairs)


def upgrade():
    connection = op.get_bind()
    connection.execute(sa.text(
        "CREATE TABLE IF NOT EXISTS id_counters ("
        " kind VARCHAR(8) NOT NULL, scope VARBINARY(16) NOT NULL, next_value BIGINT UNSIGNED NOT NULL,"
        " PRIMARY KEY (kind, scope)) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"))
    for table in _TOP_LEVEL + _BODIES + _PHENOMENA + ("facilities",):
        _reshape(connection, table)
    _number_bodies(connection, _number_top_level(connection))
