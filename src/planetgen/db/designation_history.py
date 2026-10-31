# planetgen/db/designation_history.py

"""
The writers of designation history (GEN.1138; docs/design/architecture/designation-history.md): the `designations` table holds what an object has
stopped being called, and every place that changes a stored object's name or ID writes the old one here, in the same transaction as the change.

The rules (which reasons exist, what makes a new ID acceptable) are the pure module `planetgen.galaxy.designations`; this module reads and writes
the table. The current designations stay in the objects' own `name` and `uid` columns, so an object that never changes costs no row. `object_kind`
is the object's table and `object_row_id` its `id`: the row number never changes, the `uid` is the thing that does.
"""

from planetgen.db import kind_tables
from planetgen.galaxy import designations, object_uid
from planetgen.util import generator_epoch

EMBEDDED_IN_OTHERS = ("star_systems", "planets", "asteroid_belts")
"""tuple: The tables whose ID other IDs embed (a body's ID carries its system's, a belt member's or captured body's its belt's or planet's). They
are not renumbered here: the cascade of invariant 4 is part e (GEN.1141)."""


def _epoch(epoch):
    return generator_epoch.generator_epoch() if epoch is None else epoch


def _previous_epoch(conn, kind, table, row_id):
    """The epoch the object's last retired designation of `kind` stopped in, which is when its current one began; `None` for a birth designation."""
    row = conn.execute("SELECT MAX(retired_epoch) AS epoch FROM designations WHERE designation_kind = ? AND object_kind = ? AND object_row_id = ?",
                       (kind, table, row_id)).fetchone()
    return row["epoch"] if row else None


def retire_name(conn, table, row_id, old_name, new_name, reason, epoch=None):
    """
    Writes the history row of a name an object stops having: call it in the transaction of the `UPDATE` that changes the name. Nothing is written
    for a name that was empty or does not change (a change of case only is a change).

    Returns:
        bool: Whether a row was written.

    Raises:
        DesignationError: A reason that is not on the closed list.
    """
    designations.check_reason(reason)
    if not old_name or old_name == new_name:
        return False
    conn.execute(
        "INSERT INTO designations (designation_kind, name, object_kind, object_row_id, valid_from_epoch, retired_epoch, reason)"
        " VALUES ('name', ?, ?, ?, ?, ?, ?)",
        (old_name, table, row_id, _previous_epoch(conn, designations.KIND_NAME, table, row_id), _epoch(epoch), reason))
    return True


def is_live_uid(conn, uid):
    """Whether any object holds the ID `uid` (an integer) now."""
    raw = object_uid.to_bytes(uid)
    return any(conn.execute(f"SELECT 1 AS found FROM {table} WHERE uid = ? LIMIT 1", (raw,)).fetchone() is not None
               for table in kind_tables.tables("uid"))


def is_retired_uid(conn, uid):
    """Whether `uid` (an integer) was an object's ID once."""
    return conn.execute("SELECT 1 AS found FROM designations WHERE uid = ?", (object_uid.to_bytes(uid),)).fetchone() is not None


def renumber(conn, table, row_id, new_uid, reason="admin renumber", epoch=None):
    """
    Gives the object `row_id` of `table` the ID `new_uid` (an integer): the old ID is written to `designations` and reserved for ever, the column
    changes and the object's row of `object_index` follows, all in the caller's transaction. This is the only place an ID changes.

    Raises:
        DesignationError: A reason off the list, an ID that is or ever was any object's, a table with no IDs or one whose ID others embed, or no
            such object.
    """
    designations.check_reason(reason)
    if table not in kind_tables.tables("uid"):
        raise designations.DesignationError(f"{table} has no object IDs")
    if table in EMBEDDED_IN_OTHERS:
        raise designations.DesignationError(f"the ID of a row of {table} is embedded in other IDs; its bodies renumber with it")
    row = conn.execute(f"SELECT uid FROM {table} WHERE id = ?", (row_id,)).fetchone()
    if row is None:
        raise designations.DesignationError(f"no such object: {table} {row_id}")
    old_uid = object_uid.from_bytes(row["uid"])
    designations.check_new_uid(new_uid, {new_uid} if is_live_uid(conn, new_uid) else (), {new_uid} if is_retired_uid(conn, new_uid) else ())
    conn.execute(
        "INSERT INTO designations (designation_kind, uid, object_kind, object_row_id, valid_from_epoch, retired_epoch, reason)"
        " VALUES ('uid', ?, ?, ?, ?, ?, ?)",
        (object_uid.to_bytes(old_uid), table, row_id, _previous_epoch(conn, designations.KIND_UID, table, row_id), _epoch(epoch), reason))
    conn.execute(f"UPDATE {table} SET uid = ? WHERE id = ?", (object_uid.to_bytes(new_uid), row_id))
    kind = {name: kind for kind, name in kind_tables.kind_tables("uid").items()}[table]
    conn.execute("UPDATE object_index SET uid = ? WHERE kind = ? AND uid = ?",
                 (object_uid.to_bytes(new_uid), kind, object_uid.to_bytes(old_uid)))
    return new_uid


def history(conn, table, row_id):
    """The `Retired` records of one object, newest first."""
    rows = conn.execute(
        "SELECT designation_kind, uid, name, object_kind, object_row_id, valid_from_epoch, retired_epoch, reason, gone_epoch FROM designations"
        " WHERE object_kind = ? AND object_row_id = ? ORDER BY retired_epoch DESC, id DESC", (table, row_id)).fetchall()
    return [designations.Retired(
        row["designation_kind"], object_uid.from_bytes(row["uid"]) if row["designation_kind"] == designations.KIND_UID else row["name"],
        row["object_kind"], row["object_row_id"], row["valid_from_epoch"], row["retired_epoch"], row["reason"], row["gone_epoch"])
        for row in rows]
