"""
planetgen/db/designation_rows.py

The database side of designation history (GEN.1138; docs/design/architecture/designation-history.md; the rules are `galaxy/designations.py`).
The `designations` table holds what an object stopped being called: this module is the one place that writes it. A name changes through
`set_name` (the old name is kept, then the row changes), an ID through `renumber` (the old ID is kept and reserved for ever, then the row, its
`object_index` row and the live pointers at it change). Nothing else in the package updates a `name` column of a named table or a `uid`.

`object_kind` of a row is the table that held the object and `object_row_id` its `id`, which never changes. A retired name is not reserved (the
collision bump hands a stored object's name to another); a retired ID is, by the unique index on `designation_uid` and the objects' own unique
`uid` columns.

What `renumber` leaves to the cascade of GEN.1141: the bodies of a system whose ID changes (each is renumbered by its own call), the pointers with
no index (`object_index.parent_uid`, which the next index build rewrites, `orbital_events.containing_uid`, `sectors.max_mass_id`), and the
tombstone of a deleted object. `edit_journal` and `change_plan_ops` name the ID an edit was made under and keep it; a lookup of that ID
redirects.
"""
from planetgen.db import kind_tables, storage_checksum
from planetgen.galaxy import designations, object_uid
from planetgen.util import generator_epoch

UID_TABLES = kind_tables.tables("uid")
"""tuple: The tables whose objects have an ID (`uid`)."""

NAMED_TABLES = ("sectors",) + tuple(table for table in UID_TABLES if table != "asteroid_belts")
"""tuple: The tables whose objects have a name that is kept when it changes. A sector has a name and no stored ID (its ID is its address); a
belt has an ID and no name (it prints as its ID, and its members as `<ID>.<rank>`)."""

_LIVE_POINTERS = (("constellations", "system_uid"), ("constellations", "name_id"), ("object_locks", "object_uid"))
"""tuple: `(table, column)` of the live pointers at an ID that have an index, moved with the object. Not the journals, which keep the ID they were
written under, and not the three unindexed ones (module note)."""


def _check_table(table, allowed, what):
    if table not in allowed:
        raise ValueError(f"not a table with {what}: {table!r}")


def _valid_from(conn, designation_kind, table, row_id):
    """The epoch the current designation of this kind began in: when the last one was retired, `None` if it is the birth designation."""
    row = conn.execute(
        "SELECT MAX(retired_epoch) AS epoch FROM designations WHERE designation_kind = ? AND object_kind = ? AND object_row_id = ?",
        (designation_kind, table, row_id)).fetchone()
    return row["epoch"]


def retire_name(conn, table, row_id, old_name, reason):
    """Keeps `old_name`, the name object `row_id` of `table` stops having now (a row of `designations`; the name itself is changed by the caller)."""
    _check_table(table, NAMED_TABLES, "a kept name")
    designations.check_reason(reason)
    conn.execute(
        "INSERT INTO designations (designation_kind, designation_name, object_kind, object_row_id, valid_from_epoch, retired_epoch, reason)"
        " VALUES (?, ?, ?, ?, ?, ?, ?)",
        (designations.KIND_NAME, old_name, table, row_id, _valid_from(conn, designations.KIND_NAME, table, row_id),
         generator_epoch.generator_epoch(), reason))


def set_name(conn, table, row_id, new_name, reason, old_name=None, extra_set=None):
    """
    Gives one object a new name and keeps the old one: the writer of every `name` of a named table.

    Args:
        conn (Connection): Open; the history row and the change are one transaction of the caller's.
        table (str): One of `NAMED_TABLES`.
        row_id (int): The object's `id`.
        new_name (str): The name it has from now on.
        reason (str or None): One of `designations.REASONS`; `None` when the old name was a placeholder nobody was ever called by (a body built
            again under a temporary name), which keeps nothing.
        old_name (str, optional): The name it has, when the caller already read it; read here otherwise.
        extra_set (str, optional): More assignments for the same UPDATE (`"modified_at = CURRENT_TIMESTAMP(3)"`).

    Returns:
        bool: `False` if no such object exists. The old name is kept only when it differs from `new_name` in its characters (a change of case
            is a change); an unchanged name still gets the UPDATE, as the callers' other columns need.
    """
    _check_table(table, NAMED_TABLES, "a kept name")
    if reason is not None:
        designations.check_reason(reason)
    if old_name is None:
        row = conn.execute(f"SELECT name FROM {table} WHERE id = ?", (row_id,)).fetchone()
        if row is None:
            return False
        old_name = row["name"]
    if reason is not None and old_name is not None and old_name != new_name:
        retire_name(conn, table, row_id, old_name, reason)
    assignments = "name = ?" + (f", {extra_set}" if extra_set else "")
    conn.execute(f"UPDATE {table} SET {assignments} WHERE id = ?", (new_name, row_id))
    _mark(conn, table, row_id)
    return True


def _mark(conn, table, row_id):
    """The name and the ID are columns the sector's storage checksum covers: the sector of the row learns of the change (DB.44)."""
    if table == "sectors":
        storage_checksum.mark_sector(conn, row_id)
    elif table in storage_checksum.CONTENT_TABLES:
        storage_checksum.mark_rows(conn, table, [row_id])


def _retired(row):
    return designations.Retired(
        row["designation_kind"],
        object_uid.from_bytes(row["designation_uid"]) if row["designation_kind"] == designations.KIND_UID else row["designation_name"],
        row["object_kind"], row["object_row_id"], row["valid_from_epoch"], row["retired_epoch"], row["reason"], row["gone_epoch"],
        row["retired_at"])


_COLUMNS = ("designation_kind, designation_uid, designation_name, object_kind, object_row_id, valid_from_epoch, retired_epoch, reason,"
            " gone_epoch, retired_at")


def history(conn, table, row_id):
    """The designations object `row_id` of `table` used to have, newest first, as `designations.Retired` records (IDs as integers)."""
    rows = conn.execute(
        f"SELECT {_COLUMNS} FROM designations WHERE object_kind = ? AND object_row_id = ? ORDER BY retired_epoch DESC, id DESC",
        (table, row_id)).fetchall()
    return [_retired(row) for row in rows]


def former_designations(conn, table, row_id):
    """The `former_designations` field of object `row_id` of `table`: `designations.former_json` of its history, newest first, `[]` for an object that never changed."""
    return designations.former_json(history(conn, table, row_id))


def retired_uid(conn, uid):
    """The `Retired` record of ID `uid` (an integer), or `None` when no object ever stopped having it."""
    row = conn.execute(f"SELECT {_COLUMNS} FROM designations WHERE designation_uid = ?", (object_uid.to_bytes(uid),)).fetchone()
    return None if row is None else _retired(row)


def retired_names(conn, name):
    """The `Retired` records of `name`, newest first: every object that stopped having it."""
    rows = conn.execute(
        f"SELECT {_COLUMNS} FROM designations WHERE designation_kind = ? AND designation_name = ? ORDER BY retired_epoch DESC, id DESC",
        (designations.KIND_NAME, name)).fetchall()
    return [_retired(row) for row in rows]


def live_holder(conn, uid):
    """`(table, row id)` of the object that has ID `uid` now, or `None`."""
    raw = object_uid.to_bytes(uid)
    for table in UID_TABLES:
        row = conn.execute(f"SELECT id FROM {table} WHERE uid = ? LIMIT 1", (raw,)).fetchone()
        if row is not None:
            return table, row["id"]
    return None


def renumber(conn, table, row_id, new_uid, reason):
    """
    Gives one object a new ID. The old ID is kept and reserved for ever; the object's row and its `object_index` row take the new one, and so do
    the live pointers at it that have an index (`_LIVE_POINTERS`). The bodies of a system are not touched: GEN.1141's cascade renumbers each with
    its own call.

    Args:
        conn (Connection): Open; one transaction of the caller's.
        table (str): One of `UID_TABLES`.
        row_id (int): The object's `id`.
        new_uid (int): The ID it has from now on.
        reason (str): One of `designations.REASONS`.

    Returns:
        int: The ID the object had.

    Raises:
        designations.DesignationError: The ID is any object's, or was, or the object does not exist.
    """
    _check_table(table, UID_TABLES, "an ID")
    designations.check_reason(reason)
    row = conn.execute(f"SELECT uid FROM {table} WHERE id = ?", (row_id,)).fetchone()
    if row is None:
        raise designations.DesignationError(f"no object with row {row_id} in {table}")
    if row["uid"] is None:
        raise designations.DesignationError(f"row {row_id} of {table} has no ID yet, so it has none to give up")
    old_uid = object_uid.from_bytes(row["uid"])
    designations.check_new_uid(
        new_uid, {new_uid} if live_holder(conn, new_uid) is not None else set(),
        {new_uid} if retired_uid(conn, new_uid) is not None else set())
    old_raw, new_raw = object_uid.to_bytes(old_uid), object_uid.to_bytes(new_uid)
    conn.execute(
        "INSERT INTO designations (designation_kind, designation_uid, object_kind, object_row_id, valid_from_epoch, retired_epoch, reason)"
        " VALUES (?, ?, ?, ?, ?, ?, ?)",
        (designations.KIND_UID, old_raw, table, row_id, _valid_from(conn, designations.KIND_UID, table, row_id),
         generator_epoch.generator_epoch(), reason))
    conn.execute(f"UPDATE {table} SET uid = ? WHERE id = ?", (new_raw, row_id))
    index_kind = next(kind for kind, owner in kind_tables.kind_tables("uid").items() if owner == table)
    conn.execute("UPDATE object_index SET uid = ? WHERE kind = ? AND uid = ?", (new_raw, index_kind, old_raw))
    for pointer_table, column in _LIVE_POINTERS:
        conn.execute(f"UPDATE {pointer_table} SET {column} = ? WHERE {column} = ?", (new_raw, old_raw))
    _mark(conn, table, row_id)
    return old_uid
