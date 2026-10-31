# planetgen/db/designation_history.py

"""
The writers of designation history (GEN.1138; docs/design/architecture/designation-history.md): the `designations` table holds what an object has
stopped being called, and every place that changes a stored object's name or ID writes the old one here, in the same transaction as the change.

The rules (which reasons exist, what makes a new ID acceptable) are the pure module `planetgen.galaxy.designations`; this module reads and writes
the table. The current designations stay in the objects' own `name` and `uid` columns, so an object that never changes costs no row. `object_kind`
is the object's table and `object_row_id` its `id`: the row number never changes, the `uid` is the thing that does.
"""

import contextlib

from planetgen.db import kind_tables
from planetgen.galaxy import designations, object_uid
from planetgen.util import generator_epoch

BODY_TABLES = ("stars", "planets", "moons", "asteroid_belts", "comets")
"""tuple: The tables whose IDs embed their system's (a body's ID carries its system's sector and serial). A system that changes ID renumbers them
all in the same transaction (invariant 4). A belt member's or captured body's designation (`<belt or planet ID>.<rank>`) is derived from its
parent's ID and never stored, so it follows without a row."""

REFERENCE_TABLES = (("constellations", "system_uid"),)
"""tuple: `(table, column)` of the other places that hold a system's or body's ID as a live reference and so move with it. Records of what
happened (the edit journal, change plans, orbital events) keep the ID they were written with; the history resolves it."""


@contextlib.contextmanager
def paused(conn):
    """No history is written on `conn` inside the block: for the stand-in rows a system's content is written through before they are dropped
    (`store.replace_star_system_content`), whose names were never an object's."""
    before = getattr(conn, "history_paused", False)
    conn.history_paused = True
    try:
        yield
    finally:
        conn.history_paused = before


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
    if not old_name or old_name == new_name or getattr(conn, "history_paused", False):
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

    A system takes every body with it (invariant 4, `designations.cascade`): each body's ID keeps its body number under the system's new one and
    gets its own history row, children after the parent. A body of a system takes an ID of that system (its ID embeds it). A belt's or planet's
    members need nothing: their designation is derived from the parent's.

    Returns:
        int: `new_uid`.

    Raises:
        DesignationError: A reason off the list, an ID that is or ever was any object's, a table with no IDs, a body ID in another system, an
            object that is locked for an edit, or no such object.
    """
    designations.check_reason(reason)
    if table not in kind_tables.tables("uid"):
        raise designations.DesignationError(f"{table} has no object IDs")
    row = conn.execute(f"SELECT uid FROM {table} WHERE id = ?", (row_id,)).fetchone()
    if row is None or row["uid"] is None:
        raise designations.DesignationError(f"no such object: {table} {row_id}")
    old_uid = object_uid.from_bytes(row["uid"])
    moved = {old_uid: new_uid}
    members = []
    if table == "star_systems":
        for body_table in BODY_TABLES:
            members.extend((body_table, body["id"], object_uid.from_bytes(body["uid"])) for body in conn.execute(
                f"SELECT id, uid FROM {body_table} WHERE star_system_id = ? AND uid IS NOT NULL ORDER BY id", (row_id,)).fetchall())
        moved.update(designations.cascade(old_uid, new_uid, [uid for _t, _i, uid in members]))
    elif table in BODY_TABLES:
        system = conn.execute(f"SELECT s.uid FROM star_systems s JOIN {table} t ON t.star_system_id = s.id WHERE t.id = ?", (row_id,)).fetchone()
        if system is not None and system["uid"] is not None and not object_uid.same_system(new_uid, object_uid.from_bytes(system["uid"])):
            raise designations.DesignationError(f"{object_uid.format_id(new_uid)} is not an ID of the system this {table[:-1]} is in")
    _refuse_locked(conn, moved)
    for target in moved.values():
        designations.check_new_uid(target, {target} if is_live_uid(conn, target) else (), {target} if is_retired_uid(conn, target) else ())
    if len(set(moved.values())) != len(moved):
        raise designations.DesignationError("two objects would be given one ID")
    steps = [(table, row_id, old_uid)] + [(t, i, uid) for t, i, uid in members]
    for step_table, step_id, uid in steps:
        _retire_uid(conn, step_table, step_id, uid, reason if step_table == table and step_id == row_id else "parent change", epoch)
    for step_table, step_id, uid in steps:
        conn.execute(f"UPDATE {step_table} SET uid = ? WHERE id = ?", (object_uid.to_bytes(moved[uid]), step_id))
    kinds = {name: kind for kind, name in kind_tables.kind_tables("uid").items()}
    for step_table, _step_id, uid in steps:
        conn.execute("UPDATE object_index SET uid = ? WHERE kind = ? AND uid = ?",
                     (object_uid.to_bytes(moved[uid]), kinds[step_table], object_uid.to_bytes(uid)))
    for uid, target in moved.items():
        conn.execute("UPDATE object_index SET parent_uid = ? WHERE parent_uid = ?", (object_uid.to_bytes(target), object_uid.to_bytes(uid)))
        for reference_table, column in REFERENCE_TABLES:
            conn.execute(f"UPDATE {reference_table} SET {column} = ? WHERE {column} = ?", (object_uid.to_bytes(target), object_uid.to_bytes(uid)))
    return new_uid


def _retire_uid(conn, table, row_id, uid, reason, epoch):
    conn.execute(
        "INSERT INTO designations (designation_kind, uid, object_kind, object_row_id, valid_from_epoch, retired_epoch, reason)"
        " VALUES ('uid', ?, ?, ?, ?, ?, ?)",
        (object_uid.to_bytes(uid), table, row_id, _previous_epoch(conn, designations.KIND_UID, table, row_id), _epoch(epoch), reason))


def _refuse_locked(conn, moved):
    """An object locked for an edit (`object_locks`) keeps its ID until the lock goes."""
    for uid in moved:
        if conn.execute("SELECT 1 AS held FROM object_locks WHERE object_uid = ? LIMIT 1", (object_uid.to_bytes(uid),)).fetchone() is not None:
            raise designations.DesignationError(f"{object_uid.format_id(uid)} is locked for an edit")


def history(conn, table, row_id):
    """The `Retired` records of one object, newest first."""
    rows = conn.execute(
        "SELECT designation_kind, uid, name, object_kind, object_row_id, valid_from_epoch, retired_epoch, reason, gone_epoch FROM designations"
        " WHERE object_kind = ? AND object_row_id = ? ORDER BY retired_epoch DESC, id DESC", (table, row_id)).fetchall()
    return [designations.Retired(
        row["designation_kind"], object_uid.from_bytes(row["uid"]) if row["designation_kind"] == designations.KIND_UID else row["name"],
        row["object_kind"], row["object_row_id"], row["valid_from_epoch"], row["retired_epoch"], row["reason"], row["gone_epoch"])
        for row in rows]


def former_designations(conn, table, row_ids):
    """
    `{row id: [{"designation", "kind", "reason", "retired_epoch"}, ...]}` for the objects `row_ids` of `table`, newest first, for the "Formerly" lines
    and the `former_designations` field of an answer (an object with none is not in the dict). A retired ID is printed (GEN.1140).
    """
    row_ids = list(dict.fromkeys(row_ids))
    found = {}
    for start in range(0, len(row_ids), 1000):
        chunk = row_ids[start:start + 1000]
        for row in conn.execute(
                "SELECT designation_kind, uid, name, object_row_id, valid_from_epoch, retired_epoch, reason, gone_epoch FROM designations"
                f" WHERE object_kind = ? AND object_row_id IN ({', '.join('?' * len(chunk))})", (table, *chunk)).fetchall():
            found.setdefault(row["object_row_id"], []).append(designations.Retired(
                row["designation_kind"], object_uid.from_bytes(row["uid"]) if row["designation_kind"] == designations.KIND_UID else row["name"],
                table, row["object_row_id"], row["valid_from_epoch"], row["retired_epoch"], row["reason"], row["gone_epoch"]))
    return {row_id: [{"designation": text, "kind": kind, "reason": reason, "retired_epoch": epoch}
                     for kind, text, reason, epoch in designations.former_lines(rows)]
            for row_id, rows in found.items()}


def attach_former(conn, table, rows):
    """Sets `former_designations` (a list, empty for none) on each dict of `rows`, the objects of `table` that carry their row number as `id`."""
    rows = [row for row in rows if isinstance(row, dict)]
    found = former_designations(conn, table, [row["id"] for row in rows if row.get("id") is not None])
    for row in rows:
        row["former_designations"] = found.get(row.get("id"), [])
    return rows


def tombstone(conn, table, where="1 = 1", params=(), epoch=None):
    """
    Marks the rows of `table` that match `where` (SQL over the alias `t`) as deleted: call it in the transaction, before the `DELETE`. The current
    ID of each is retired with reason `deleted` and the epoch (so it is never given again and a link to it answers 410), and the rows the object
    already had get `gone_epoch`.
    """
    epoch = _epoch(epoch)
    if table in kind_tables.tables("uid"):
        conn.execute(
            "INSERT INTO designations (designation_kind, uid, object_kind, object_row_id, retired_epoch, reason, gone_epoch)"
            f" SELECT 'uid', t.uid, '{table}', t.id, ?, 'deleted', ? FROM {table} t WHERE t.uid IS NOT NULL AND ({where})", (epoch, epoch, *params))
    conn.execute(
        f"UPDATE designations d JOIN {table} t ON d.object_kind = '{table}' AND d.object_row_id = t.id SET d.gone_epoch = ?"
        f" WHERE d.gone_epoch IS NULL AND ({where})", (epoch, *params))


def tombstone_systems(conn, column, value, epoch=None):
    """`tombstone` for the systems with `column` (`id` or `sector_id`) equal to `value`, and every body and facility in them."""
    if column not in ("id", "sector_id"):
        raise ValueError(f"systems are chosen by id or sector_id, not {column!r}")
    tombstone(conn, "star_systems", f"t.{column} = ?", (value,), epoch)
    inside = f"t.star_system_id IN (SELECT s.id FROM star_systems s WHERE s.{column} = ?)"
    for table in (*BODY_TABLES, "facilities"):
        tombstone(conn, table, inside, (value,), epoch)


def snapshot_bodies(conn, system_id):
    """`{table: {row id: uid bytes}}` of a system's bodies, taken before `replace_system_content` deletes and writes them again under the same IDs."""
    return {table: {row["id"]: bytes(row["uid"]) for row in conn.execute(
        f"SELECT id, uid FROM {table} WHERE star_system_id = ? AND uid IS NOT NULL", (system_id,)).fetchall()} for table in BODY_TABLES}


def carry_over(conn, snapshot, epoch=None):
    """
    After a system's bodies were written again (`snapshot_bodies` before): the history of each body follows it to its new row (matched by the ID it
    kept), and a body whose ID did not come back is deleted: its ID is retired with the tombstone.
    """
    epoch = _epoch(epoch)
    for table, rows in snapshot.items():
        for old_id, uid in rows.items():
            new = conn.execute(f"SELECT id FROM {table} WHERE uid = ?", (uid,)).fetchone()
            if new is not None:
                conn.execute("UPDATE designations SET object_row_id = ? WHERE object_kind = ? AND object_row_id = ?", (new["id"], table, old_id))
                continue
            conn.execute(
                "INSERT INTO designations (designation_kind, uid, object_kind, object_row_id, retired_epoch, reason, gone_epoch)"
                " VALUES ('uid', ?, ?, ?, ?, 'deleted', ?)", (uid, table, old_id, epoch, epoch))
            conn.execute("UPDATE designations SET gone_epoch = ? WHERE object_kind = ? AND object_row_id = ? AND gone_epoch IS NULL",
                         (epoch, table, old_id))


def resolve_member(conn, text):
    """
    What a belt member's or captured body's designation (`<parent ID>.<rank hex>`) names now, however old the parent ID in it is
    (`designations.resolve_member`): `(Resolution of the parent, current designation or None)`. The rank is kept; a deleted parent gives `None`.
    """
    def parent(uid):
        raw = object_uid.to_bytes(uid)
        for table in kind_tables.tables("uid"):
            found = conn.execute(f"SELECT id FROM {table} WHERE uid = ?", (raw,)).fetchone()
            if found is not None:
                return designations.Resolution(designations.CURRENT, table, found["id"])
        row = conn.execute("SELECT object_kind, object_row_id, retired_epoch, gone_epoch FROM designations WHERE uid = ?", (raw,)).fetchone()
        held = [] if row is None else [designations.Retired(designations.KIND_UID, uid, row["object_kind"], row["object_row_id"], None,
                                                            row["retired_epoch"], "rename", row["gone_epoch"])]
        found = designations.resolve_uid(uid, lambda _uid: None, held)
        if found.status == designations.REDIRECT and conn.execute(
                f"SELECT 1 AS here FROM {found.object_kind} WHERE id = ?", (found.object_row_id,)).fetchone() is None:
            return designations.Resolution(designations.GONE, found.object_kind, found.object_row_id, found.formerly)
        return found

    def current_uid(table, row_id):
        return object_uid.from_bytes(conn.execute(f"SELECT uid FROM {table} WHERE id = ?", (row_id,)).fetchone()["uid"])

    return designations.resolve_member(text, parent, current_uid)
