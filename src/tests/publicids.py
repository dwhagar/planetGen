# tests/publicids.py

"""
The printed object IDs (API.23) of rows a test seeded by row id: the API, the pages and their URLs speak
IDs, while the seeding helpers still return row ids.

`pid(kind, row_id)` is the printed ID of a row in the current test's database (`mysql_config` registers it),
so a test writes `f"/api/systems/{pid('system', system_id)}"`. A value that is already printed (text with a
letter or a hyphen) is returned as it is, so a helper can take either.
"""

from planetgen.api import ids
from planetgen.db import store

_config = None


def use(config):
    """Registers (or clears) the database `pid` reads; `mysql_config` calls this."""
    global _config
    _config = config


def _printed_text(value):
    return isinstance(value, str) and not value.isdigit() and ids.looks_printed(value)


def _printed_or_issue(conn, kind, row_id):
    """The row's printed ID; a row a test inserted by hand has no `uid` yet, so it is issued one first (the way
    `store.assign_uids` does for any row that arrives outside a sector save)."""
    if kind != ids.SECTOR:
        table = ids.OBJECT_TABLES[kind]
        row = conn.execute(f"SELECT uid FROM {table} WHERE id = ?", (row_id,)).fetchone()
        if row is not None and row["uid"] is None:
            if kind == "system":
                store.assign_uids(conn, system_ids=[row_id])
            else:
                store.assign_uids(conn, phenomenon=(table, row_id))
            conn.commit()
    return ids.printed(conn, kind, row_id)


def pid(kind, row_id, config=None):
    """The printed ID of one row of `kind` (`None` stays `None`)."""
    if row_id is None or _printed_text(row_id):
        return row_id
    config = config or _config
    if hasattr(config, "execute"):
        return _printed_or_issue(config, kind, int(row_id))
    conn = store.get_connection(config)
    try:
        return _printed_or_issue(conn, kind, int(row_id))
    finally:
        conn.close()


def pids(kind, row_ids, config=None):
    """The printed IDs of several rows, in order."""
    return [pid(kind, row_id, config) for row_id in row_ids]


def psys(row_id, config=None):
    """`pid("system", row_id)`, for a spot where quotes are scarce (inside a single-quoted f-string)."""
    return pid("system", row_id, config)


def psec(row_id, config=None):
    """`pid("sector", row_id)`, for the same spots."""
    return pid("sector", row_id, config)
