"""
planetgen/api/lookup.py

`GET /api/lookup?text=` (GEN.1140): what a printed Object ID or a name is, now. An ID or a name an object used to have answers with the object
that has it now (`docs/design/architecture/designation-history.md` section 5); a retired ID of a deleted object answers 410. The rules are
`galaxy/designations.resolve_uid` and `resolve_name`; this module reads the stored holders and the retired rows for them.
"""
from planetgen.db import designation_rows
from planetgen.galaxy import designations, object_uid
from . import ids

KIND_OF_TABLE = {table: kind for kind, table in ids.OBJECT_TABLES.items()}
KIND_OF_TABLE["sectors"] = ids.SECTOR
"""dict: The kind (the word before the colon of a reference) of each table that holds named objects."""

HOLDERS_PER_TABLE = 3
"""int: The most current holders of one name read from one table; a name is not unique among bodies."""

ALSO_LIMIT = 10
"""int: The most other past holders of a name listed in `also`."""


def _holders_of_name(conn, name):
    """`[(table, row id)]` of the objects named exactly `name` now, sectors and systems first."""
    holders = []
    for table in designation_rows.NAMED_TABLES:
        rows = conn.execute(f"SELECT id FROM {table} WHERE name = ? ORDER BY id LIMIT {HOLDERS_PER_TABLE}", (name,)).fetchall()
        holders.extend((table, row["id"]) for row in rows)
    return holders


def _object_answer(conn, table, row_id):
    """`{kind, ref, id, name, former_designations}` of one existing object, or `None` when the row is gone."""
    row = conn.execute(f"SELECT name FROM {table} WHERE id = ?", (row_id,)).fetchone()
    if row is None:
        return None
    kind = KIND_OF_TABLE[table]
    printed = ids.printed(conn, kind, row_id)
    return {"kind": kind, "ref": f"{kind}:{printed}", "id": printed, "name": row["name"],
            "former_designations": designation_rows.former_designations(conn, table, row_id)}


def lookup(conn, text):
    """
    Resolves a printed ID or a name.

    Returns:
        tuple: `(http status, body)`. 200 with `status` `current` (the object has the ID or name now) or `redirect` (it used to, and the body holds
            the object that has it now, with `retired` naming what was asked for); 410 `gone` for a designation of a deleted object; 404
            `unknown` when nobody ever had it. `also` lists other objects that used to have the same name, newest first.

    Raises:
        ValueError: Empty text.
    """
    text = str(text or "").strip()
    if not text:
        raise ValueError("lookup needs text: an Object ID or a name")
    try:
        uid = object_uid.parse_id(text.upper())
    except ValueError:
        uid = None
    if uid is not None:
        return _by_id(conn, text, uid)
    return _by_name(conn, text)


def _by_id(conn, text, uid):
    retired = designation_rows.retired_uid(conn, uid)
    found = designations.resolve_uid(uid, lambda wanted: designation_rows.live_holder(conn, wanted), [] if retired is None else [retired])
    return _answer(conn, text, found, [])


def _by_name(conn, text):
    retired = designation_rows.retired_names(conn, text)
    name = retired[0].designation if retired else text
    found = designations.resolve_name(name, _holders_of_name(conn, name), retired)
    return _answer(conn, text, found, retired)


def _answer(conn, text, found, retired):
    if found.status == designations.UNKNOWN:
        return 404, {"query": text, "status": designations.UNKNOWN, "error": f"nothing has, or ever had, the ID or name {text}"}
    if found.status == designations.GONE:
        return designations.STATUS_GONE, {"query": text, "status": designations.GONE,
                                          "error": f"the object that had {text} was deleted", "gone_epoch": found.formerly[0].gone_epoch}
    body = _object_answer(conn, found.object_kind, found.object_row_id)
    if body is None:     # the row went without its tombstone: the same answer as a deleted object
        return designations.STATUS_GONE, {"query": text, "status": designations.GONE, "error": f"the object that had {text} was deleted"}
    body = {"query": text, "status": found.status, **body}
    if found.status == designations.REDIRECT:
        body["retired"] = text
    also = []
    for row in retired:
        if (row.object_kind, row.object_row_id) != (found.object_kind, found.object_row_id) and row.gone_epoch is None and len(also) < ALSO_LIMIT:
            other = _object_answer(conn, row.object_kind, row.object_row_id)
            if other is not None and other["ref"] != body["ref"]:
                also.append({"kind": other["kind"], "ref": other["ref"], "name": other["name"], "retired_epoch": row.retired_epoch})
    if also:
        body["also"] = also
    return 200, body
