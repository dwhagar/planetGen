# planetgen/population/facility_attrs.py

"""
What a facility carries beyond its place (POP.10): a type with the values of
that type's fields, an affiliation and five ratings.

`resolve` checks what a caller sent against the database and answers a plan;
`write` stores the plan. Both leave the transaction to the caller and raise
`FacilityTypeError` (its `status` is the HTTP status) for anything refused.
"""

from planetgen.population import facility_ratings as ratings
from planetgen.population.facility_types import FacilityTypeError, definition_fields, get_type, parse_value, value_column

MAX_NOTE_LENGTH = 4000
"""int: The longest `ratings_note`."""


def _resolve_type(conn, kind, attrs, current):
    """The type a facility ends up with: `(type_id or None, changed)`; checks it exists, is active and fits `kind`."""
    if "type_id" not in attrs:
        return (current or {}).get("type_id"), False
    type_id = attrs["type_id"]
    if type_id is None:
        return None, (current or {}).get("type_id") is not None
    found = get_type(conn, type_id)
    if found is None:
        raise FacilityTypeError(f"no such facility type: {type_id}", status=404)
    if found["kind"] != kind:
        raise FacilityTypeError(f"{found['name']!r} is a {found['kind']} type; this facility is a {kind}")
    changed = type_id != (current or {}).get("type_id")
    if changed and not found["active"]:
        raise FacilityTypeError(f"{found['name']!r} is retired")
    return type_id, changed


def _resolve_affiliation(conn, attrs):
    """The affiliation columns for the `affiliation_*` keys sent."""
    polity_id, name = attrs.get("affiliation_polity_id"), attrs.get("affiliation_name")
    if polity_id is not None and name is not None:
        raise FacilityTypeError("send 'affiliation_polity_id' or 'affiliation_name', not both")
    if polity_id is not None:
        row = conn.execute("SELECT name FROM polities WHERE id = ?", (polity_id,)).fetchone()
        if row is None:
            raise FacilityTypeError(f"no such polity: {polity_id}", status=404)
        return {"affiliation_polity_id": polity_id, "affiliation_name": row["name"]}
    if name is not None:
        name = name.strip()
        if not name or len(name) > 255:
            raise FacilityTypeError("'affiliation_name' must be 1 to 255 characters")
        return {"affiliation_polity_id": None, "affiliation_name": name}
    return {"affiliation_polity_id": None, "affiliation_name": None}


def _resolve_values(conn, type_id, type_changed, given, current_values):
    """
    The field values to store: `(upserts, deletes)`, where `upserts` is `{field_id: (column, value)}` and
    `deletes` the field ids to clear. A new type starts clean: its required fields must be given or have a
    default.
    """
    fields = definition_fields(conn, type_id) if type_id is not None else []
    by_key = {field["key"]: field for field in fields}
    unknown = sorted(set(given) - set(by_key))
    if unknown:
        raise FacilityTypeError(f"this facility's type has no field {', '.join(map(repr, unknown))}")
    upserts, deletes = {}, set()
    held = {} if type_changed else dict(current_values)
    for key, raw in given.items():
        field = by_key[key]
        if raw is None:
            if field["required"]:
                raise FacilityTypeError(f"field {key!r} is required")
            deletes.add(field["id"])
            held.pop(field["id"], None)
        else:
            upserts[field["id"]] = (value_column(field["datatype"]), parse_value(field, raw))
            held[field["id"]] = True
    for field in fields:
        if field["id"] in held or field["id"] in upserts:
            continue
        if field["default_text"] is not None:
            upserts[field["id"]] = (value_column(field["datatype"]), parse_value(field, field["default_text"]))
        elif field["required"]:
            raise FacilityTypeError(f"field {field['key']!r} is required")
    if type_changed:
        deletes |= set(current_values) - set(upserts)
    return upserts, sorted(deletes - set(upserts))


def current_attrs(conn, facility_id):
    """What a facility holds now: `{"type_id", "values": {field_id: True}}`, or `None` for an unknown facility."""
    row = conn.execute("SELECT type_id FROM facilities WHERE id = ?", (facility_id,)).fetchone()
    if row is None:
        return None
    held = {r["field_id"]: True for r in conn.execute(
        "SELECT field_id FROM facility_field_values WHERE facility_id = ?", (facility_id,)).fetchall()}
    return {"type_id": row["type_id"], "values": held}


def resolve(conn, kind, attrs, current=None):
    """
    Checks the type, affiliation, ratings and field values a caller sent.

    Args:
        conn (Connection): An open connection.
        kind (str): The facility's kind (its type must be of that kind).
        attrs (dict): Any of `type_id`, `affiliation_polity_id`, `affiliation_name`, `ratings` (`{dimension:
            colour or None}`), `ratings_note` and `fields` (`{key: value or None}`); only the keys present
            change anything.
        current (dict, optional): `current_attrs` of a facility being changed; `None` for a new one.

    Returns:
        dict: `columns` (the `facilities` columns to set), `upserts` and `deletes` (see `write`).
    """
    columns = {}
    type_id, type_changed = _resolve_type(conn, kind, attrs, current)
    if "type_id" in attrs:
        columns["type_id"] = type_id
    if "affiliation_polity_id" in attrs or "affiliation_name" in attrs:
        columns.update(_resolve_affiliation(conn, attrs))
    if attrs.get("ratings"):
        for name, color in attrs["ratings"].items():
            try:
                columns[f"rating_{name}"] = None if color is None else ratings.value_of(color)
            except ValueError as exc:
                raise FacilityTypeError(f"rating {name!r}: {exc}") from None
        columns["ratings_source"] = "admin"
        columns["ratings_set_at"] = "NOW"
    if "ratings_note" in attrs:
        note = attrs["ratings_note"]
        if note is not None and len(note) > MAX_NOTE_LENGTH:
            raise FacilityTypeError(f"'ratings_note' must be at most {MAX_NOTE_LENGTH} characters")
        columns["ratings_note"] = note or None
    if "fields" in attrs and attrs["fields"] and type_id is None:
        raise FacilityTypeError("'fields' needs a 'type_id'")
    held = (current or {}).get("values", {})
    if type_id is None:
        upserts, deletes = {}, sorted(held)
    elif current is not None and not type_changed and not attrs.get("fields"):
        upserts, deletes = {}, []
    else:
        upserts, deletes = _resolve_values(conn, type_id, type_changed or current is None, attrs.get("fields") or {},
                                           held)
    return {"columns": columns, "upserts": upserts, "deletes": deletes}


def write(conn, facility_id, plan):
    """Stores a `resolve` plan for one facility: the columns, then the field values cleared and set."""
    columns = dict(plan["columns"])
    set_now = columns.pop("ratings_set_at", None) == "NOW"
    clauses = [f"`{name}` = ?" for name in columns] + (["ratings_set_at = CURRENT_TIMESTAMP"] if set_now else [])
    if clauses:
        conn.execute(f"UPDATE facilities SET {', '.join(clauses)} WHERE id = ?", (*columns.values(), facility_id))
    for field_id in plan["deletes"]:
        conn.execute("DELETE FROM facility_field_values WHERE facility_id = ? AND field_id = ?",
                     (facility_id, field_id))
    for field_id, (column, value) in plan["upserts"].items():
        conn.execute(
            f"INSERT INTO facility_field_values (facility_id, field_id, {column}) VALUES (?, ?, ?) "
            f"ON DUPLICATE KEY UPDATE {column} = ?", (facility_id, field_id, value, value))
