# planetgen/population/facility_types.py

"""
Facility types and their typed fields (POP.10, schema v89).

A **type** refines a facility kind ("Agricultural colony" is a colony), so
`kind` keeps driving `tuning.FACILITY_RULES`. An admin defines types, and on
each type up to `MAX_FIELDS` typed **fields** (`int`, `decimal`, `text`,
`bool` or `choice`); a facility stores one value row per field it has a
value for. Everything here takes an open connection and the caller owns the
transaction. Rules broken raise `FacilityTypeError`, whose `status` is the
HTTP status the API answers with.

See `docs/design/habitability/population-and-politics.md`, "Facility types and ratings".
"""

import math
import re

from planetgen import tuning
from planetgen.population.facility_icons import ICON_NAMES

DATATYPES = ("int", "decimal", "text", "bool", "choice")
"""tuple: Every `facility_type_fields.datatype`."""

MAX_FIELDS = 24
"""int: The most fields one type may have."""

MAX_CHOICES = 32
"""int: The most options one `choice` field may list."""

MAX_CHOICE_LENGTH = 64
"""int: The longest option of a `choice` field, in characters."""

MAX_TEXT_LENGTH = 255
"""int: The longest `text` value, `default_text`, label or name (`VARCHAR(255)`)."""

MAX_UNIT_LENGTH = 32
"""int: `facility_type_fields.unit` is `VARCHAR(32)`."""

MAX_DESCRIPTION_LENGTH = 4000
"""int: The longest type description."""

KEY_PATTERN = re.compile(r"^[a-z][a-z0-9_]{0,31}$")
"""re.Pattern: What a field key looks like; the key names the field in the API."""

_RANGE_TYPES = ("int", "decimal")
_TYPE_COLUMNS = ("id", "name", "kind", "icon", "description", "sort_order", "active")
_FIELD_COLUMNS = ("id", "type_id", "key", "label", "datatype", "unit", "required", "min_value", "max_value",
                  "choices", "default_text", "sort_order")
_TYPE_EDITABLE = ("name", "kind", "icon", "description", "sort_order", "active")
_FIELD_EDITABLE = ("label", "unit", "required", "min_value", "max_value", "choices", "default_text", "sort_order",
                   "datatype")


class FacilityTypeError(ValueError):
    """A facility type or field that breaks a rule. `status` is 400 (bad value), 404 (no such type or field) or
    409 (the change clashes with facilities or another type)."""

    def __init__(self, message, status=400):
        super().__init__(message)
        self.status = status


def split_choices(text):
    """The options of a stored `choices` text, in order."""
    return [line for line in (text or "").split("\n") if line]


def _options(field):
    choices = field.get("choices")
    return list(choices) if isinstance(choices, (list, tuple)) else split_choices(choices)


def parse_value(field, raw):
    """
    One value of a field, checked against its definition.

    Args:
        field (dict): A field as `get_type` or `definition_fields` returns it.
        raw: The value as sent (a JSON value) or `default_text` (a string).

    Returns:
        The value as Python stores it: `int`, `float`, `str` or `bool`.

    Raises:
        FacilityTypeError: 400 naming the field when the value does not fit.
    """
    label, datatype = field["key"], field["datatype"]
    if datatype == "bool":
        if isinstance(raw, str) and raw.strip().lower() in ("true", "false"):
            return raw.strip().lower() == "true"
        if isinstance(raw, bool):
            return raw
        raise FacilityTypeError(f"field {label!r} must be true or false")
    if datatype in _RANGE_TYPES:
        if isinstance(raw, bool) or not isinstance(raw, (int, float, str)):
            raise FacilityTypeError(f"field {label!r} must be a number")
        try:
            number = float(raw)
        except ValueError:
            raise FacilityTypeError(f"field {label!r} must be a number") from None
        if not math.isfinite(number):
            raise FacilityTypeError(f"field {label!r} must be a finite number")
        if datatype == "int":
            if number != int(number):
                raise FacilityTypeError(f"field {label!r} must be a whole number")
            value = raw if isinstance(raw, int) else int(number)
        else:
            value = number
        if field.get("min_value") is not None and value < field["min_value"]:
            raise FacilityTypeError(f"field {label!r} must be at least {field['min_value']:g}")
        if field.get("max_value") is not None and value > field["max_value"]:
            raise FacilityTypeError(f"field {label!r} must be at most {field['max_value']:g}")
        return value
    if not isinstance(raw, str):
        raise FacilityTypeError(f"field {label!r} must be text")
    if datatype == "choice":
        if raw not in _options(field):
            raise FacilityTypeError(f"field {label!r} must be one of: {', '.join(_options(field))}")
        return raw
    if len(raw) > MAX_TEXT_LENGTH:
        raise FacilityTypeError(f"field {label!r} must be at most {MAX_TEXT_LENGTH} characters")
    return raw


def _clean_field(spec, existing=None):
    """
    A field spec checked and normalised: the columns of a `facility_type_fields` row (less `id` and `type_id`).
    `existing` is the stored field a change applies to, so only the keys in `spec` are checked anew.
    """
    row = dict(existing or {})
    row.update(spec)
    key, datatype = row.get("key"), row.get("datatype")
    if not isinstance(key, str) or not KEY_PATTERN.match(key):
        raise FacilityTypeError("a field 'key' starts with a lowercase letter and has up to 32 lowercase letters, "
                                "digits or underscores")
    if datatype not in DATATYPES:
        raise FacilityTypeError(f"field {key!r}: 'datatype' must be one of: {', '.join(DATATYPES)}")
    label = (row.get("label") or "").strip()
    if not label or len(label) > MAX_TEXT_LENGTH:
        raise FacilityTypeError(f"field {key!r}: 'label' must be 1 to {MAX_TEXT_LENGTH} characters")
    unit = (row.get("unit") or "").strip() or None
    if unit is not None and (len(unit) > MAX_UNIT_LENGTH or datatype not in _RANGE_TYPES):
        raise FacilityTypeError(f"field {key!r}: only a number takes a 'unit', up to {MAX_UNIT_LENGTH} characters")
    low, high = row.get("min_value"), row.get("max_value")
    if (low is not None or high is not None) and datatype not in _RANGE_TYPES:
        raise FacilityTypeError(f"field {key!r}: only a number takes 'min_value' and 'max_value'")
    if low is not None and high is not None and low > high:
        raise FacilityTypeError(f"field {key!r}: 'min_value' is above 'max_value'")
    choices = row.get("choices")
    if datatype == "choice":
        options = split_choices(choices) if isinstance(choices, str) else list(choices or [])
        if not 2 <= len(options) <= MAX_CHOICES:
            raise FacilityTypeError(f"field {key!r}: a choice lists 2 to {MAX_CHOICES} options")
        if len(set(options)) != len(options):
            raise FacilityTypeError(f"field {key!r}: the options must all differ")
        for option in options:
            if not option.strip() or option != option.strip() or len(option) > MAX_CHOICE_LENGTH or "\n" in option:
                raise FacilityTypeError(f"field {key!r}: each option is 1 to {MAX_CHOICE_LENGTH} characters, "
                                        "with no spaces at either end")
        choices = "\n".join(options)
    elif choices:
        raise FacilityTypeError(f"field {key!r}: only a choice takes 'choices'")
    else:
        choices = None
    cleaned = {
        "key": key, "label": label, "datatype": datatype, "unit": unit, "required": 1 if row.get("required") else 0,
        "min_value": low, "max_value": high, "choices": choices, "default_text": None,
        "sort_order": int(row.get("sort_order") or 0),
    }
    default = row.get("default_text")
    if default not in (None, ""):
        value = parse_value(cleaned, default)
        cleaned["default_text"] = ("true" if value else "false") if datatype == "bool" else str(value)
    return cleaned


def _clean_type(spec, existing=None):
    row = dict(existing or {})
    row.update(spec)
    name = (row.get("name") or "").strip()
    if not name or len(name) > MAX_TEXT_LENGTH:
        raise FacilityTypeError(f"'name' must be 1 to {MAX_TEXT_LENGTH} characters")
    if row.get("kind") not in tuning.FACILITY_KINDS:
        raise FacilityTypeError(f"'kind' must be one of: {', '.join(sorted(tuning.FACILITY_KINDS))}")
    icon = row.get("icon") or None
    if icon is not None and icon not in ICON_NAMES:
        raise FacilityTypeError(f"'icon' must be one of the listed icons ({len(ICON_NAMES)} names)")
    description = row.get("description") or None
    if description is not None and len(description) > MAX_DESCRIPTION_LENGTH:
        raise FacilityTypeError(f"'description' must be at most {MAX_DESCRIPTION_LENGTH} characters")
    return {"name": name, "kind": row["kind"], "icon": icon, "description": description,
            "sort_order": int(row.get("sort_order") or 0), "active": 0 if row.get("active") is False or row.get("active") == 0 else 1}


def _field_dict(row):
    field = {column: row[column] for column in _FIELD_COLUMNS}
    field["required"] = bool(field["required"])
    field["choices"] = split_choices(field["choices"]) if field["datatype"] == "choice" else None
    return field


def _type_dict(row, fields, facility_count):
    found = {column: row[column] for column in _TYPE_COLUMNS}
    found["active"] = bool(found["active"])
    found["fields"] = fields
    found["facilities"] = facility_count
    return found


def _raw_fields(conn, type_id):
    return conn.execute(
        f"SELECT {', '.join('`' + c + '`' for c in _FIELD_COLUMNS)} FROM facility_type_fields WHERE type_id = ? "
        "ORDER BY sort_order, id", (type_id,)).fetchall()


def definition_fields(conn, type_id):
    """The fields of a type as the checks in `parse_value` read them (`choices` still the stored text), in sort
    order. Empty for an unknown type."""
    return [dict(row) for row in _raw_fields(conn, type_id)]


def get_type(conn, type_id):
    """One type with its fields and how many facilities use it, or `None`."""
    row = conn.execute(f"SELECT {', '.join(_TYPE_COLUMNS)} FROM facility_types WHERE id = ?", (type_id,)).fetchone()
    if row is None:
        return None
    count = conn.execute("SELECT COUNT(*) AS n FROM facilities WHERE type_id = ?", (type_id,)).fetchone()["n"]
    return _type_dict(row, [_field_dict(f) for f in _raw_fields(conn, type_id)], count)


def list_types(conn, kind=None, active_only=False):
    """Every type (with its fields and facility count) in sort order, optionally one kind or only the active."""
    where, params = [], []
    if kind is not None:
        where.append("kind = ?")
        params.append(kind)
    if active_only:
        where.append("active = 1")
    rows = conn.execute(
        f"SELECT {', '.join(_TYPE_COLUMNS)} FROM facility_types"
        + (f" WHERE {' AND '.join(where)}" if where else "") + " ORDER BY sort_order, name, id", params).fetchall()
    fields = {}
    for row in conn.execute(
            f"SELECT {', '.join('`' + c + '`' for c in _FIELD_COLUMNS)} FROM facility_type_fields "
            "ORDER BY sort_order, id").fetchall():
        fields.setdefault(row["type_id"], []).append(_field_dict(row))
    counts = {row["type_id"]: row["n"] for row in conn.execute(
        "SELECT type_id, COUNT(*) AS n FROM facilities WHERE type_id IS NOT NULL GROUP BY type_id").fetchall()}
    return [_type_dict(row, fields.get(row["id"], []), counts.get(row["id"], 0)) for row in rows]


def _require_type(conn, type_id):
    found = get_type(conn, type_id)
    if found is None:
        raise FacilityTypeError(f"no such facility type: {type_id}", status=404)
    return found


def _name_taken(conn, name, except_id=None):
    row = conn.execute("SELECT id FROM facility_types WHERE name = ?", (name,)).fetchone()
    if row is not None and row["id"] != except_id:
        raise FacilityTypeError(f"a facility type named {name!r} already exists", status=409)


def create_type(conn, spec, fields=()):
    """
    Adds a type and its fields.

    Args:
        spec (dict): `name`, `kind` (required) and `icon`, `description`, `sort_order`, `active`.
        fields (iterable of dict): Field specs: `key`, `label`, `datatype` (required) and `unit`, `required`,
            `min_value`, `max_value`, `choices` (a list), `default_text`, `sort_order`.

    Returns:
        int: The new type's id.
    """
    cleaned = _clean_type(spec)
    fields = [_clean_field(field) for field in fields]
    if len(fields) > MAX_FIELDS:
        raise FacilityTypeError(f"a type has at most {MAX_FIELDS} fields")
    keys = [field["key"] for field in fields]
    if len(set(keys)) != len(keys):
        raise FacilityTypeError("each field needs its own 'key'")
    _name_taken(conn, cleaned["name"])
    type_id = conn.execute(
        "INSERT INTO facility_types (name, kind, icon, description, sort_order, active) VALUES (?, ?, ?, ?, ?, ?)",
        (cleaned["name"], cleaned["kind"], cleaned["icon"], cleaned["description"], cleaned["sort_order"],
         cleaned["active"])).lastrowid
    for field in fields:
        _insert_field(conn, type_id, field)
    return type_id


def update_type(conn, type_id, changes):
    """
    Changes the named columns of a type (`name`, `kind`, `icon`, `description`, `sort_order`, `active`).

    Raises:
        FacilityTypeError: 404 for an unknown type, 409 for a name already taken or a new kind while facilities
            use the type (a facility's kind and its type must agree).
    """
    found = _require_type(conn, type_id)
    cleaned = _clean_type(changes, existing=found)
    if cleaned["kind"] != found["kind"] and found["facilities"]:
        raise FacilityTypeError(f"{found['facilities']} facilities use this type, so its kind stays "
                                f"{found['kind']!r}; add a new type instead", status=409)
    if cleaned["name"] != found["name"]:
        _name_taken(conn, cleaned["name"], except_id=type_id)
    columns = {column: cleaned[column] for column in _TYPE_EDITABLE if column in changes}
    if columns:
        conn.execute(f"UPDATE facility_types SET {', '.join(f'{column} = ?' for column in columns)} WHERE id = ?",
                     (*columns.values(), type_id))


def delete_type(conn, type_id):
    """
    Deletes a type and its fields.

    Raises:
        FacilityTypeError: 404 for an unknown type; 409 while facilities use it (make it inactive instead).
    """
    found = _require_type(conn, type_id)
    if found["facilities"]:
        raise FacilityTypeError(f"{found['facilities']} facilities use this type; set it inactive instead of "
                                "deleting it", status=409)
    conn.execute("DELETE FROM facility_types WHERE id = ?", (type_id,))


def _insert_field(conn, type_id, field):
    return conn.execute(
        "INSERT INTO facility_type_fields (type_id, `key`, label, datatype, unit, required, min_value, max_value, "
        "choices, default_text, sort_order) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (type_id, field["key"], field["label"], field["datatype"], field["unit"], field["required"],
         field["min_value"], field["max_value"], field["choices"], field["default_text"], field["sort_order"]),
    ).lastrowid


def _values_exist(conn, field_id):
    return conn.execute("SELECT COUNT(*) AS n FROM facility_field_values WHERE field_id = ?",
                        (field_id,)).fetchone()["n"]


def add_field(conn, type_id, spec):
    """Adds one field to a type; returns its id. 409 for a key already on the type or a 25th field."""
    found = _require_type(conn, type_id)
    field = _clean_field(spec)
    if len(found["fields"]) >= MAX_FIELDS:
        raise FacilityTypeError(f"a type has at most {MAX_FIELDS} fields", status=409)
    if any(other["key"] == field["key"] for other in found["fields"]):
        raise FacilityTypeError(f"this type already has a field {field['key']!r}", status=409)
    if field["required"] and found["facilities"] and field["default_text"] is None:
        raise FacilityTypeError(f"{found['facilities']} facilities use this type, so a required field needs a "
                                "'default_text' they can take", status=409)
    field_id = _insert_field(conn, type_id, field)
    _fill_default(conn, type_id, field_id, field)
    return field_id


def _fill_default(conn, type_id, field_id, field):
    """Gives every facility of the type the field's default as its value (a field added after facilities exist)."""
    if field["default_text"] is None:
        return
    value = parse_value(field, field["default_text"])
    column = value_column(field["datatype"])
    conn.execute(
        f"INSERT INTO facility_field_values (facility_id, field_id, {column}) "
        f"SELECT id, ?, ? FROM facilities WHERE type_id = ?", (field_id, value, type_id))


def value_column(datatype):
    """The `facility_field_values` column a datatype is stored in."""
    return {"int": "value_int", "decimal": "value_dec", "bool": "value_bool"}.get(datatype, "value_text")


def update_field(conn, type_id, field_id, changes):
    """
    Changes a field (`label`, `unit`, `required`, `min_value`, `max_value`, `choices`, `default_text`,
    `sort_order`, `datatype`). The key never changes. While facilities hold values for the field, its datatype
    stays, no option in use is dropped and no stored number falls outside a narrowed range (409).
    """
    found = _require_type(conn, type_id)
    old = next((field for field in found["fields"] if field["id"] == field_id), None)
    if old is None:
        raise FacilityTypeError(f"no such field {field_id} on facility type {type_id}", status=404)
    stored = dict(old, choices="\n".join(old["choices"]) if old["choices"] else None)
    if "choices" in changes and isinstance(changes["choices"], list):
        changes = dict(changes, choices="\n".join(changes["choices"]))
    cleaned = _clean_field(changes, existing=stored)
    in_use = _values_exist(conn, field_id)
    if in_use:
        _check_in_use(conn, field_id, old, cleaned, in_use)
    columns = {column: cleaned[column] for column in _FIELD_EDITABLE if column in changes}
    if columns:
        conn.execute(f"UPDATE facility_type_fields SET {', '.join(f'`{c}` = ?' for c in columns)} WHERE id = ?",
                     (*columns.values(), field_id))


def _check_in_use(conn, field_id, old, cleaned, in_use):
    if cleaned["datatype"] != old["datatype"]:
        raise FacilityTypeError(f"{in_use} facilities hold a value for this field, so its datatype stays "
                                f"{old['datatype']!r}", status=409)
    if cleaned["datatype"] == "choice":
        kept = split_choices(cleaned["choices"])
        used = [row["value_text"] for row in conn.execute(
            "SELECT DISTINCT value_text FROM facility_field_values WHERE field_id = ?", (field_id,)).fetchall()]
        dropped = [option for option in used if option not in kept]
        if dropped:
            raise FacilityTypeError(f"facilities still use the options {', '.join(map(repr, dropped))}", status=409)
    if cleaned["datatype"] in _RANGE_TYPES:
        column = value_column(cleaned["datatype"])
        row = conn.execute(f"SELECT MIN({column}) AS low, MAX({column}) AS high FROM facility_field_values "
                           "WHERE field_id = ?", (field_id,)).fetchone()
        if cleaned["min_value"] is not None and row["low"] is not None and row["low"] < cleaned["min_value"]:
            raise FacilityTypeError(f"a facility holds {row['low']:g}, below the new 'min_value'", status=409)
        if cleaned["max_value"] is not None and row["high"] is not None and row["high"] > cleaned["max_value"]:
            raise FacilityTypeError(f"a facility holds {row['high']:g}, above the new 'max_value'", status=409)


def delete_field(conn, type_id, field_id):
    """Deletes a field nobody holds a value for; 409 otherwise, 404 for an unknown field."""
    row = conn.execute("SELECT id FROM facility_type_fields WHERE id = ? AND type_id = ?",
                       (field_id, type_id)).fetchone()
    if row is None:
        raise FacilityTypeError(f"no such field {field_id} on facility type {type_id}", status=404)
    in_use = _values_exist(conn, field_id)
    if in_use:
        raise FacilityTypeError(f"{in_use} facilities hold a value for this field; it cannot be deleted", status=409)
    conn.execute("DELETE FROM facility_type_fields WHERE id = ?", (field_id,))


EXAMPLE_TYPES = (
    {
        "name": "Research outpost", "kind": "outpost", "icon": "flask", "sort_order": 10,
        "description": "A small crewed post that studies the world or system it watches.",
        "fields": (
            {"key": "staff", "label": "Staff", "datatype": "int", "min_value": 1, "max_value": 10000,
             "required": True, "default_text": "12", "sort_order": 1},
            {"key": "focus", "label": "Focus", "datatype": "choice", "sort_order": 2,
             "choices": ["Survey", "Laboratory", "Relay"], "default_text": "Survey"},
            {"key": "armed", "label": "Armed", "datatype": "bool", "sort_order": 3},
        ),
    },
)
"""tuple: The example types `install_examples` adds: the mechanism is there to show what a type can do."""


def install_examples(conn):
    """Adds each of `EXAMPLE_TYPES` that is not there yet (by name). Returns the names added."""
    added = []
    for spec in EXAMPLE_TYPES:
        if conn.execute("SELECT 1 FROM facility_types WHERE name = ?", (spec["name"],)).fetchone() is None:
            create_type(conn, {key: value for key, value in spec.items() if key != "fields"}, spec["fields"])
            added.append(spec["name"])
    return added
