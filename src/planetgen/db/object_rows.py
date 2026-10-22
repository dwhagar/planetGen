# planetgen/db/object_rows.py

"""
One insert for every kind of object (DB.74): a row is written from a list of fields, not from a hand-written `INSERT`.

A `Field` names the column or columns it fills and how to read their values from the object and from the call's
context (the owning system, the sector, the placement). `store.py` declares one list of fields per kind table; `insert_row`
turns the list into the `INSERT`, so a column is spelled once (in its kind's list) instead of in the column list, the
placeholders and the value tuple of a function of its own. The constructors build the common fields:

* `attr(column, name)`: `getattr(obj, name)`, the column's own name by default.
* `optional(column, name)`: the same, `None` when the object lacks it (an object restored from older data).
* `context(column, key)`: the value the caller passed (`star_system_id`, `sector_id`).
* `computed(columns, function)`: `function(obj, context)`, one value for one column or a tuple for several.

The same lists drive the load: `load_dict(fields, row)` turns a stored row back into the dict an object's `from_dict`
takes, each field naming the keys it fills and inverting its own conversion (`scaled` divides the kilometres back to
AU). In the 8.4 revision the lists also name which columns go to the shared core table and which to a profile table.
"""

import dataclasses


@dataclasses.dataclass(frozen=True)
class Field:
    """
    The column or columns one source fills.

    Attributes:
        columns (tuple): The column names, in the order `get` returns their values.
        get (callable): `get(obj, context)`: the value of a one-column field, or a tuple of `len(columns)` values.
        load (callable or None): `load(row)`: the dict of `from_dict` keys the stored columns fill; `None` for a
            column that is not part of the object (the owning ids).
    """

    columns: tuple
    get: object
    load: object = None


def attr(column, name=None):
    """A column filled from the object's attribute `name` (the column's own name when not given)."""
    attribute = name or column
    return Field((column,), lambda obj, context: getattr(obj, attribute), lambda row: {attribute: row[column]})


def optional(column, name=None):
    """As `attr`, `None` for an object that lacks the attribute."""
    attribute = name or column
    return Field((column,), lambda obj, context: getattr(obj, attribute, None), lambda row: {attribute: row[column]})


def context(column, key=None):
    """A column filled from the context the caller passed (`context[key]`, the column's name by default)."""
    name = key or column
    return Field((column,), lambda obj, context: context[name])


def scaled(column, name, factor, null=None, nullable=False):
    """
    A column holding the attribute `name` times `factor` (AU stored as km): the load divides it back.

    With `nullable` a `None` attribute is stored as NULL and a NULL loads as `null` (the object's own default).
    """
    def get(obj, context):
        value = getattr(obj, name)
        return None if nullable and value is None else value * factor

    def load(row):
        value = row[column]
        return {name: null if nullable and value is None else value / factor}

    return Field((column,), get, load)


def each(names):
    """Columns filled from the attributes of the same names, `None` for one the object lacks."""
    return Field(tuple(names), lambda obj, context: tuple(getattr(obj, name, None) for name in names),
                 lambda row: {name: row[name] for name in names})


def computed(columns, function, load=None):
    """
    Columns filled by `function(obj, context)`: one value for a column name, a tuple for a tuple of them.

    `load(row)` returns the `from_dict` keys they fill; without it the field is write-only.
    """
    if isinstance(columns, str):
        return Field((columns,), function, load)
    return Field(tuple(columns), function, load)


def columns_of(fields):
    """Every column the list fills, in order."""
    return tuple(column for field in fields for column in field.columns)


def values_of(fields, obj, context=None):
    """The values of the list for `obj`, in the order of `columns_of`."""
    context = context or {}
    out = []
    for field in fields:
        value = field.get(obj, context)
        if len(field.columns) == 1:
            out.append(value)
        else:
            value = tuple(value)
            if len(value) != len(field.columns):
                raise ValueError(f"{field.columns} got {len(value)} values")
            out.extend(value)
    return out


def load_dict(fields, row, only=None):
    """
    The `from_dict` dict `fields` read back from the stored `row` (fields with no `load` add nothing).

    `only` limits it to the fields whose first column is named: the satellite row of an anchored remnant holds a
    few columns of its own, the rest of the object being on its `stars` row.
    """
    data = {}
    for field in fields:
        if field.load is not None and (only is None or field.columns[0] in only):
            data.update(field.load(row))
    return data


def insert_row(conn, table, fields, obj, context=None):
    """
    Writes `obj` to `table` as one row filled by `fields`; returns the new row's id.

    The statement is a plain `INSERT INTO <table> (...) VALUES (...)` naming no `id`, so a table that takes its ids
    from the id blocks (`store.ID_BLOCK_TABLES`) still does.
    """
    columns = columns_of(fields)
    cur = conn.execute(f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({', '.join('?' * len(columns))})",
                       tuple(values_of(fields, obj, context)))
    return cur.lastrowid
