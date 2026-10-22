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

The same lists drive the load (DB.74, next slice) and, in the 8.4 revision, name which columns go to the shared core
table and which to a profile table.
"""

import dataclasses


@dataclasses.dataclass(frozen=True)
class Field:
    """
    The column or columns one source fills.

    Attributes:
        columns (tuple): The column names, in the order `get` returns their values.
        get (callable): `get(obj, context)`: the value of a one-column field, or a tuple of `len(columns)` values.
    """

    columns: tuple
    get: object


def attr(column, name=None):
    """A column filled from the object's attribute `name` (the column's own name when not given)."""
    attribute = name or column
    return Field((column,), lambda obj, context: getattr(obj, attribute))


def optional(column, name=None):
    """As `attr`, `None` for an object that lacks the attribute."""
    attribute = name or column
    return Field((column,), lambda obj, context: getattr(obj, attribute, None))


def context(column, key=None):
    """A column filled from the context the caller passed (`context[key]`, the column's name by default)."""
    name = key or column
    return Field((column,), lambda obj, context: context[name])


def computed(columns, function):
    """Columns filled by `function(obj, context)`: one value for a column name, a tuple for a tuple of them."""
    if isinstance(columns, str):
        return Field((columns,), function)
    return Field(tuple(columns), function)


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
