# planetgen/db/region.py

"""
Everything inside a sphere, read by the cell key (NAV.60)
=========================================================

`region_rows` answers "which rows of this table lie within `radius` parsecs of a point?" from the table's `cell_key`
column and its covering index `(cell_key, x, y, z)` (design note `docs/design/galaxy/spherical-region-queries.md`, sections 6.3
and 6.4), for any table with a position, whether or not the row has a sector:

    SELECT p.id, ..., d2
    FROM (VALUES ROW(lo1, hi1), ROW(lo2, hi2), ...) AS rg(lo, hi)        -- at most 256 runs of keys (GEN.276)
    STRAIGHT_JOIN table p FORCE INDEX (idx_table_cell) ON p.cell_key BETWEEN rg.lo AND rg.hi
    WHERE (x - cx)^2 + (y - cy)^2 + (z - cz)^2 <= r^2
    UNION ALL
    SELECT ... FROM table p FORCE INDEX (idx_table_cell) WHERE p.cell_key IS NULL AND <the same test>     -- a row outside the box has no key

The runs are written into the statement as integers; the point and the radius as finite numbers. The SQL test is a
little loose (a billionth), and the caller measures the exact distance with `detmath`, so the answer does not depend on
how MySQL rounds a sum.

The `sectors` table is never read.
"""

import math

from planetgen.galaxy import cellkey

SLACK = 1e-9
"""float: How much fatter than `radius` the SQL test is, relative and absolute (parsecs), so a row on the surface is
never lost to a different rounding of the sum of squares; the caller does the exact test."""


class RegionError(ValueError):
    """A point, radius or table that cannot be searched."""


def index_name(table):
    """The covering index every table with a cell key carries: `idx_<table>_cell`."""
    return f"idx_{table}_cell"


def _number(value):
    """A finite number as SQL text (a float, never a string from outside)."""
    value = float(value)
    if not math.isfinite(value):
        raise RegionError("the point and the radius are finite numbers")
    return repr(value)


def values_list(runs):
    """The `VALUES ROW(lo, hi), ...` text of a list of key runs, as integers."""
    return "VALUES " + ", ".join(f"ROW({int(low)}, {int(high)})" for low, high in runs)


def region_statement(table, columns, position, center, radius, where="", runs=None):
    """
    The statement and its parameters for the rows of `table` within `radius` of `center`.

    Args:
        table (str): A table with a `cell_key` column and `idx_<table>_cell`.
        columns (str): What to select, in terms of the alias `p` (`"p.id, p.name"`); the squared distance follows as
            the last column `d2`.
        position (tuple): The three SQL expressions of the position in parsecs, in terms of `p`
            (`("p.center_x_pc", "p.center_y_pc", "p.center_z_pc")`, or `"p.position_x_mpc / 1000.0"`).
        center (sequence): The point `(x, y, z)`, parsecs, galaxy frame.
        radius (float): Parsecs, zero or more.
        where (str): More conditions on `p`, with no `?` (`region_rows` takes `?` values through `where_params`; the
            text stands in both branches).
        runs (list, optional): The key runs; `sphere_runs(center, radius)` by default.

    Returns:
        tuple: `(sql, runs)`; the SQL holds the point and the radius as numbers and `?` only where `where` does. A
        sphere wholly outside the box has no runs, and the statement is the NULL-key branch alone.
    """
    if runs is None:
        runs = cellkey.sphere_runs(center, radius)
    reach = float(radius) * (1.0 + SLACK) + SLACK
    squares = " + ".join(f"({axis} - {_number(c)}) * ({axis} - {_number(c)})" for axis, c in zip(position, center))
    test = f"{squares} <= {_number(reach * reach)}"
    extra = f" AND ({where})" if where else ""
    null_branch = (
        f"SELECT {columns}, {squares} AS d2 FROM {table} p FORCE INDEX ({index_name(table)}) "
        f"WHERE p.cell_key IS NULL AND {position[0]} IS NOT NULL AND {test}{extra}"
    )
    if not runs:        # the sphere lies wholly outside the box: only the rows without a key can be in it
        return null_branch, runs
    sql = (
        f"SELECT {columns}, {squares} AS d2 "
        f"FROM ({values_list(runs)}) AS rg(lo, hi) "
        f"STRAIGHT_JOIN {table} p FORCE INDEX ({index_name(table)}) ON p.cell_key BETWEEN rg.lo AND rg.hi "
        f"WHERE {test}{extra} "
        f"UNION ALL {null_branch}"
    )
    return sql, runs


def region_rows(conn, table, columns, position, center, radius, where="", where_params=(), runs=None):
    """
    The rows of `table` within `radius` of `center`, as `region_statement` reads them.

    Returns:
        list: The rows (mappings with the selected columns and `d2`), unordered. The SQL test is a little loose: the
        caller measures the exact distance and keeps the rows within the radius.
    """
    sql, runs = region_statement(table, columns, position, center, radius, where, runs)
    return conn.execute(sql, tuple(where_params) * (2 if runs else 1)).fetchall()
