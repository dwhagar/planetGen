# planetgen/db/region.py

"""
Everything inside a sphere, read by the cell key (NAV.60)
=========================================================

`region_rows` answers "which rows of this table lie within `radius` parsecs of a point?" from the table's `cell_key`
column and its covering index `(cell_key, x, y, z)` (design note `docs/design/galaxy/spherical-region-queries.md`, sections 6.3
and 6.4), for any table with a position, whether or not the row has a sector:

    SELECT p.id, ..., POW(x - cx, 2) + POW(y - cy, 2) + POW(z - cz, 2) AS d2
    FROM (VALUES ROW(lo1, hi1), ROW(lo2, hi2), ...) AS rg(lo, hi)        -- at most 256 runs of keys (GEN.276)
    STRAIGHT_JOIN table p FORCE INDEX (idx_table_cell) ON p.cell_key BETWEEN rg.lo AND rg.hi
    HAVING d2 <= r^2
    UNION ALL
    SELECT ... FROM table p FORCE INDEX (idx_table_cell) WHERE p.cell_key IS NULL AND ... HAVING d2 <= r^2     -- a row outside the box has no key

The squared distance is worked out once per candidate row, as `POW(axis - c, 2)` (PERF.216): MySQL does not share a repeated
sub-expression, so `(axis - c) * (axis - c)` is worked out twice and the test and the select list worked it out again;
at 50 pc, where about 90,000 candidate rows are tested, that was most of the statement's time (the index read is 0.5 us a row, the
test was 2.7).

The runs are written into the statement as integers; the point and the radius as finite numbers. The SQL test is a
little loose (a billionth), and the caller measures the exact distance with `detmath`, so the answer does not depend on
how MySQL rounds a sum.

`capped_count` answers "how many, up to a cap?" from the same statement, stopped early (NAV.61): the count of a sphere holding a
million rows costs about what the count of one holding the cap does.

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


def region_statement(table, columns, position, center, radius, where="", runs=None, with_distance=True):
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
        with_distance (bool): Whether the squared distance is selected as the last column `d2` (the default) or only tested, in
            `WHERE`, for a statement that counts rows and needs no distance.

    Returns:
        tuple: `(sql, runs)`; the SQL holds the point and the radius as numbers and `?` only where `where` does. A
        sphere wholly outside the box has no runs, and the statement is the NULL-key branch alone.
    """
    if runs is None:
        runs = cellkey.sphere_runs(center, radius)
    reach = float(radius) * (1.0 + SLACK) + SLACK
    squares = " + ".join(f"POW({axis} - {_number(c)}, 2)" for axis, c in zip(position, center))
    limit = _number(reach * reach)
    selected = f"{columns}, {squares} AS d2" if with_distance else columns
    extra = f"({where})" if where else ""
    if with_distance:
        join_where, join_tail = (f" WHERE {extra}" if extra else ""), f" HAVING d2 <= {limit}"
        null_where, null_tail = f" AND {extra}" if extra else "", join_tail
    else:
        join_where, join_tail = f" WHERE {extra + ' AND ' if extra else ''}{squares} <= {limit}", ""
        null_where, null_tail = f"{' AND ' + extra if extra else ''} AND {squares} <= {limit}", ""
    null_branch = (
        f"SELECT {selected} FROM {table} p FORCE INDEX ({index_name(table)}) "
        f"WHERE p.cell_key IS NULL AND {position[0]} IS NOT NULL{null_where}{null_tail}"
    )
    if not runs:        # the sphere lies wholly outside the box: only the rows without a key can be in it
        return null_branch, runs
    sql = (
        f"SELECT {selected} "
        f"FROM ({values_list(runs)}) AS rg(lo, hi) "
        f"STRAIGHT_JOIN {table} p FORCE INDEX ({index_name(table)}) ON p.cell_key BETWEEN rg.lo AND rg.hi"
        f"{join_where}{join_tail} "
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


def capped_count(conn, table, position, center, radius, cap, where="", where_params=()):
    """
    How many rows of `table` lie within `radius` of `center`, counting no further than `cap` (NAV.61).

    The statement is `region_statement` stopped by `LIMIT cap` inside a `COUNT(*)`, so a sphere holding a million rows is
    counted as fast as one holding `cap`: 11 ms at most on the design note's measurements. The SQL test is the loose one
    (a billionth), so a row on the surface may count; a count is a guard and a label, not a list.

    Args:
        cap (int): The most to count, 1 or more.

    Returns:
        int: The count, at most `cap`; `cap` means "`cap` or more".
    """
    cap = int(cap)
    if cap < 1:
        raise RegionError("a count is capped at 1 or more")
    sql, runs = region_statement(table, "1 AS one", position, center, radius, where, with_distance=False)
    row = conn.execute(f"SELECT COUNT(*) AS n FROM ({sql} LIMIT {cap}) AS counted",
                       tuple(where_params) * (2 if runs else 1)).fetchone()
    return int(row["n"])
