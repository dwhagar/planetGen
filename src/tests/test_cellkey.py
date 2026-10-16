# tests/test_cellkey.py

"""The cell key (region queries): MySQL's generated column equals the Python key on every point, edges and outside included."""

import math
import random

import numpy as np
import pytest

from planetgen.db import store
from planetgen.galaxy import cellkey


def _points(count=20000):
    rng = random.Random(7)
    reach, cell = cellkey.REACH_PC, cellkey.CELL_PC
    points = [(rng.uniform(-reach, reach), rng.uniform(-reach, reach), rng.uniform(-reach, reach)) for _ in range(count)]
    # Near the Sun, where a real galaxy's cubes are small compared with its sectors.
    points += [(8200 + rng.uniform(-60, 60), rng.uniform(-60, 60), rng.uniform(-30, 30)) for _ in range(count // 4)]
    edges = (0.0, cell, -cell, cell / 2, 1e-9, -1e-9, reach - cell, reach - 1e-6, reach - cell / 2, -reach + 1e-6,
             -reach + cell, -reach)
    points += [(a, b, c) for a in edges for b in (0.0, edges[2], reach - 1e-6) for c in edges[:6]]
    # On and just outside the box.
    points += [(reach, 0.0, 0.0), (0.0, -reach, 0.0), (0.0, 0.0, reach + 1.0), (1e9, 0.0, 0.0), (-reach - 1e-6, 1.0, 1.0)]
    return points


def test_the_constants_make_one_cube_per_integer_of_21_bits():
    assert 2 * cellkey.REACH_PC / cellkey.CELL_PC == 2 ** cellkey.BITS
    assert cellkey.KEY_BITS == 63
    assert cellkey.cell_key(-cellkey.REACH_PC + 1e-6, -cellkey.REACH_PC + 1e-6, -cellkey.REACH_PC + 1e-6) == 0
    assert cellkey.cell_key(cellkey.REACH_PC - 1e-6, cellkey.REACH_PC - 1e-6, cellkey.REACH_PC - 1e-6) == 2 ** 63 - 1
    assert cellkey.cell_key(cellkey.REACH_PC, 0.0, 0.0) is None and cellkey.cell_key(0.0, 0.0, -1e9) is None


def test_spread_puts_bit_i_at_bit_3i():
    for bit in range(cellkey.BITS):
        assert cellkey.spread(1 << bit) == 1 << (3 * bit)
    assert cellkey.spread(2 ** 21 - 1) == int("1" + "001" * 20, 2)


def test_the_vector_key_equals_the_scalar_key():
    points = _points()
    keys, valid = cellkey.cell_keys(points)
    for point, key, ok in zip(points, keys, valid):
        expected = cellkey.cell_key(*point)
        assert (expected is not None) == bool(ok)
        assert expected is None or expected == int(key)


def test_nearby_cubes_share_a_run():
    """The 8 cubes of an aligned 2x2x2 block have consecutive keys: what lets a sphere become a list of runs."""
    origin = (10.0, -20.0, 30.0)
    half = cellkey.CELL_PC
    base = cellkey.cell_index(origin[0]), cellkey.cell_index(origin[1]), cellkey.cell_index(origin[2])
    even = tuple(c - c % 2 for c in base)
    keys = sorted(cellkey.cell_key(*(((e + d) * cellkey.CELL_PC - cellkey.REACH_PC + half / 2) for e, d in zip(even, bits)))
                  for bits in ((i, j, k) for i in (0, 1) for j in (0, 1) for k in (0, 1)))
    assert keys == list(range(keys[0], keys[0] + 8)) and keys[0] % 8 == 0


@pytest.mark.parametrize("unit", ["pc", "mpc"])
def test_mysqls_generated_column_equals_the_python_key(mysql_config, unit):
    """The same expression the schema uses, on a scratch table, against `cell_key` for the whole point list."""
    points = _points()
    columns = ("x", "y", "z") if unit == "pc" else ("x_mpc", "y_mpc", "z_mpc")
    expression = cellkey.sql_cell_key(*(c if unit == "pc" else f"{c} / 1000.0" for c in columns))
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DROP TABLE IF EXISTS tmp_cell_keys")
        conn.execute(f"CREATE TABLE tmp_cell_keys (id INT PRIMARY KEY, {columns[0]} DOUBLE NOT NULL, {columns[1]} DOUBLE NOT NULL,"
                     f" {columns[2]} DOUBLE NOT NULL, cell_key BIGINT UNSIGNED GENERATED ALWAYS AS ({expression}) STORED,"
                     f" KEY idx_cell (cell_key, {columns[0]}, {columns[1]}, {columns[2]}))")
        scale = 1.0 if unit == "pc" else 1000.0
        conn.executemany(f"INSERT INTO tmp_cell_keys (id, {', '.join(columns)}) VALUES (?, ?, ?, ?)",
                         [(i, *(c * scale for c in p)) for i, p in enumerate(points)])
        rows = conn.execute(f"SELECT id, {', '.join(columns)}, cell_key FROM tmp_cell_keys ORDER BY id").fetchall()
        conn.execute("DROP TABLE tmp_cell_keys")
    finally:
        conn.close()
    assert len(rows) == len(points)
    for row in rows:
        stored = tuple(row[c] / 1000.0 if unit == "mpc" else row[c] for c in columns)
        expected = cellkey.cell_key(*stored)
        assert row["cell_key"] == expected, (row["id"], stored)
