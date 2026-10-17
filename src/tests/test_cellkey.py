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


# --- sphere_runs (GEN.276) -----------------------------------------------------------------------------------

def _in_runs(runs, keys):
    lows = np.array([low for low, high in runs], dtype=np.uint64)
    highs = np.array([high for low, high in runs], dtype=np.uint64)
    where = np.searchsorted(lows, keys, side="right") - 1
    return (where >= 0) & (keys <= highs[np.maximum(where, 0)])


def test_runs_are_ascending_merged_and_within_the_budget():
    rng = random.Random(11)
    for _ in range(60):
        center = (rng.uniform(-9000, 9000), rng.uniform(-60, 60), rng.uniform(-30, 30))
        r = rng.choice((0.0, 0.2, 1.0, 3.7, 25.0, 50.0, 400.0))
        budget = rng.choice((1, 4, 16, 64, 256))
        runs = cellkey.sphere_runs(center, r, budget)
        assert 1 <= len(runs) <= budget
        assert all(low <= high < 2 ** cellkey.KEY_BITS for low, high in runs)
        assert all(a[1] + 1 < b[0] for a, b in zip(runs, runs[1:])), "ascending and not touching"


def test_every_point_inside_the_sphere_has_its_key_in_a_run():
    rng = np.random.default_rng(5)
    for center, r in (((8200.0, 12.0, 3.0), 50.0), ((8200.0, 12.0, 3.0), 0.3), ((0.0, 0.0, 0.0), 7.0),
                      ((-1234.5, 300.25, -77.0), 1.0), ((262000.0, 0.0, 100.0), 500.0), ((1e5, -2e5, 1.5e5), 4000.0)):
        runs = cellkey.sphere_runs(center, r)
        direction = rng.normal(size=(20000, 3))
        direction /= np.linalg.norm(direction, axis=1)[:, None]
        radius = r * rng.random(20000) ** (1 / 3)
        points = np.array(center) + direction * radius[:, None]
        points = np.vstack([points, np.array(center) + direction[:2000] * r])           # on the surface too
        points = points[np.all(np.abs(points) < cellkey.REACH_PC, axis=1)]
        assert len(points) > 0
        assert _in_runs(runs, cellkey.cell_keys(points)[0]).all()


def test_the_runs_are_tight_not_the_whole_box():
    runs = cellkey.sphere_runs((8200.0, 0.0, 0.0), 50.0)
    cells = sum(high - low + 1 for low, high in runs)
    ball = 4 / 3 * math.pi * (50.0 / cellkey.CELL_PC) ** 3
    assert ball <= cells < 3 * ball


def test_a_smaller_budget_gives_fewer_and_larger_runs_that_still_cover():
    center, r = (8200.0, 12.0, 3.0), 50.0
    fine, coarse = cellkey.sphere_runs(center, r, 256), cellkey.sphere_runs(center, r, 8)
    assert len(coarse) <= 8 < len(fine)
    assert sum(h - l for l, h in coarse) >= sum(h - l for l, h in fine)
    points = np.array(center) + np.random.default_rng(3).uniform(-r / 2, r / 2, size=(5000, 3))
    assert _in_runs(coarse, cellkey.cell_keys(points)[0]).all()


def test_the_runs_of_a_sphere_beyond_the_box_are_only_the_part_inside():
    assert cellkey.sphere_runs((cellkey.REACH_PC + 100.0, 0.0, 0.0), 10.0) == []
    assert cellkey.sphere_runs((1e9, 1e9, 1e9), 5.0) == []
    runs = cellkey.sphere_runs((cellkey.REACH_PC - 1.0, 0.0, 0.0), 10.0)
    inside = np.array([[cellkey.REACH_PC - 5.0, 0.0, 0.0], [cellkey.REACH_PC - 0.5, 3.0, 3.0]])
    assert _in_runs(runs, cellkey.cell_keys(inside)[0]).all()
    whole = cellkey.sphere_runs((0.0, 0.0, 0.0), 1e7)
    assert whole == [(0, 2 ** cellkey.KEY_BITS - 1)]


def test_a_radius_of_zero_gives_the_cell_of_the_point():
    center = (8200.13, -4.4, 2.01)
    key = cellkey.cell_key(*center)
    runs = cellkey.sphere_runs(center, 0.0)
    assert any(low <= key <= high for low, high in runs) and sum(h - l + 1 for l, h in runs) <= 8


def test_another_reach_scales_the_cell_with_it():
    reach = 1024.0
    center, r = (100.0, -20.0, 5.0), 3.0
    runs = cellkey.sphere_runs(center, r, reach_pc=reach)
    rng = np.random.default_rng(9)
    points = np.array(center) + rng.uniform(-r / 2, r / 2, size=(2000, 3))
    cell = 2 * reach / 2 ** cellkey.BITS
    index = np.floor((points + reach) / cell).astype(np.uint64)
    keys = (cellkey._weave(index[:, 0], index[:, 1], index[:, 2]))
    assert _in_runs(runs, keys).all()


@pytest.mark.parametrize("bad", [((0, 0), 1.0), ((0, 0, 0), -1.0), ((0, 0, float("nan")), 1.0), ((0, 0, 0), float("inf"))])
def test_nonsense_is_refused(bad):
    with pytest.raises(ValueError):
        cellkey.sphere_runs(*bad)
    with pytest.raises(ValueError):
        cellkey.sphere_runs((0, 0, 0), 1.0, budget=0)
