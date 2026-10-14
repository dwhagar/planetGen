"""
DB.23, GEN.201: `phenomenon_scatter.map_luminosity_sol`, the luminosity a star of a row's mass would have on the
Galaxy Map (GEN.200). A scatter writes it with each row, and `plan --fill-map-luminosity` fills it for a galaxy
scattered before it, to the same values.
"""

import math

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.generation import run_plan
from planetgen.generation import phenomenon_scatter as scatter
from planetgen.physics import stellar_evolution

from tests.test_phenomenon_scatter import _plan_args, _seed_galaxy

pytestmark = pytest.mark.db


def _rows(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute("SELECT id, kind, mass_solar, map_luminosity_sol FROM phenomenon_scatter ORDER BY id"
                            ).fetchall()
    finally:
        conn.close()


def _agree(stored, mass):
    return math.isclose(stored, stellar_evolution.map_luminosity_sol(mass), rel_tol=1e-5)


def test_the_columns_include_the_map_luminosity():
    assert scatter.PHENOMENON_SCATTER_COLUMNS[-1] == "map_luminosity_sol" == store.PHENOMENON_SCATTER_COLUMNS[-1]


def test_a_scatter_row_carries_the_luminosity_of_its_mass():
    rows = [dict(zip(scatter.PHENOMENON_SCATTER_COLUMNS, row)) for row in scatter.special_rows(
        [(0, 6)], 4.0, 5)] + [dict(zip(scatter.PHENOMENON_SCATTER_COLUMNS, scatter.nucleus_row(seed)))
                              for seed in range(30)]
    assert rows
    for row in rows:
        if row["mass_solar"] is None:
            assert row["map_luminosity_sol"] is None
        else:
            assert row["map_luminosity_sol"] == stellar_evolution.map_luminosity_sol(row["mass_solar"])


def test_the_sql_expression_gives_the_python_function_across_every_piece():
    pieces = [piece["max_mass_sol"] for piece in tuning.MAP_MASS_LUMINOSITY_PIECES if piece["max_mass_sol"]]
    masses = [0.1, 0.43, 1.0, 1.999, 2.0, 5.0, 19.99, 20.0, 1e3, 4e6, 1e8] + [m * (1 - 1e-9) for m in pieces]
    expression = run_plan.map_luminosity_sql("m")
    for mass in masses:
        assert math.isclose(_sql_value(expression, mass), stellar_evolution.map_luminosity_sol(mass),
                            rel_tol=1e-9), mass


@pytest.fixture
def sql_conn(mysql_config):
    conn = store.get_connection(mysql_config)
    yield conn
    conn.close()


def _sql_value(expression, mass, _cache={}):
    conn = _cache["conn"]
    return conn.execute(f"SELECT {expression} AS v FROM (SELECT {mass!r} AS m) t").fetchone()["v"]


@pytest.fixture(autouse=True)
def _sql_connection(sql_conn):
    _sql_value.__defaults__[0]["conn"] = sql_conn


def test_a_plan_writes_it_and_the_fill_gives_an_old_galaxy_the_same_values(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_plan_args(mysql_config))
    written = _rows(mysql_config)
    massive = [row for row in written if row["mass_solar"] is not None]
    assert len(massive) > 50
    assert all(_agree(row["map_luminosity_sol"], row["mass_solar"]) for row in massive)
    assert all(row["map_luminosity_sol"] is None for row in written if row["mass_solar"] is None)

    # A galaxy scattered before v84: no values. The fill restores them, in small batches, and only once.
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE phenomenon_scatter SET map_luminosity_sol = NULL")
        conn.commit()
    finally:
        conn.close()
    args = _plan_args(mysql_config)
    args.fill_map_luminosity = True
    original = run_plan.MAP_LUMINOSITY_BATCH_IDS
    run_plan.MAP_LUMINOSITY_BATCH_IDS = 1000
    try:
        assert run_plan.fill_map_luminosity(args) == len(massive)
        assert run_plan.fill_map_luminosity(args) == 0
    finally:
        run_plan.MAP_LUMINOSITY_BATCH_IDS = original
    again = _rows(mysql_config)
    for before, after in zip(written, again):
        if before["mass_solar"] is None:
            assert after["map_luminosity_sol"] is None
        else:
            assert math.isclose(after["map_luminosity_sol"], before["map_luminosity_sol"], rel_tol=1e-5)


def test_a_stopped_fill_carries_on_where_it_was(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_plan_args(mysql_config))
    conn = store.get_connection(mysql_config)
    try:
        half = conn.execute("SELECT MAX(id) DIV 2 AS n FROM phenomenon_scatter").fetchone()["n"]
        conn.execute("UPDATE phenomenon_scatter SET map_luminosity_sol = NULL WHERE id > ?", (half,))
        conn.commit()
    finally:
        conn.close()
    args = _plan_args(mysql_config)
    args.fill_map_luminosity = True
    filled = run_plan.fill_map_luminosity(args)
    assert 0 < filled
    assert all(row["map_luminosity_sol"] is not None for row in _rows(mysql_config) if row["mass_solar"] is not None)
