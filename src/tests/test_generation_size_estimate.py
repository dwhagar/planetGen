# tests/test_generation_size_estimate.py

"""
PERF.26 The size estimate is off (docs/TODO.md): what a real fill adds to
the galaxy database, table by table, against the size per system the
estimate uses before anything was measured and the one it measures after
a run.

Needs a MySQL server (the `mysql_config` fixture); skips without one.
"""

import pytest

from planetgen.db import store
from planetgen.generation import stats as generationStats
from tests.test_galaxy_gen import _mysql_argv, _plan_wide_galaxy, _run_cli

MARGIN = 1.5
"""float: The estimate may be off by at most this factor either way."""


def _table_bytes(conn):
    """Every table's data and index bytes, with fresh statistics."""
    tables = [row["t"] for row in conn.execute(
        "SELECT table_name AS t FROM information_schema.tables"
        " WHERE table_schema = DATABASE() AND table_type = 'BASE TABLE'").fetchall()]
    conn.execute("ANALYZE TABLE " + ", ".join(tables)).fetchall()
    try:
        conn.execute("SET SESSION information_schema_stats_expiry = 0")
    except Exception:  # noqa: BLE001 -- MariaDB reads them live
        pass
    rows = conn.execute(
        "SELECT table_name AS t, data_length + index_length AS b FROM information_schema.tables"
        " WHERE table_schema = DATABASE() AND table_type = 'BASE TABLE'").fetchall()
    return {row["t"]: int(row["b"] or 0) for row in rows}


@pytest.fixture
def filled(mysql_config):
    """A planned galaxy before and after a fill of 12 sectors of 40
    systems each."""
    _plan_wide_galaxy(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        before = _table_bytes(conn)
        _run_cli(["--ring", "40", "--layer", "0", "--limit", "12", "--num-systems", "40", "--yes"]
                 + _mysql_argv(mysql_config))
        after = _table_bytes(conn)
        systems = conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"]
        measured = generationStats.GenerationStats(None).measure_size(conn, mysql_config.database)
    finally:
        conn.close()
    growth = sum(after.values()) - sum(before.values())
    return {"growth_per_system": growth / systems, "systems": systems, "after": after, "measured": measured}


def test_the_default_size_per_system_matches_a_real_fill(filled):
    assert filled["systems"] == 480
    estimate = generationStats.DEFAULT_BYTES_PER_SYSTEM * (1 + generationStats.SIZE_MARGIN)
    ratio = estimate / filled["growth_per_system"]
    assert 1 / MARGIN <= ratio <= MARGIN, (estimate, filled["growth_per_system"])


def test_the_measured_size_per_system_matches_what_the_fill_added(filled):
    measured = filled["measured"]
    assert measured["systems"] == filled["systems"]
    ratio = measured["bytes_per_system"] / filled["growth_per_system"]
    assert 1 / MARGIN <= ratio <= MARGIN, (measured, filled["growth_per_system"])


def test_the_galaxy_wide_tables_are_left_out(filled):
    kept = sum(size for table, size in filled["after"].items()
               if table not in generationStats.GALAXY_WIDE_TABLES)
    assert filled["measured"]["total_bytes"] == kept
    assert generationStats.GALAXY_WIDE_TABLES <= set(filled["after"])
