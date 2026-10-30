# tests/test_scatter_member_clear.py

"""
PERF.180: emptying a scatter removes its members and nodes without one statement over millions of rows: the tables are truncated when
the other scatter has none, else the rows go in batches that each commit, and the other scatter's rows are never touched.
"""

import pytest

from planetgen.db import scatter_systems, store


def _put(conn, source, count):
    for number in range(count):
        conn.execute("INSERT INTO scatter_members (source, source_id, member_no, parent_node_no, slot, kind, mass_solar, initial_mass_solar, seed)"
                     " VALUES (?, ?, 0, 0, 0, 'main-sequence', 1.0, 1.0, 1)", (source, number + 1))
        conn.execute("INSERT INTO scatter_nodes (source, source_id, node_no, mass_solar, separation_km, eccentricity, inclination_deg,"
                     " ascending_node_deg, arg_periapsis_deg, mean_anomaly_deg, epoch_unix, period_years, kind)"
                     " VALUES (?, ?, 0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 'close')", (source, number + 1))
    conn.commit()


def _held(conn):
    return {(table, row["source"]): row["n"] for table in ("scatter_members", "scatter_nodes")
            for row in conn.execute(f"SELECT source, COUNT(*) AS n FROM {table} GROUP BY source").fetchall()}


@pytest.fixture
def conn(mysql_config):
    connection = store.get_connection(mysql_config)
    for table in ("scatter_members", "scatter_nodes"):
        connection.execute(f"TRUNCATE TABLE {table}")
    yield connection
    connection.close()


def test_the_other_scatters_members_are_kept_and_the_rows_go_in_batches(conn, monkeypatch):
    monkeypatch.setattr(scatter_systems, "DELETE_BATCH", 2)
    _put(conn, 1, 5)
    _put(conn, 2, 3)
    scatter_systems.delete_all_member_trees(conn, "bright_stars")
    assert _held(conn) == {("scatter_members", 2): 3, ("scatter_nodes", 2): 3}


def test_with_no_other_scatter_the_tables_are_emptied(conn):
    _put(conn, 1, 4)
    scatter_systems.delete_all_member_trees(conn, "bright_stars")
    assert _held(conn) == {}
    _put(conn, 2, 2)
    scatter_systems.delete_all_member_trees(conn, "phenomenon_scatter")
    assert _held(conn) == {}
    scatter_systems.delete_all_member_trees(conn, "phenomenon_scatter")     # nothing to remove is fine
