# tests/test_migration_galaxy_containers.py

"""
DB.85: migration 0130 adds the empty `galaxy_containers` table; its CHECKs and keys refuse a malformed row, and a
cluster's row goes with its component.
"""

import pytest

from planetgen.db import store
from planetgen.generation import clusters
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

_UID = "UNHEX('000000000000000000000001')"


def _insert(conn, kind="superbubble", extent=100.0, uid=_UID, bubble="NULL, NULL, NULL"):
    conn.execute(
        "INSERT INTO galaxy_containers (uid, kind, bubble_col, bubble_row, bubble_n, center_x_pc, center_y_pc,"
        f" center_z_pc, extent_pc, shell_index, band_index, wedge_index) VALUES ({uid}, '{kind}', {bubble}, 0, 0, 0, {extent},"
        " 0, 0, 0)")


def test_a_v129_galaxy_gets_an_empty_containers_table(mysql_config):
    load_old_schema(mysql_config, 129)
    assert store.schema_status(mysql_config)[0] == 129
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM galaxy_containers").fetchone()["n"] == 0
    finally:
        conn.close()


def test_the_checks_and_keys_refuse_a_malformed_row(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        _insert(conn)
        for bad in (lambda: _insert(conn, kind="galaxy", uid="UNHEX('000000000000000000000002')"),
                    lambda: _insert(conn, extent=0.5, uid="UNHEX('000000000000000000000003')"),
                    lambda: _insert(conn)):
            with pytest.raises(Exception):
                bad()
        _insert(conn, uid="UNHEX('000000000000000000000004')", bubble="1, 2, 0")
        with pytest.raises(Exception):
            _insert(conn, uid="UNHEX('000000000000000000000005')", bubble="1, 2, 0")
    finally:
        conn.rollback()
        conn.close()


def test_a_cluster_row_goes_with_its_component(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, clusters.generate_clusters(7, count=2), kind="cluster")
        component = store.get_galaxy_components(conn, kind="cluster")[0].id
        conn.execute(
            "INSERT INTO galaxy_containers (uid, kind, name, component_id, center_x_pc, center_y_pc, center_z_pc,"
            " extent_pc, shell_index, band_index, wedge_index) VALUES (UNHEX('000000000000000000000009'),"
            " 'globular_cluster', 'NGC', ?, 0, 0, 0, 50, 0, 0, 0)", (component,))
        conn.execute("DELETE FROM galaxy_components WHERE id = ?", (component,))
        assert conn.execute("SELECT COUNT(*) AS n FROM galaxy_containers").fetchone()["n"] == 0
    finally:
        conn.rollback()
        conn.close()
