# tests/test_migration_state_vectors.py

"""
GEN.289: migration 0109 drops `orbital_phase_deg` from `planets` and `moons`, points a circumbinary planet at its
system, and can be run again; new rows carry their parent.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema
from tests.test_db_persistence import _make_system_with_moons_and_belt
from tests.test_migration_point_columns import _columns, _insert

pytestmark = pytest.mark.fresh_database


def test_a_v108_database_loses_the_phase_columns_and_keeps_its_rows(mysql_config):
    load_old_schema(mysql_config, 108)
    assert store.schema_status(mysql_config)[0] == 108
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert "orbital_phase_deg" in _columns(conn, "planets") | _columns(conn, "moons")
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        _insert(conn, "planets", body_type="t", star_system_id=7, star_id=None, orbital_phase_deg=33.0, parent_type="star", parent_id=None)
        _insert(conn, "planets", body_type="t", star_system_id=7, star_id=4, orbital_phase_deg=44.0, parent_type="star", parent_id=4)
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert "orbital_phase_deg" not in _columns(conn, "planets")
        assert "orbital_phase_deg" not in _columns(conn, "moons")
        rows = conn.execute("SELECT star_id, parent_type, parent_id FROM planets ORDER BY id").fetchall()
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()
    assert (rows[0]["parent_type"], rows[0]["parent_id"]) == ("system", 7)  # a circumbinary planet goes round its system
    assert (rows[1]["parent_type"], rows[1]["parent_id"]) == ("star", 4)


def test_a_migration_run_again_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 108)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0108')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION


def test_a_saved_planet_and_moon_carry_their_parent(mysql_config):
    system, cfg = _make_system_with_moons_and_belt()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        planets = conn.execute("SELECT id, star_id, star_system_id, parent_type, parent_id FROM planets"
                               " WHERE star_system_id = ?", (system_id,)).fetchall()
        moons = conn.execute("SELECT planet_id, parent_type, parent_id FROM moons WHERE star_system_id = ?",
                             (system_id,)).fetchall()
    finally:
        conn.close()
    assert planets and moons
    for row in planets:
        if row["star_id"] is None:
            assert (row["parent_type"], row["parent_id"]) == ("system", row["star_system_id"])
        else:
            assert (row["parent_type"], row["parent_id"]) == ("star", row["star_id"])
    for row in moons:
        assert (row["parent_type"], row["parent_id"]) == ("planet", row["planet_id"])
