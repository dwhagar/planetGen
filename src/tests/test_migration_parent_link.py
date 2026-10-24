# tests/test_migration_parent_link.py

"""
DB.73: migration 0140 gives stars, belts, comets, black holes, neutron stars, rogue planets and interstellar comets the
parent link planets and moons already had, backfilled from `parent_node_id` and `star_id`, widens `parent_type` and its
CHECK, lets a moon have no planet or system, and can be run again.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import insert_filled_row as _insert, load_old_schema, table_columns as _columns

pytestmark = pytest.mark.fresh_database

NEW_TABLES = ("stars", "asteroid_belts", "comets", "black_holes", "neutron_stars", "rogue_planets",
              "interstellar_comets")


def undo_parent_link(conn):
    """Takes a current database back to before v140 (the columns gone, the rows kept), for the tests that rebuild an older one."""
    conn.execute("ALTER TABLE moons DROP FOREIGN KEY fk_moons_rogue_planet, DROP KEY idx_moons_rogue_planet_id,"
                 " DROP COLUMN rogue_planet_id, ALGORITHM=COPY")
    for table in NEW_TABLES:
        conn.execute(f"ALTER TABLE {table} DROP CHECK chk_{table}_parent_type, DROP KEY idx_{table}_parent,"
                     " DROP COLUMN parent_type, DROP COLUMN parent_id, ALGORITHM=COPY")


def _old_database(mysql_config):
    load_old_schema(mysql_config, 139)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        _insert(conn, "star_systems", id=1, binary_configuration="wide", is_binary=1)
        _insert(conn, "star_systems", id=2)
        _insert(conn, "stars", id=1, star_system_id=1, role="primary", parent_node_id=0)
        _insert(conn, "stars", id=2, star_system_id=1, role="secondary", parent_node_id=0)
        _insert(conn, "stars", id=3, star_system_id=2, role="single")
        _insert(conn, "asteroid_belts", id=1, star_system_id=1, star_id=2, orbital_index=0, density="sparse")
        _insert(conn, "asteroid_belts", id=2, star_system_id=2, orbital_index=0, density="sparse")
        _insert(conn, "comets", id=1, star_system_id=1, star_id=1, orbit_type="parabolic")
        _insert(conn, "comets", id=2, star_system_id=2, orbit_type="parabolic")
        _insert(conn, "black_holes", id=1, star_id=3, has_accretion_disk=0)
        _insert(conn, "black_holes", id=2, has_accretion_disk=0)
        _insert(conn, "neutron_stars", id=1, star_id=1, pulsar_type="young")
        _insert(conn, "rogue_planets", id=1, planet_type="t", has_internal_heat=0, has_moons=0)
        _insert(conn, "interstellar_comets", id=1, is_active=1)
        _insert(conn, "planets", id=1, star_system_id=2, star_id=3, orbital_index=0, body_type="t", parent_type="star", parent_id=3)
        _insert(conn, "moons", id=1, planet_id=1, star_system_id=2, orbital_index=0, body_type="t", parent_type="planet", parent_id=1)
        conn.execute("SET FOREIGN_KEY_CHECKS = 1")
        conn.commit()
    finally:
        conn.close()


def _links(conn, table):
    return {row["id"]: (row["parent_type"], row["parent_id"]) for row in conn.execute(
        f"SELECT id, parent_type, parent_id FROM {table}").fetchall()}


def test_a_v139_database_gets_the_link_backfilled_from_what_said_it_before(mysql_config):
    _old_database(mysql_config)
    assert store.schema_status(mysql_config)[0] == 139
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
        for table in NEW_TABLES:
            assert {"parent_type", "parent_id"} <= _columns(conn, table), table
        assert _links(conn, "stars") == {1: ("node", 0), 2: ("node", 0), 3: ("system", 2)}
        assert _links(conn, "asteroid_belts") == {1: ("star", 2), 2: ("system", 2)}
        assert _links(conn, "comets") == {1: ("star", 1), 2: ("system", 2)}
        assert _links(conn, "black_holes") == {1: ("star", 3), 2: ("free", None)}
        assert _links(conn, "neutron_stars") == {1: ("star", 1)}
        assert _links(conn, "rogue_planets") == {1: ("free", None)}
        assert _links(conn, "interstellar_comets") == {1: ("free", None)}
        assert _links(conn, "planets") == {1: ("star", 3)} and _links(conn, "moons") == {1: ("planet", 1)}
    finally:
        conn.close()


def test_the_check_names_six_parents_and_a_moon_may_orbit_a_rogue_planet(mysql_config):
    _old_database(mysql_config)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("stars", "planets", "moons") + NEW_TABLES[1:]:
            with pytest.raises(Exception, match=f"chk_{table}_parent_type"):
                conn.execute(f"UPDATE {table} SET parent_type = 'galaxy' WHERE id = 1")
            conn.rollback()
        conn.execute("UPDATE planets SET parent_type = 'rogue_planet', parent_id = 1 WHERE id = 1")
        _insert(conn, "moons", planet_id=None, star_system_id=None, rogue_planet_id=1, orbital_index=0, body_type="t",
                parent_type="rogue_planet", parent_id=1)
        conn.commit()
        conn.execute("DELETE FROM rogue_planets WHERE id = 1")
        assert conn.execute("SELECT COUNT(*) AS n FROM moons WHERE rogue_planet_id IS NOT NULL").fetchone()["n"] == 0
        conn.rollback()
    finally:
        conn.close()


def test_a_migration_run_again_changes_nothing(mysql_config):
    _old_database(mysql_config)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        before = {table: _links(conn, table) for table in NEW_TABLES}
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0139')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert {table: _links(conn, table) for table in NEW_TABLES} == before
    finally:
        conn.close()
