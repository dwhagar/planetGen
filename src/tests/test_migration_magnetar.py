# tests/test_migration_magnetar.py

"""
GEN.153: migration 0096 lets `neutron_stars.pulsar_type` hold `magnetar`, replacing the CHECK that named three types.
Nothing else changes: the quadrant CHECK stays, and the rows already there are untouched.
"""

import importlib

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

migration = importlib.import_module("planetgen.db.migrations.versions.0096_magnetar_pulsar_type")


def _insert(conn, pulsar_type):
    conn.execute(
        "INSERT INTO neutron_stars (name, mass_solar, radius_km, spin_period_ms, magnetic_field_gauss, pulsar_type,"
        " surface_temperature_k, luminosity_w, age_gy) VALUES ('Core', 1.4, 12, 5000, 1e14, ?, 1e6, 1e25, 1e-6)",
        (pulsar_type,))


def _clauses(conn):
    return {row["n"]: row["c"] for row in conn.execute(
        "SELECT cc.constraint_name AS n, cc.check_clause AS c FROM information_schema.check_constraints cc"
        " JOIN information_schema.table_constraints tc ON tc.constraint_schema = cc.constraint_schema"
        " AND tc.constraint_name = cc.constraint_name"
        " WHERE cc.constraint_schema = DATABASE() AND tc.table_name = 'neutron_stars'").fetchall()}


def test_a_v95_galaxy_refuses_a_magnetar_until_it_is_migrated(mysql_config):
    load_old_schema(mysql_config, 95)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _insert(conn, "young")
        conn.commit()
        with pytest.raises(Exception, match="neutron_stars_chk_1"):
            _insert(conn, "magnetar")
    finally:
        conn.close()
    assert store.schema_status(mysql_config)[0] == 95
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for pulsar_type in migration.PULSAR_TYPES:
            _insert(conn, pulsar_type)
        with pytest.raises(Exception, match="neutron_stars_chk_1"):
            _insert(conn, "pulsar")
        conn.rollback()
        assert conn.execute("SELECT COUNT(*) AS n FROM neutron_stars WHERE pulsar_type = 'young'").fetchone()["n"] == 1
    finally:
        conn.close()


def test_the_check_names_exactly_the_four_types_and_the_quadrant_check_stays(mysql_config):
    load_old_schema(mysql_config, 95)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        clauses = _clauses(conn)
    finally:
        conn.close()
    assert all(name in clauses["neutron_stars_chk_1"] for name in migration.PULSAR_TYPES)
    assert clauses["neutron_stars_chk_1"].count("_utf8mb4") == len(migration.PULSAR_TYPES)
    assert "quadrant" in clauses["neutron_stars_chk_2"]


def test_a_stopped_migration_starts_again(mysql_config):
    """The old CHECK already dropped (the run stopped between the two halves of the statement is not possible, but a
    galaxy migrated by hand may lack it): the new one is added and a second run does nothing."""
    load_old_schema(mysql_config, 95)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("ALTER TABLE neutron_stars DROP CHECK neutron_stars_chk_1")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _insert(conn, "magnetar")
    finally:
        conn.close()
