# tests/test_migration_orbit_nodes.py

"""
DB.53: migration 0104 adds `orbit_nodes` and the star and planet columns of a built system's orbit tree, converts a
stored binary into one node, and can be run again. The new tables and columns keep the names of the scatter's.
"""

import importlib
import math

import pytest

from planetgen.db import models, store
from planetgen.physics import constants
from tests.db_schema_support import load_old_schema
from tests.test_migration_point_columns import _columns, _insert

pytestmark = pytest.mark.fresh_database

MIGRATION = importlib.import_module("planetgen.db.migrations.versions.0104_orbit_nodes")
TABLES = models.metadata.tables
AU_KM = constants.AU_M / 1000.0


def _names(table):
    return {column.name for column in TABLES[table].columns}


def test_a_node_holds_the_columns_of_a_scatter_node_and_a_star_those_of_a_scatter_member():
    node, scatter = _names("orbit_nodes"), _names("scatter_nodes")
    assert {"star_system_id", "node_no", "parent_node_no", "disk_kind", "disk_mass_solar", "disk_inner_au",
            "disk_outer_au"} <= node
    assert scatter - {"source", "source_id", "id"} <= node, "a scatter node's column is missing from orbit_nodes"
    members = _names("scatter_members")
    for name in ("state", "flags", "offset_x_au", "offset_y_au", "offset_z_au", "velocity_x_kms", "velocity_y_kms",
                 "velocity_z_kms", "decay_age_yr", "epoch_unix"):
        assert name in members and name in _names("stars"), name
    assert {"parent_node_id", "disk_kind", "disk_mass_solar", "disk_inner_au", "disk_outer_au"} <= _names("stars")
    assert "host_node_id" in _names("planets")


def test_the_migration_lists_what_the_models_hold():
    stars = {name for name, _ in MIGRATION.STAR_COLUMNS}
    assert stars <= _names("stars")
    assert {name for name, _ in MIGRATION.PLANET_COLUMNS} <= _names("planets")


def test_a_v103_database_gets_the_table_and_the_columns(mysql_config):
    load_old_schema(mysql_config, 103)
    assert store.schema_status(mysql_config)[0] == 103
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert _names("orbit_nodes") == _columns(conn, "orbit_nodes")
        assert {name for name, _ in MIGRATION.STAR_COLUMNS} <= _columns(conn, "stars")
        assert "host_node_id" in _columns(conn, "planets")
        assert store.detect_schema_version(conn) == 104
    finally:
        conn.close()


def test_a_migration_run_again_changes_nothing(mysql_config):
    load_old_schema(mysql_config, 103)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0103')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION


def _system(conn, system_id, configuration, **given):
    conn.execute("SET FOREIGN_KEY_CHECKS = 0")
    _insert(conn, "star_systems", id=system_id, binary_configuration=configuration, is_binary=int(bool(configuration)),
            **given)
    conn.execute("SET FOREIGN_KEY_CHECKS = 1")


def test_a_wide_pair_becomes_one_node_and_every_star_gets_its_offset(mysql_config):
    load_old_schema(mysql_config, 103)
    solar = constants.SOLAR_MASS_TO_KG
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _system(conn, 1, "wide", binary_separation_km=2.0 * AU_KM, binary_eccentricity=0.3,
                binary_mutual_orbital_inclination_deg=12.0, binary_mutual_orbital_ascending_node_deg=40.0,
                binary_mutual_orbital_phase_deg=77.0, binary_epoch_unix=5.0, binary_mutual_orbital_period_years=3.5,
                binary_primary_position_x_km=0.5 * AU_KM, binary_primary_position_y_km=0.0,
                binary_primary_position_z_km=0.0, binary_secondary_position_x_km=-1.5 * AU_KM,
                binary_secondary_position_y_km=0.0, binary_secondary_position_z_km=0.0)
        _system(conn, 2, None)
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        _insert(conn, "stars", star_system_id=1, role="primary", mass_kg=1.0 * solar)
        _insert(conn, "stars", star_system_id=1, role="secondary", mass_kg=0.5 * solar)
        _insert(conn, "stars", star_system_id=2, role="single", mass_kg=0.8 * solar)
        conn.execute("SET FOREIGN_KEY_CHECKS = 1")
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        nodes = conn.execute("SELECT * FROM orbit_nodes").fetchall()
        assert len(nodes) == 1
        node = nodes[0]
        assert node["star_system_id"] == 1 and node["node_no"] == 0 and node["parent_node_no"] is None
        assert node["kind"] == "wide" and node["eccentricity"] == 0.3 and node["inclination_deg"] == 12.0
        assert node["mean_anomaly_deg"] == 77.0 and node["epoch_unix"] == 5.0 and node["period_years"] == 3.5
        assert math.isclose(node["mass_solar"], 1.5, rel_tol=1e-9)
        assert math.isclose(node["separation_km"], 2.0 * AU_KM)
        stars = {row["role"]: row for row in conn.execute("SELECT * FROM stars").fetchall()}
        assert stars["primary"]["parent_node_id"] == 0 and stars["secondary"]["parent_node_id"] == 0
        assert stars["single"]["parent_node_id"] is None
        assert math.isclose(stars["primary"]["offset_x_au"], 0.5) and math.isclose(stars["secondary"]["offset_x_au"], -1.5)
        assert stars["single"]["offset_x_au"] == 0 and stars["single"]["offset_y_au"] == 0
        assert stars["single"]["offset_z_au"] == 0
        assert stars["primary"]["velocity_x_kms"] is None and stars["primary"]["epoch_unix"] is None
    finally:
        conn.close()


def test_deleting_a_system_deletes_its_nodes(mysql_config):
    load_old_schema(mysql_config, 103)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        _system(conn, 1, "close", binary_separation_km=0.0)
        conn.commit()
    finally:
        conn.close()
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM orbit_nodes").fetchone()["n"] == 1
        conn.execute("DELETE FROM star_systems WHERE id = 1")
        conn.commit()
        assert conn.execute("SELECT COUNT(*) AS n FROM orbit_nodes").fetchone()["n"] == 0
    finally:
        conn.close()
