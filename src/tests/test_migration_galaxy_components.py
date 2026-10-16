# tests/test_migration_galaxy_components.py

"""
GEN.159: migration 0097 adds the empty `galaxy_components` table; the store writes and reads its rows, the table's
CHECKs refuse a malformed row, and a component's `cell_key` is the one `galaxy/cellkey.py` computes.
"""

import dataclasses

import pytest

from planetgen.db import store
from planetgen.galaxy import cellkey, components
from planetgen.generation import clusters
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _tables(conn):
    return {row["TABLE_NAME"] for row in conn.execute(
        "SELECT TABLE_NAME FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()").fetchall()}


def test_a_v96_galaxy_gets_an_empty_components_table(mysql_config):
    load_old_schema(mysql_config, 96)
    assert store.schema_status(mysql_config)[0] == 96
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert "galaxy_components" in _tables(conn)
        assert conn.execute("SELECT COUNT(*) AS n FROM galaxy_components").fetchone()["n"] == 0
        assert store.get_galaxy_components(conn) == []
    finally:
        conn.close()


def test_a_migration_run_again_leaves_the_table_and_its_rows(mysql_config):
    load_old_schema(mysql_config, 96)
    store.migrate_database(mysql_config)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        store.replace_galaxy_components(conn, clusters.generate_clusters(5, count=3))
        conn.commit()
        conn.execute("DELETE FROM alembic_version")
        conn.execute("INSERT INTO alembic_version (version_num) VALUES ('0096')")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert len(store.get_galaxy_components(conn)) == 3
    finally:
        conn.close()


def test_clusters_round_trip_through_the_store(mysql_config):
    found = clusters.generate_clusters(b"round trip", count=20, epoch_unix=1.5e9)
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, found)
        conn.commit()
        back = store.get_galaxy_components(conn)
        assert [c.id for c in back] == sorted(c.id for c in back) and all(c.id for c in back)
        assert [_without_id(c) for c in back] == found
        assert store.get_galaxy_components(conn, kind="satellite") == []
        assert len(store.get_galaxy_components(conn, kind="cluster")) == 20
        store.replace_galaxy_components(conn, found[:5], kind="cluster")
        conn.commit()
        assert len(store.get_galaxy_components(conn)) == 5
    finally:
        conn.close()


def _without_id(component):
    return dataclasses.replace(component, id=None)


def test_replacing_a_kind_leaves_the_others(mysql_config):
    cluster = clusters.generate_clusters(1, count=2)
    satellite = components.Component(
        kind="satellite", name="Large Cloud", source="catalogue", catalogue_id="LMC", center_pc=(-1000.0, -40000.0, -27000.0),
        velocity_kms=(-57.0, -226.0, 221.0), epoch_unix=0.0, mass_solar=1.4e10, half_mass_radius_pc=1500.0,
        profile="exponential", age_gy=9.0, metallicity_fe_h=-0.5, seed=3)
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, cluster + [satellite])
        store.replace_galaxy_components(conn, clusters.generate_clusters(2, count=4), kind="cluster")
        conn.commit()
        assert len(store.get_galaxy_components(conn, kind="cluster")) == 4
        stored = store.get_galaxy_components(conn, kind="satellite")
        assert [_without_id(c) for c in stored] == [satellite]
        with pytest.raises(ValueError):
            store.replace_galaxy_components(conn, [satellite], kind="cluster")
    finally:
        conn.rollback()
        conn.close()


def test_a_cluster_inside_a_satellite_names_its_parent_and_loses_it_when_the_parent_goes(mysql_config):
    parent = components.Component(
        kind="satellite", name="Dwarf", source="generated", center_pc=(1.0, 2.0, 3.0), velocity_kms=(0.0, 0.0, 0.0),
        epoch_unix=0.0, mass_solar=1e8, half_mass_radius_pc=300.0, profile="plummer", age_gy=10.0, metallicity_fe_h=-1.8,
        seed=9)
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, [parent])
        parent_id = store.get_galaxy_components(conn)[0].id
        child = clusters.generate_clusters(4, count=1)[0]
        child = dataclasses.replace(child, parent_component_id=parent_id)
        conn.execute("DELETE FROM galaxy_components WHERE kind = 'cluster'")
        store.replace_galaxy_components(conn, [child], kind="cluster")
        assert store.get_galaxy_components(conn, kind="cluster")[0].parent_component_id == parent_id
        conn.execute("DELETE FROM galaxy_components WHERE kind = 'satellite'")
        assert store.get_galaxy_components(conn, kind="cluster")[0].parent_component_id is None
    finally:
        conn.rollback()
        conn.close()


def test_the_cell_key_is_the_one_the_cell_key_module_computes(mysql_config):
    found = clusters.generate_clusters(b"keys", count=10)
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, found)
        conn.commit()
        rows = conn.execute("SELECT name, cell_key FROM galaxy_components").fetchall()
        keys = {row["name"]: row["cell_key"] for row in rows}
        for c in found:
            assert keys[c.name] == cellkey.cell_key(*c.center_pc), c
    finally:
        conn.close()


@pytest.mark.parametrize("change", (
    "kind = 'galaxy'", "source = 'guess'", "profile = 'nfw'", "core_collapsed = 2", "king_w0 = NULL", "tidal_radius_pc = NULL"))
def test_the_table_refuses_a_malformed_row(mysql_config, change):
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, clusters.generate_clusters(1, count=1))
        conn.commit()
        with pytest.raises(Exception) as caught:
            conn.execute(f"UPDATE galaxy_components SET {change}")
        assert "chk_galaxy_components" in str(caught.value) or "constraint" in str(caught.value).lower()
    finally:
        conn.rollback()
        conn.close()


def test_a_kind_and_name_are_unique(mysql_config):
    one = clusters.generate_clusters(1, count=1)
    conn = store.get_connection(mysql_config)
    try:
        with pytest.raises(Exception) as caught:
            store.replace_galaxy_components(conn, one + one)
        assert "uq_galaxy_components_name" in str(caught.value) or "Duplicate" in str(caught.value)
    finally:
        conn.rollback()
        conn.close()
