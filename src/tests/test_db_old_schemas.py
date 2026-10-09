# tests/test_db_old_schemas.py

"""
TEST.8: every checked-in old galaxy schema (`fixtures/old_schemas/`, the
Alembic baseline v61 and later) migrates to the current version and ends up
exactly the shape `schema.sql` gives a new database: the same tables,
columns (type, nullability, default), indexes, foreign keys and CHECKs.
Starting from the real old schema, not a new database with columns
dropped, also covers drift between `schema.sql` and the revisions.
"""

import pytest

from planetgen.db import alembic_runner, store
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem
from tests.db_schema_support import (
    load_old_schema,
    old_schema_versions,
    schema_differences,
    schema_snapshot,
    scratch_database,
)

_fresh_snapshots = {}


def fresh_snapshot(beside):
    """The shape of a new database on `beside`'s server (once per server)."""
    server = (beside.host, beside.port)
    if server not in _fresh_snapshots:
        with scratch_database(beside) as fresh:
            store.get_connection(fresh).close()
            _fresh_snapshots[server] = schema_snapshot(fresh)
    return _fresh_snapshots[server]


def test_every_version_before_the_current_one_has_a_fixture():
    versions = old_schema_versions()
    assert versions[0] == alembic_runner.BASELINE_VERSION
    assert set(range(alembic_runner.BASELINE_VERSION, store.SCHEMA_VERSION)) <= set(versions)


@pytest.mark.parametrize("version", old_schema_versions())
def test_old_schema_migrates_to_the_fresh_shape(mysql_config, version):
    load_old_schema(mysql_config, version)
    assert store.schema_status(mysql_config)[0] == version

    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION

    differences = schema_differences(schema_snapshot(mysql_config), fresh_snapshot(mysql_config))
    assert not differences, f"v{version} migrated differs from a new database:\n" + "\n".join(differences)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        versions = [row["version"] for row in conn.execute("SELECT version FROM schema_migrations").fetchall()]
    finally:
        conn.close()
    assert sorted(versions) == [version, *range(version + 1, store.SCHEMA_VERSION + 1)]


def test_a_baseline_database_saves_and_loads_a_sector(mysql_config):
    load_old_schema(mysql_config, alembic_runner.BASELINE_VERSION)
    store.migrate_database(mysql_config)

    cfg = SystemConfig()
    cfg.STAR_TYPE = "K2V"
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    sector = SpaceSector("Old Schema Sector", edge_ly=11.5)
    sector.add_system(system, position=(1.0, 2.0, 3.0), system_config=cfg)
    sector_id = store.save_sector(sector, config=mysql_config)

    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        reloaded = store.load_sector(conn, sector_id)
    finally:
        conn.close()
    assert reloaded.name == sector.name
    assert [entry.star_system.name for entry in reloaded.entries] == [system.name]
    assert len(reloaded.entries[0].star_system.planets) == len(system.planets)
