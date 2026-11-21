# tests/test_db_old_schemas.py

"""
TEST.8, TEST.729: the Alembic baseline (v130) and the newest `NEWEST` old galaxy schemas (`fixtures/old_schemas/`) migrate to the current version and end up
exactly the shape `schema.sql` gives a new database: the same tables,
columns (type, nullability, default), indexes, foreign keys and CHECKs.
Starting from the real old schema, not a new database with columns
dropped, also covers drift between `schema.sql` and the revisions.
The run from the baseline passes through every revision, and the fixtures hold schema only (no rows), so a start point between adds nothing to it; the newest revisions
are the ones Boss's server has not run yet, so they are also started from the schema just before each. Any other fixture a test loads by number is kept for that test.
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

# A fresh database with no schema, which these tests build and migrate themselves (TEST.139 stage 3).
pytestmark = pytest.mark.fresh_database

NEWEST = 3
"""int: How many of the newest old schemas are started from besides the baseline (TEST.729)."""

_fresh_snapshots = {}


def checked_versions():
    """The old schemas the shape test starts from: the baseline and the `NEWEST` versions before the current one."""
    return [alembic_runner.BASELINE_VERSION, *range(max(alembic_runner.BASELINE_VERSION + 1, store.SCHEMA_VERSION - NEWEST), store.SCHEMA_VERSION)]


def fresh_snapshot(beside):
    """The shape of a new database on `beside`'s server (once per server)."""
    server = (beside.host, beside.port)
    if server not in _fresh_snapshots:
        with scratch_database(beside) as fresh:
            store.get_connection(fresh).close()
            _fresh_snapshots[server] = schema_snapshot(fresh)
    return _fresh_snapshots[server]


def test_the_baseline_and_the_newest_versions_before_the_current_one_have_a_fixture():
    versions = old_schema_versions()
    assert versions[0] == alembic_runner.BASELINE_VERSION
    assert set(checked_versions()) <= set(versions)
    assert len(checked_versions()) == NEWEST + 1


@pytest.mark.parametrize("version", checked_versions())
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
    assert len(reloaded.entries[0].star_system.primary.orbiters) == len(system.primary.orbiters)
