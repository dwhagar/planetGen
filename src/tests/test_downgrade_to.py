# tests/test_downgrade_to.py

"""
TEST.328: `db_schema_support.downgrade_to` takes a database at the current schema to the shape of an older release, derived
from the checked-in schema of that release, so a migration test does not keep a list of what each newer migration added.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import (
    downgrade_to,
    load_old_schema,
    schema_differences,
    schema_snapshot,
    scratch_database,
)

pytestmark = pytest.mark.fresh_database


@pytest.mark.parametrize("version", [131, 140, 150, 165, 173])
def test_a_downgraded_database_has_the_shape_of_the_old_release(mysql_config, version):
    store.get_connection(mysql_config).close()          # the current schema
    downgrade_to(mysql_config, version)
    with scratch_database(mysql_config) as old:
        load_old_schema(old, version)
        differences = schema_differences(schema_snapshot(mysql_config), schema_snapshot(old))
    assert not differences, f"v{version}:\n" + "\n".join(differences)
    assert store.schema_status(mysql_config)[0] == version


def test_a_downgrade_keeps_the_rows_and_the_migrations_bring_it_back(mysql_config):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.util import draw

    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = True
    cfg.BINARY_SYSTEM = False
    with draw.bound(5):
        system = StarSystem(system_config=cfg)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            store.insert_star_system(conn, system, cfg)
        planets = conn.execute("SELECT COUNT(*) AS n FROM planets").fetchone()["n"]
    finally:
        conn.close()
    assert planets > 0
    downgrade_to(mysql_config, 150)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM planets").fetchone()["n"] == planets
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    with scratch_database(mysql_config) as fresh:
        store.get_connection(fresh).close()
        differences = schema_differences(schema_snapshot(mysql_config), schema_snapshot(fresh))
    assert not differences, "\n".join(differences)
