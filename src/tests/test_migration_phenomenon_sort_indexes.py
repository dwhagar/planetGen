# tests/test_migration_phenomenon_sort_indexes.py

"""
PERF.211: migration 0186 gives each of the eight phenomenon tables the indexes the Phenomena table's placed, descriptor and radius sorts read in order
(`idx_<table>_placed_name`, `idx_<table>_descriptor_name`, and `idx_<table>_radius_name` for the three tables with a radius). Some key parts are
expressions. Nothing else changes (`test_db_old_schemas` compares the whole shape with a new database's, expressions included).
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema, migration_number

pytestmark = pytest.mark.fresh_database

TABLES = ("nebulae", "asteroid_fields", "black_holes", "neutron_stars", "supernova_remnants", "rogue_planets", "interstellar_comets", "quasars")
RADIUS = ("nebulae", "asteroid_fields", "supernova_remnants")


def _indexes(conn):
    found = {}
    for row in conn.execute("SELECT TABLE_NAME AS t, INDEX_NAME AS i, SEQ_IN_INDEX AS s, COLUMN_NAME AS c, EXPRESSION AS e "
                            "FROM information_schema.STATISTICS WHERE TABLE_SCHEMA = DATABASE() ORDER BY TABLE_NAME, INDEX_NAME, SEQ_IN_INDEX").fetchall():
        if row["i"].endswith("_name"):
            found.setdefault((row["t"], row["i"]), []).append(row["c"] if row["c"] else f"({row['e']})")
    return found


def test_the_previous_schemas_tables_get_the_sort_indexes(mysql_config):
    number = migration_number("phenomenon_sort_indexes")
    load_old_schema(mysql_config, number - 1)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        before = _indexes(conn)
        assert ("nebulae", "idx_nebulae_placed_name") not in before and ("nebulae", "idx_nebulae_name") in before
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        after = _indexes(conn)
        for table in TABLES:
            for sort in ("placed", "descriptor") + (("radius",) if table in RADIUS else ()):
                assert (table, f"idx_{table}_{sort}_name") in after, (table, sort)
            assert after[(table, f"idx_{table}_placed_name")][-1] == "name"
            placed = after[(table, f"idx_{table}_placed_name")]
            assert (placed[0] == "star_id") == (table in ("black_holes", "neutron_stars")), table
            assert "is not null" in " ".join(placed).lower(), placed
        assert not any(name.startswith("idx_quasars_radius") or name.startswith("idx_rogue_planets_radius") for _, name in after)
        assert "case" in " ".join(after[("rogue_planets", "idx_rogue_planets_descriptor_name")]).lower()
        assert after[("nebulae", "idx_nebulae_descriptor_name")] == ["nebula_type", "name"]
    finally:
        conn.close()
