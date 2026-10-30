# tests/test_migration_sector_sort_keys.py

"""
PERF.208: `sector_system_counts` holds a count, a density and a Quadrant for every sector. Migration 0181 gives the table of the schema before it the
two columns and their indexes and fills it again, one row for each sector, the empty ones and the ones with no galaxy position included.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import insert_filled_row, load_old_schema, migration_number

pytestmark = pytest.mark.fresh_database

SECTORS = (       # name, centre (x, y) or None, systems
    ("North east", (10.0, 20.0), 2), ("North west", (-10.0, 20.0), 1), ("South west", (-10.0, -20.0), 0),
    ("South east", (10.0, -20.0), 3), ("Nowhere", None, 1),
)


def test_the_previous_schemas_table_is_filled_again_with_a_row_for_every_sector(mysql_config):
    number = migration_number("sector_sort_keys")
    load_old_schema(mysql_config, number - 1)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not conn.execute("SHOW COLUMNS FROM sector_system_counts LIKE 'quadrant'").fetchall()
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        for name, centre, count in SECTORS:
            if centre is None:
                sector_id = conn.execute("INSERT INTO sectors (name, edge_mpc) VALUES (?, 4000)", (name,)).lastrowid
            else:
                sector_id = conn.execute(
                    "INSERT INTO sectors (name, edge_mpc, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc) VALUES (?, 4000, ?, ?, 0, 22)",
                    (name, centre[0], centre[1])).lastrowid
            for index in range(count):
                insert_filled_row(conn, "star_systems", name=f"{name} {index}", sector_id=sector_id)
            if count:       # the table of the schema before has a row only for a sector with systems
                conn.execute("INSERT INTO sector_system_counts (sector_id, system_count) VALUES (?, ?)", (sector_id, count))
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        rows = {row["name"]: row for row in conn.execute(
            "SELECT sec.name, c.system_count, c.system_density, c.quadrant FROM sectors sec "
            "JOIN sector_system_counts c ON c.sector_id = sec.id").fetchall()}
        assert set(rows) == {name for name, _, _ in SECTORS}
        for name, centre, count in SECTORS:
            assert rows[name]["system_count"] == count, name
            assert (rows[name]["system_density"] is None) == (centre is None), name
        assert rows["North east"]["system_density"] == pytest.approx(2 / 4000 ** 3)
        assert {name: rows[name]["quadrant"] for name in rows} == {
            "North east": "I", "North west": "II", "South west": "III", "South east": "IV", "Nowhere": None}
        names = {row["Key_name"] for row in conn.execute("SHOW INDEX FROM sector_system_counts").fetchall()}
        assert {"idx_sector_system_counts_density", "idx_sector_system_counts_quadrant"} <= names
    finally:
        conn.close()
