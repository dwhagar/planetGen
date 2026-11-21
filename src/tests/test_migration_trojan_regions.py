# tests/test_migration_trojan_regions.py

"""
GEN.1150: migration 0271 adds the five columns of a giant's trojan region to `asteroid_belts` (`lagrange_point`, `libration_amplitude_deg`, `libration_phase_deg`, `mean_eccentricity`,
`mean_inclination_deg`). The rows made before keep their figures and get NULL, and the migration can be run over a table that already has a column (a stopped run is started again).
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

COLUMNS = ("lagrange_point", "libration_amplitude_deg", "libration_phase_deg", "mean_eccentricity", "mean_inclination_deg")


def _has_column(conn, name):
    return bool(conn.execute(f"SHOW COLUMNS FROM asteroid_belts LIKE '{name}'").fetchall())


def test_the_columns_are_added_empty_and_the_rows_keep_their_figures(mysql_config):
    load_old_schema(mysql_config, 270)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert not any(_has_column(conn, name) for name in COLUMNS)
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        conn.execute("INSERT INTO asteroid_belts (id, star_system_id, orbital_index, distance_km, lower_limit_km, upper_limit_km, density, composition_summary, parent_type, parent_id) "
                     "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (1, 1, 1, 1.0, 0.5, 1.5, "typical", "rock", "star", 1))
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert all(_has_column(conn, name) for name in COLUMNS)
        row = conn.execute(f"SELECT distance_km, {', '.join(COLUMNS)} FROM asteroid_belts WHERE id = 1").fetchone()
        assert row["distance_km"] == 1.0 and all(row[name] is None for name in COLUMNS)
    finally:
        conn.close()
