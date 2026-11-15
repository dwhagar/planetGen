# tests/test_migration_sky_draw_parts.py

"""
OPS.253: migration 0260 gives `sky_tier_draws` the `part_index` and `part_count` columns and puts them in the unique key. A draw from before v260 keeps its row as a whole sphere (0 and 0).
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

INSERT = ("INSERT INTO sky_tier_draws (kind, band_low, band_high, center_x_pc, center_y_pc, center_z_pc, radius_pc, limit_mag, complete)"
          " VALUES ('luminosity', 100, 200, 1, 2, 3, 60, 7.5, 1)")


def test_a_v259_ledger_keeps_its_spheres_as_whole_ones_and_may_hold_parts_beside_them(mysql_config):
    load_old_schema(mysql_config, 258)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute(INSERT)
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT part_index, part_count FROM sky_tier_draws").fetchall() == [{"part_index": 0, "part_count": 0}]
        key = [row["COLUMN_NAME"] for row in conn.execute(
            "SELECT COLUMN_NAME FROM information_schema.STATISTICS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'sky_tier_draws'"
            " AND INDEX_NAME = 'uq_sky_tier_draws_cell' ORDER BY SEQ_IN_INDEX").fetchall()]
        assert key[-2:] == ["part_count", "part_index"]
        with pytest.raises(Exception):
            conn.execute(INSERT)                                        # the same whole sphere twice is still refused
        conn.rollback()
        conn.execute(INSERT.replace("complete)", "complete, part_index, part_count)").replace("7.5, 1)", "7.5, 0, 1, 4)"))
        conn.commit()
        assert conn.execute("SELECT COUNT(*) AS n FROM sky_tier_draws").fetchone()["n"] == 2
    finally:
        conn.close()
