# tests/test_facility_velocity_migration.py

"""GEN.125: revision 0064 adds the velocity columns to `facilities` and works them out for the
stand-alone facilities already stored. Needs MySQL."""

import math

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.db


def test_migrating_from_v63_gives_stand_alone_facilities_the_rotation_curves_velocity(mysql_config):
    load_old_schema(mysql_config, 63)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("INSERT INTO sectors (name, edge_mpc, center_x_pc, center_y_pc, center_z_pc,"
                     " galactic_radius_pc, ring_index, layer_index, ring_slot_index)"
                     " VALUES ('Harbor', 4000, 0, 8000, 5, 8000, 100, 0, 1)")
        sector_id = conn.execute("SELECT MAX(id) AS id FROM sectors").fetchone()["id"]
        conn.execute("INSERT INTO facilities (name, kind, placement, host_type, sector_id, center_x_pc, center_y_pc,"
                     " center_z_pc, galactic_radius_pc) VALUES ('Waypoint', 'station', 'standalone', 'space', ?,"
                     " 0, 8000, 5, 8000)", (sector_id,))
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT velocity_x_kms, velocity_y_kms, velocity_z_kms FROM facilities").fetchone()
    finally:
        conn.close()
    # At (0, +8 kpc) the counterclockwise tangent points along -x, at the rotation speed there.
    assert row["velocity_y_kms"] == pytest.approx(0.0, abs=1e-9) and row["velocity_z_kms"] == 0.0
    assert row["velocity_x_kms"] < -100.0 and abs(row["velocity_x_kms"]) < 300.0
