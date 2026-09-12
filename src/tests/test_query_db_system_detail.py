# tests/test_query_db_system_detail.py

"""
Coverage for `queryDb.system_detail`'s `binary_mutual_position_x/y/z_km`
keys -- added so `html/lib/systemmap.py`'s System Map could place a
binary pair's two stars at their real mass-weighted offsets from the
system's own barycenter, instead of the old map's fixed schematic offset.
Before this, `system_detail` queried the underlying `star_systems` row
(which already carried these columns, reused unchanged from a 'close'
pair's own mutual-orbit tracking -- see `schema.sql`'s "v14"/"v15" notes)
but never copied them into its own returned dict.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable,
same as every other database-backed test in this suite.
"""

import pytest

import queryDb
from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.systemData import StarSystem


def _insert(mysql_config, system, cfg):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
        return system_id
    finally:
        conn.close()


def test_system_detail_exposes_binary_mutual_position_for_a_close_binary(mysql_config):
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = True
    cfg.WIDE_BINARY = False
    system = StarSystem(system_config=cfg)
    assert system.binary_type == "close"

    system_id = _insert(mysql_config, system, cfg)

    conn = _db.get_connection(mysql_config)
    try:
        raw_row = conn.execute("SELECT * FROM star_systems WHERE id = ?", (system_id,)).fetchone()
        detail = queryDb.system_detail(conn, system_id)
    finally:
        conn.close()

    assert detail["binary_configuration"] == "close"
    for axis in ("x", "y", "z"):
        key = f"binary_mutual_position_{axis}_km"
        assert detail[key] is not None
        assert detail[key] == pytest.approx(raw_row[key])


def test_system_detail_binary_mutual_position_is_none_for_a_single_star(mysql_config):
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)

    system_id = _insert(mysql_config, system, cfg)

    conn = _db.get_connection(mysql_config)
    try:
        detail = queryDb.system_detail(conn, system_id)
    finally:
        conn.close()

    assert detail["binary_configuration"] is None
    assert detail["binary_mutual_position_x_km"] is None
    assert detail["binary_mutual_position_y_km"] is None
    assert detail["binary_mutual_position_z_km"] is None
