# tests/test_query_db_list_order.py

"""
Coverage for the order `queryDb.list_sectors` (the browse page's sector
table, `GET /api/sectors`) and `queryDb.sector_detail` (the sector page's
systems table, `GET /api/sectors/<id>`) hand rows back in: sectors
nearest the galactic core first, with never-placed sectors after them by
name, and a sector's systems nearest its own center first.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable,
same as every other database-backed test in this suite.
"""

import math

import pytest

import queryDb
from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem


def _small_system_config():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = False
    return cfg


def _save_sector(mysql_config, name, center_pc=None, positions=((0.0, 0.0, 0.0),)):
    sector = SpaceSector(name, edge_ly=10.0)
    for position in positions:
        cfg = _small_system_config()
        sector.add_system(StarSystem(system_config=cfg), position=position, system_config=cfg)
    galaxy_position = None
    if center_pc is not None:
        galaxy_position = {
            "center_x_pc": center_pc[0], "center_y_pc": center_pc[1], "center_z_pc": center_pc[2],
            "galactic_radius_pc": math.dist(center_pc, (0.0, 0.0, 0.0)),
        }
    return _db.save_sector(sector, config=mysql_config, galaxy_position=galaxy_position)


def test_list_sectors_orders_by_distance_from_core_then_unplaced_by_name(mysql_config):
    # Names chosen so alphabetical order is the opposite of distance order.
    far_id = _save_sector(mysql_config, "Aardvark", center_pc=(900.0, 0.0, 0.0))
    near_id = _save_sector(mysql_config, "Zebra", center_pc=(0.0, 100.0, 0.0))
    middle_id = _save_sector(mysql_config, "Mongoose", center_pc=(0.0, 0.0, -400.0))
    unplaced_b_id = _save_sector(mysql_config, "Bravo Unplaced")
    unplaced_a_id = _save_sector(mysql_config, "Alpha Unplaced")

    conn = queryDb.open_readonly(mysql_config)
    try:
        everything = queryDb.list_sectors(conn)
        second_page = queryDb.list_sectors(conn, limit=2, offset=2)
    finally:
        conn.close()

    assert [row["id"] for row in everything] == [near_id, middle_id, far_id, unplaced_a_id, unplaced_b_id]
    assert [row["galactic_radius_ly"] is None for row in everything] == [False, False, False, True, True]
    assert everything[0]["galactic_radius_ly"] == pytest.approx(100.0 * 3.2616, rel=1e-3)
    # Paging walks the same order rather than re-sorting each page.
    assert [row["id"] for row in second_page] == [far_id, unplaced_a_id]


def test_sector_detail_orders_systems_by_distance_from_sector_center(mysql_config):
    positions = [(4.0, 0.0, 0.0), (0.0, -1.0, 0.0), (1.5, 1.5, 1.5), (0.0, 0.0, 3.0)]
    sector_id = _save_sector(mysql_config, "Ordered Sector", center_pc=(50.0, 0.0, 0.0), positions=positions)

    conn = queryDb.open_readonly(mysql_config)
    try:
        detail = queryDb.sector_detail(conn, sector_id)
    finally:
        conn.close()

    distances = [system["center_distance_ly"] for system in detail["systems"]]
    expected = sorted(math.dist(p, (0.0, 0.0, 0.0)) for p in positions)
    assert distances == pytest.approx(expected, rel=1e-6)
