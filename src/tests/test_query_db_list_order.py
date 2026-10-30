# tests/test_query_db_list_order.py

"""
Coverage for the order `query.list_sectors` (the browse page's sector
table, `GET /api/sectors`) and `query.sector_detail` (the sector page's
systems table, `GET /api/sectors/<id>`) hand rows back in: sectors
nearest the galactic core first, with never-placed sectors after them by
name, and a sector's systems nearest its own center first.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable,
same as every other database-backed test in this suite.
"""

import math

import pytest

from planetgen.db import query
from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem


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
    return store.save_sector(sector, config=mysql_config, galaxy_position=galaxy_position)


def test_list_sectors_orders_by_distance_from_core_then_unplaced_by_name(mysql_config):
    # Names chosen so alphabetical order is the opposite of distance order.
    far_id = _save_sector(mysql_config, "Aardvark", center_pc=(900.0, 0.0, 0.0))
    near_id = _save_sector(mysql_config, "Zebra", center_pc=(0.0, 100.0, 0.0))
    middle_id = _save_sector(mysql_config, "Mongoose", center_pc=(0.0, 0.0, -400.0))
    unplaced_b_id = _save_sector(mysql_config, "Bravo Unplaced")
    unplaced_a_id = _save_sector(mysql_config, "Alpha Unplaced")

    conn = query.open_readonly(mysql_config)
    try:
        everything = query.list_sectors(conn)
        second_page = query.list_sectors(conn, limit=2, offset=2)
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

    conn = query.open_readonly(mysql_config)
    try:
        detail = query.sector_detail(conn, sector_id)
    finally:
        conn.close()

    distances = [system["center_distance_ly"] for system in detail["systems"]]
    expected = sorted(math.dist(p, (0.0, 0.0, 0.0)) for p in positions)
    assert distances == pytest.approx(expected, rel=1e-6)


# --- paging by key (PERF.137) ---------------------------------------------------------

def _galaxy_of_sectors(mysql_config):
    """Seven sectors: two at the same distance from the core, three more placed, two unplaced; the ids are in
    the order they are saved, which is not the order of any sort."""
    spots = {"Quail": (100.0, 0.0, 0.0), "Lynx": (-100.0, 0.0, 0.0), "Wolf": (0.0, -300.0, 0.0),
             "Ibis": (0.0, 50.0, 0.0), "Bison": (0.0, 0.0, 700.0)}
    made = {name: _save_sector(mysql_config, name, center_pc=spot) for name, spot in spots.items()}
    made["Yak Unplaced"] = _save_sector(mysql_config, "Yak Unplaced")
    made["Ant Unplaced"] = _save_sector(mysql_config, "Ant Unplaced")
    return made


def _ids_by_key(conn, **kwargs):
    """Every sector id from pages of two, each started `after` the one before."""
    seen, after = [], None
    while True:
        page = query.list_sectors(conn, limit=2, after=after, **kwargs)
        if not page:
            return seen
        seen += [row["id"] for row in page]
        after = page[-1]["id"]
        assert len(seen) < 50


@pytest.mark.parametrize("sort, descending", [(None, False), ("name", False), ("name", True),
                                              ("distance", False), ("distance", True)])
@pytest.mark.parametrize("quadrants", [(), ("I", "II"), ("III", "unplaced"), ("unplaced",), ("IV",)])
def test_a_page_after_a_sector_continues_the_whole_list_across_placed_and_unplaced(
        mysql_config, sort, descending, quadrants):
    _galaxy_of_sectors(mysql_config)
    conn = query.open_readonly(mysql_config)
    try:
        whole = [row["id"] for row in query.list_sectors(conn, sort=sort, descending=descending, quadrants=quadrants)]
        by_key = _ids_by_key(conn, sort=sort, descending=descending, quadrants=quadrants)
        by_offset = [row["id"] for offset in range(0, len(whole) + 2, 2)
                     for row in query.list_sectors(conn, limit=2, offset=offset, sort=sort, descending=descending,
                                                   quadrants=quadrants)]
    finally:
        conn.close()
    assert len(set(whole)) == len(whole)
    assert by_key == whole == by_offset


def test_the_distance_sort_and_the_name_sort_put_sectors_in_their_order(mysql_config):
    made = _galaxy_of_sectors(mysql_config)
    conn = query.open_readonly(mysql_config)
    try:
        def names(**kwargs):
            return [row["name"] for row in query.list_sectors(conn, **kwargs)]

        nearest = names()
        assert nearest == ["Ibis", "Quail", "Lynx", "Wolf", "Bison", "Ant Unplaced", "Yak Unplaced"]
        assert names(sort="distance") == nearest
        assert names(sort="distance", descending=True) == ["Bison", "Wolf", "Lynx", "Quail", "Ibis",
                                                           "Ant Unplaced", "Yak Unplaced"]
        assert names(sort="name") == sorted(made)
        assert names(sort="name", descending=True) == sorted(made, reverse=True)
        # An offset past the placed sectors lands among the unplaced ones, and past the end there is nothing.
        assert names(limit=3, offset=4) == ["Bison", "Ant Unplaced", "Yak Unplaced"]
        assert names(limit=3, offset=6) == ["Yak Unplaced"]
        assert names(limit=3, offset=7) == []
    finally:
        conn.close()



# --- PERF.207: planets and moons page by key -----------------------------------------------------------------------

@pytest.mark.parametrize("chunk", [1, 2, 1000])
def test_planets_and_moons_page_by_key_in_the_order_one_long_read_gives(mysql_config, monkeypatch, chunk):
    """`after` (the row id the previous page ended on) starts the next page in the same order, for planets and for
    moons, with and without a filter; an id that names nothing starts at the top and `offset` is not used with it.
    Unfiltered pages (PERF.215) read the bodies of `chunk` systems at a time, so a page spans chunks when it is small."""
    from tests.test_api import _save_wide_binary_with_moons

    monkeypatch.setattr(query, "BODY_CHUNK_SYSTEMS", chunk)

    for _ in range(3):
        _save_wide_binary_with_moons(mysql_config)
    conn = query.open_readonly(mysql_config)
    try:
        for read, kind in ((query.list_planets, "planets"), (query.list_moons, "moons")):
            for filters in ({}, {"min_radius_km": 1.0}):
                everything = [row["id"] for row in read(conn, **filters)]
                assert len(everything) >= 5, (kind, "the seeded galaxy has several")
                paged, after = [], None
                for _ in range(len(everything) + 1):
                    page = read(conn, limit=2, offset=7 if after else 0, after=after, **filters)
                    if not page:
                        break
                    paged += [row["id"] for row in page]
                    after = page[-1]["id"]
                assert paged == everything, (kind, filters)
            assert [row["id"] for row in read(conn, limit=3, after=10 ** 9)] == everything[:3]
    finally:
        conn.close()
