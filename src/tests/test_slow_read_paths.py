# tests/test_slow_read_paths.py

"""
The pages that timed out on a big galaxy (PERF.68, PERF.70, PERF.74):

- the Systems list sorted by Sector, Octant or Binary is read a group at a time (no sort of the table) and
  returns the same page the full sort would, for every offset, direction and filter;
- the Sectors list reads each sector's system count from `sector_system_counts`;
- the coarse tile's scattered-point query has an index of its own.
"""

import pytest

from planetgen.db import query as queryDb
from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem


def _sector(config, name, count, binary_every=2):
    sector = SpaceSector(name, edge_ly=11.5)
    for index in range(count):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "K1V"
        cfg.BINARY_SYSTEM = index % binary_every == 0
        cfg.MOONS = False
        system = StarSystem(system_config=cfg)
        sector.add_system(system, position=(float(index % 3) - 1, float(index % 5) - 2, float(index % 2) - 1),
                          system_config=cfg)
    return store.save_sector(sector, config=config)


@pytest.fixture
def galaxy(mysql_config):
    for name, count in (("Beta Sector", 5), ("Alpha Sector", 4), ("Empty Sector", 0), ("Gamma Sector", 6)):
        _sector(mysql_config, name, count)
    conn = store.get_connection(mysql_config)
    try:
        # Two standalone systems: no sector, no octant.
        conn.execute("UPDATE star_systems SET sector_id = NULL, quadrant = NULL, position_x_mpc = NULL,"
                     " position_y_mpc = NULL, position_z_mpc = NULL, location = NULL ORDER BY id LIMIT 2")
        conn.commit()
    finally:
        conn.close()
    conn = queryDb.open_readonly(mysql_config)
    try:
        yield conn
    finally:
        conn.close()


def _ids(rows):
    return [row["id"] for row in rows]


@pytest.mark.parametrize("sort", ["sector", "octant", "binary"])
@pytest.mark.parametrize("descending", [False, True])
@pytest.mark.parametrize("filters", [{}, {"binary": True}, {"octants": ("I", "V")}])
def test_a_page_sorted_by_group_is_the_page_the_full_sort_gives(galaxy, sort, descending, filters):
    everything = _ids(queryDb.list_systems(galaxy, sort=sort, descending=descending, **filters))
    assert everything
    size = 4
    for offset in range(0, len(everything) + size, 3):
        page = _ids(queryDb.list_systems(galaxy, limit=size, offset=offset, sort=sort, descending=descending,
                                         **filters))
        assert page == everything[offset:offset + size], (sort, descending, filters, offset)


def test_the_group_walk_reads_no_more_than_a_page_of_rows(galaxy, monkeypatch):
    sql = []
    real = store.Connection._run
    monkeypatch.setattr(store.Connection, "_run", lambda conn, text, params: sql.append(text) or real(conn, text, params))
    queryDb.list_systems(galaxy, limit=3, sort="octant")
    assert not any("ORDER BY ss.quadrant" in text or "ss.quadrant IS NULL," in text for text in sql)


def test_the_sectors_list_gives_each_sectors_system_count(galaxy):
    counts = {row["name"]: row["system_count"] for row in queryDb.list_sectors(galaxy)}
    assert counts["Alpha Sector"] + counts["Beta Sector"] + counts["Gamma Sector"] == 13
    assert counts["Empty Sector"] == 0
    by_systems = [row["system_count"] for row in queryDb.list_sectors(galaxy, sort="systems", descending=True)]
    assert by_systems == sorted(by_systems, reverse=True)
    stored = {row["sector_id"]: row["system_count"] for row in galaxy.execute(
        "SELECT sector_id, system_count FROM sector_system_counts").fetchall()}
    assert sum(stored.values()) == 13 and len(stored) == 3


def test_the_scattered_point_query_has_an_index_of_its_own(galaxy):
    names = {row["Key_name"] for row in galaxy.execute("SHOW INDEX FROM phenomenon_scatter").fetchall()}
    assert "idx_phenomenon_scatter_class" in names
