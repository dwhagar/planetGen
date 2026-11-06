# tests/test_slow_read_paths.py

"""
The pages that timed out on a big galaxy (PERF.68, PERF.70, PERF.74):

- the Systems list sorted by Sector, Octant or Binary is read a group at a time (no sort of the table) and
  returns the same page the full sort would, for every offset, direction and filter;
- the Sectors list reads each sector's system count from `sector_system_counts`;
- the scattered-point query of a coarse node has an index of its own.
"""

import pytest

from planetgen.db import query as queryDb
from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem


def _sector(config, name, count, binary_every=2, star_types=("K1V",)):
    sector = SpaceSector(name, edge_ly=11.5)
    for index in range(count):
        cfg = SystemConfig()
        cfg.STAR_TYPE = star_types[index % len(star_types)]
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


@pytest.fixture
def mixed_galaxy(mysql_config):
    """`galaxy`, with G, K and M stars alternating, so a star type filter keeps some of each group (PERF.219)."""
    for name, count in (("Beta Sector", 6), ("Alpha Sector", 5), ("Empty Sector", 0), ("Gamma Sector", 7)):
        _sector(mysql_config, name, count, star_types=("G2V", "K1V", "M3V"))
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE star_systems SET sector_id = NULL, quadrant = NULL, position_x_mpc = NULL,"
                     " position_y_mpc = NULL, position_z_mpc = NULL, location = NULL ORDER BY id LIMIT 3")
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


@pytest.mark.parametrize("sort", ["sector", "octant"])
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
    assert sum(stored.values()) == 13 and len(stored) == 4       # every sector has a row (PERF.208), the empty one too


def test_the_scattered_point_query_has_an_index_of_its_own(galaxy):
    names = {row["Key_name"] for row in galaxy.execute("SHOW INDEX FROM phenomenon_scatter").fetchall()}
    assert "idx_phenomenon_scatter_class" in names


# --- PERF.76: a node piece that runs past its limit leaves the node incomplete ---------------------------

SLOW = "SELECT COUNT(*) AS n FROM (SELECT SLEEP(3)) slow"


def test_a_statement_limit_lowers_the_limit_for_the_block_and_restores_it(mysql_config):
    import pymysql

    conn = queryDb.open_readonly(mysql_config)  # no limit
    try:
        with conn.statement_limit(0.3):
            with pytest.raises(pymysql.err.OperationalError) as caught:
                conn.execute(SLOW).fetchone()
            assert caught.value.args[0] in store.STATEMENT_TIMEOUT_ERRORS
        assert conn.execute("SELECT SLEEP(0.6) AS s").fetchone()["s"] == 0
    finally:
        conn.close()
    limited = queryDb.open_readonly(mysql_config, statement_timeout_s=0.3)
    try:
        with limited.statement_limit(5):  # never raises a limit that is already lower
            with pytest.raises(pymysql.err.OperationalError):
                limited.execute(SLOW).fetchone()
    finally:
        limited.close()


def test_a_node_whose_slowest_piece_times_out_is_served_marked_incomplete(galaxy, monkeypatch):
    from planetgen.db import nodes
    monkeypatch.setattr(queryDb, "TILE_PIECE_SECONDS", 0.4)
    real = queryDb.galaxy_scattered_points_in_box

    def slow(conn, *args, **kwargs):
        conn.execute(SLOW).fetchone()
        return real(conn, *args, **kwargs)

    keys = ["243.1.1.1", "81.2.1.1"]
    monkeypatch.setattr(queryDb, "galaxy_scattered_points_in_box", slow)
    found = nodes.galaxy_nodes(galaxy, keys)["nodes"]
    assert all(node["incomplete"] == ["scattered"] for node in found.values())
    monkeypatch.setattr(queryDb, "galaxy_scattered_points_in_box", real)
    assert "incomplete" not in nodes.galaxy_nodes(galaxy, keys[:1])["nodes"][keys[0]]


# --- PERF.77: a count nobody has stored stops at the cap -------------------------------------------------

def test_a_count_with_no_stored_answer_stops_at_the_cap(galaxy, monkeypatch):
    monkeypatch.setattr(queryDb, "_counted", lambda conn, key, compute, fallback: fallback(conn))
    monkeypatch.setattr(queryDb, "COUNT_FALLBACK_CAP", 5)
    capped = queryDb.count_systems(galaxy, star_type_prefix="K")
    assert capped == 4 and queryDb.is_capped(capped)
    monkeypatch.setattr(queryDb, "COUNT_FALLBACK_CAP", 500)
    whole = queryDb.count_systems(galaxy, star_type_prefix="K")
    assert whole == 15 and not queryDb.is_capped(whole)
    assert queryDb.is_capped(queryDb.CappedCount(10000)) and not queryDb.is_capped(7)


# --- PERF.75: a page after the previous page's last system -----------------------------------------------

from tests.test_api import client  # noqa: E402,F401


@pytest.mark.parametrize("descending", [False, True])
def test_a_keyset_page_is_the_page_offset_paging_gives(galaxy, descending):
    everything = queryDb.list_systems(galaxy, sort="name", descending=descending)
    for index in range(0, len(everything), 4):
        after = everything[index - 1]["id"] if index else None
        page = queryDb.list_systems(galaxy, limit=4, offset=999, sort="name", descending=descending, after=after) \
            if after else queryDb.list_systems(galaxy, limit=4, sort="name", descending=descending)
        assert _ids(page) == _ids(everything[index:index + 4])
    # A system that is gone, or another sort, falls back to the offset.
    assert _ids(queryDb.list_systems(galaxy, limit=3, offset=2, sort="name", after=10 ** 9)) == \
        _ids(queryDb.list_systems(galaxy, limit=3, offset=2, sort="name"))
    assert _ids(queryDb.list_systems(galaxy, limit=3, offset=2, sort="octant", after=10 ** 9)) == \
        _ids(queryDb.list_systems(galaxy, limit=3, offset=2, sort="octant"))
    # a list filtered by sector reads by offset in every sort (it is not a group walk)
    sector = next(row["sector_id"] for row in everything if row["sector_id"] is not None)
    assert _ids(queryDb.list_systems(galaxy, limit=3, offset=2, sort="octant", sector_id=sector,
                                     after=everything[0]["id"])) == \
        _ids(queryDb.list_systems(galaxy, limit=3, offset=2, sort="octant", sector_id=sector))


@pytest.mark.parametrize("sort", ["sector", "octant"])
@pytest.mark.parametrize("descending", [False, True])
@pytest.mark.parametrize("filters", [{}, {"binary": True}, {"octants": ("I", "V")}])
def test_a_page_sorted_by_group_after_a_system_is_the_page_the_full_sort_gives(galaxy, sort, descending, filters):
    """PERF.206: the group walks (sector, octant) start after the system the previous page ended on, in any group
    and in the systems with no sector or octant at the end, and give the pages of one long read; `offset` is ignored."""
    everything = _ids(queryDb.list_systems(galaxy, sort=sort, descending=descending, **filters))
    assert len(everything) > 4
    size = 3
    for index in range(0, len(everything), size):
        page = _ids(queryDb.list_systems(galaxy, limit=size, offset=999 if index else 0, sort=sort,
                                         descending=descending, after=everything[index - 1] if index else None, **filters))
        assert page == everything[index:index + size], (sort, descending, filters, index)


STAR_TYPE_FILTERS = [{"star_type_prefix": "G"}, {"star_type_prefix": "K1", "binary": True},
                     {"star_type_prefix": "M", "octants": ("I", "V")}, {"star_type_prefix": "Z"}]


def test_a_star_type_filter_pages_by_key_and_by_offset_as_the_full_sort_does(mixed_galaxy):
    """PERF.219: a star type filter is one more condition inside each group of the walk; offset pages and keyset pages both give the
    pages of the unlimited (joined, sorted) list, a system counted once whichever of its stars matches. (One galaxy for all the
    combinations: building one costs more than reading it.)"""
    size = 3
    for sort in ("sector", "octant"):
        for descending in (False, True):
            for filters in STAR_TYPE_FILTERS:
                where = (sort, descending, filters)
                everything = _ids(queryDb.list_systems(mixed_galaxy, sort=sort, descending=descending, **filters))
                assert bool(everything) == (filters["star_type_prefix"] != "Z"), where
                assert len(set(everything)) == len(everything), where
                for offset in range(0, len(everything) + size, 2):
                    page = _ids(queryDb.list_systems(mixed_galaxy, limit=size, offset=offset, sort=sort, descending=descending, **filters))
                    assert page == everything[offset:offset + size], (*where, "offset", offset)
                for index in range(0, len(everything), size):
                    page = _ids(queryDb.list_systems(mixed_galaxy, limit=size, offset=999 if index else 0, sort=sort,
                                                     descending=descending, after=everything[index - 1] if index else None, **filters))
                    assert page == everything[index:index + size], (*where, "after", index)


def test_a_keyset_page_with_a_star_type_filter_may_start_after_a_system_the_filter_leaves_out(mixed_galaxy):
    """The key is a position in the sort, not a row of the filtered list: a link made before the filter changed still lands."""
    for sort in ("sector", "octant"):
        every = _ids(queryDb.list_systems(mixed_galaxy, sort=sort))
        kept = set(_ids(queryDb.list_systems(mixed_galaxy, sort=sort, star_type_prefix="G")))
        left_out = [number for number, system in enumerate(every) if system not in kept and number < len(every) - 6]
        assert left_out
        for number in left_out[::3]:
            expected = [system for system in every[number + 1:] if system in kept][:4]
            assert _ids(queryDb.list_systems(mixed_galaxy, limit=4, sort=sort, star_type_prefix="G", after=every[number])) == expected, \
                (sort, number)


def test_the_star_type_walk_never_joins_or_sorts_the_table(mixed_galaxy, monkeypatch):
    sql = []
    real = store.Connection._run
    monkeypatch.setattr(store.Connection, "_run", lambda conn, text, params: sql.append(text) or real(conn, text, params))
    queryDb.list_systems(mixed_galaxy, limit=3, sort="octant", star_type_prefix="G")
    walked = [text for text in sql if "FROM star_systems ss" in text]
    assert walked and not any("JOIN stars" in text or "DISTINCT" in text for text in walked)
    assert any("EXISTS (SELECT 1 FROM stars s" in text for text in walked)


def test_an_octant_offset_page_with_a_star_type_reads_as_before_and_a_sector_one_walks(mixed_galaxy, monkeypatch):
    """Skipping an octant by counting its matching systems costs a pass over it, so only a page by key walks the octant sort with a star
    type; the sector sort counts a sector at a time and walks either way (measured, PERF.219)."""
    sql = []
    real = store.Connection._run
    monkeypatch.setattr(store.Connection, "_run", lambda conn, text, params: sql.append(text) or real(conn, text, params))
    queryDb.list_systems(mixed_galaxy, limit=3, offset=4, sort="octant", star_type_prefix="G")
    assert any("JOIN stars" in text for text in sql)
    sql.clear()
    queryDb.list_systems(mixed_galaxy, limit=3, offset=4, sort="sector", star_type_prefix="G")
    assert not any("JOIN stars" in text for text in sql)


def test_the_systems_route_pages_by_key(client, galaxy):
    first = client.get("/api/systems?limit=4").get_json()["items"]
    second = client.get("/api/systems?limit=4&offset=4").get_json()["items"]
    keyed = client.get(f"/api/systems?limit=4&offset=999&after={first[-1]['id']}").get_json()["items"]
    assert [item["id"] for item in keyed] == [item["id"] for item in second]
    assert client.get("/api/systems?after=not-an-id").status_code == 400


# --- PERF.69: the Phenomena table's counts are stored like the others' --------------------------------------

def test_the_phenomena_counts_go_through_the_stored_count_cache(galaxy, monkeypatch):
    asked = []

    def counted(conn, key, compute, fallback, stamp=None):
        asked.append(key[0])
        return fallback(conn)

    monkeypatch.setattr(queryDb, "_counted", counted)
    assert queryDb.count_phenomena(galaxy) >= 0
    assert queryDb.phenomena_facets(galaxy) == {"type": [], "descriptor": []}
    assert asked == ["count_phenomena", "phenomena_facets"]
