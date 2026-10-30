# tests/test_sector_sort_keys.py

"""
PERF.208: the Sectors list's sorts by systems, density and position page by key. `sector_system_counts` holds a row for
every sector with its system count, density and Quadrant, and each sort walks an index of it, the sector id breaking
ties in the sort's direction. A page by `after` is the page `offset` gives, and both are the order worked out here from
the rows. Sectors with no galaxy position come last in the density and position sorts; a sector made since the last
refresh follows the others, by id.
"""

import pytest

from planetgen.db import query, store

from tests.test_query_db_list_order import _save_sector

SECTORS = (       # name, where it is in the galaxy (pc) or None, how many systems it holds
    ("A", (500.0, 400.0, 0.0), 3), ("B", (-500.0, 300.0, 0.0), 1), ("C", (-300.0, -200.0, 0.0), 3),
    ("D", (200.0, -700.0, 0.0), 2), ("E", (800.0, 100.0, 0.0), 0), ("F", None, 2), ("G", None, 3),
    ("H", (100.0, 50.0, 0.0), 1), ("I", (-100.0, 900.0, 0.0), 3),
)
QUADRANT_OF = {"A": "I", "B": "II", "C": "III", "D": "IV", "E": "I", "H": "I", "I": "II"}


@pytest.fixture
def sectors(mysql_config, monkeypatch):
    """The sectors above saved, two of them twice the edge (so the density sort differs from the systems sort), the
    counts refreshed; the background refresh is switched off so a test controls when the rows are made."""
    monkeypatch.setattr(query, "_counted", lambda *args, **kwargs: None)
    for name, center, count in SECTORS:
        _save_sector(mysql_config, name, center_pc=center, positions=[(float(i), 0.0, 0.0) for i in range(count)])
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE sectors SET edge_mpc = edge_mpc * 2 WHERE name IN ('A', 'C', 'F')")
        conn.commit()
        query.refresh_sector_system_counts(conn)
    finally:
        conn.close()
    reader = query.open_readonly(mysql_config)
    try:
        yield reader
    finally:
        reader.close()


def _expected(rows, sort, descending, quadrants=()):
    """The ids of `rows` in the order the sort promises, worked out from the rows themselves."""
    def quadrant(row):
        if not row["placed"]:
            return None
        return ("I" if row["center_x_pc"] >= 0 else "II") if row["center_y_pc"] >= 0 else ("IV" if row["center_x_pc"] >= 0 else "III")

    key = {"systems": lambda row: row["system_count"], "density": lambda row: row["system_count"] / row["edge_mpc"] ** 3,
           "position": quadrant}[sort]
    wanted = [row for row in rows if not quadrants or (quadrant(row) or "unplaced") in quadrants]
    with_value = [row for row in wanted if sort == "systems" or row["placed"]]
    without = [row for row in wanted if sort != "systems" and not row["placed"]]
    ordered = sorted(with_value, key=lambda row: (key(row), row["id"]), reverse=descending)
    return [row["id"] for row in ordered + sorted(without, key=lambda row: row["id"], reverse=descending)]


SORTS = [(sort, descending, quadrants) for sort in ("systems", "density", "position") for descending in (False, True)
         for quadrants in ((), ("II",), ("unplaced",), ("I", "unplaced"))]


@pytest.mark.parametrize("sort, descending, quadrants", SORTS)
def test_every_page_by_key_or_by_offset_is_the_order_worked_out_from_the_rows(sectors, sort, descending, quadrants):
    rows = query.list_sectors(sectors)
    assert len(rows) == len(SECTORS)
    everything = [row["id"] for row in query.list_sectors(sectors, sort=sort, descending=descending, quadrants=quadrants)]
    assert everything == _expected(rows, sort, descending, quadrants)
    for size in (1, 2, 4):
        by_offset = []
        for offset in range(0, len(everything) + size, size):
            by_offset += [row["id"] for row in query.list_sectors(
                sectors, limit=size, offset=offset, sort=sort, descending=descending, quadrants=quadrants)]
        assert by_offset == everything, (size, "offset")
        by_key, after = [], None
        for _ in range(len(everything) + 2):
            page = query.list_sectors(sectors, limit=size, offset=5 if after else 0, after=after, sort=sort,
                                      descending=descending, quadrants=quadrants)
            if not page:
                break
            by_key += [row["id"] for row in page]
            after = page[-1]["id"]
        assert by_key == everything, (size, "key")


def test_the_refresh_gives_every_sector_a_row_with_its_count_density_and_quadrant(sectors, mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        query.refresh_sector_system_counts(conn)
        stored = {row["name"]: row for row in conn.execute(
            "SELECT sec.name, c.system_count, c.system_density, c.quadrant, sec.edge_mpc FROM sector_system_counts c "
            "JOIN sectors sec ON sec.id = c.sector_id").fetchall()}
    finally:
        conn.close()
    assert set(stored) == {name for name, _, _ in SECTORS}
    for name, center, count in SECTORS:
        row = stored[name]
        assert row["system_count"] == count and row["quadrant"] == QUADRANT_OF.get(name)
        if center is None:
            assert row["system_density"] is None
        else:
            assert row["system_density"] == pytest.approx(count / row["edge_mpc"] ** 3)


@pytest.mark.parametrize("sort", ["systems", "density", "position"])
def test_a_sector_made_since_the_last_refresh_follows_the_others_and_can_be_paged_to(sectors, mysql_config, sort):
    before = [row["id"] for row in query.list_sectors(sectors, sort=sort)]
    new_id = _save_sector(mysql_config, "New", center_pc=(1.0, 1.0, 0.0), positions=[(0.0, 0.0, 0.0)])
    reader = query.open_readonly(mysql_config)
    try:
        everything = [row["id"] for row in query.list_sectors(reader, sort=sort)]
        assert everything == before + [new_id]
        for size in (1, 3):
            paged, after = [], None
            for _ in range(len(everything) + 2):
                page = query.list_sectors(reader, limit=size, after=after, sort=sort)
                if not page:
                    break
                paged += [row["id"] for row in page]
                after = page[-1]["id"]
            assert paged == everything
        assert [row["id"] for row in query.list_sectors(reader, limit=2, offset=len(before) - 1, sort=sort)] == before[-1:] + [new_id]
    finally:
        reader.close()


@pytest.mark.parametrize("sort", ["systems", "density", "position"])
def test_a_deleted_sector_with_a_stored_row_is_not_listed_and_does_not_break_the_pages(sectors, mysql_config, sort):
    victim = query.list_sectors(sectors, sort=sort)[3]["id"]
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM sectors WHERE id = ?", (victim,))
        conn.commit()
    finally:
        conn.close()
    reader = query.open_readonly(mysql_config)
    try:
        everything = [row["id"] for row in query.list_sectors(reader, sort=sort)]
        assert victim not in everything and len(everything) == len(SECTORS) - 1
        paged, after = [], None
        for _ in range(len(everything) + 2):
            page = query.list_sectors(reader, limit=2, after=after, sort=sort)
            if not page:
                break
            paged += [row["id"] for row in page]
            after = page[-1]["id"]
        assert paged == everything
        # a page asked to start after the deleted sector starts where it stood, while its stored row is there
        assert [row["id"] for row in query.list_sectors(reader, limit=2, after=victim, sort=sort)] == everything[3:5]
    finally:
        reader.close()
