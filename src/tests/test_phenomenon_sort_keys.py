# tests/test_phenomenon_sort_keys.py

"""
PERF.211: the Phenomena table's sorts by type, descriptor, radius, placed and sector read the built phenomena in index order. A page
from the top, and each page after a built row (`after=(type, id)`), is the page of the one long list the whole union gives sorted
(`limit=None`), ties settling by name, type and id in the direction of the sort, a phenomenon with no sector after those with one.
Each phenomenon table has the indexes the walk needs, and they serve the walk's query without a sort.
"""

import os
import random

import pytest

from planetgen.db import query, store

from tests.db_schema_support import insert_filled_row

pytestmark = pytest.mark.fresh_database

NAMES = ("Alpha", "Beta", "Gamma", "Delta", "Epsilon")
TABLES = {            # table -> (descriptor column and its values, or None for the ones the query works out)
    "nebulae": ("nebula_type", ("emission", "reflection", "dark")),
    "asteroid_fields": ("density", ("dense", "sparse", "typical")),
    "black_holes": (None, None),
    "neutron_stars": ("pulsar_type", ("young", "millisecond", "magnetar", "non-pulsing")),
    "supernova_remnants": ("morphology", ("shell", "plerion", "composite")),
    "rogue_planets": (None, None),
    "interstellar_comets": (None, None),
    "quasars": (None, None),
}
REQUIRED = {          # the columns the table's CHECKs refuse to take a neutral value for
    "supernova_remnants": {"progenitor_type": "Type Ia"},
    "black_holes": {"parent_type": "free", "mass_class": "stellar"},
    "neutron_stars": {"parent_type": "free"},
    "rogue_planets": {"parent_type": "free"},
    "interstellar_comets": {"parent_type": "free"},
}
SECTOR_NAMES = ("Sector B", "Sector D", "Sector A", "Sector C")


def _row(table, rng, sector_ids, serial):
    """The columns of a phenomenon row that the sorts read: names that repeat across the tables, a few descriptors and radii (so
    that keys tie), a sector or none, a position or none."""
    row = {"name": f"{rng.choice(NAMES)} {rng.randrange(3)}", "uid": os.urandom(12),
           "sector_id": rng.choice(sector_ids + [None, None]),
           **rng.choice([{"center_x_pc": 1.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 1.0},
                         {"center_x_pc": None, "center_y_pc": None, "center_z_pc": None, "galactic_radius_pc": None}])}
    column, values = TABLES[table]
    if column:
        row[column] = rng.choice(values)
    if table in ("nebulae", "asteroid_fields", "supernova_remnants"):
        row["radius_ly"] = rng.choice([0.5, 1.0, 1.0, 2.0, 7.5])
    if table == "black_holes":
        row["has_accretion_disk"] = rng.randrange(2)
    if table == "rogue_planets":
        row["mass_bin"] = rng.choice(["terrestrial", "brown-dwarf", "jupiter"])
        row["planet_type"] = rng.choice(["t", "g"])
    if table == "interstellar_comets":
        row["is_active"] = rng.randrange(2)
    if table == "quasars":
        row["is_radio_loud"] = rng.randrange(2)
    if table in ("black_holes", "neutron_stars"):
        row["star_id"] = rng.choice([None, None, None, 1])      # one in four belongs to a star system and is left out
    return row


@pytest.fixture
def phenomena(fresh_mysql_config):
    conn = store.get_connection(fresh_mysql_config)
    try:
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        rng = random.Random(211)
        sector_ids = [conn.execute("INSERT INTO sectors (name, edge_mpc) VALUES (?, 4000)", (name,)).lastrowid
                      for name in SECTOR_NAMES]
        for table in TABLES:
            for serial in range(14):
                insert_filled_row(conn, table, **REQUIRED.get(table, {}), **_row(table, rng, sector_ids, serial))
        conn.commit()
    finally:
        conn.close()
    reader = query.open_readonly(fresh_mysql_config)
    try:
        yield reader
    finally:
        reader.close()


def _ids(rows):
    return [(row["type"], row["id"]) for row in rows]


def _paged(conn, size, **kwargs):
    """Every built row, `size` at a time: the first page from the top, each later one after the last row of the one before."""
    rows, after = [], None
    for _ in range(500):
        page = query.list_phenomena(conn, limit=size, after=after, **kwargs)
        rows += page
        if len(page) < size:
            return rows
        after = (page[-1]["type"], page[-1]["id"])
    raise AssertionError("the pages never ended")


FILTERS = (
    {}, {"placed": True}, {"placed": False}, {"types": ["nebula", "rogue_planet", "quasar"]},
    {"descriptors": ["active", "dense", "shell", "terrestrial"]}, {"types": ["black_hole"], "placed": True},
)


@pytest.mark.parametrize("sort", ["name", "type", "descriptor", "radius", "placed", "sector"])
@pytest.mark.parametrize("descending", [False, True])
def test_every_page_by_key_is_the_page_of_the_whole_sorted_union(phenomena, sort, descending):
    for filters in FILTERS:
        whole = query.list_phenomena(phenomena, sort=sort, descending=descending, **filters)
        assert whole, filters
        for size in (1, 3, 7, 100):
            assert _ids(_paged(phenomena, size, sort=sort, descending=descending, **filters)) == _ids(whole), (filters, size)
            assert _ids(query.list_phenomena(phenomena, limit=size, sort=sort, descending=descending, **filters)) == \
                _ids(whole[:size]), (filters, size, "first page")


@pytest.mark.parametrize("walk", [(0, 1), (0, 2), (10 ** 9, 5000)])
@pytest.mark.parametrize("descending", [False, True])
def test_the_sector_sort_walks_the_sectors_or_sorts_them_whole_and_agrees(phenomena, monkeypatch, walk, descending):
    """The walk reads the sectors a few at a time (one, two) and the whole sort is used under the row floor; either is the order of
    the sort."""
    minimum, chunk = walk
    whole = query.list_phenomena(phenomena, sort="sector", descending=descending)
    monkeypatch.setattr(query, "SECTOR_WALK_MIN_ROWS", minimum)
    monkeypatch.setattr(query, "SECTOR_CHUNK", chunk)
    for filters in FILTERS:
        wanted = query.list_phenomena(phenomena, sort="sector", descending=descending, **filters)
        for size in (1, 4, 9):
            assert _ids(_paged(phenomena, size, sort="sector", descending=descending, **filters)) == _ids(wanted), (filters, size)
    assert _ids(_paged(phenomena, 5, sort="sector", descending=descending)) == _ids(whole)


def test_the_order_is_the_one_the_sort_promises(phenomena):
    """The keys themselves: radius ascending puts the point-like types (0 ly) first, placed puts the unplaced first, the sector sort puts
    the rows that have no sector last either way, and a tie follows the direction of the sort."""
    for descending in (False, True):
        radius = query.list_phenomena(phenomena, sort="radius", descending=descending)
        values = [row["radius_ly"] for row in radius]
        assert values == sorted(values, reverse=descending)
        placed = query.list_phenomena(phenomena, sort="placed", descending=descending)
        flags = [row["placed"] for row in placed]
        assert flags == sorted(flags, reverse=descending)
        sectors = [row["sector_name"] for row in query.list_phenomena(phenomena, sort="sector", descending=descending)]
        named = [name for name in sectors if name is not None]
        assert sectors[:len(named)] == named and named == sorted(named, reverse=descending)
        assert len(named) < len(sectors)
        descriptors = [row["descriptor"] for row in query.list_phenomena(phenomena, sort="descriptor", descending=descending)]
        assert descriptors == sorted(descriptors, key=str.lower, reverse=descending)
    names = [(row["name"], row["type"], row["id"]) for row in query.list_phenomena(phenomena, sort="type")
             if row["type"] == "nebula"]
    assert names == sorted(names)
    names = [(row["name"], row["type"], row["id"]) for row in query.list_phenomena(phenomena, sort="type", descending=True)
             if row["type"] == "nebula"]
    assert names == sorted(names, reverse=True)


def test_a_cursor_on_a_row_that_is_gone_starts_at_the_top(phenomena):
    rows = query.list_phenomena(phenomena, limit=3, sort="descriptor")
    assert _ids(query.list_phenomena(phenomena, limit=3, sort="descriptor", after=("nebula", 10 ** 9))) == _ids(rows)
    assert _ids(query.list_phenomena(phenomena, limit=3, sort="sector", after=("quasar", 10 ** 9))) == \
        _ids(query.list_phenomena(phenomena, limit=3, sort="sector"))


@pytest.mark.parametrize("sort", ["placed", "descriptor", "radius"])
def test_each_tables_index_serves_the_walk_without_a_sort(fresh_mysql_config, phenomena, sort):
    """`EXPLAIN` of the query one table runs for a page: it reads `idx_<table>_<sort>_name` and has no filesort (a table with the
    same key for every row, the type and a radius of 0, walks the name index and is not asked for)."""
    conn = store.get_connection(fresh_mysql_config)
    try:
        for table, type_label, descriptor_expr, radius_expr in query._PHENOMENON_TABLES:
            if sort == "radius" and radius_expr == "0":
                continue
            key = query._built_sort_key(sort, type_label, descriptor_expr, radius_expr)[0]
            lead = "t.star_id IS NULL AND " if table in ("black_holes", "neutron_stars") else ""
            index = f"idx_{table}_{sort}_name"
            for direction in ("ASC", "DESC"):       # forced: a table this small is scanned whatever the indexes are, and a key that is not
                plan = conn.execute(               # the index's expression cannot be read in its order, so the plan sorts
                    f"EXPLAIN SELECT t.id FROM {table} t FORCE INDEX ({index}) WHERE {lead}{key} IS NOT NULL "
                    f"ORDER BY {key} {direction}, t.name {direction}, t.id {direction} LIMIT 10").fetchall()
                assert plan[0]["key"] == index, (table, sort, plan[0])
                assert "filesort" not in (plan[0]["Extra"] or ""), (table, sort, direction, plan[0])
    finally:
        conn.close()
