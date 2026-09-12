# tests/test_query_db_planets_moons.py

"""
Coverage for `queryDb.list_planets`/`list_moons` -- the `planets`/`moons`
CLI subcommands that close the gap `docs/TODO.md`'s "Open items" > "Search"
flagged: `queryDb.py`'s `systems` subcommand only ever filtered by star
type/sector, with no equivalent at all for a planet's/moon's own class or
size the way the web/API faceted search (`queryDb.search`) already
supports.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable,
same as every other database-backed test in this suite.
"""

import pytest

from queryDb import list_moons, list_planets
from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem


def _make_system_with_moons():
    """Retries generation (bounded) until a system with at least one moon
    and at least two distinct planet classes comes out -- needed so the
    class-filter tests below have both a matching and a non-matching
    planet to distinguish."""
    for _ in range(30):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.MOONS = True
        cfg.MAX_PLANETS = True
        # Pinned False for the same reason test_db_persistence.py's own
        # helper pins it: a binary system's secondary star can add planets
        # this test's own bookkeeping (picking two planets by name/class
        # off `system.planets` alone) doesn't account for.
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        planets = [obj for obj in system.planets if obj.body_type != "a"]
        classes = {p.planet_class for p in planets}
        if len(classes) >= 2 and any(p.moons for p in planets):
            return system, cfg
    pytest.fail("could not generate a system with 2+ planet classes and at least one moon")


def _insert(mysql_config, system, cfg):
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            system_id = _db.insert_star_system(conn, system, cfg)
    finally:
        conn.close()
    return system_id


def test_list_planets_returns_every_planet_with_no_filter(mysql_config):
    system, cfg = _make_system_with_moons()
    system_id = _insert(mysql_config, system, cfg)
    expected_names = sorted(p.name for p in system.planets if p.body_type != "a")

    conn = _db.get_connection(mysql_config)
    try:
        rows = list_planets(conn, system_id=system_id)
    finally:
        conn.close()

    assert sorted(r["name"] for r in rows) == expected_names
    assert all(r["system_name"] == system.name for r in rows)


def test_list_planets_filters_by_class(mysql_config):
    system, cfg = _make_system_with_moons()
    system_id = _insert(mysql_config, system, cfg)
    planets = [p for p in system.planets if p.body_type != "a"]
    target_class = planets[0].planet_class
    expected_names = sorted(p.name for p in planets if p.planet_class == target_class)

    conn = _db.get_connection(mysql_config)
    try:
        rows = list_planets(conn, planet_class=target_class, system_id=system_id)
    finally:
        conn.close()

    assert sorted(r["name"] for r in rows) == expected_names
    assert all(r["planet_class"] == target_class for r in rows)


def test_list_planets_filters_by_radius_range(mysql_config):
    system, cfg = _make_system_with_moons()
    system_id = _insert(mysql_config, system, cfg)
    planets = sorted((p for p in system.planets if p.body_type != "a"), key=lambda p: p.radius_km)
    # A midpoint bound that necessarily excludes the smallest planet (and
    # keeps at least the largest), so the filter is provably doing
    # something rather than vacuously matching everything.
    min_radius = (planets[0].radius_km + planets[-1].radius_km) / 2
    expected_names = sorted(p.name for p in planets if p.radius_km >= min_radius)
    assert expected_names, "test fixture needs at least one planet above the midpoint"
    assert len(expected_names) < len(planets), "test fixture needs at least one planet below the midpoint"

    conn = _db.get_connection(mysql_config)
    try:
        rows = list_planets(conn, min_radius_km=min_radius, system_id=system_id)
    finally:
        conn.close()

    assert sorted(r["name"] for r in rows) == expected_names
    assert all(r["radius_km"] >= min_radius for r in rows)


def test_list_planets_filters_by_sector_id(mysql_config):
    sector = SpaceSector("Query Test Sector", edge_ly=20.0)
    system_in, cfg_in = _make_system_with_moons()
    sector.add_system(system_in, position=(1.0, 2.0, 3.0), system_config=cfg_in)

    cfg_out = SystemConfig()
    cfg_out.STAR_TYPE = "M5V"
    system_out = StarSystem(system_config=cfg_out)

    sector_id = _db.save_sector(sector, config=mysql_config)
    _insert(mysql_config, system_out, cfg_out)

    conn = _db.get_connection(mysql_config)
    try:
        rows = list_planets(conn, sector_id=sector_id)
    finally:
        conn.close()

    expected_names = sorted(p.name for p in system_in.planets if p.body_type != "a")
    assert sorted(r["name"] for r in rows) == expected_names


def test_list_moons_filters_by_class_and_radius(mysql_config):
    system, cfg = _make_system_with_moons()
    system_id = _insert(mysql_config, system, cfg)
    all_moons = [m for p in system.planets for m in p.moons]
    assert all_moons, "test fixture must actually contain moons"
    target_class = all_moons[0].planet_class
    expected_names = sorted(m.name for m in all_moons if m.planet_class == target_class)

    conn = _db.get_connection(mysql_config)
    try:
        rows = list_moons(conn, planet_class=target_class, system_id=system_id)
    finally:
        conn.close()

    assert sorted(r["name"] for r in rows) == expected_names
    assert all(r["planet_class"] == target_class for r in rows)
    # Each moon row also reports its parent planet's name, for display.
    assert all(r["planet_name"] for r in rows)
