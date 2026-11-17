# tests/test_sector_stats_updates.py

"""
PERF.275: `store.set_sector_bright_levels` and `set_sector_bright_masses` write the sectors that take one value in one statement (a sky share wrote one statement for every cell it
visited). What they store is what one statement a sector stored; and `spherical.sector_position_pc` is cached without changing what it gives.
"""

from planetgen.db import store
from planetgen.galaxy import spherical
from tests.test_bright_star_scatter import EDGE_PC
from tests.test_sky_charting import SPHERE_CENTER, galaxy  # noqa: F401 (the fixture)


def _addresses(conn, count):
    """`count` cells of a sphere, with their `sector_stats` rows made (`lock_sector_stats` makes the row of a cell that has none)."""
    cells = [tuple(row[:3]) for row in spherical.enumerate_sectors_within_radius(SPHERE_CENTER, 3.0 * EDGE_PC, EDGE_PC)][:count]
    store.lock_sector_stats(conn, [(address, 1.0, 1.0) for address in cells])
    conn.commit()
    return cells


def _stored(conn, addresses):
    out = {}
    for address in addresses:
        row = conn.execute("SELECT bright_level_sol, bright_mass_sol FROM sector_stats WHERE shell_index = ? AND band_index = ? AND wedge_index = ?", address).fetchone()
        out[address] = (row["bright_level_sol"], row["bright_mass_sol"])
    return out


def test_levels_are_stored_a_sector_each_whatever_the_statements_are_cut_in(galaxy, monkeypatch):
    conn = galaxy
    addresses = _addresses(conn, 9)
    assert len(addresses) == 9
    monkeypatch.setattr(store, "SECTOR_STATS_UPDATE_ROWS", 2)           # five statements for the sectors that take one level
    levels = {address: (500.0 if index % 3 else 1000.0) for index, address in enumerate(addresses)}
    store.set_sector_bright_levels(conn, levels)
    conn.commit()
    assert {address: level for address, (level, _mass) in _stored(conn, addresses).items()} == levels


def test_masses_are_stored_with_their_levels_and_a_none_level_leaves_the_stored_one(galaxy, monkeypatch):
    conn = galaxy
    addresses = _addresses(conn, 7)
    monkeypatch.setattr(store, "SECTOR_STATS_UPDATE_ROWS", 2)
    store.set_sector_bright_levels(conn, {address: 123.0 for address in addresses})
    states = {address: ((5.0, 250.0) if index % 2 else (8.0, None)) for index, address in enumerate(addresses)}
    store.set_sector_bright_masses(conn, states)
    conn.commit()
    stored = _stored(conn, addresses)
    for address, (mass, level) in states.items():
        assert stored[address] == (123.0 if level is None else level, mass)


def test_the_cached_sector_center_is_the_center_it_was_and_a_bad_address_still_raises():
    spherical.sector_position_pc.cache_clear()
    first = spherical.sector_position_pc(4, 0, 5, EDGE_PC)
    assert spherical.sector_position_pc(4, 0, 5, EDGE_PC) == first == spherical.sector_position_pc.__wrapped__(4, 0, 5, EDGE_PC)
    for _ in range(2):
        try:
            spherical.sector_position_pc(4, 0, 9999, EDGE_PC)
        except ValueError:
            continue
        raise AssertionError("an address that is not a cell is not a cell the second time either")
