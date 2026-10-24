# tests/test_tile_reads.py

"""
PERF.157: the two pieces of a Galaxy Map tile that read the most rows when the warm-up built the tiles
(`galaxy_scattered_points_in_box` and `galaxy_bright_stars_in_box`) answer the same tile whichever way they read it, and
a box picks the way that reads fewer rows.
"""

import pytest

from planetgen.db import query, store
from planetgen.galaxy.spherical import sector_address_at
from tests.test_db_persistence import _bright_galaxy

pytestmark = pytest.mark.fresh_database

EDGE_PC = 4.0
BOXES = [
    ((0.0, 0.0, -1000.0), (100.0, 100.0, 1000.0)),
    ((-5000.0, -5000.0, -1000.0), (5000.0, 5000.0, 1000.0)),
    ((0.0, 0.0, 300.0), (100.0, 100.0, 1000.0)),
    ((-256.0, -256.0, -256.0), (0.0, 0.0, 0.0)),
    ((200.0, -400.0, -50.0), (700.0, -100.0, 50.0)),
    ((-30.0, -30.0, -30.0), (30.0, 30.0, 30.0)),
]


def _scatter_rows(conn):
    rows = []
    for index in range(240):
        position = (30.0 * (index % 24) - 330.0, 17.0 * (index % 13) - 100.0, 9.0 * (index % 7) - 30.0)
        if index % 5 == 0:
            position = (position[0] * 10, position[1] * 10, position[2] * 10)
        shell, band, wedge = sector_address_at(position, EDGE_PC)
        kind, subtype = (("black-hole", "stellar"), ("neutron-star", None), ("black-hole", "supermassive"),
                         ("black-hole", "intermediate"))[index % 4]
        rows.append((shell, band, wedge, kind, subtype, *(round(v * 1000) for v in position), index + 1,
                     1.4 + index * 0.37))
    with conn:
        conn.execute("DELETE FROM phenomenon_scatter")
        conn.executemany(
            "INSERT INTO phenomenon_scatter (shell_index, band_index, wedge_index, kind, subtype, position_x_mpc,"
            " position_y_mpc, position_z_mpc, seed, mass_solar) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", rows)


def test_scattered_points_are_the_same_with_one_range_per_shell_or_one_for_all(mysql_config, monkeypatch):
    conn = store.get_connection(mysql_config)
    try:
        _scatter_rows(conn)
        for coarse in (False, True):
            for lo, hi in BOXES:
                for limit in (3, 50):
                    answers = []
                    for ranges in (0, query.SCATTERED_MAX_SHELL_RANGES):
                        monkeypatch.setattr(query, "SCATTERED_MAX_SHELL_RANGES", ranges)
                        answers.append(query.galaxy_scattered_points_in_box(conn, lo, hi, EDGE_PC, limit, coarse))
                    assert answers[0] == answers[1]
        points = query.galaxy_scattered_points_in_box(conn, *BOXES[1], EDGE_PC, 50, False)
        assert len(points) > 10, "the boxes hold points, so the comparison means something"
    finally:
        conn.close()


def test_bright_stars_are_the_same_walking_the_luminosity_index_or_reading_the_addresses(mysql_config, monkeypatch):
    """Where both ways see the same candidates (a box holding fewer stars than a pool), the choice does not show."""
    conn = store.get_connection(mysql_config)
    try:
        _bright_galaxy(conn)
        for lo, hi in BOXES:
            for limit in (4, 12, 1000):
                answers = []
                for walk in (True, False):
                    monkeypatch.setattr(query, "_luminosity_walk_is_cheaper", lambda *args, walk=walk: walk)
                    answers.append(query.galaxy_bright_stars_in_box(conn, lo, hi, EDGE_PC, limit=limit))
                assert answers[0] == answers[1]
        assert query.galaxy_bright_stars_in_box(conn, *BOXES[1], EDGE_PC, limit=30)
    finally:
        conn.close()


def test_the_walk_is_chosen_for_a_box_with_many_stars_and_the_addresses_for_a_sliver():
    """A box holding a good share of the stars finds a list's worth quickly; one holding a sliver reads its own
    addresses, however high its area share looks (PERF.157: a box a kiloparsec above the plane)."""
    assert query._luminosity_walk_is_cheaper(800_000, 1.0, 1.0, 40)
    assert not query._luminosity_walk_is_cheaper(800_000, 0.016, 0.125, 50)
