"""
MAP.179: a node of the Galaxy Map lists the scattered planetary nebulae and supernova remnants no sector has built (`phenomenon_scatter`
rows of those kinds with `built_at` empty) as marks; a built row is a stored nebula, listed in `clouds`, so it is never listed twice.
"""

import numpy as np
import pytest

from planetgen import tuning
from planetgen.db import nodes, query, store
from planetgen.galaxy import nebula_blocks
from planetgen.galaxy.drill import DrillBlock, format_drill_key
from planetgen.galaxy.spherical import sector_address_at

pytestmark = pytest.mark.db

EDGE_PC = float(tuning.DEFAULT_SECTOR_EDGE_PC)
BOX = ((-400.0, -400.0, -100.0), (400.0, 400.0, 100.0))


def _position(index):
    return (40.0 * (index % 10) - 180.0, 25.0 * (index % 7) - 80.0, 6.0 * (index % 5) - 12.0)


def _row(index, kind, subtype, built):
    position = _position(index)
    shell, band, wedge = sector_address_at(position, EDGE_PC)
    return (shell, band, wedge, kind, subtype, *(round(v * 1000) for v in position), index + 1, built)


@pytest.fixture
def scattered(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        rows = [_row(0, "planetary-nebula", None, None), _row(1, "supernova-remnant", None, None), _row(2, "planetary-nebula", None, "2026-01-01 00:00:00"),
                _row(3, "black-hole", "stellar", None), _row(4, "supernova-remnant", None, None)]
        with conn:
            conn.execute("DELETE FROM phenomenon_scatter")
            conn.executemany(
                "INSERT INTO phenomenon_scatter (shell_index, band_index, wedge_index, kind, subtype, position_x_mpc, position_y_mpc,"
                " position_z_mpc, seed, built_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", rows)
    finally:
        conn.close()
    return mysql_config


def test_the_unbuilt_scattered_nebulae_of_a_box_are_marks_and_a_built_one_or_a_black_hole_is_not(scattered):
    conn = store.get_connection(scattered)
    try:
        found = query.galaxy_scattered_nebulae_in_box(conn, *BOX, EDGE_PC, 50)
        assert sorted((p["type"], p["name"]) for p in found) == [
            ("planetary_nebula", "Planetary nebula (uncharted)"), ("supernova_remnant", "Supernova remnant (uncharted)"),
            ("supernova_remnant", "Supernova remnant (uncharted)")]
        assert all(p["scattered"] and not p["charted"] and p["id"].startswith("s") and p["luminosity_sol"] == 0.0 for p in found)
        assert len({p["id"] for p in found}) == 3
        assert len(query.galaxy_scattered_nebulae_in_box(conn, *BOX, EDGE_PC, 2)) == 2
    finally:
        conn.close()


def test_a_node_lists_them_among_its_points_and_a_coarse_node_does_not(scattered, monkeypatch):
    monkeypatch.setattr(nodes, "_outer_shell", lambda conn: 1000)
    conn = store.get_connection(scattered)
    try:
        position = _position(0)
        shell, band, wedge = nebula_blocks.locate(np.array([position[0]]), np.array([position[1]]), np.array([position[2]]), 81, EDGE_PC)
        small = DrillBlock(81, int(shell[0]), int(band[0]), int(wedge[0]))
        listed = nodes.galaxy_nodes(conn, [format_drill_key(small)])["nodes"][format_drill_key(small)]["points"]
        kinds = {p["type"] for p in listed}
        assert "planetary_nebula" in kinds
        assert all(p["id"].startswith("s") for p in listed if p["type"] in ("planetary_nebula", "supernova_remnant"))
        coarse = DrillBlock(729, 1, 0, 0)
        points = nodes.galaxy_nodes(conn, [format_drill_key(coarse)])["nodes"][format_drill_key(coarse)]["points"]
        assert not [p for p in points if p["type"] in ("planetary_nebula", "supernova_remnant")]
    finally:
        conn.close()
