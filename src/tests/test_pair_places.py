# tests/test_pair_places.py

"""
GEN.640 slice a: the system's own pair is read from its orbit tree. `orbit_trees.pair_places_km` gives the secondary's position relative to the primary and each star's place about the pair's center
of mass from the pair node, and it equals what the `star_systems.binary_mutual_*`, `binary_primary_position_*` and `binary_secondary_position_*` columns hold for a close pair of any tree; for a wide
pair, whose old column placed the secondary on a circle of its separation, the node's ellipse is the position (the orbit the scene's tree already draws). The orbit update turns the node, so the
position follows it and no pair column is written.
"""

import math

import pytest

from planetgen.db import orbit_trees, store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

AU_KM = 1.495978707e8


def _store(conn, system, config):
    with conn:
        sector_id = conn.execute(
            "INSERT INTO sectors (name, edge_mpc, center_x_pc, center_y_pc, center_z_pc, galactic_radius_pc)"
            " VALUES (?, ?, ?, ?, ?, ?)", (f"Pair {system.name}", 11500, 8000.0, 0.0, 0.0, 8000.0)).lastrowid
        return store.insert_star_system(conn, system, config, sector_id=sector_id, position=(1.0, 2.0, 3.0), location="x")


def _system(seed, stars, wide):
    config = SystemConfig()
    config.STARS = stars
    config.PLANETS = False
    config.WIDE_BINARY = wide
    with draw.bound(seed):
        return StarSystem(config), config


def _columns(conn, system_id):
    row = conn.execute("SELECT * FROM star_systems WHERE id = ?", (system_id,)).fetchone()
    return {name: [row[f"binary_{name}_position_{axis}_km"] for axis in "xyz"] for name in ("mutual", "primary", "secondary")}


@pytest.mark.parametrize("stars", [2, 3, 4])
def test_a_close_pairs_places_are_what_the_columns_held(mysql_config, stars):
    conn = store.get_connection(mysql_config)
    try:
        checked = 0
        for seed in range(30):
            system, config = _system(seed, stars, wide=False)
            if system.binary_type != "close":
                continue
            system_id = _store(conn, system, config)
            held, read = _columns(conn, system_id), orbit_trees.pair_places_km(conn, system_id)
            for key, name in (("mutual", "mutual_km"), ("primary", "primary_km"), ("secondary", "secondary_km")):
                assert read[name] == pytest.approx(held[key], rel=1e-6, abs=1.0), (seed, key)
            checked += 1
        assert checked >= 3
    finally:
        conn.close()


def test_a_wide_pair_is_placed_on_its_ellipse_with_the_columns_circle_as_its_mean_distance(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        checked = 0
        for seed in range(30):
            system, config = _system(seed, 2, wide=True)
            if system.binary_type != "wide":
                continue
            system_id = _store(conn, system, config)
            held, read = _columns(conn, system_id), orbit_trees.pair_places_km(conn, system_id)
            node = read["node"]
            radius = math.dist(read["mutual_km"], (0, 0, 0))
            assert node.separation_km * (1.0 - node.eccentricity) * 0.999999 <= radius <= node.separation_km * (1.0 + node.eccentricity) * 1.000001
            if node.eccentricity < 1.0e-9:
                assert read["mutual_km"] == pytest.approx(held["mutual"], rel=1e-6, abs=1.0)
            checked += 1
        assert checked >= 3
    finally:
        conn.close()


def test_a_single_star_has_no_pair(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        for seed in range(40):
            system, config = _system(seed, 1, wide=None)
            if system.binary_type is None:
                assert orbit_trees.pair_places_km(conn, _store(conn, system, config)) is None
                return
        pytest.fail("no single star drawn")
    finally:
        conn.close()


def test_the_orbit_update_moves_the_pair_the_node_places(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        for seed in range(40):
            system, config = _system(seed, 2, wide=False)
            if system.binary_type == "close":
                break
        system_id = _store(conn, system, config)
        before = orbit_trees.pair_places_km(conn, system_id)
        clock = store.orbit_clock(conn, 3.25 * before["node"].period_years)
        store.advance_orbital_phases(conn, clock)
        conn.commit()
        after = orbit_trees.pair_places_km(conn, system_id)
        assert after["node"].mean_anomaly_deg != before["node"].mean_anomaly_deg
        assert after["mutual_km"] != before["mutual_km"]
        assert math.dist(after["mutual_km"], (0, 0, 0)) == pytest.approx(math.dist(before["mutual_km"], (0, 0, 0)), rel=1e-9)
    finally:
        conn.close()
