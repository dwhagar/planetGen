# tests/test_gen_pair_positions.py

"""
GEN.709: the places stored for the stars of a pair (the columns the pair classes carry, `binary_*_position_*` for a close pair,
`*_position_*_au` for a wide one) are the places the system's tree of orbits (`orbit_nodes`) puts them at its epoch. They
once differed by the node's argument of periapsis (a circular pair), by the node's ellipse (a wide pair), and for a pair
inside a larger tree by the whole tree (the scene took the pair's own frame for the system's).
"""

import math

import pytest

from planetgen.db import orbit_trees
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, system_tree
from planetgen.util import draw


def _system(seed, stars=None):
    cfg = SystemConfig()
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = True
    if stars is not None:
        cfg.STARS = stars
    with draw.bound(seed):
        return StarSystem(system_config=cfg), cfg


def _tree_offsets(system):
    slots = system.star_slots or ((0, 0), (0, 1))
    return system_tree.star_offsets_au([star.mass / constants.SOLAR_MASS_TO_KG for star in system.stars], slots,
                                       system.orbit_nodes, system.orbit_nodes[0].epoch_unix)


@pytest.mark.parametrize("stars", [2, 3, 5])
def test_the_stored_star_places_are_the_places_of_the_tree(stars):
    checked = {"close": 0, "wide": 0}
    for seed in range(60):
        system, _cfg = _system(seed, stars)
        if len(system.stars) != stars:
            continue
        places = orbit_trees.star_places(system)
        want = _tree_offsets(system)
        for (_node, got, _slot), expected in zip(places, want):
            assert got == pytest.approx(expected, abs=1e-6 * max(1.0, max(abs(c) for c in expected)))
        checked[system.binary_type] += 1
    assert sum(checked.values()) >= 5, checked


def test_a_close_pair_is_tracked_by_the_phase_its_node_puts_it_at():
    for seed in range(80):
        system, _cfg = _system(seed, 2)
        if system.binary_type != "close":
            continue
        (node,) = system.orbit_nodes
        assert node.eccentricity == 0.0 and node.arg_periapsis_deg == 0.0
        assert system.star.binary_mutual_orbital_phase_deg == pytest.approx(node.mean_anomaly_deg)
        return
    pytest.fail("no close pair in 80 draws")


def _scene_of(system, cfg, mysql_config):
    from planetgen.db import store
    from planetgen.galaxy.sector import SpaceSector
    from planetgen.web.maps.systemscene import build_scene

    sector = SpaceSector("Pair Positions", edge_ly=11.5)
    sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        return build_scene(conn, conn.execute("SELECT id FROM star_systems ORDER BY id DESC LIMIT 1").fetchone()["id"])
    finally:
        conn.close()


@pytest.mark.parametrize("stars", [2, 3, 5])
def test_the_scene_puts_a_close_pairs_stars_where_the_tree_does(mysql_config, stars):
    from planetgen.physics.body_positions import positions_at

    for seed in range(120):
        system, cfg = _system(seed, stars)
        if len(system.stars) == stars and system.binary_type == "close":
            break
    else:
        pytest.fail(f"no close pair in a system of {stars} stars in 120 draws")
    scene = _scene_of(system, cfg, mysql_config)
    here = positions_at(scene, 0.0)
    close = [star for star in scene["stars"] if star["role"] in ("primary", "secondary")]
    assert len(close) == 2
    for star in close:
        assert here[star["ref"]] == pytest.approx(star["position_km"], rel=1e-6, abs=1e4)
