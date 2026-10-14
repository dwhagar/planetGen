"""MAP.186: a binary pair sits in the same place on the Galaxy Map's tiles and on the Sector Map."""

import math

import pytest

from planetgen.db import query
from planetgen.db import store as _db
from planetgen.galaxy.geometry import BINARY_PAIR_DIRECTION, binary_companion_offset_pc
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.web import sector_page
from planetgen.web.maps import starmap

CENTER_PC = (300.0, 400.0, 20.0)  # off the X axis, so the sector's own axes are turned
SYSTEM_POSITION_LY = (1.0, 0.5, -0.4)


def _link(name, **params):
    return f"/{name}"


def _binary_system():
    for _ in range(400):
        cfg = SystemConfig()
        cfg.BINARY_SYSTEM = True
        cfg.PLANETS = False
        system = StarSystem(system_config=cfg)
        if system.binary_type == "close":
            return system, cfg
    pytest.skip("no close binary generated")


@pytest.fixture
def binary_sector(mysql_config):
    system, cfg = _binary_system()
    sector = SpaceSector("Pair Sector", edge_ly=13.0)
    sector.add_system(system, position=SYSTEM_POSITION_LY, system_config=cfg)
    position = {
        "center_x_pc": CENTER_PC[0], "center_y_pc": CENTER_PC[1], "center_z_pc": CENTER_PC[2],
        "galactic_radius_pc": math.dist(CENTER_PC, (0.0, 0.0, 0.0)),
        "ring_index": 80, "layer_index": 5, "ring_slot_index": 3,
    }
    sector_id = _db.save_sector(sector, config=mysql_config, galaxy_position=position)
    conn = _db.get_connection(mysql_config)
    try:
        yield conn, sector_id
    finally:
        conn.close()


def _tile_stars(conn):
    lo = tuple(c - 10.0 for c in CENTER_PC)
    hi = tuple(c + 10.0 for c in CENTER_PC)
    stars = query.galaxy_generated_stars_in_box(conn, lo, hi, 0.0, 1)
    assert len(stars) == 2
    return stars


def _scene_stars_pc(conn, sector_id, compress):
    """The scene's two dots as galaxy-frame parsecs, via the scene's own scale."""
    detail = query.sector_detail(conn, sector_id)
    systems = [sector_page._map_system({**row, "stars": row["stars"]}) for row in detail["systems"]]
    scene = starmap.map_scene_data(
        _link, detail["edge_mpc"], (80, 5, 3), CENTER_PC, systems, compress_binaries=compress)
    half = scene["halfEdgePc"]
    scale = half / scene["sceneHalfPx"]
    return scene, [
        (CENTER_PC[0] + star["x"] * scale, CENTER_PC[1] - star["y"] * scale, CENTER_PC[2] + star["z"] * scale)
        for star in scene["stars"]
    ]


def test_companion_offset_is_the_real_separation_along_one_direction():
    offset = binary_companion_offset_pc(1.5e8)
    assert math.hypot(*offset) == pytest.approx(1.5e8 / 3.0856775814913673e13, rel=1e-6)
    assert offset == pytest.approx(tuple(math.hypot(*offset) * c for c in BINARY_PAIR_DIRECTION))
    assert binary_companion_offset_pc(None) == (0.0, 0.0, 0.0)


def test_tile_stars_sit_where_the_sector_map_puts_the_system(binary_sector):
    conn, sector_id = binary_sector
    tile = {star["name"]: (star["x"], star["y"], star["z"]) for star in _tile_stars(conn)}
    _scene, real = _scene_stars_pc(conn, sector_id, compress=False)
    primary = real[0]
    assert any(point == pytest.approx(primary, abs=0.01) for point in tile.values())
    # Both stars, in either view, are within the pair's real gap (a tiny fraction of a parsec) of the system.
    for point in list(tile.values()) + real:
        assert math.dist(point, primary) < 0.01


def test_the_system_is_rotated_into_the_galaxy_frame_in_the_tiles(binary_sector):
    conn, sector_id = binary_sector
    # Local +X points away from the galactic axis: from (300, 400) that is (0.6, 0.8), not (1, 0).
    row = conn.execute("SELECT position_x_mpc, position_y_mpc, position_z_mpc FROM star_systems"
                       " WHERE sector_id = ?", (sector_id,)).fetchone()
    local = (row["position_x_mpc"] / 1000.0, row["position_y_mpc"] / 1000.0, row["position_z_mpc"] / 1000.0)
    ax, ay = (0.6, 0.8), (-0.8, 0.6)
    expected = (
        CENTER_PC[0] + local[0] * ax[0] + local[1] * ay[0],
        CENTER_PC[1] + local[0] * ax[1] + local[1] * ay[1],
        CENTER_PC[2] + local[2],
    )
    star = _tile_stars(conn)[0]
    assert (star["x"], star["y"], star["z"]) == pytest.approx(expected, abs=0.01)


def test_the_compressed_pair_is_drawn_apart_and_the_scene_says_so(binary_sector):
    conn, sector_id = binary_sector
    compressed, _points = _scene_stars_pc(conn, sector_id, compress=True)
    real, _points = _scene_stars_pc(conn, sector_id, compress=False)
    assert compressed["binariesCompressed"] is True and real["binariesCompressed"] is False
    gap_px = math.dist((compressed["stars"][0]["x"], compressed["stars"][0]["y"]),
                       (compressed["stars"][1]["x"], compressed["stars"][1]["y"]))
    real_gap_px = math.dist((real["stars"][0]["x"], real["stars"][0]["y"]),
                            (real["stars"][1]["x"], real["stars"][1]["y"]))
    assert gap_px > 0.5
    assert real_gap_px < 0.01
