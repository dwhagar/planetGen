# tests/test_sector_paths_db.py

"""GEN.123: sector paths are computed from the stored sector, saved as spline knots and read back. Needs MySQL."""

import math
import sys

import pytest

from planetgen.db import sector_paths, store
from planetgen.galaxy import geometry
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.compact_remnant import BlackHole
from planetgen.generation.system import StarSystem
from planetgen.physics import constants

ADDRESS = (1500, 2, 700)
EDGE_LY = 11.5
EDGE_PC = EDGE_LY * constants.LIGHTYEAR_M / constants.PARSEC_M


def _save(mysql_config, with_black_hole=True, systems=2):
    center_pc = geometry.sector_position_pc(*ADDRESS, EDGE_PC)
    sector = SpaceSector("Path Round Trip", edge_ly=EDGE_LY)
    for index in range(systems):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        system.name = f"Pathstar {index}"
        sector.add_system(system, position=(float(index) - 0.5, 0.0, 0.0), system_config=cfg)
    if with_black_hole:
        sector.add_phenomenon(BlackHole(SystemConfig()), "black-hole", position=(0.0, 2.0, 0.0))
    sector.place_in_galaxy(tuple(c * constants.PARSEC_M / constants.LIGHTYEAR_M for c in center_pc))
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": center_pc[0], "center_y_pc": center_pc[1], "center_z_pc": center_pc[2],
        "galactic_radius_pc": math.sqrt(sum(c * c for c in center_pc)),
        "ring_index": ADDRESS[0], "layer_index": ADDRESS[1], "ring_slot_index": ADDRESS[2],
    })
    return sector_id, center_pc


def test_bodies_with_nothing_near_them_get_straight_two_knot_paths(mysql_config):
    sector_id, center_pc = _save(mysql_config, with_black_hole=False)
    conn = store.get_connection(mysql_config)
    try:
        assert sector_paths.compute_sector_paths(conn, sector_id) == 2
        conn.commit()
        paths = sector_paths.load_sector_paths(conn, sector_id)
        systems = conn.execute("SELECT id, position_x_mpc, position_y_mpc, position_z_mpc FROM star_systems"
                               " ORDER BY id").fetchall()
    finally:
        conn.close()
    assert sorted(paths) == [("star_systems", row["id"]) for row in systems]
    for row in systems:
        path = paths[("star_systems", row["id"])]
        assert len(path.knots) == 2 and path.duration_s > 0
        start = geometry.local_to_galaxy_pc(center_pc, tuple(row[f"position_{a}_mpc"] / 1000.0 for a in "xyz"))
        assert path.start.position == pytest.approx(tuple(c * constants.PARSEC_M for c in start), rel=1e-9)
        speed = math.sqrt(sum(v * v for v in path.start.velocity))
        assert speed > 1.0e5  # a star on the galaxy's rotation
        assert path.end.velocity == pytest.approx(path.start.velocity, rel=1e-9)  # no mass bent it


def test_a_heavy_black_hole_bends_the_path_of_a_slow_neighbour_and_recomputing_replaces_the_rows(mysql_config):
    sector_id, center_pc = _save(mysql_config, systems=1)
    ex, _ey, _ez = geometry.sector_orientation(center_pc)
    conn = store.get_connection(mysql_config)
    try:
        # A stellar black hole made heavy, at the sector's center; the star 0.8 pc from it drifts at 20 km/s.
        conn.execute("UPDATE black_holes SET mass_solar = 1.0e7, center_x_pc = ?, center_y_pc = ?, center_z_pc = ?",
                     tuple(center_pc))
        conn.execute("UPDATE star_systems SET position_x_mpc = -300, position_y_mpc = 800, position_z_mpc = 0,"
                     " velocity_x_kms = ?, velocity_y_kms = ?, velocity_z_kms = ?",
                     tuple(20.0 * c for c in ex))
        conn.commit()
        assert sector_paths.compute_sector_paths(conn, sector_id) == 1
        conn.commit()
        path = sector_paths.load_sector_paths(conn, sector_id).popitem()[1]
        again = sector_paths.compute_sector_paths(conn, sector_id)
        conn.commit()
        counts = (conn.execute("SELECT COUNT(*) AS n FROM sector_paths").fetchone()["n"],
                  conn.execute("SELECT COUNT(*) AS n FROM sector_path_knots").fetchone()["n"])
    finally:
        conn.close()
    assert len(path.knots) > 2
    cos = sum(a * b for a, b in zip(path.start.velocity, path.end.velocity)) / (
        math.sqrt(sum(v * v for v in path.start.velocity)) * math.sqrt(sum(v * v for v in path.end.velocity)))
    assert math.degrees(math.acos(max(-1.0, min(1.0, cos)))) > 10.0
    assert again == 1 and counts == (1, len(path.knots))


def test_the_saved_spline_passes_through_its_knots_and_a_missing_sector_is_an_error(mysql_config):
    sector_id, _center = _save(mysql_config, with_black_hole=False, systems=1)
    conn = store.get_connection(mysql_config)
    try:
        sector_paths.compute_sector_paths(conn, sector_id)
        conn.commit()
        path = sector_paths.load_sector_paths(conn, sector_id).popitem()[1]
        with pytest.raises(ValueError):
            sector_paths.compute_sector_paths(conn, sector_id + 1000)
    finally:
        conn.close()
    for knot in path.knots:
        assert path.position_at(knot.t_s) == pytest.approx(knot.position, rel=1e-12)
    assert path.exited


def test_deleting_a_sector_deletes_its_paths(mysql_config):
    sector_id, _center = _save(mysql_config, with_black_hole=False, systems=1)
    conn = store.get_connection(mysql_config)
    try:
        sector_paths.compute_sector_paths(conn, sector_id)
        conn.execute("DELETE FROM sectors WHERE id = ?", (sector_id,))
        conn.commit()
        assert conn.execute("SELECT COUNT(*) AS n FROM sector_paths").fetchone()["n"] == 0
        assert conn.execute("SELECT COUNT(*) AS n FROM sector_path_knots").fetchone()["n"] == 0
    finally:
        conn.close()


def test_settling_covers_a_sector_and_its_neighbours_and_nothing_far(mysql_config, monkeypatch):
    near_id, _center = _save(mysql_config, systems=1)
    monkeypatch.setattr(sys.modules[__name__], "ADDRESS", (1500, 2, 701))
    neighbour_id, _center = _save(mysql_config, with_black_hole=False, systems=1)
    monkeypatch.setattr(sys.modules[__name__], "ADDRESS", (1500, 2, 760))
    far_id, _center = _save(mysql_config, with_black_hole=False, systems=1)
    conn = store.get_connection(mysql_config)
    try:
        assert sector_paths.sectors_to_settle(conn, [near_id]) == sorted([near_id, neighbour_id])
        progress = []
        saved = sector_paths.settle_sectors(conn, [near_id], lambda done, total: progress.append((done, total)))
        assert saved == 2 and progress == [(1, 2), (2, 2)]  # the sector's system and the neighbour's
        assert sorted(sector_paths.load_sector_paths(conn, neighbour_id)) != []
        assert sector_paths.load_sector_paths(conn, far_id) == {}
        assert sector_paths.sector_ids_since(conn, "2999-01-01") == []
        assert far_id in sector_paths.sector_ids_since(conn, "2000-01-01")
    finally:
        conn.close()
