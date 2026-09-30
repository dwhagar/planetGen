# tests/test_galactic_motion.py

"""The correlative update's galactic motion (TODO item 32): everything
moves along its galactic orbit and is refiled into the sector it drifts
into, with octant, location, containment and nearest systems following."""

import math

import pytest

from stellarObjects import _db
from stellarObjects.config import SystemConfig
from stellarObjects.galaxyGeometry import (
    galaxy_to_local_pc, local_to_galaxy_pc, sector_position_pc, slot_angle_bounds,
)
from stellarObjects.spaceSector import SpaceSector
from stellarObjects.systemData import StarSystem
from stellarObjects.utils import pc_to_ly

EDGE_PC = 4.0
RING = 10


def test_rotation_follows_galactic_phase():
    x, y, z = _db._rotate_about_axis((10.0, 0.0, 1.0), math.pi / 2)
    assert (x, y, z) == pytest.approx((0.0, 10.0, 1.0))
    assert _db._galactic_turn(250e6, 0.25) == pytest.approx(2 * math.pi)
    assert _db._galactic_turn(1.0, None) == 0.0


def _place(slot, positions_ly, name):
    center = sector_position_pc(RING, 0, slot, EDGE_PC)
    sector = SpaceSector(name, edge_ly=pc_to_ly(EDGE_PC))
    for position in positions_ly:
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = False
        cfg.BINARY_SYSTEM = False
        sector.add_system(StarSystem(system_config=cfg), position=position, system_config=cfg)
    return sector, {"center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
                    "galactic_radius_pc": math.hypot(*center), "ring_index": RING, "layer_index": 0,
                    "ring_slot_index": slot}


def test_a_system_that_drifts_across_a_boundary_is_refiled(mysql_config):
    # A system just short of slot 0's leading edge (local +Y is the
    # direction of travel), and a second sector ahead of it in slot 1.
    edge_ly = pc_to_ly(EDGE_PC)
    sector, position = _place(0, [(0.0, edge_ly / 2 - 0.3, 0.0), (0.0, -edge_ly / 2 + 0.5, 0.0)], "Trailing")
    first = _db.save_sector(sector, config=mysql_config, galaxy_position=position)
    sector, position = _place(1, [(0.0, 0.0, 0.0)], "Leading")
    second = _db.save_sector(sector, config=mysql_config, galaxy_position=position)

    conn = _db.get_connection(mysql_config)
    try:
        mover, stay = [r["id"] for r in conn.execute(
            "SELECT id FROM star_systems WHERE sector_id = ? ORDER BY position_y_mpc DESC", (first,)).fetchall()]
        with conn:
            conn.execute("UPDATE stars SET galactic_orbital_period_gy = 1.0")
        # Turn far enough to carry the mover about a light-year forward.
        radius_pc = math.hypot(*sector_position_pc(RING, 0, 0, EDGE_PC))
        elapsed = 1e9 * (1.0 / 3.2616) / (2 * math.pi * radius_pc)
        with conn:
            motion = _db.advance_galactic_positions(conn, elapsed)
            rewritten = _db.refresh_after_motion(conn, motion["sectors"])
        assert motion["refiled"] >= 1 and {first, second} <= motion["sectors"]
        row = conn.execute("SELECT * FROM star_systems WHERE id = ?", (mover,)).fetchone()
        assert row["sector_id"] == second
        assert row["location"].startswith("Leading")
        assert rewritten >= 1
        # It sits where the turn put it, in its new sector's frame.
        center = sector_position_pc(RING, 0, 1, EDGE_PC)
        point = local_to_galaxy_pc(center, (row["position_x_mpc"] / 1000, row["position_y_mpc"] / 1000,
                                            row["position_z_mpc"] / 1000))
        low, high = slot_angle_bounds(RING, 1)
        assert low <= math.atan2(point[1], point[0]) <= high
        assert conn.execute("SELECT sector_id FROM star_systems WHERE id = ?", (stay,)).fetchone()["sector_id"] == first
        near = conn.execute("SELECT sector_id FROM nearest_systems WHERE object_table = 'star_systems'"
                            " AND object_id = ? LIMIT 1", (mover,)).fetchone()
        assert near["sector_id"] == second
    finally:
        conn.close()


def test_phenomena_and_facilities_move_and_the_quasar_stays(mysql_config):
    sector, position = _place(3, [(0.0, 0.0, 0.0)], "Drift")
    sector_id = _db.save_sector(sector, config=mysql_config, galaxy_position=position)
    conn = _db.get_connection(mysql_config)
    try:
        from stellarObjects.nebulaData import Nebula
        center = (position["center_x_pc"], position["center_y_pc"], position["center_z_pc"])
        with conn:
            nebula = Nebula(SystemConfig())
            nebula.galactic_orbital_period_gy = 1.0
            nebula_id = _db.insert_nebula(conn, nebula, sector_id=sector_id, placement={
                "center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
                "galactic_radius_pc": math.hypot(*center)})
            facility_id = _db.add_facility(conn, "Drifter", "station", "standalone", "space", sector_id)
            motion = _db.advance_galactic_positions(conn, 1e6)
        assert motion["moved"] >= 3
        row = conn.execute("SELECT center_x_pc, center_y_pc FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
        angle = math.atan2(row["center_y_pc"], row["center_x_pc"]) - math.atan2(center[1], center[0])
        assert angle == pytest.approx(2 * math.pi * 1e6 / 1e9)
        moved = conn.execute("SELECT center_x_pc FROM facilities WHERE id = ?", (facility_id,)).fetchone()
        assert moved["center_x_pc"] != pytest.approx(center[0], abs=1e-12)
    finally:
        conn.close()


def test_orbital_facilities_advance_their_phase(mysql_config):
    from tests.test_facilities import _ids, _system_with_worlds
    system_id = _system_with_worlds(mysql_config)
    conn = _db.get_connection(mysql_config)
    try:
        ids = _ids(conn, system_id)
        with conn:
            facility_id = _db.add_facility(conn, "Ring", "station", "orbital", "planet", ids["giant"], phase_deg=0.0)
        period = conn.execute("SELECT orbit_period_years FROM facilities WHERE id = ?", (facility_id,)).fetchone()[
            "orbit_period_years"]
        with conn:
            assert _db.advance_facility_orbits(conn, period / 4) == 1
        phase = conn.execute("SELECT orbit_phase_deg FROM facilities WHERE id = ?", (facility_id,)).fetchone()
        assert phase["orbit_phase_deg"] == pytest.approx(90.0)
    finally:
        conn.close()


def test_update_orbits_runs_end_to_end(mysql_config, monkeypatch, capsys):
    """The whole correlative update against a real current-schema
    database: a generated sector with planets and phenomena, run twice
    (the first run only sets the clock)."""
    import updateOrbits

    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    sector = SpaceSector("Clockwork", edge_ly=pc_to_ly(EDGE_PC))
    sector.add_system(StarSystem(system_config=cfg), position=(0.0, 0.0, 0.0), system_config=cfg)
    _sector, position = _place(5, [], "unused")
    _db.save_sector(sector, config=mysql_config, galaxy_position=position)

    argv = ["updateOrbits.py", "--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
            "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password or "",
            "--mysql-database", mysql_config.database]
    monkeypatch.setattr("sys.argv", argv)
    updateOrbits.main()
    conn = _db.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("UPDATE orbit_simulation_state SET last_updated_at = NOW() - INTERVAL 400 DAY")
    finally:
        conn.close()
    updateOrbits.main()
    out = capsys.readouterr().out
    assert "Moved 1 object(s) along their galactic orbits; 0 changed sector." in out
