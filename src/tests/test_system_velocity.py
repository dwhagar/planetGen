# tests/test_system_velocity.py

"""GEN.121: a star system's galactic velocity (rotation plus runaway motion) is set at placement, stored on
`star_systems`, turned with the position by the galactic orbit update and put back on a loaded sector."""

import math
import random

import pytest

from planetgen.db import store
from planetgen.galaxy import geometry
from planetgen.galaxy.sector import SpaceSector
from planetgen.galaxy.system_position import galactic_velocity_ms, peculiar_velocity_ms, random_unit_vector
from planetgen.generation import run_sector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants

ADDRESS = (1500, 2, 700)
EDGE_LY = 11.5
EDGE_PC = EDGE_LY * constants.LIGHTYEAR_M / constants.PARSEC_M


def _norm(v):
    return math.sqrt(sum(c * c for c in v))


def _runaway_system(speed_kms=400.0, direction=(0.0, 0.0, 1.0)):
    for _ in range(300):
        cfg = SystemConfig()
        cfg.PLANETS = True
        system = StarSystem(system_config=cfg)
        if system.planets:
            break
    system.runaway_class = "hypervelocity"
    system.runaway_speed_kms = speed_kms
    system.runaway_direction = direction
    return system, cfg


def _sector(system, cfg):
    center_pc = geometry.sector_position_pc(*ADDRESS, EDGE_PC)
    sector = SpaceSector("Velocity Round Trip", edge_ly=EDGE_LY)
    sector.add_system(system, position=(1.0, -0.5, 0.5), system_config=cfg)
    sector.place_in_galaxy(tuple(c * constants.PARSEC_M / constants.LIGHTYEAR_M for c in center_pc))
    return sector, center_pc


def _save(sector, center_pc, mysql_config):
    return store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": center_pc[0], "center_y_pc": center_pc[1], "center_z_pc": center_pc[2],
        "galactic_radius_pc": _norm(center_pc),
        "ring_index": ADDRESS[0], "layer_index": ADDRESS[1], "ring_slot_index": ADDRESS[2],
    })


def test_a_random_unit_vector_is_a_direction_and_follows_the_generator_it_is_given():
    vectors = [random_unit_vector() for _ in range(50)]
    assert all(_norm(v) == pytest.approx(1.0) for v in vectors)
    assert random_unit_vector(random.Random(7)) == random_unit_vector(random.Random(7))
    assert peculiar_velocity_ms(object()) == (0.0, 0.0, 0.0)


def test_flagged_fast_stars_get_a_direction_and_it_survives_a_dict_round_trip():
    system, cfg = _runaway_system()
    system.runaway_class = system.runaway_speed_kms = system.runaway_direction = None
    sector = SpaceSector("Fast", edge_ly=EDGE_LY)
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
    for _ in range(2000):
        run_sector.flag_fast_stars(sector)
        if system.runaway_class:
            break
    assert system.runaway_class
    assert _norm(system.runaway_direction) == pytest.approx(1.0)
    again = StarSystem.from_dict(system.to_dict())
    assert again.runaway_direction == pytest.approx(system.runaway_direction)


def test_a_runaway_system_moves_beyond_the_rotation_curve_and_its_stars_and_bodies_share_it():
    system, cfg = _runaway_system(speed_kms=400.0, direction=(0.0, 0.0, 1.0))
    sector, _center = _sector(system, cfg)
    entry = sector.entries[0]
    rotation = galactic_velocity_ms(entry.spatial.get_coordinates("galactic", "cartesian"),
                                    system.star.galactic_orbital_speed_kms)
    velocity = entry.spatial.get_velocity_vector("galactic")
    assert velocity == pytest.approx((rotation[0], rotation[1], rotation[2] + 400.0e3), rel=1e-9, abs=1e-6)
    for star in system.stars:
        assert star.spatial.get_velocity_vector("galactic")[2] == pytest.approx(400.0e3, rel=1e-3)
    planet = next(p for p in system.planets if getattr(p, "spatial", None))
    # A planet's galactic velocity is its star's plus its own orbital speed, so the runaway motion is in it.
    assert planet.spatial.get_velocity_vector("galactic")[2] == pytest.approx(400.0e3, abs=planet.orbital_speed_kms * 1e3 * 1.01)


def test_a_saved_runaway_system_loads_back_with_the_velocity_it_had(mysql_config):
    system, cfg = _runaway_system(speed_kms=250.0, direction=(0.6, 0.0, 0.8))
    sector, center_pc = _sector(system, cfg)
    sector_id = _save(sector, center_pc, mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        row = conn.execute("SELECT velocity_x_kms, velocity_y_kms, velocity_z_kms FROM star_systems LIMIT 1").fetchone()
        loaded = store.load_sector(conn, sector_id)
    finally:
        conn.close()
    before = sector.entries[0].spatial.get_velocity_vector("galactic")
    assert (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"]) == pytest.approx(
        tuple(v / 1000.0 for v in before), rel=1e-9, abs=1e-9)
    after = loaded.entries[0]
    assert after.spatial.get_velocity_vector("galactic") == pytest.approx(before, rel=1e-6, abs=1e-3)
    for star, original in zip(after.star_system.stars, sector.entries[0].star_system.stars):
        assert star.spatial.get_velocity_vector("galactic") == pytest.approx(
            original.spatial.get_velocity_vector("galactic"), rel=1e-6, abs=1e-3)


def test_the_galactic_orbit_update_turns_the_velocity_with_the_position(mysql_config):
    system, cfg = _runaway_system(speed_kms=300.0, direction=(0.0, 0.0, 1.0))
    sector, center_pc = _sector(system, cfg)
    sector_id = _save(sector, center_pc, mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        def velocity():
            row = conn.execute("SELECT velocity_x_kms, velocity_y_kms, velocity_z_kms FROM star_systems").fetchone()
            return (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"])

        start = velocity()
        period_gy = conn.execute(
            "SELECT s.galactic_orbital_period_gy AS p FROM stars s JOIN star_systems ss ON ss.id = s.star_system_id"
            " WHERE s.role IN ('single', 'primary')").fetchone()["p"]
        store.advance_galactic_positions(conn, elapsed_years=period_gy * 1e9 / 4.0)
        turned = velocity()
        assert _norm(turned) == pytest.approx(_norm(start), rel=1e-9)
        assert turned[2] == pytest.approx(start[2], rel=1e-9)  # the axis is unchanged
        assert turned[0] == pytest.approx(-start[1], rel=1e-6, abs=1e-6)  # a quarter turn counterclockwise
        assert turned[1] == pytest.approx(start[0], rel=1e-6, abs=1e-6)
        assert store.load_sector(conn, sector_id).entries[0].spatial.get_velocity_vector("galactic") == pytest.approx(
            tuple(v * 1000.0 for v in turned), rel=1e-6, abs=1e-3)
    finally:
        conn.close()


def test_migrate_v60_to_v61_works_out_the_velocity_of_systems_saved_before(mysql_config):
    ordinary = _runaway_system()
    ordinary[0].runaway_class = ordinary[0].runaway_speed_kms = ordinary[0].runaway_direction = None
    sector, center_pc = _sector(*ordinary)
    _save(sector, center_pc, mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        before = conn.execute("SELECT velocity_x_kms AS x, velocity_y_kms AS y, velocity_z_kms AS z"
                              " FROM star_systems").fetchone()
        conn.execute("ALTER TABLE star_systems DROP COLUMN velocity_x_kms, DROP COLUMN velocity_y_kms,"
                     " DROP COLUMN velocity_z_kms")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (60)")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        after = conn.execute("SELECT velocity_x_kms AS x, velocity_y_kms AS y, velocity_z_kms AS z"
                             " FROM star_systems").fetchone()
    finally:
        conn.close()
    assert (after["x"], after["y"], after["z"]) == pytest.approx((before["x"], before["y"], before["z"]),
                                                                rel=1e-6, abs=1e-9)
    assert _norm((after["x"], after["y"], after["z"])) > 100.0  # a star on the rotation curve, ~220 km/s
