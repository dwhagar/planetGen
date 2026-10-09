# tests/test_spatial_position_db.py

"""GEN.74 part 3: the stored columns map to the position objects.

A sector saved with its place in the galaxy loads back with every entry and
body at the same place, and the columns hold exactly what the objects hold
(planet offsets in km, system places in mpc). Needs MySQL (`mysql_config`).
"""

import pytest

from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants

CENTER_PC = (8000.0, -1200.0, 40.0)


def _sector_with_moons():
    sector = SpaceSector("Spatial Round Trip", edge_ly=11.5)
    for _ in range(300):
        cfg = SystemConfig()
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        system = StarSystem(system_config=cfg)
        if any(getattr(p, "moons", None) for p in system.planets):
            break
    sector.add_system(system, position=(1.25, -0.5, 0.75), system_config=cfg)
    sector.place_in_galaxy(tuple(c * constants.PARSEC_M / constants.LIGHTYEAR_M for c in CENTER_PC))
    return sector


def test_a_saved_sector_loads_back_with_every_entry_and_body_in_place(mysql_config):
    sector = _sector_with_moons()
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": CENTER_PC[0], "center_y_pc": CENTER_PC[1], "center_z_pc": CENTER_PC[2],
        "galactic_radius_pc": 8100.0, "ring_index": 3, "layer_index": 0, "ring_slot_index": 5,
    })
    conn = store.get_connection(mysql_config)
    try:
        loaded = store.load_sector(conn, sector_id)
        row = conn.execute("SELECT * FROM planets WHERE star_system_id = (SELECT id FROM star_systems LIMIT 1)"
                           " ORDER BY id LIMIT 1").fetchone()
    finally:
        conn.close()
    before, after = sector.entries[0], loaded.entries[0]
    assert after.position == pytest.approx(before.position, rel=1e-6)
    assert after.spatial.get_coordinates("galactic", "cartesian") == pytest.approx(
        before.spatial.get_coordinates("galactic", "cartesian"), rel=1e-9)
    assert after.spatial.mass_kg == pytest.approx(before.spatial.mass_kg, rel=1e-9)
    planets = [p for p in after.star_system.planets if getattr(p, "spatial", None)]
    original = [p for p in before.star_system.planets if getattr(p, "spatial", None)]
    assert len(planets) == len(original) and planets
    for loaded_planet, planet in zip(planets, original):
        assert loaded_planet.spatial.get_coordinates("galactic", "cartesian") == pytest.approx(
            planet.spatial.get_coordinates("galactic", "cartesian"), rel=1e-9)
        for loaded_moon, moon in zip(loaded_planet.moons, planet.moons):
            assert loaded_moon.spatial.get_coordinates("galactic", "cartesian") == pytest.approx(
                moon.spatial.get_coordinates("galactic", "cartesian"), rel=1e-9)
    # Mass and mu travel with the object, stars included.
    for loaded_planet, planet in zip(planets, original):
        assert loaded_planet.spatial.mass_kg == pytest.approx(planet.spatial.mass_kg, rel=1e-9)
        assert loaded_planet.spatial.mu == pytest.approx(planet.spatial.mu, rel=1e-9)
    for loaded_star, star in zip(after.star_system.stars, before.star_system.stars):
        assert loaded_star.spatial.get_coordinates("galactic", "cartesian") == pytest.approx(
            star.spatial.get_coordinates("galactic", "cartesian"), rel=1e-9)
        assert loaded_star.spatial.mass_kg == pytest.approx(star.spatial.mass_kg, rel=1e-9)
    # The columns are the objects' own numbers, in km.
    first = original[0]
    assert row["position_x_km"] == pytest.approx(first.spatial.get_coordinates("system", "cartesian")[0]
                                                 * constants.AU_TO_KM, rel=1e-9)


def test_a_loaded_sector_knows_every_objects_cell_and_velocity(mysql_config):
    """GEN.121, GEN.124: the entries and the bodies of their systems come back with their sector address and
    the velocity they had when saved."""
    from planetgen.galaxy import geometry

    address = (1500, 2, 700)
    edge_pc = 11.5 * constants.LIGHTYEAR_M / constants.PARSEC_M
    center_pc = geometry.sector_position_pc(*address, edge_pc)
    sector = SpaceSector("Addressed Round Trip", edge_ly=11.5)
    for _ in range(300):
        cfg = SystemConfig()
        cfg.PLANETS = True
        system = StarSystem(system_config=cfg)
        if system.planets:
            break
    sector.add_system(system, position=(1.0, -0.5, 0.5), system_config=cfg)
    sector.place_in_galaxy(tuple(c * constants.PARSEC_M / constants.LIGHTYEAR_M for c in center_pc))
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": center_pc[0], "center_y_pc": center_pc[1], "center_z_pc": center_pc[2],
        "galactic_radius_pc": float((center_pc[0] ** 2 + center_pc[1] ** 2 + center_pc[2] ** 2) ** 0.5),
        "ring_index": address[0], "layer_index": address[1], "ring_slot_index": address[2],
    })
    conn = store.get_connection(mysql_config)
    try:
        loaded = store.load_sector(conn, sector_id)
    finally:
        conn.close()
    before, after = sector.entries[0], loaded.entries[0]
    assert before.spatial.sector_address == after.spatial.sector_address == address
    assert after.spatial.get_velocity_vector("galactic") == pytest.approx(
        before.spatial.get_velocity_vector("galactic"), rel=1e-6)
    assert max(abs(v) for v in after.spatial.get_velocity_vector("galactic")) > 1.0e4  # a star moving ~220 km/s
    for star, original in zip(after.star_system.stars, before.star_system.stars):
        assert star.spatial.sector_address == address
        assert star.spatial.get_velocity_vector("galactic") == pytest.approx(
            original.spatial.get_velocity_vector("galactic"), rel=1e-6)
    planets = [p for p in after.star_system.planets if getattr(p, "spatial", None)]
    originals = [p for p in before.star_system.planets if getattr(p, "spatial", None)]
    assert planets
    for planet, original in zip(planets, originals):
        assert planet.spatial.sector_address == address
        assert planet.spatial.get_velocity_vector("galactic") == pytest.approx(
            original.spatial.get_velocity_vector("galactic"), rel=1e-6, abs=1e-3)


def test_a_loaded_sector_carries_the_time_its_orbits_were_last_advanced(mysql_config):
    """GEN.121: the epoch of every stored position and velocity is when the orbit update last ran."""
    sector = _sector_with_moons()
    sector_id = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        assert store.get_orbit_epoch_unix(conn) is None
        assert store.load_sector(conn, sector_id).entries[0].spatial.epoch_unix is None
        clock = store.orbit_clock(conn, 0.001)
        store.advance_orbital_phases(conn, clock)
        store.finish_orbit_update(conn, clock)
        epoch = store.get_orbit_epoch_unix(conn)
        assert epoch is not None and epoch > 1.7e9
        loaded = store.load_sector(conn, sector_id)
    finally:
        conn.close()
    entry = loaded.entries[0]
    assert entry.spatial.epoch_unix == epoch
    system = entry.star_system
    assert all(star.spatial.epoch_unix == epoch for star in system.stars)
    bodies = [p for p in system.planets if getattr(p, "spatial", None)]
    assert bodies and all(p.spatial.epoch_unix == epoch for p in bodies)
    assert all(m.spatial.epoch_unix == epoch for p in bodies for m in p.moons)
