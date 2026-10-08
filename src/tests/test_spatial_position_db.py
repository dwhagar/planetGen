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
    # The columns are the objects' own numbers, in km.
    first = original[0]
    assert row["position_x_km"] == pytest.approx(first.spatial.get_coordinates("system", "cartesian")[0]
                                                 * constants.AU_TO_KM, rel=1e-9)
