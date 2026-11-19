# tests/test_sun_eclipse_rows.py

"""GEN.1434: the eclipses of a planet's other suns by its moons, and of a moon's star by its planet, on the stored sheet (`describe/clock_rows`), against a stored multiple system."""

import math

import pytest

from planetgen.db import store
from planetgen.describe import clock_rows, sheet
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, orbit_elements, sky_clock, sky_suns, sky_view, spin, vec3
from planetgen.util import draw


def _expected(planet):
    """The `(label number, moon name, Eclipse)` a planet's moons give its other suns by the physics alone."""
    light_au = orbit_elements.flux_distance_au(planet.distance, planet.eccentricity or 0.0)
    sky = sky_view.bodies(planet.star, light_au)[1:clock_rows.sky_size_rows.MAX_OTHER_ROWS + 1]
    normal = spin.orbit_normal(planet.orbital_inclination_deg, planet.orbital_ascending_node_deg)
    _, _, placed = sky_suns.placed(planet.star, (planet.spin_axis_x, planet.spin_axis_y, planet.spin_axis_z), normal, light_au)
    found = []
    for body in sky:
        source = next((item for item in placed if item.star is body.star), None)
        if source is None:
            continue
        latitude = math.asin(max(-1.0, min(1.0, vec3.unit(source.offset_au)[2])))
        for moon in planet.moons:
            eclipse = sky_clock.companion_eclipse(math.radians(body.angular_deg) / 2.0, latitude, moon.radius, moon.distance * constants.AU_TO_KM, moon.period, moon.orbital_inclination_deg or 0.0)
            if eclipse:
                found.append((moon.name, eclipse))
    return found


def _binary_with_an_eclipsed_sun():
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = True
    cfg.MOONS = True
    for seed in range(300):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if len(system.stars) < 2:
            continue
        for planet in system.planets:
            if getattr(planet, "body_type", None) == "a" or getattr(planet, "star", None) is None or not getattr(planet, "moons", None) or planet.spin_axis_x is None:
                continue
            if _expected(planet):
                return system, cfg, planet
    pytest.fail("no multiple system with a planet whose moon eclipses another sun in 300 seeds")


@pytest.mark.fresh_database
def test_a_stored_planet_lists_the_eclipses_of_its_other_suns_by_its_moons(mysql_config):
    system, cfg, planet = _binary_with_an_eclipsed_sun()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        row = conn.execute("SELECT id FROM planets WHERE star_system_id = ? AND name = ?", (system_id, planet.name)).fetchone()
        found = [item for item in clock_rows.rows(conn, "planet", row["id"]) if item[1].startswith(clock_rows.SUN_ECLIPSE_KEY_PREFIX)]
        expected = sorted(_expected(planet), key=lambda pair: pair[1].interval_yr)[:clock_rows.MAX_ECLIPSE_ROWS]
        assert [item[4] for item in found] == pytest.approx([eclipse.interval_yr for _, eclipse in expected], rel=1e-3)
        assert all(" by " in item[2] and item[2].startswith("Eclipses of its ") for item in found)
        sheet_keys = [item["key"] for group in sheet.sheet(conn, "planet", row["id"])["groups"] for item in group["rows"]]
        assert any(key.startswith(clock_rows.SUN_ECLIPSE_KEY_PREFIX) for key in sheet_keys)
        moon = conn.execute("SELECT id FROM moons WHERE parent_type = 'planet' AND parent_id = ? ORDER BY id LIMIT 1", (row["id"],)).fetchone()
        by_planet = [item for item in clock_rows.rows(conn, "moon", moon["id"]) if item[1] == f"{clock_rows.ECLIPSE_KEY_PREFIX}1"]
        assert len(by_planet) <= 1 and all(item[2] == "Eclipses of its star by its planet" for item in by_planet)
    finally:
        conn.rollback()
        conn.close()
