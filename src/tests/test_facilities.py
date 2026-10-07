# tests/test_facilities.py

"""Facilities (schema v42): the placement rules, the orbit
math, storage and the API."""

import pytest

from planetgen.db import store
from planetgen.population import facilities
from planetgen.physics import constants
from planetgen.generation.config import SystemConfig
from planetgen.physics.planets import calculate_orbital_period_years
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem
from planetgen.physics.orbits import circular_orbital_speed_kms
from planetgen.db import query


@pytest.mark.parametrize("kind, placement, host_type, body_type, allowed", [
    ("colony", "terrestrial", "planet", "t", True),
    ("colony", "terrestrial", "planet", "g", False),   # gas giants: orbital only
    ("starbase", "orbital", "planet", "g", True),
    ("colony", "terrestrial", "moon", "t", True),
    ("station", "orbital", "moon", "t", True),
    ("outpost", "orbital", "star", None, True),
    ("colony", "orbital", "star", None, False),
    ("mining-colony", "asteroid", "asteroid_belt", None, True),
    ("mining-colony", "asteroid", "asteroid_field", None, False),
    ("outpost", "asteroid", "asteroid_field", None, True),
    ("starbase", "standalone", "space", None, True),
    ("colony", "standalone", "space", None, False),
    ("outpost", "terrestrial", "star", None, False),
])
def test_placement_rules(kind, placement, host_type, body_type, allowed):
    assert (facilities.check_facility(kind, placement, host_type, body_type) is None) == allowed


def test_orbit_matches_how_moons_orbit():
    earth_mass, earth_radius = 5.972e24, 6371.0
    orbit = facilities.orbit_for(earth_mass, earth_radius, 42164.0)
    distance_au = 42164.0 / constants.AU_TO_KM
    period = calculate_orbital_period_years(distance_au, earth_mass)
    assert orbit["period_years"] == pytest.approx(period)
    assert orbit["orbital_speed_kms"] == pytest.approx(circular_orbital_speed_kms(distance_au, period))
    # A geostationary orbit: about a day, about 3.07 km/s.
    assert orbit["period_years"] * 365.25 == pytest.approx(1.0, rel=0.01)
    assert orbit["orbital_speed_kms"] == pytest.approx(3.07, rel=0.01)
    assert facilities.orbit_for(earth_mass, earth_radius)["distance_km"] == pytest.approx(3 * earth_radius)
    with pytest.raises(ValueError):
        facilities.orbit_for(earth_mass, earth_radius, 1000.0)


def test_host_placements_follow_the_rules():
    assert facilities.host_placements("planet", "t") == ["terrestrial", "orbital"]
    assert facilities.host_placements("planet", "g") == ["orbital"]
    assert facilities.host_placements("moon", "t") == ["terrestrial", "orbital"]
    assert facilities.host_placements("star") == ["orbital"]
    assert facilities.host_placements("asteroid_belt") == ["asteroid"]


def test_orbit_slider_is_logarithmic_from_surface_to_sphere():
    lowest, highest = facilities.orbit_limits(6371.0, 1.5e6)
    assert lowest == pytest.approx(6371.0 * 1.01) and highest == 1.5e6
    assert facilities.distance_from_step(0, lowest, highest) == pytest.approx(lowest)
    assert facilities.distance_from_step(1000, lowest, highest) == pytest.approx(highest)
    # Halfway along is the geometric mean, and each step the same factor.
    assert facilities.distance_from_step(500, lowest, highest) == pytest.approx((lowest * highest) ** 0.5)
    assert facilities.distance_from_step(-5, lowest, highest) == pytest.approx(lowest)
    for step in (0, 1, 250, 999, 1000):
        distance = facilities.distance_from_step(step, lowest, highest)
        assert facilities.step_for_distance(distance, lowest, highest) == step
    # No stored sphere: a fixed number of host radii.
    assert facilities.orbit_limits(6371.0, None)[1] == pytest.approx(6371.0 * 1000)
    assert facilities.orbit_limits(6371.0, 100.0)[1] == pytest.approx(6371.0 * 1000)


def test_orbit_must_stay_inside_the_sphere_of_influence():
    with pytest.raises(ValueError, match="sphere of influence"):
        facilities.orbit_for(5.972e24, 6371.0, 2.0e6, highest_km=1.5e6)
    # The default orbit comes in to the sphere's edge for a tiny sphere.
    assert facilities.orbit_for(5.972e24, 6371.0, highest_km=10_000.0)["distance_km"] == pytest.approx(10_000.0)


def test_star_host_orbits_a_close_pair_as_one():
    stars = [{"id": 1, "mass_kg": 2.0e30, "radius_km": 7.0e5, "heliosphere_radius_km": 1.0e10},
             {"id": 2, "mass_kg": 1.0e30, "radius_km": 5.0e5, "heliosphere_radius_km": 8.0e9}]
    assert facilities.star_host(stars, 2, "wide") == (1.0e30, 5.0e5, 8.0e9)
    assert facilities.star_host(stars, 1, "close", 3.0e6, 2.0e10) == pytest.approx((3.0e30, 3.7e6, 2.0e10))


def _system_with_worlds(mysql_config):
    """A saved single-star system with a terrestrial planet, a gas giant
    and a belt -- retried, since generation is random."""
    for _ in range(200):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        types = {p.body_type for p in system.planets}
        if {"t", "g", "a"} <= types:
            return store.save_system(system, cfg, config=mysql_config)
    pytest.fail("could not generate a system with a terrestrial planet, a gas giant and a belt")


def _ids(conn, system_id):
    one = lambda sql: conn.execute(sql, (system_id,)).fetchone()["id"]  # noqa: E731
    return {
        "star": one("SELECT id FROM stars WHERE star_system_id = ?"),
        "terrestrial": one("SELECT id FROM planets WHERE star_system_id = ? AND body_type = 't' LIMIT 1"),
        "giant": one("SELECT id FROM planets WHERE star_system_id = ? AND body_type = 'g' LIMIT 1"),
        "belt": one("SELECT id FROM asteroid_belts WHERE star_system_id = ? LIMIT 1"),
    }


def test_facilities_are_stored_on_their_hosts(mysql_config):
    system_id = _system_with_worlds(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        ids = _ids(conn, system_id)
        with conn:
            colony = store.add_facility(conn, "New Hope", "colony", "terrestrial", "planet", ids["terrestrial"])
            yard = store.add_facility(conn, "High Yard", "starbase", "orbital", "planet", ids["giant"],
                                    distance_km=5.0e5, phase_deg=370.0)
            store.add_facility(conn, "Sunwatch", "outpost", "orbital", "star", ids["star"],
                             distance_km=constants.AU_TO_KM * 0.5)
            store.add_facility(conn, "Rockpile", "mining-colony", "asteroid", "asteroid_belt", ids["belt"])
        with pytest.raises(store.FacilityError):
            store.add_facility(conn, "Floaters", "colony", "terrestrial", "planet", ids["giant"])
        with pytest.raises(store.FacilityError, match="sphere of influence"):
            store.add_facility(conn, "Runaway", "station", "orbital", "planet", ids["terrestrial"],
                             distance_km=constants.AU_TO_KM)
        with pytest.raises(store.FacilityError):
            store.add_facility(conn, "Pinned", "outpost", "asteroid", "asteroid_belt", ids["belt"], distance_km=1.0e8)
        with pytest.raises(store.FacilityError) as missing:
            store.add_facility(conn, "Nowhere", "colony", "terrestrial", "planet", 10 ** 12)
        assert missing.value.not_found
        conn.rollback()

        listed = {f["name"]: f for f in query.facilities_for_system(conn, system_id)}
        assert set(listed) == {"New Hope", "High Yard", "Sunwatch", "Rockpile"}
        assert listed["High Yard"]["orbit_phase_deg"] == pytest.approx(10.0)
        assert listed["High Yard"]["orbit_period_years"] > 0
        assert listed["New Hope"]["orbit_distance_km"] is None
        # A belt facility gets a spot in the belt and an orbit around the
        # star from there (ADM.9).
        belt = conn.execute("SELECT lower_limit_km, upper_limit_km FROM asteroid_belts WHERE id = ?",
                            (ids["belt"],)).fetchone()
        rockpile = listed["Rockpile"]
        assert belt["lower_limit_km"] <= rockpile["orbit_distance_km"] <= belt["upper_limit_km"]
        assert rockpile["orbit_period_years"] > 0 and rockpile["orbital_speed_kms"] > 0
        assert 0 <= rockpile["orbit_phase_deg"] < 360
        assert listed["New Hope"]["host_id"] == ids["terrestrial"]
        assert query.colonized_body_ids(conn, system_id)["planets"] == {ids["terrestrial"]}
        assert query.facility_detail(conn, yard)["kind"] == "starbase"

        # Position updates move orbital and belt facilities alike.
        quarter = rockpile["orbit_period_years"] / 4
        with conn:
            assert store.advance_facility_orbits(conn, quarter) == 3
        moved = conn.execute("SELECT orbit_phase_deg FROM facilities WHERE name = 'Rockpile'").fetchone()
        assert (moved["orbit_phase_deg"] - rockpile["orbit_phase_deg"]) % 360.0 == pytest.approx(90.0)

        orbit = store.facility_orbit(conn, "planet", ids["terrestrial"])
        hill = conn.execute("SELECT radius_km, hill_radius_km FROM planets WHERE id = ?",
                            (ids["terrestrial"],)).fetchone()
        assert orbit["min_distance_km"] == pytest.approx(hill["radius_km"] * 1.01)
        assert orbit["max_distance_km"] == pytest.approx(hill["hill_radius_km"])

        with conn:
            assert store.delete_facility(conn, colony)
            conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,))
        assert conn.execute("SELECT COUNT(*) AS n FROM facilities").fetchone()["n"] == 0
    finally:
        conn.close()


def test_standalone_facilities_park_in_their_sector(mysql_config):
    sector = SpaceSector("Harbor", edge_ly=13.0)
    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position={
        "center_x_pc": 0.0, "center_y_pc": 50.0, "center_z_pc": 0.0, "galactic_radius_pc": 50.0,
    })
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            store.add_facility(conn, "Waypoint", "station", "standalone", "space", sector_id, offset_ly=(1.0, 0.0, 0.0))
        with pytest.raises(store.FacilityError):
            store.add_facility(conn, "Too Far", "station", "standalone", "space", sector_id, offset_ly=(9.0, 0.0, 0.0))
        conn.rollback()
        (waypoint,) = query.facilities_in_sector(conn, sector_id)
        # Local +X points away from the galactic axis (+y here).
        assert waypoint["center_y_pc"] == pytest.approx(50.0 + 1.0 / 3.2616, rel=1e-3)
    finally:
        conn.close()


def test_migrate_v41_to_v42_creates_facilities(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DROP TABLE facilities")
        conn.execute("DELETE FROM schema_migrations WHERE version IN (42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (41)")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM facilities").fetchone()["n"] == 0
    finally:
        conn.close()
