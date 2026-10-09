# tests/test_orbit_velocity.py

"""GEN.121: every body's velocity, set at generation, stored, kept in step with its position and loaded back."""

import math

import pytest

from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.galaxy.system_position import galactic_velocity_ms
from planetgen.generation.comet import Comet
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.system import StarSystem
from planetgen.physics import constants as pc
from planetgen.physics import state_vectors as sv

KM_S_PER_AU_YR = pc.AU_TO_KM / pc.SECONDS_PER_YEAR


def _speed(body):
    return math.sqrt(body.velocity_x_kms ** 2 + body.velocity_y_kms ** 2 + body.velocity_z_kms ** 2)


def _system(predicate):
    for _ in range(300):
        system = StarSystem(system_config=SystemConfig())
        if predicate(system):
            return system
    pytest.fail("no system matching")


def _with_planets_and_comets(system):
    return any(isinstance(p, Planet) for p in system.planets) and bool(system.comets)


def test_a_planet_moon_and_comet_move_at_their_orbital_speed_along_their_orbit():
    system = _system(lambda s: _with_planets_and_comets(s) and any(getattr(p, "moons", None) for p in s.planets))
    planet = next(p for p in system.planets if isinstance(p, Planet) and p.moons)
    for body in (planet, planet.moons[0], system.comets[0]):
        assert _speed(body) == pytest.approx(body.orbital_speed_kms, rel=1e-9)
    for body in (planet, planet.moons[0]):  # on a circle the velocity is perpendicular to the radius
        radius = (body.position_x, body.position_y, body.position_z)
        velocity = (body.velocity_x_kms, body.velocity_y_kms, body.velocity_z_kms)
        assert sum(r * v for r, v in zip(radius, velocity)) == pytest.approx(0.0, abs=1e-6 * body.distance * _speed(body))
    assert planet.spatial.get_speed("system") == pytest.approx(planet.orbital_speed_kms * 1000.0, rel=1e-9)


def test_a_comet_velocity_conserves_the_orbits_energy():
    system = _system(_with_planets_and_comets)
    comet = system.comets[0]
    mu = 4.0 * math.pi ** 2 * comet.primary_mass_solar
    position = (comet.position_x_au, comet.position_y_au, comet.position_z_au)
    velocity = tuple(v / KM_S_PER_AU_YR for v in (comet.velocity_x_kms, comet.velocity_y_kms, comet.velocity_z_kms))
    elements = sv.elements_from_state(position, velocity, mu)
    assert elements["periapsis_distance"] == pytest.approx(comet.perihelion_distance_au, rel=1e-6)
    assert elements["inclination"] == pytest.approx(math.radians(comet.inclination_deg), rel=1e-6, abs=1e-9)


def test_the_velocity_survives_a_dict_round_trip_and_staging():
    system = _system(_with_planets_and_comets)
    planet = next(p for p in system.planets if isinstance(p, Planet))
    clone = Planet.from_dict(planet.to_dict(), system.star, system.system_config)
    assert (clone.velocity_x_kms, clone.velocity_y_kms, clone.velocity_z_kms) == pytest.approx(
        (planet.velocity_x_kms, planet.velocity_y_kms, planet.velocity_z_kms))
    comet = system.comets[0]
    again = Comet.from_dict(comet.to_dict(), system.system_config)
    assert (again.velocity_x_kms, again.velocity_y_kms, again.velocity_z_kms) == pytest.approx(
        (comet.velocity_x_kms, comet.velocity_y_kms, comet.velocity_z_kms))
    assert again.spatial.get_speed("system") == pytest.approx(comet.spatial.get_speed("system"))


def test_a_body_with_no_position_stages_its_velocity_until_it_has_one():
    comet = Comet(SystemConfig(), primary_mass_solar=1.0, orbit_type="elliptical")
    comet.spatial = None
    comet.velocity_x_kms, comet.velocity_y_kms, comet.velocity_z_kms = 1.0, 2.0, 3.0
    assert comet.spatial is None and comet.velocity_y_kms == 2.0
    comet.set_position_au(1.0, 0.0, 0.0)
    assert comet.spatial.get_velocity_vector("system") == pytest.approx((1000.0, 2000.0, 3000.0))
    assert comet.velocity_z_kms == pytest.approx(3.0)


def test_stars_move_on_the_rotation_curve_and_their_bodies_carry_it():
    speed = 220.0
    assert galactic_velocity_ms((8000.0, 0.0, 5.0), speed) == pytest.approx((0.0, 220_000.0, 0.0))
    assert galactic_velocity_ms((0.0, -3.0, 0.0), speed) == pytest.approx((220_000.0, 0.0, 0.0))
    assert galactic_velocity_ms((0.0, 0.0, 9.0), speed) == (0.0, 0.0, 0.0)
    assert galactic_velocity_ms((1.0, 0.0, 0.0), None) == (0.0, 0.0, 0.0)
    assert galactic_velocity_ms((1.0, 0.0, 0.0), 0.0) == (0.0, 0.0, 0.0)

    system = _system(_with_planets_and_comets)
    sector = SpaceSector("Moving", edge_ly=13.05)
    sector.place_in_galaxy((26000.0, 0.0, 0.0))
    entry = sector.add_system(system, position=(1.0, 0.0, 0.0))
    assert system.star.galactic_orbital_speed_kms > 0
    galactic = entry.spatial.get_velocity_vector("galactic")
    assert math.sqrt(sum(v * v for v in galactic)) == pytest.approx(system.star.galactic_orbital_speed_kms * 1000.0, rel=1e-6)
    planet = next(p for p in system.planets if isinstance(p, Planet))
    # A planet's galactic velocity is its star's plus its own round it; its system velocity is its own.
    own = tuple(c * 1000.0 for c in (planet.velocity_x_kms, planet.velocity_y_kms, planet.velocity_z_kms))
    assert planet.spatial.get_velocity_vector("system") == pytest.approx(own)
    anchor = entry.spatial.get_velocity_vector("galactic") if not getattr(system, "binary_type", None) == "wide" else None
    if anchor is not None:
        assert planet.spatial.get_velocity_vector("galactic") == pytest.approx(
            tuple(a + o for a, o in zip(anchor, own)), abs=1e-6)
    for moon in getattr(planet, "moons", [])[:1]:
        moon_own = tuple(c * 1000.0 for c in (moon.velocity_x_kms, moon.velocity_y_kms, moon.velocity_z_kms))
        assert moon.spatial.get_velocity_vector("galactic") == pytest.approx(
            tuple(a + o for a, o in zip(planet.spatial.get_velocity_vector("galactic"), moon_own)), abs=1e-6)


def test_a_phenomenon_entry_moves_with_its_own_galactic_speed():
    from planetgen.generation.phenomena.compact_remnant import BlackHole

    sector = SpaceSector("Hole", edge_ly=13.05)
    hole = BlackHole(SystemConfig())
    sector.place_in_galaxy((0.0, 26000.0, 0.0))
    entry = sector.add_phenomenon(hole, "black-hole", position=(0.0, 0.0, 0.0))
    assert hole.galactic_orbital_speed_kms > 0
    velocity = entry.spatial.get_velocity_vector("galactic")
    assert velocity[0] < 0.0 and velocity[1] == pytest.approx(0.0, abs=1e-6)
    assert math.sqrt(sum(v * v for v in velocity)) == pytest.approx(hole.galactic_orbital_speed_kms * 1000.0, rel=1e-6)


# --- The database ------------------------------------------------------------------------

def _saved_system(mysql_config, with_moons=True):
    system = _system(lambda s: _with_planets_and_comets(s) and (not with_moons or any(
        getattr(p, "moons", None) for p in s.planets)))
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, system.system_config)
    finally:
        conn.close()
    return system, system_id


def test_a_saved_system_keeps_every_bodys_velocity(mysql_config):
    system, system_id = _saved_system(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        planets = conn.execute("SELECT * FROM planets WHERE star_system_id = ? ORDER BY orbital_index",
                               (system_id,)).fetchall()
        moons = conn.execute("SELECT * FROM moons WHERE star_system_id = ?", (system_id,)).fetchall()
        comets = conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchall()
        for row in list(planets) + list(moons) + list(comets):
            speed = math.sqrt(row["velocity_x_kms"] ** 2 + row["velocity_y_kms"] ** 2 + row["velocity_z_kms"] ** 2)
            assert speed == pytest.approx(row["orbital_speed_kms"], rel=1e-9)
        loaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()
    original = next(p for p in system.planets if isinstance(p, Planet) and p.moons)
    reloaded = next(p for p in loaded.planets if isinstance(p, Planet) and p.name == original.name)
    assert (reloaded.velocity_x_kms, reloaded.velocity_y_kms, reloaded.velocity_z_kms) == pytest.approx(
        (original.velocity_x_kms, original.velocity_y_kms, original.velocity_z_kms))
    assert reloaded.spatial.get_speed("system") == pytest.approx(original.spatial.get_speed("system"))
    again = loaded.comets[0]
    assert again.spatial.get_speed("system") == pytest.approx(system.comets[0].spatial.get_speed("system"))


def test_advancing_the_orbits_moves_the_velocity_with_the_position(mysql_config):
    from planetgen.physics.orbits import circular_orbital_velocity_au_per_year

    _system, system_id = _saved_system(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        store.advance_orbital_phases(conn, elapsed_years=0.37)
        store.advance_comet_orbits(conn, elapsed_years=0.37)
        for table in ("planets", "moons"):
            for row in conn.execute(f"SELECT * FROM {table} WHERE star_system_id = ?", (system_id,)).fetchall():
                expected = circular_orbital_velocity_au_per_year(
                    row["distance_km"] / pc.AU_TO_KM, row["orbital_inclination_deg"], row["orbital_ascending_node_deg"],
                    row["orbital_phase_deg"], row["period_years"])
                got = (row["velocity_x_kms"], row["velocity_y_kms"], row["velocity_z_kms"])
                assert got == pytest.approx(tuple(v * KM_S_PER_AU_YR for v in expected), rel=1e-9, abs=1e-9)
        for row in conn.execute("SELECT * FROM comets WHERE star_system_id = ?", (system_id,)).fetchall():
            speed = math.sqrt(row["velocity_x_kms"] ** 2 + row["velocity_y_kms"] ** 2 + row["velocity_z_kms"] ** 2)
            assert speed == pytest.approx(row["orbital_speed_kms"], rel=1e-9)
            assert speed > 0.0
    finally:
        conn.close()


