# tests/test_trojan_planet.py

"""
GEN.1151 slice A: a giant whose 1:1 zone grew a body of planet mass holds it as a planet of its own at L4 or L5 (`generation/trojan_planet.py`): the giant's orbit, 60 degrees ahead or behind it in mean anomaly,
stored with the point and the giant's row.
"""

import math

import pytest

from planetgen.db import store
from planetgen.generation import trojan_planet
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import trojan_bodies, trojans
from planetgen.util import draw


def _system(seed=3):
    config = SystemConfig()
    config.PLANETS = True
    config.MAX_PLANETS = True
    config.BINARY_SYSTEM = False
    config.FORCE_CLASS = "J"
    with draw.bound(seed):
        return StarSystem(system_config=config), config


@pytest.fixture(scope="module")
def made():
    return _system()


def _giant(system):
    return next(planet for planet in system.planets if planet.body_type == "g")


def _force(monkeypatch, ceiling=2.0, survival=1.0):
    monkeypatch.setattr(trojan_planet, "ceiling_of", lambda system, star, giant, giants: (ceiling, 10.0 * ceiling))
    monkeypatch.setattr(trojan_planet, "survival_of", lambda star, giant, giants: survival)


def test_a_planet_at_a_point_shares_its_giants_orbit_and_stands_60_degrees_from_it_in_mean_anomaly(made, monkeypatch):
    system, _config = made
    giant = _giant(system)
    _force(monkeypatch)
    found = trojan_planet.build(system, system.star, system.planets)
    assert len(found) == 1
    planet = found[0]
    assert planet.lagrange_host is giant and planet.lagrange_point in ("L4", "L5") and planet.body_type != "a"
    assert planet.distance == giant.distance and planet.period == giant.period and planet.eccentricity == giant.eccentricity
    assert planet.orbital_inclination_deg == giant.orbital_inclination_deg and planet.orbital_ascending_node_deg == giant.orbital_ascending_node_deg
    host_orbit, own_orbit = giant.orbit_from_vector(), planet.orbit_from_vector()
    assert own_orbit["arg_periapsis"] == pytest.approx(host_orbit["arg_periapsis"], abs=1e-6)
    ahead = (trojan_bodies.mean_anomaly_rad(own_orbit["true_anomaly"], planet.eccentricity) - trojan_bodies.mean_anomaly_rad(host_orbit["true_anomaly"], giant.eccentricity)) % (2.0 * math.pi)
    assert math.degrees(ahead) == pytest.approx(60.0 if planet.lagrange_point == "L4" else 300.0, abs=1e-4)
    assert planet.mass > 0.0 and planet.planet_class is not None


def test_the_mass_is_between_the_lightest_planet_and_the_ceiling_and_the_draws_are_the_giants_own(made, monkeypatch):
    system, _config = made
    _force(monkeypatch, ceiling=1.5)
    first = trojan_planet.build(system, system.star, system.planets)
    second = trojan_planet.build(system, system.star, system.planets)
    mass_earth = first[0].mass / 5.972e24
    assert 0.1 <= mass_earth <= 1.5 + 1e-9
    assert [(p.mass, p.lagrange_point, p.position_x, p.position_y) for p in first] == [(p.mass, p.lagrange_point, p.position_x, p.position_y) for p in second]


def test_a_giant_whose_zone_grew_less_than_a_planet_or_that_cannot_keep_it_holds_none(made, monkeypatch):
    system, _config = made
    _force(monkeypatch, ceiling=0.05)
    assert trojan_planet.build(system, system.star, system.planets) == []
    _force(monkeypatch, ceiling=2.0, survival=0.0)
    assert trojan_planet.build(system, system.star, system.planets) == []


def test_the_zone_is_lost_when_the_giant_moved_from_where_it_formed():
    assert trojan_bodies.zone_kept(5.2, 5.2) and trojan_bodies.zone_kept(5.2, 4.5)
    assert not trojan_bodies.zone_kept(0.05, 5.0)


def test_the_ceiling_is_the_gate_of_the_zone_the_giant_formed_in(made):
    system, _config = made
    giant = _giant(system)
    ceiling, solids = trojan_planet.ceiling_of(system, system.star, giant, [giant])
    assert ceiling >= 0.0 and solids >= 0.0 and ceiling <= trojan_bodies.FORMATION_EFFICIENCY * solids + 1e-12
    giant.formation_distance_au = giant.distance * 3.0           # it moved in: the zone is gone
    try:
        assert trojan_planet.ceiling_of(system, system.star, giant, [giant]) == (0.0, 0.0)
    finally:
        giant.formation_distance_au = None


@pytest.mark.fresh_database
def test_a_planet_at_a_point_is_saved_with_its_giants_row_and_read_back_beside_it(mysql_config, monkeypatch):
    system, config = _system(4)
    giant = _giant(system)
    _force(monkeypatch)
    found = trojan_planet.build(system, system.star, system.planets)
    assert found
    system.planets[:] = sorted([*system.planets, *found], key=lambda body: body.distance)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, config)
        rows = conn.execute("SELECT id, orbital_index, lagrange_point, lagrange_host_id FROM planets WHERE star_system_id = ? ORDER BY orbital_index", (system_id,)).fetchall()
        giant_row = rows[system.planets.index(giant)]
        trojan_rows = [row for row in rows if row["lagrange_point"] is not None]
        assert len(trojan_rows) == len(found) and all(row["lagrange_host_id"] == giant_row["id"] for row in trojan_rows)
        assert all(row["lagrange_point"] is None and row["lagrange_host_id"] is None for row in rows if row not in trojan_rows)
        loaded = store.load_star_system(conn, system_id)
        back = [body for body in loaded.planets if getattr(body, "lagrange_point", None) and body.body_type != "a"]
        assert [body.lagrange_point for body in back] == [found[0].lagrange_point]
        assert all(body.lagrange_host is not None and body.lagrange_host.db_id == giant_row["id"] and body.lagrange_host.body_type == "g" for body in back)
    finally:
        conn.rollback()
        conn.close()


def test_the_point_is_a_lagrange_point_of_the_pair_and_the_room_holds_the_body(made, monkeypatch):
    system, _config = made
    giant = _giant(system)
    _force(monkeypatch)
    planet = trojan_planet.build(system, system.star, system.planets)[0]
    assert planet.lagrange_point in trojans.POINTS
    assert trojans.stable(system.star.mass, giant.mass, planet.mass)
