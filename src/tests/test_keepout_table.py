# tests/test_keepout_table.py

"""
NAV.54: the one keep-out table (docs/design/orbits/course-avoidance.md section 2): an asteroid field keeps out its stored radius, a
supermassive hole its sphere of influence, a quasar the larger of that and its light, a moon its Hill radius against its planet, a star
inside its own system its radiation radius, and a nebula, remnant or belt is passed through with a note.
"""

import math

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import keepout
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.compact_remnant import BlackHole
from planetgen.generation.phenomena.quasar import Quasar
from planetgen.physics import constants
from planetgen.physics.orbits import calculate_hill_sphere
from tests.publicids import pid
from tests.test_api import _save_wide_binary_with_moons, client  # noqa: F401

KM_PER_PC = constants.PARSEC_M / 1000.0
KM_PER_LY = 9.4607304725808e12
AU_KM = 1.495978707e8


def _keep_out(client, kind, object_id):
    return client.get(f"/api/objects/{kind}:{pid(kind, object_id)}").get_json()


def test_an_asteroid_field_keeps_out_its_stored_radius():
    found = keepout.asteroid_field_keep_out(0.5)
    assert found.basis == "radius" and found.note is None
    assert found.radius_km == pytest.approx(0.5 * KM_PER_LY, rel=1e-6)


def test_a_supermassive_hole_keeps_out_its_sphere_of_influence_not_its_horizon():
    found = keepout.supermassive_keep_out(4.0e6)
    horizon_km = 2 * constants.G * 4.0e6 * constants.SOLAR_MASS_TO_KG / constants.SPEED_OF_LIGHT_M_S ** 2 / 1000.0
    assert found.basis == "influence" and found.radius_km > 1e4 * horizon_km
    # Parsecs, not astronomical units: G M / sigma^2 for the bulge dispersion its mass sets (the M-sigma relation).
    assert 1.0 < found.radius_km / KM_PER_PC < 10.0
    # A heavier hole reaches farther.
    assert keepout.sphere_of_influence_km(1e8) > keepout.sphere_of_influence_km(4e6)


def test_a_quasar_keeps_out_the_larger_of_its_influence_and_its_light():
    influence = keepout.sphere_of_influence_km(1.0e8)
    dim = keepout.quasar_keep_out(1.0e8, 1.0e30)
    assert (dim.basis, dim.radius_km) == ("influence", pytest.approx(influence))
    bright = keepout.quasar_keep_out(1.0e8, 1.0e41)
    light = math.sqrt(1.0e41 / (4 * math.pi * tuning.KEEP_OUT_QUASAR_FLUX_W_M2)) / 1000.0
    assert bright.basis == "radiation" and bright.radius_km == pytest.approx(light) and bright.radius_km > influence


def test_a_star_inside_its_own_system_keeps_out_its_radiation_radius():
    # The Sun: 10 R_sun is 0.047 AU, 50 kW/m^2 is 0.17 AU, so the light sets it.
    sun = keepout.star_in_system_keep_out(6.957e5, 3.828e26)
    assert sun.basis == "radiation" and sun.radius_km / AU_KM == pytest.approx(0.17, abs=0.01)
    # A dim star's floor is ten of its radii.
    dwarf = keepout.star_in_system_keep_out(1.0e5, 1.0e22)
    assert dwarf.basis == "radius" and dwarf.radius_km == pytest.approx(10 * 1.0e5)
    # An O3 star (1e6 solar luminosities) keeps out about 165 AU.
    o3 = keepout.star_in_system_keep_out(1.2e7, 1.0e6 * 3.828e26)
    assert o3.radius_km / AU_KM == pytest.approx(165, rel=0.05)


def test_the_notes_of_a_cloud_and_a_belt_say_why_a_course_passes_through():
    assert "thin" in keepout.PASS_THROUGH_NOTE and "No mass" not in keepout.PASS_THROUGH_NOTE
    assert "thin" in keepout.BELT_NOTE


def test_the_rows_of_the_table_are_read_from_the_stored_objects(client, mysql_config):
    placement = {"center_x_pc": 5.0, "center_y_pc": 1.0, "center_z_pc": 2.0, "galactic_radius_pc": 5.5}
    center = {"center_x_pc": 0.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 0.0}
    conn = store.get_connection(mysql_config)
    try:
        field = AsteroidField(SystemConfig())
        field_id = store.insert_asteroid_field(conn, field, placement=placement)
        hole = BlackHole(SystemConfig(), mass_class="supermassive")
        hole_id = store.insert_black_hole(conn, hole, placement=center)
        quasar = Quasar(SystemConfig())
        quasar_id = store.insert_quasar(conn, quasar, placement=center)
        conn.commit()
    finally:
        conn.close()
    got = _keep_out(client, "asteroid_field", field_id)["keep_out"]
    assert got["basis"] == "radius" and got["radius_km"] == pytest.approx(field.radius_ly * KM_PER_LY, rel=1e-6)
    got = _keep_out(client, "black_hole", hole_id)["keep_out"]
    assert got["basis"] == "influence" and got["radius_km"] == pytest.approx(keepout.sphere_of_influence_km(hole.mass_solar))
    got = _keep_out(client, "quasar", quasar_id)["keep_out"]
    assert got["basis"] in ("influence", "radiation") and got["radius_km"] >= keepout.sphere_of_influence_km(quasar.black_hole_mass_solar)


def test_a_star_says_its_radius_inside_its_own_system_and_a_moon_its_hill_radius_against_its_planet(client, mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        star = conn.execute("SELECT id, radius_km, luminosity_w FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchone()
        moon = conn.execute(
            "SELECT m.id, m.distance_km, m.mass_kg, p.mass_kg AS planet_mass_kg, m.hill_radius_km FROM moons m "
            "JOIN planets p ON p.id = m.planet_id WHERE m.star_system_id = ? ORDER BY m.id", (system_id,)).fetchone()
    finally:
        conn.close()
    found = _keep_out(client, "star", star["id"])
    inside = found["keep_out_in_system"]
    expected = keepout.star_in_system_keep_out(star["radius_km"], star["luminosity_w"])
    assert (inside["basis"], inside["radius_km"]) == (expected.basis, pytest.approx(expected.radius_km))
    assert found["keep_out"]["basis"] == "perimeter"      # between systems it is still the galactic Hill radius
    got = _keep_out(client, "moon", moon["id"])["keep_out"]
    against_planet = calculate_hill_sphere(moon["distance_km"] * 1000.0, moon["mass_kg"], moon["planet_mass_kg"]) / 1000.0
    assert got["basis"] == "hill" and got["radius_km"] == pytest.approx(moon["hill_radius_km"])
    assert got["radius_km"] == pytest.approx(against_planet, rel=1e-6)
