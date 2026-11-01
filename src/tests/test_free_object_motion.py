# tests/test_free_object_motion.py

"""
GEN.336: a free object (black hole, neutron star, rogue planet, interstellar comet, asteroid field, nebula, supernova
remnant) is stored with a movement vector and a mass instead of a scalar speed. The vector is the tangent of the galaxy's
rotation at its place, at the speed it was generated with; the orbit update turns it with the position; a star reads its
system's.
"""

import math

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.generation.system import StarSystem
from planetgen.physics import bulk_mass

PLACE = {"center_x_pc": 3000.0, "center_y_pc": 4000.0, "center_z_pc": 10.0, "galactic_radius_pc": 5000.0}

KINDS = (
    ("black_holes", lambda: BlackHole(SystemConfig(), mass_class="stellar"), store.insert_black_hole),
    ("neutron_stars", lambda: NeutronStar(SystemConfig()), store.insert_neutron_star),
    ("rogue_planets", lambda: RoguePlanet(SystemConfig()), store.insert_rogue_planet),
    ("interstellar_comets", lambda: InterstellarComet(SystemConfig()), store.insert_interstellar_comet),
    ("asteroid_fields", lambda: AsteroidField(SystemConfig()), store.insert_asteroid_field),
    ("nebulae", lambda: Nebula(SystemConfig()), store.insert_nebula),
    ("supernova_remnants", lambda: SupernovaRemnant(SystemConfig()), store.insert_supernova_remnant),
)

pytestmark = pytest.mark.fresh_database


def _vector(row):
    return tuple(row[f"velocity_{axis}_kms"] for axis in "xyz")


@pytest.mark.parametrize("table,make,insert", KINDS, ids=[kind[0] for kind in KINDS])
def test_a_placed_object_stores_the_tangent_of_the_rotation_at_its_speed(mysql_config, table, make, insert):
    obj = make()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            object_id = insert(conn, obj, placement=dict(PLACE))
        row = conn.execute(f"SELECT * FROM {table} WHERE id = ?", (object_id,)).fetchone()
    finally:
        conn.close()
    vx, vy, vz = _vector(row)
    assert math.sqrt(vx * vx + vy * vy + vz * vz) == pytest.approx(obj.galactic_orbital_speed_kms, rel=1e-9)
    assert vx * PLACE["center_x_pc"] + vy * PLACE["center_y_pc"] == pytest.approx(0.0, abs=1e-6)  # along the tangent
    assert vz == 0.0
    assert (PLACE["center_x_pc"] * vy - PLACE["center_y_pc"] * vx) > 0  # counterclockwise from galactic north
    assert "galactic_orbital_speed_kms" not in row


@pytest.mark.parametrize("table,make,insert", KINDS, ids=[kind[0] for kind in KINDS])
def test_an_object_with_no_place_keeps_its_speed_along_the_tangent_at_the_default_place(mysql_config, table, make, insert):
    obj = make()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            object_id = insert(conn, obj)
        row = conn.execute(f"SELECT * FROM {table} WHERE id = ?", (object_id,)).fetchone()
    finally:
        conn.close()
    assert _vector(row) == pytest.approx((0.0, obj.galactic_orbital_speed_kms, 0.0))
    assert store._speed_kms(row) == pytest.approx(obj.galactic_orbital_speed_kms)


def test_the_bulk_mass_of_a_cloud_or_a_comet_is_derived_and_a_field_has_a_stated_zero(mysql_config):
    nebula, remnant = Nebula(SystemConfig()), SupernovaRemnant(SystemConfig())
    comet, field = InterstellarComet(SystemConfig()), AsteroidField(SystemConfig())
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            ids = (store.insert_nebula(conn, nebula, placement=dict(PLACE)),
                   store.insert_supernova_remnant(conn, remnant, placement=dict(PLACE)),
                   store.insert_interstellar_comet(conn, comet, placement=dict(PLACE)),
                   store.insert_asteroid_field(conn, field, placement=dict(PLACE)))
        masses = [conn.execute(f"SELECT mass_kg FROM {table} WHERE id = ?", (object_id,)).fetchone()["mass_kg"]
                  for table, object_id in zip(("nebulae", "supernova_remnants", "interstellar_comets", "asteroid_fields"),
                                              ids)]
    finally:
        conn.close()
    assert masses[0] == pytest.approx(bulk_mass.gas_cloud_mass_kg(nebula.density_cm3, nebula.radius_ly))
    assert masses[1] == pytest.approx(bulk_mass.gas_cloud_mass_kg(remnant.density_cm3, remnant.radius_ly))
    assert masses[2] == pytest.approx(comet.makeup.mass_kg(comet.nucleus_diameter_km))
    assert masses[3] == 0.0


def test_the_orbit_update_turns_the_vector_with_the_position(mysql_config):
    nebula = Nebula(SystemConfig())
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            nebula_id = store.insert_nebula(conn, nebula, placement=dict(PLACE))
            conn.execute("UPDATE nebulae SET galactic_orbital_period_gy = 0.25")
        before = conn.execute("SELECT * FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
        with conn:
            store.advance_galactic_positions(conn, store.orbit_clock(conn, 5.0e7))  # a fifth of the orbit
        after = conn.execute("SELECT * FROM nebulae WHERE id = ?", (nebula_id,)).fetchone()
    finally:
        conn.close()
    assert (after["center_x_pc"], after["center_y_pc"]) != pytest.approx((before["center_x_pc"], before["center_y_pc"]))
    # The stepper (GEN.340) carries the vector with the position through the galaxy's potential: the speed follows the energy along the path, so it stays
    # within a few percent of the start over a fifth of an orbit (it was a rotation, exactly the same length, before), and the path is still the tangent's way round.
    assert store._speed_kms(after) == pytest.approx(store._speed_kms(before), rel=0.1)
    vx, vy, _vz = _vector(after)
    radial = (vx * after["center_x_pc"] + vy * after["center_y_pc"]) / (store._speed_kms(after) * math.hypot(after["center_x_pc"], after["center_y_pc"]))
    assert abs(radial) < 0.5                                                      # mostly along the tangent
    assert after["center_x_pc"] * vy - after["center_y_pc"] * vx > 0              # counterclockwise


def test_a_system_and_its_stars_carry_the_speed_of_the_systems_vector(mysql_config):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = False
    cfg.BINARY_SYSTEM = True
    system = StarSystem(system_config=cfg)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        row = conn.execute("SELECT * FROM star_systems WHERE id = ?", (system_id,)).fetchone()
        assert "binary_galactic_orbital_speed_kms" not in row
        reloaded = store.load_star_system(conn, system_id)
    finally:
        conn.close()
    assert reloaded.star.galactic_orbital_speed_kms == pytest.approx(store._speed_kms(row))
    assert store._speed_kms(row) == pytest.approx(system.star.galactic_orbital_speed_kms)
