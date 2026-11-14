# tests/test_rogue_body_storage.py

"""GEN.1111 part 6: a planet built with no star is stored in `rogue_planets` through the same field list a bound planet's `planets` row is, with its
child rows and its moons, and read back."""

import pytest

from planetgen.db import store
from planetgen.generation import ejection_events, rogue_landing
from planetgen.generation.phenomena import starless_host
from planetgen.galaxy.sector import SpaceSector
from planetgen.physics import constants
from planetgen.util import draw
from tests.test_gen_builder_internals import make_system

pytestmark = pytest.mark.fresh_database

EARTH_KG = constants.EARTH_MASS_TO_KG
DONOR_UID = bytes(range(1, 13))


def _event():
    return ejection_events.EjectionEvent(DONOR_UID.hex(), 0, 3.0 * EARTH_KG, 2.0 * EARTH_KG, 5.0, 5.0, 3.0, 3000.0 * EARTH_KG, True, 1.0, 2.0, 1.5,
                                         (8.0, -3.0, 2.0), (4.0, 5.0, -1.0), 2.0e6, "thrown")


def _rogue(seed, mass_earth, moons):
    with draw.bound(f"rogue-storage-donor-{seed}"):    # the donor star is drawn from a stream of its own, so the moon below does not depend on the run seed
        system = make_system()
    donor = system.star
    with draw.bound(f"rogue-storage-{seed}"):
        planet = starless_host.build_ejected_planet(system.system_config, donor, 5.0, mass_earth * EARTH_KG, moon_count=moons)
    landing = rogue_landing.landing((8000.0, 100.0, 20.0), _event(), 10.0, 220.0)
    rogue_landing.apply_motion(planet, landing)
    return planet, landing


@pytest.fixture
def world(mysql_config, redis_server):
    sector_id = store.save_sector(SpaceSector("Rogue Sector", edge_ly=10.0), config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        yield conn, sector_id
    finally:
        conn.close()


def _insert(conn, sector_id, planet, landing, index=0):
    with conn:
        row_id = store.insert_rogue_body(
            conn, planet, sector_id=sector_id, placement=rogue_landing.placement_of(landing), donor_uid=DONOR_UID, event_index=index,
            ejection_delay_gy=1.0, since_ejection_gy=1.5, ejection_speed_kms=9.0, rogue_note="Thrown out of a system 1.5 Gy ago.")
        store.assign_uids(conn, phenomenon=("rogue_planets", row_id))
    return row_id


def test_a_rogue_planet_is_stored_with_every_planet_column_and_its_ejection(world):
    conn, sector_id = world
    planet, landing = _rogue(0, 3.0, 1)
    row_id = _insert(conn, sector_id, planet, landing)
    row = conn.execute("SELECT * FROM rogue_planets WHERE id = ?", (row_id,)).fetchone()
    assert row["donor_uid"] == DONOR_UID and row["event_index"] == 0 and row["rogue_note"].startswith("Thrown out")
    assert row["mass_kg"] == pytest.approx(planet.mass) and row["body_type"] == planet.body_type
    assert row["center_x_pc"] == pytest.approx(landing.position_pc[0]) and row["galactic_radius_pc"] > 7000.0
    assert row["velocity_x_kms"] == pytest.approx(landing.velocity_km_s[0])
    assert row["galactic_orbital_period_gy"] > 0.0 and row["planet_type"] == planet.planet_type
    assert row["surface_temperature_k"] == pytest.approx(planet.surface_temperature)
    moons = conn.execute("SELECT parent_type, parent_id, rogue_planet_id FROM moons WHERE rogue_planet_id = ?", (row_id,)).fetchall()
    assert len(moons) == len(planet.moons) == 1 and moons[0]["parent_type"] == "rogue_planet" and moons[0]["parent_id"] == row_id


def test_the_child_rows_of_a_rogue_planet_are_the_ones_a_planet_keeps(world):
    conn, sector_id = world
    planet, landing = _rogue(1, 80.0, 0)
    row_id = _insert(conn, sector_id, planet, landing)
    layers = conn.execute("SELECT COUNT(*) AS n FROM rogue_planet_interior_layers WHERE rogue_planet_id = ?", (row_id,)).fetchone()["n"]
    assert layers == len(planet.interior_layers or [])
    gases = conn.execute("SELECT gas, partial_pa FROM rogue_planet_gases WHERE rogue_planet_id = ?", (row_id,)).fetchall()
    assert {g["gas"]: g["partial_pa"] for g in gases} == pytest.approx({k: v for k, v in store.atmosphere.gas_partials_pa(planet).items() if v > 0.0})
    spectrum = conn.execute("SELECT COUNT(*) AS n FROM rogue_planet_reflection_spectrum WHERE rogue_planet_id = ?", (row_id,)).fetchone()["n"]
    assert spectrum == len(planet.reflection_spectrum_visible or []) + len(planet.reflection_spectrum_non_visible or [])


def test_a_stored_rogue_planet_reads_back_as_the_planet_it_was(world):
    conn, sector_id = world
    planet, landing = _rogue(2, 3.0, 2)
    row_id = _insert(conn, sector_id, planet, landing)
    row = conn.execute("SELECT * FROM rogue_planets WHERE id = ?", (row_id,)).fetchone()
    back = store.load_rogue_body(conn, row, planet.system_config)
    assert back.star is None and back.db_id == row_id
    for name in ("mass", "radius", "surface_temperature", "planet_class", "body_type", "gravity", "atmosphere", "density", "internal_heat_flux_w_m2",
                 "rotation_period_hours", "bond_albedo", "hydrosphere", "life_class"):
        assert getattr(back, name) == pytest.approx(getattr(planet, name)) if isinstance(getattr(planet, name), float) else getattr(back, name) == getattr(planet, name), name
    assert store.atmosphere.gas_partials_pa(back) == pytest.approx(store.atmosphere.gas_partials_pa(planet))
    assert len(back.moons) == 2 and [m.mass for m in back.moons] == pytest.approx([m.mass for m in planet.moons])


def test_saving_a_donor_again_replaces_its_rogues_and_never_duplicates_them(world):
    conn, sector_id = world
    planet, landing = _rogue(3, 3.0, 1)
    _insert(conn, sector_id, planet, landing, index=0)
    with pytest.raises(Exception):
        _insert(conn, sector_id, planet, landing, index=0)       # the same donor and event: refused
    conn.rollback()
    _insert(conn, sector_id, planet, landing, index=1)
    assert conn.execute("SELECT COUNT(*) AS n FROM rogue_planets WHERE donor_uid = ?", (DONOR_UID,)).fetchone()["n"] == 2
    with conn:
        assert store.delete_rogues_of_donor(conn, DONOR_UID) == 2
    assert conn.execute("SELECT COUNT(*) AS n FROM rogue_planets WHERE donor_uid = ?", (DONOR_UID,)).fetchone()["n"] == 0
    assert conn.execute("SELECT COUNT(*) AS n FROM moons WHERE rogue_planet_id IS NOT NULL").fetchone()["n"] == 0
    assert conn.execute("SELECT COUNT(*) AS n FROM rogue_planet_gases").fetchone()["n"] == 0
