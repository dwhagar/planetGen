# tests/test_influx.py

"""GEN.1440, GEN.1453 slice A: the belt dust and comet influx of a planet (`physics/influx.py`) and its derived row on the sheet (`describe/influx_row.py`)."""

import math

import pytest

from planetgen.db import store
from planetgen.describe import influx_row, sheet
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import influx
from planetgen.util import draw

BELT = influx.Reservoir("belt", 1.0e21, 2.2, 3.2, 4.6e9)
EARTH_KM = 6371.0
EARTH_KG = 5.972e24
MAIN_BELT = influx.Reservoir("belt", 2.4e21, 2.2, 3.2, 4.6e9)
JUPITER = influx.Giant(317.8, 11.2, 1.8)


def flux(planet_au, reservoirs=(BELT,), giants=()):
    return influx.planet_influx_kg_per_yr(EARTH_KM, EARTH_KG, planet_au, 1.0, 1.0, list(reservoirs), list(giants))[0]


def test_a_belt_sends_its_collisional_loss_to_the_planets_inside_it():
    assert influx.supply_kg_per_yr(BELT) == pytest.approx(influx.BELT_DUST_SHARE * 1.0e21 / 4.6e9)
    assert flux(1.0) > 0.0
    assert flux(10.0) == 0.0, "nothing is sent outward"


def test_the_influx_grows_with_the_belt_and_the_planets_size_and_falls_with_distance():
    big = influx.Reservoir("belt", 2.0e21, 2.2, 3.2, 4.6e9)
    assert flux(1.0, [big]) == pytest.approx(2.0 * flux(1.0))
    assert influx.planet_influx_kg_per_yr(2.0 * EARTH_KM, 8.0 * EARTH_KG, 1.0, 1.0, 1.0, [BELT], [])[0] > 4.0 * flux(1.0), "a bigger planet takes more than its cross-section"
    assert flux(0.5) > flux(1.0) > flux(1.5) > 0.0, "the stream thins outward"


def test_a_sun_with_the_main_belt_gives_the_earth_about_the_measured_influx():
    # GEN.1543: the isotropic sphere gave 45 kg a year against the 1e7 measured (28 tonnes a day, [R1397] recalled). A Sun and a main-belt sized reservoir are the test input, not a content of the generator.
    earth = influx.planet_influx_kg_per_yr(EARTH_KM, EARTH_KG, 1.0, 1.0, 1.0, [MAIN_BELT], [])[0]
    assert 3.0e6 < earth < 3.0e7
    assert influx.planet_influx_kg_per_yr(EARTH_KM, EARTH_KG, 1.0, 1.0, 1.0, [MAIN_BELT], [influx.Giant(317.8, 11.2, 5.2)])[0] == earth, "a Jupiter beyond the belt shields nothing of it"


def test_the_grains_drift_in_at_the_poynting_robertson_speed_the_star_gives_them():
    beta = influx.radiation_ratio(1.0, 1.0)
    assert 0.005 < beta < 0.012    # 3 L / (16 pi G M c rho s) for a 30 micrometer, 2500 kg per cubic meter grain
    assert influx.drift_speed_m_s(1.0, 1.0, 1.0) == pytest.approx(2.0 * beta * 1.32712e20 / (2.99792458e8 * 1.495978707e11), rel=1e-3)
    assert influx.drift_speed_m_s(2.0, 1.0, 1.0) == pytest.approx(influx.drift_speed_m_s(1.0, 1.0, 1.0) / 2.0), "the drift falls as 1 / a"
    assert influx.radiation_ratio(1.0, 4.0) == pytest.approx(4.0 * beta), "a brighter star pushes harder"


def test_a_planet_with_a_stronger_pull_takes_a_larger_share_by_focusing():
    light = influx.capture_share(EARTH_KM, 0.1 * EARTH_KG, 1.0, 1.0, 1.0)
    heavy = influx.capture_share(EARTH_KM, 10.0 * EARTH_KG, 1.0, 1.0, 1.0)
    assert 0.0 < light < heavy
    assert influx.capture_share(EARTH_KM, EARTH_KG, 1.0, 1.0, 0.0) == 0.0, "a star with no light drives no drift to read"
    assert influx.capture_share(EARTH_KM, EARTH_KG, 1.0, 1.0, 1.0) <= 1.0


def test_a_planet_inside_the_ring_gets_the_part_of_the_ring_outside_its_orbit():
    assert influx.reach(BELT, 2.7) == pytest.approx(0.5)
    assert influx.reach(BELT, 2.2) == 1.0 and influx.reach(BELT, 3.2) == 0.0


def test_a_giant_between_the_belt_and_the_planet_shields_it():
    plain = flux(1.0)
    shielded = flux(1.0, giants=[JUPITER])
    theta = influx.scatter.safronov_number(317.8, 11.2, 1.8, 1.0)
    assert shielded == pytest.approx(plain * math.exp(-theta / influx.scatter.THETA_EJECT))
    assert 0.0 < shielded < plain
    assert flux(1.0, giants=[influx.Giant(317.8, 11.2, 0.5)]) == plain, "a giant inside the planet's orbit is not in the way"
    assert flux(1.0, giants=[influx.Giant(317.8, 11.2, 5.2)]) == plain, "nor one beyond the belt (a Jupiter outside an asteroid belt shields the inner planets from nothing of it)"


def test_a_comet_reaches_a_planet_its_perihelion_passes_inside():
    comet = influx.Reservoir("comet", 1.0e14, 0.5, 0.5)
    assert flux(1.0, [comet]) > 0.0 and flux(0.3, [comet]) == 0.0
    assert influx.supply_kg_per_yr(comet) == pytest.approx(1.0e14 / influx.COMET_LIFETIME_YR)


def test_a_system_with_nothing_to_send_gives_nothing():
    assert influx.planet_influx_kg_per_yr(EARTH_KM, EARTH_KG, 1.0, 1.0, 1.0, [], [JUPITER])[0] == 0.0


def _saved_system(mysql_config):
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    config.PLANETS = True
    config.MOONS = False
    config.BINARY_SYSTEM = False
    for seed in range(1, 60):
        draw.set_run_seed(seed)
        system = StarSystem(system_config=config)
        belts = [b for b in system.primary.orbiters if b.body_type == 'a' and b.mass_kg]
        # a small planet the belt's dust reaches: one inside its outer edge (a planet beyond it gets nothing from it)
        if any(b.body_type not in ('a', 'g') and b.distance < belt.upper_limit for b in system.primary.orbiters for belt in belts):
            return store.save_system(system, config, config=mysql_config)
    pytest.fail("no system with a belt that reaches a small planet in 59 seeds")


def test_a_small_planet_of_a_stored_system_with_a_belt_has_the_row_on_its_sheet(mysql_config):
    saved = _saved_system(mysql_config)
    with store.get_connection(mysql_config) as conn:
        rows = conn.execute("SELECT id, distance_km FROM planets WHERE star_system_id = ? AND body_type NOT IN ('g') ORDER BY distance_km", (saved,)).fetchall()
        values = {row["id"]: influx_row.influx_kg_per_yr(conn, row["id"]) for row in rows}
        assert any(value for value in values.values()), values
        best = next(row["id"] for row in rows if values[row["id"]])
        page = sheet.sheet(conn, "planet", best)
    shown = [row for group in page["groups"] for row in group["rows"] if row["key"] == influx_row.INFLUX_KEY]
    assert len(shown) == 1 and shown[0]["raw"] == pytest.approx(values[best]) and shown[0]["raw"] > 0.0


def test_a_giant_has_no_influx_row(mysql_config):
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    config.MOONS = False
    config.BINARY_SYSTEM = False
    for seed in range(1, 60):
        draw.set_run_seed(seed)
        system = StarSystem(system_config=config)
        if any(b.body_type == 'g' for b in system.primary.orbiters):
            saved = store.save_system(system, config, config=mysql_config)
            with store.get_connection(mysql_config) as conn:
                giant = conn.execute("SELECT id FROM planets WHERE star_system_id = ? AND body_type = 'g'", (saved,)).fetchone()["id"]
                assert influx_row.influx_kg_per_yr(conn, giant) is None
            return
    pytest.fail("no giant in 59 seeds")
