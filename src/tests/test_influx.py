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
JUPITER = influx.Giant(317.8, 11.2, 1.8)


def flux(planet_au, reservoirs=(BELT,), giants=()):
    return influx.planet_influx_kg_per_yr(EARTH_KM, planet_au, 1.0, list(reservoirs), list(giants))[0]


def test_a_belt_sends_its_collisional_loss_to_the_planets_inside_it():
    assert influx.supply_kg_per_yr(BELT) == pytest.approx(influx.BELT_DUST_SHARE * 1.0e21 / 4.6e9)
    assert flux(1.0) > 0.0
    assert flux(10.0) == 0.0, "nothing is sent outward"


def test_the_influx_grows_with_the_planet_and_the_belt_and_is_a_share_of_the_belts_dust():
    big = influx.Reservoir("belt", 2.0e21, 2.2, 3.2, 4.6e9)
    assert flux(1.0, [big]) == pytest.approx(2.0 * flux(1.0))
    assert influx.planet_influx_kg_per_yr(2.0 * EARTH_KM, 1.0, 1.0, [BELT], [])[0] > 3.5 * flux(1.0), "the radius squared, less the share that saturates"
    assert flux(1.0) < influx.supply_kg_per_yr(BELT), "a planet never takes more than the belt sends"
    assert flux(0.5) > flux(1.0), "the nearer the star the larger the share of the passes"


def test_a_grain_crosses_a_planets_orbit_tens_of_thousands_of_times_on_its_way_in():
    assert influx.crossings(1.0, 1.0) == pytest.approx(3.2e4, rel=0.05)
    assert influx.crossings(4.0, 1.0) == pytest.approx(2.0 * influx.crossings(1.0, 1.0)), "N goes as the square root of the distance"
    assert influx.crossings(1.0, 1.0, luminosity_solar=2.0) == pytest.approx(0.5 * influx.crossings(1.0, 1.0)), "a brighter star blows the grains faster"
    assert influx.grain_beta(1.0) == influx.GRAIN_BETA
    assert influx.crossings(0.0, 1.0) == 0.0


def test_the_earth_takes_from_a_solar_system_belt_within_a_factor_of_twenty_of_the_measured_influx():
    """The asteroid belt alone (2.4e21 kg, 4.6 Gyr) is not all of the Sun's dust (the Jupiter-family comets are the larger share), the measured influx is 1e7 kg a year: it was 45 kg a year before the crossings."""
    asteroid_belt = influx.Reservoir("belt", 2.4e21, 2.2, 3.2, 4.6e9)
    earth = influx.planet_influx_kg_per_yr(EARTH_KM, 1.0, 1.0, [asteroid_belt], [], 1.0)[0]
    assert 1.0e7 / 20.0 < earth < 1.0e7 * 20.0


def test_a_comet_nucleus_is_not_multiplied_by_the_crossings():
    comet = influx.Reservoir("comet", 1.0e14, 0.5, 0.5)
    expected = influx.supply_kg_per_yr(comet) * influx.cross_section_share(EARTH_KM, 1.0)
    assert flux(1.0, [comet]) == pytest.approx(expected)


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
    assert influx.planet_influx_kg_per_yr(EARTH_KM, 1.0, 1.0, [], [JUPITER])[0] == 0.0


def _saved_system(mysql_config):
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    config.PLANETS = True
    config.MOONS = False
    config.BINARY_SYSTEM = False
    for seed in range(1, 60):
        draw.set_run_seed(seed)
        system = StarSystem(system_config=config)
        if any(b.body_type == 'a' for b in system.planets) and any(b.body_type not in ('a', 'g') for b in system.planets):
            return store.save_system(system, config, config=mysql_config)
    pytest.fail("no system with a belt and a small planet in 59 seeds")


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
        if any(b.body_type == 'g' for b in system.planets):
            saved = store.save_system(system, config, config=mysql_config)
            with store.get_connection(mysql_config) as conn:
                giant = conn.execute("SELECT id FROM planets WHERE star_system_id = ? AND body_type = 'g'", (saved,)).fetchone()["id"]
                assert influx_row.influx_kg_per_yr(conn, giant) is None
            return
    pytest.fail("no giant in 59 seeds")
