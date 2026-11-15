# tests/test_sky_size_rows.py

"""GEN.1434: how wide the sources look in a body's sky (`describe/sky_size_rows`), against the Sun, the Earth and the Moon, a close white dwarf, and a stored multiple system."""

import math

import pytest

from planetgen.db import store
from planetgen.describe import sheet, sky_size_rows
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants, orbit_elements, sky_view
from planetgen.util import draw

KM = constants.AU_TO_KM
SUN_R = 695_700.0


class _Result:
    def __init__(self, rows):
        self._rows = rows

    def fetchone(self):
        return self._rows[0] if self._rows else None

    def fetchall(self):
        return list(self._rows)


class _Conn:
    """Stored rows by the table named in a query (`FROM <table> WHERE id|star_id = ?`)."""

    def __init__(self, tables):
        self.tables = tables

    def execute(self, sql, params=()):
        table = sql.split(" FROM ")[1].split()[0]
        rows = self.tables.get(table, [])
        if "WHERE id = ?" in sql:
            rows = [row for row in rows if row["id"] == params[0]]
        return _Result(rows)


def _sun_earth_moon():
    star = {"id": 1, "radius_km": SUN_R}
    earth = {"id": 5, "star_system_id": 1, "star_id": 1, "radius_km": 6371.0, "distance_km": KM, "eccentricity": 0.0, "sky_suns": None}
    moon = {"id": 9, "parent_type": "planet", "parent_id": 5, "distance_km": 384_400.0, "sky_suns": None}
    return _Conn({"stars": [star], "planets": [earth], "moons": [moon]})


def test_the_sun_is_half_a_degree_wide_in_the_earths_sky():
    rows = {row[1]: row for row in sky_size_rows.rows(_sun_earth_moon(), "planet", 5)}
    assert rows[sky_size_rows.SUN_KEY][4] == pytest.approx(0.533, abs=0.002)
    assert rows[sky_size_rows.SUN_KEY][3] == ["0.533 degrees"]
    assert set(rows) == {sky_size_rows.SUN_KEY}                    # a lone star's planet: no other sun, no planet in the sky


def test_the_earth_is_nearly_four_times_as_wide_in_the_moons_sky_as_the_moon_in_the_earths():
    rows = {row[1]: row for row in sky_size_rows.rows(_sun_earth_moon(), "moon", 9)}
    assert rows[sky_size_rows.PLANET_KEY][4] == pytest.approx(1.90, abs=0.01)
    moon_from_earth = sky_size_rows.width_deg(1737.4, 384_400.0)
    assert rows[sky_size_rows.PLANET_KEY][4] / moon_from_earth == pytest.approx(3.67, abs=0.03)
    assert rows[sky_size_rows.SUN_KEY][4] == pytest.approx(0.533, abs=0.002)      # the Moon shares the Earth's distance from the Sun


def test_a_sun_is_wider_for_a_nearer_world_and_the_mean_light_distance_is_used():
    assert sky_size_rows.width_deg(SUN_R, 0.387 * KM) == pytest.approx(1.38, abs=0.01)         # the Sun from Mercury (1.40 at its mean-light distance)
    assert sky_size_rows.width_deg(SUN_R, 30.07 * KM) == pytest.approx(0.0177, abs=0.0005)     # and from Neptune
    eccentric = _sun_earth_moon()
    eccentric.tables["planets"][0]["eccentricity"] = 0.5
    seen = sky_size_rows.rows(eccentric, "planet", 5)[0][4]
    expected = sky_size_rows.width_deg(SUN_R, orbit_elements.flux_distance_au(1.0, 0.5) * KM)
    assert seen == pytest.approx(expected, rel=1e-12) and seen > 0.533      # the circle of the same mean light is nearer than the axis


def test_a_close_white_dwarf_is_a_small_hard_disk():
    degrees = sky_size_rows.width_deg(6000.0, 0.01 * KM)
    assert degrees == pytest.approx(0.46, abs=0.01)
    assert degrees == pytest.approx(math.degrees(2.0 * math.asin(6000.0 / (0.01 * KM))), rel=1e-12)


def test_what_is_missing_gives_no_row():
    assert sky_size_rows.width_deg(None, KM) == 0.0 and sky_size_rows.width_deg(SUN_R, 0.0) == 0.0
    nothing = _sun_earth_moon()
    nothing.tables["planets"][0]["distance_km"] = 0.0
    assert sky_size_rows.rows(nothing, "planet", 5) == []
    assert sky_size_rows.rows(nothing, "star", 1) == [] and sky_size_rows.rows(nothing, "belt", 1) == []
    assert sky_size_rows.rows(_sun_earth_moon(), "planet", 404) == []


def test_the_derived_keys_are_told_from_the_stored_ones():
    assert sky_size_rows.is_derived_key("sky_width_sun") and sky_size_rows.is_derived_key("sky_width_planet") and sky_size_rows.is_derived_key("sky_width_other_2")
    assert not sky_size_rows.is_derived_key("sky_suns") and not sky_size_rows.is_derived_key("clock_lag_s_per_yr")


def _binary_with_other_suns():
    """A seeded multiple system with a planet whose sky holds another source that shows, and the planet's index."""
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = True
    cfg.MOONS = False
    for seed in range(200):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if len(system.stars) < 2:
            continue
        for planet in system.planets:
            if getattr(planet, "body_type", None) == "a" or getattr(planet, "star", None) is None:
                continue
            sky = sky_view.bodies(planet.star, orbit_elements.flux_distance_au(planet.distance, planet.eccentricity or 0.0))
            if len(sky) > 1:
                return system, cfg, planet.name
    pytest.fail("no multiple system with another source in a planet's sky in 200 seeds")


@pytest.mark.fresh_database
def test_a_stored_multiple_system_lists_the_other_suns_with_the_widths_of_the_sky_text(mysql_config):
    system, cfg, name = _binary_with_other_suns()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        row = conn.execute("SELECT id FROM planets WHERE star_system_id = ? AND name = ?", (system_id, name)).fetchone()
        assert row is not None
        planet = next(p for p in system.planets if getattr(p, "name", None) == name)
        sky = sky_view.bodies(planet.star, orbit_elements.flux_distance_au(planet.distance, planet.eccentricity or 0.0))
        found = {item[1]: item for item in sky_size_rows.rows(conn, "planet", row["id"])}
        assert sky_size_rows.SUN_KEY in found
        others = [found[key] for key in found if key.startswith(sky_size_rows.OTHER_KEY_PREFIX)]
        assert 1 <= len(others) <= min(len(sky) - 1, sky_size_rows.MAX_OTHER_ROWS)
        assert others[0][4] == pytest.approx(sky[1].angular_deg, rel=1e-3)           # the brightest other source first, the width the text gives it
        sheet_keys = [item["key"] for group in sheet.sheet(conn, "planet", row["id"])["groups"] for item in group["rows"]]
        assert sky_size_rows.SUN_KEY in sheet_keys and any(key.startswith(sky_size_rows.OTHER_KEY_PREFIX) for key in sheet_keys)
    finally:
        conn.rollback()
        conn.close()
