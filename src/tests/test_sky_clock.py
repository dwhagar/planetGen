# tests/test_sky_clock.py

"""GEN.1434: the clock, the calendar of eclipses and the gravitational lens of a body (`physics/sky_clock`, `describe/clock_rows`), against the Sun, the Earth, the Moon and Venus."""

import math

import pytest

from planetgen.describe import clock_rows
from planetgen.physics import constants, eclipses, sky_clock

KM = constants.AU_TO_KM
SUN_M, SUN_R = constants.SOLAR_MASS_TO_KG, 695_700.0
EARTH_M, EARTH_R = constants.EARTH_MASS_TO_KG, 6371.0
YEAR = constants.SECONDS_PER_YEAR


def test_a_clock_on_the_earth_runs_half_a_second_a_year_behind_a_distant_one():
    earth = sky_clock.lag_per_year_s(sky_clock.clock_lag_fraction(EARTH_M, EARTH_R, [(SUN_M, KM)]))
    assert earth == pytest.approx(0.489, rel=0.02)                       # 1.5 G M / (a c^2) of the Sun at 1 AU
    own = sky_clock.lag_per_year_s(sky_clock.clock_lag_fraction(EARTH_M, EARTH_R))
    assert own == pytest.approx(0.0220, rel=0.02)                        # the Earth's own surface: 6.96e-10


def test_a_moon_adds_its_planet_and_the_planets_star_to_its_own_surface():
    moon = sky_clock.clock_lag_fraction(7.342e22, 1737.4, [(EARTH_M, 384_400.0), (SUN_M, KM)])
    parts = sky_clock.clock_lag_fraction(7.342e22, 1737.4) + sky_clock.clock_lag_fraction(0.0, 0.0, [(EARTH_M, 384_400.0)]) + sky_clock.clock_lag_fraction(0.0, 0.0, [(SUN_M, KM)])
    assert moon == pytest.approx(parts, rel=1e-12)
    assert moon > sky_clock.clock_lag_fraction(0.0, 0.0, [(SUN_M, KM)])                  # the Moon's surface and the Earth's field at its distance add to the Sun's term
    assert moon < sky_clock.clock_lag_fraction(EARTH_M, EARTH_R, [(SUN_M, KM)])         # but the Moon's small surface potential is far below the Earth's


def test_a_close_white_dwarf_slows_a_surface_clock_by_thousands_of_seconds_a_year():
    sun = sky_clock.lag_per_year_s(sky_clock.clock_lag_fraction(SUN_M, SUN_R))
    dwarf = sky_clock.lag_per_year_s(sky_clock.clock_lag_fraction(0.6 * SUN_M, 6000.0))
    assert sun == pytest.approx(67.0, rel=0.02)
    assert dwarf > 3000.0 and dwarf / sun > 50.0
    assert sky_clock.clock_lag_fraction(None, None, [(None, 5.0), (SUN_M, None)]) == 0.0      # what is missing adds nothing


def test_the_suns_lens_focuses_the_light_that_grazes_it_548_au_out():
    assert sky_clock.focal_distance_m(SUN_M, SUN_R) / constants.AU_M == pytest.approx(548.0, rel=0.01)
    assert sky_clock.focal_distance_m(0.6 * SUN_M, 6000.0) < 0.1 * constants.AU_M                      # a white dwarf focuses just outside itself
    assert sky_clock.focal_distance_m(10.0 * SUN_M, 30.0) is None                                    # inside four Schwarzschild radii there is no weak-field focus
    assert sky_clock.focal_distance_m(0.0, 10.0) is None and sky_clock.focal_distance_m(SUN_M, None) is None


def test_the_moon_eclipses_the_sun_about_every_year_and_a_fifth_for_two_hours():
    moon = sky_clock.moon_eclipse(SUN_R, KM, 1.0, 1737.4, 384_400.0, 27.3216 / 365.25, 5.145)
    assert 1.0 < moon.interval_yr < 1.6
    assert 1.5 < moon.longest_contact_h < 3.0
    assert moon.longest_total_h == 0.0 and 0.8 < moon.deepest_cover < 1.0      # the Moon is the slightly smaller disk of the two: an annular eclipse at best
    flat = sky_clock.moon_eclipse(SUN_R, KM, 1.0, 1737.4, 384_400.0, 27.3216 / 365.25, 0.0)
    assert flat.interval_yr == pytest.approx(1.0 / (365.25 / 27.3216 - 1.0), rel=1e-6)      # in the plane: every synodic period, 29.5 days
    big = sky_clock.moon_eclipse(SUN_R, KM, 1.0, 2000.0, 384_400.0, 27.3216 / 365.25, 5.145)
    assert big.longest_total_h > 0.0 and big.deepest_cover == 1.0


def test_a_tilted_orbit_eclipses_less_often_and_a_missing_moon_never():
    flat = sky_clock.moon_eclipse(SUN_R, KM, 1.0, 1737.4, 384_400.0, 0.0748, 0.0)
    tilted = sky_clock.moon_eclipse(SUN_R, KM, 1.0, 1737.4, 384_400.0, 0.0748, 20.0)
    assert tilted.interval_yr > 3.0 * flat.interval_yr
    assert sky_clock.moon_eclipse(SUN_R, KM, 1.0, 0.0, 384_400.0, 0.0748, 5.0) is None
    assert sky_clock.moon_eclipse(SUN_R, KM, 1.0, 1737.4, 0.0, 0.0748, 5.0) is None
    assert sky_clock.touching_share(0.0, 0.01, 0.1) == 0.0 and sky_clock.touching_share(0.01, 0.01, 0.0) == 1.0 and sky_clock.touching_share(0.01, 0.01, 0.015) == 1.0


def test_the_mean_eclipse_is_shorter_than_the_longest_and_agrees_with_the_orbit_averaged_share():
    """`eclipses.orbiting_the_observer` gives the share of the time the star is touched, the calendar gives how often: their ratio is the mean length, under the central one."""
    args = (SUN_R, 1.0, 1.0, 1737.4, 384_400.0, 27.3216 / 365.25, 5.145, KM)
    _, touched = eclipses.orbiting_the_observer(*args)
    moon = sky_clock.moon_eclipse(SUN_R, KM, 1.0, 1737.4, 384_400.0, 27.3216 / 365.25, 5.145)
    mean_h = touched * moon.interval_yr * YEAR / 3600.0
    assert 0.5 * moon.longest_contact_h < mean_h < moon.longest_contact_h


def test_venus_crosses_the_sun_for_the_earth_in_about_eight_hours_every_few_decades():
    venus = sky_clock.planet_eclipse(SUN_R, KM, 1.0, 6051.8, 0.723 * KM, 0.6152, 3.39)
    assert 6.0 < venus.longest_contact_h < 10.0
    assert 50.0 < venus.interval_yr < 100.0                              # the real cycle averages 61 years, the circular orbits 81
    assert venus.deepest_cover < 0.01
    _, touched = eclipses.orbiting_the_star(SUN_R, 1.0, 1.0, 6051.8, 0.723, 0.6152, 3.39, KM)
    assert touched * venus.interval_yr * YEAR / 3600.0 < venus.longest_contact_h
    assert sky_clock.planet_eclipse(SUN_R, KM, 1.0, 6371.0, 1.5 * KM, 1.8, 1.0) is None          # a blocker outside the observer's orbit


def test_the_direction_to_the_star_turns_with_the_observer_as_well_as_the_blocker():
    rate = eclipses.transit_rate(1.0, 1.0, 0.723, 0.6152)
    assert rate == pytest.approx(3.9 / 60.0 * math.pi / 180.0 * 24.0 * 365.25, rel=0.1)      # about 3.9 arc minutes an hour at a Venus transit
    assert rate > abs(2.0 * math.pi * (0.723 / 0.6152 - 1.0)) / (1.0 - 0.723)


class _Result:
    def __init__(self, rows):
        self._rows = rows

    def fetchone(self):
        return self._rows[0] if self._rows else None

    def fetchall(self):
        return list(self._rows)


class _Conn:
    """Stored rows by the table named in a query (`FROM <table> WHERE`), whatever the filter."""

    def __init__(self, tables):
        self.tables = tables

    def execute(self, sql, params=()):
        table = sql.split(" FROM ")[1].split()[0]
        rows = self.tables.get(table, [])
        if table == "planets" and sql.endswith("WHERE id = ?"):
            rows = [row for row in rows if row["id"] == params[0]]
        elif table == "planets":
            rows = [row for row in rows if row["id"] != params[0] and row["distance_km"] < params[1]]
        return _Result(rows)


def _system():
    star = {"id": 1, "mass_kg": SUN_M, "radius_km": SUN_R}
    earth = {"id": 5, "star_system_id": 1, "star_id": 1, "name": "Earth", "mass_kg": EARTH_M, "radius_km": EARTH_R, "distance_km": KM, "period_years": 1.0,
             "orbital_inclination_deg": 0.0}
    venus = {"id": 4, "star_system_id": 1, "star_id": 1, "name": "Venus", "mass_kg": 4.867e24, "radius_km": 6051.8, "distance_km": 0.723 * KM, "period_years": 0.6152,
             "orbital_inclination_deg": 3.39}
    moon = {"id": 9, "parent_type": "planet", "parent_id": 5, "name": "Moon", "mass_kg": 7.342e22, "radius_km": 1737.4, "distance_km": 384_400.0,
            "period_years": 27.3216 / 365.25, "orbital_inclination_deg": 5.145}
    return _Conn({"stars": [star], "planets": [earth, venus], "moons": [moon]})


def test_the_sheet_rows_of_a_star_a_planet_and_a_moon():
    conn = _system()
    star = {row[1]: row for row in clock_rows.rows(conn, "star", 1)}
    assert star[clock_rows.LENS_KEY][4] / constants.AU_M == pytest.approx(548.0, rel=0.01)
    assert star[clock_rows.CLOCK_KEY][4] == pytest.approx(67.0, rel=0.02)
    planet = {row[1]: row for row in clock_rows.rows(conn, "planet", 5)}
    assert planet[clock_rows.CLOCK_KEY][4] == pytest.approx(0.489, rel=0.02)
    eclipses_of_star = [row for key, row in planet.items() if clock_rows.is_derived_key(key) and key.startswith(clock_rows.ECLIPSE_KEY_PREFIX)]
    assert [row[2] for row in eclipses_of_star] == ["Eclipses of its star by Moon", "Eclipses of its star by Venus"]      # the shortest interval first
    moon = {row[1]: row for row in clock_rows.rows(conn, "moon", 9)}
    assert 0.46 < moon[clock_rows.CLOCK_KEY][4] < planet[clock_rows.CLOCK_KEY][4]      # a clock on the Moon is above the Earth's surface potential
    assert clock_rows.rows(conn, "belt", 1) == []
