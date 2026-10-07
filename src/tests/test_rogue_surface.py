"""
Rogue planet surface conditions (`rogueSurface`, schema v48): with no
star, a rogue's surface is set by its own heat. Checks Boss's worked
numbers (Earth-like rogue ~35 K, ~12 km ice lid, Jupiter-mass rogue
~99 K / 130-170 K at 1 bar, Stevenson's warm envelope base), every
regime, storage and the v48 back-fill.
"""

import math
import random

import pytest

from planetgen.db import store
from planetgen.physics import constants
from planetgen import tuning
from planetgen.physics import rogue_surface as rs
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet

import web  # noqa: F401 -- puts src/html/lib on sys.path
from api.app import create_app  # noqa: E402
from api.config import Config  # noqa: E402
from web import system_pages  # noqa: E402

EARTH = constants.EARTH_MASS_TO_KG
JUPITER = constants.JUPITER_MASS_TO_KG


def test_earth_like_rogue_matches_boss_numbers():
    flux = rs.radiogenic_flux_w_m2(EARTH, 6371, 4.5) / tuning.ROGUE_UREY_RATIO
    assert 0.08 <= flux <= 0.1
    t_eff = rs.effective_temperature_k(flux)
    assert t_eff == pytest.approx(35.5, abs=1.0)
    # D = (A / F) ln(273 / 35) ~ 11.6 km at 0.1 W/m^2.
    assert rs.ice_shell_thickness_km(0.1, 35.0) == pytest.approx(11.6, abs=0.2)


def test_no_heat_means_the_cmb_temperature():
    assert rs.effective_temperature_k(0.0) == pytest.approx(constants.COSMIC_BACKGROUND_TEMPERATURE_K)


def test_jupiter_mass_rogue_matches_boss_numbers():
    flux = rs.giant_internal_flux_w_m2(JUPITER, constants.JUPITER_RADIUS_KM, 4.5)
    assert flux == pytest.approx(5.4, rel=1e-6)
    t_eff = rs.effective_temperature_k(flux)
    assert t_eff == pytest.approx(99, abs=1)
    gravity = rs.surface_gravity_m_s2(JUPITER, constants.JUPITER_RADIUS_KM)
    assert 130 <= rs.giant_temperature_at_k(t_eff, gravity, 1e5) <= 170


def test_giants_cool_with_age_and_heavier_ones_are_hotter():
    for mass in (0.1, 1.0, 5.0, 13.0, 40.0, 80.0):
        temps = [rs.effective_temperature_k(rs.giant_internal_flux_w_m2(mass * JUPITER, 70000, t))
                 for t in (0.5, 2.0, 8.0)]
        assert temps == sorted(temps, reverse=True)
    by_mass = [rs.effective_temperature_k(rs.giant_internal_flux_w_m2(m * JUPITER, 70000, 4.5))
               for m in (0.06, 0.3, 1.0, 5.0, 13.0, 30.0, 80.0)]
    assert by_mass == sorted(by_mass)
    hottest = rs.effective_temperature_k(rs.giant_internal_flux_w_m2(80 * JUPITER, 55000, 0.1))
    assert hottest <= tuning.ROGUE_MAX_EFFECTIVE_TEMPERATURE_K + 1e-6


def test_stevenson_envelope_warms_the_ground():
    """30 K at 0.1 bar, down a gamma = 1.4 adiabat to 1,000 bar: ~417 K."""
    assert rs.adiabat_temperature_k(30.0, 1e4, 1e8) == pytest.approx(417, abs=2)
    assert rs.adiabat_temperature_k(30.0, 1e4, 1e7) > 200


def test_nitrogen_frost_gives_a_pluto_like_trace():
    assert 0.5 < rs.nitrogen_vapor_pressure_pa(40.0) < 20
    assert rs.nitrogen_vapor_pressure_pa(20.0) < 1e-6


def test_every_regime_occurs_and_conditions_are_consistent():
    seen = set()
    cases = [("terrestrial", 't', 0.15 * EARTH, 3200), ("terrestrial", 't', 1.0 * EARTH, 6371),
             ("sub-neptune", 't', 8.0 * EARTH, 12000), ("saturn", 'g', 0.3 * JUPITER, 60000),
             ("brown-dwarf", 'g', 40 * JUPITER, 60000)]
    for seed in range(400):
        rng = random.Random(seed)
        for mass_bin, kind, mass, radius in cases:
            c = rs.rogue_surface_conditions(mass, radius, kind, mass_bin, seed % 3 == 0, rng)
            seen.add(c["surface_regime"])
            assert c["surface_regime"] in rs.SURFACE_REGIMES
            lo, hi = tuning.ROGUE_PLANET_AGE_RANGE_GY
            assert lo <= c["age_gy"] <= hi
            assert c["internal_heat_flux_w_m2"] > 0
            assert c["surface_temperature_k"] >= constants.COSMIC_BACKGROUND_TEMPERATURE_K
            assert math.isfinite(c["surface_temperature_k"]) and c["surface_pressure_pa"] >= 0
            if kind == 'g':
                assert c["surface_pressure_pa"] == 1e5 and c["has_internal_heat"]
            if c["has_liquid_water"]:
                assert c["ocean_depth_km"] > 0
            if c["surface_regime"] == "ice-shell-ocean":
                assert c["ice_shell_thickness_km"] > 0 and c["has_liquid_water"]
            if c["surface_regime"] in ("bare-rock", "frozen-atmosphere", "ice-shell-ocean", "ice-world"):
                assert c["surface_temperature_k"] == pytest.approx(c["effective_temperature_k"])
            if c["surface_regime"] == "hydrogen-envelope":
                assert c["surface_temperature_k"] > c["effective_temperature_k"]
    assert seen == set(rs.SURFACE_REGIMES)


def test_seeded_conditions_repeat():
    a = rs.rogue_surface_conditions(EARTH, 6371, 't', "terrestrial", True, random.Random("Wanderer"))
    b = rs.rogue_surface_conditions(EARTH, 6371, 't', "terrestrial", True, random.Random("Wanderer"))
    assert a == b


def test_rogue_planet_carries_and_describes_its_surface():
    for mass_bin in tuning.ROGUE_PLANET_MASS_BIN_CHOICES:
        planet = RoguePlanet(SystemConfig(), mass_bin=mass_bin)
        text = planet.to_paragraph_list()[1]
        assert "only heat is its own" in text
        data = planet.to_dict()
        for field in rs.ROGUE_SURFACE_FIELDS:
            assert field in data
        old = dict(data)
        for field in rs.ROGUE_SURFACE_FIELDS:
            old.pop(field)
        first = RoguePlanet.from_dict(old, SystemConfig())
        again = RoguePlanet.from_dict(old, SystemConfig())
        assert first.surface_regime == again.surface_regime
        assert first.age_gy == again.age_gy


def test_phenomenon_page_rows():
    class _Config(Config):
        TESTING = True
        SECRET_KEY = "test"

    detail = {"planet_type": "t", "age_gy": 4.5, "internal_heat_flux_w_m2": 0.087, "effective_temperature_k": 35.1,
              "surface_regime": "ice-shell-ocean", "surface_temperature_k": 35.1, "surface_pressure_pa": 0.3,
              "ice_shell_thickness_km": 13.4, "ocean_depth_km": 120.0, "has_internal_heat": 1}
    with create_app(_Config).test_request_context("/"):
        fields = {label: text for label, text, _url in system_pages.phenomenon_fields("rogue_planet", detail)}
    assert fields["Surface"] == "Ice shell over a liquid ocean"
    assert fields["Surface Temperature"] == "35 K (-238 °C)"
    assert fields["Surface Pressure"] == "0.3 Pa (a trace)"
    assert fields["Ice Thickness"] == "13.4 km"
    assert fields["Liquid Ocean Depth"] == "120 km"
    assert fields["Internal Heat Flow"] == "0.087 W/m²"
    assert fields["Geologically Active"] == "Yes"


def test_columns_follow_the_fields():
    assert tuple(column for column, _kind in store.ROGUE_SURFACE_COLUMNS) == rs.ROGUE_SURFACE_FIELDS


def test_stored_and_backfilled_by_the_v48_migration(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            planet = RoguePlanet(SystemConfig(), mass_bin="sub-neptune")
            planet_id = store.insert_rogue_planet(conn, planet)
        row = conn.execute("SELECT * FROM rogue_planets WHERE id = ?", (planet_id,)).fetchone()
        assert row["surface_regime"] == planet.surface_regime
        assert row["surface_temperature_k"] == pytest.approx(planet.surface_temperature_k)
        assert row["has_liquid_water"] == int(planet.has_liquid_water)
        name = row["name"]
        conn.execute("ALTER TABLE rogue_planets DROP COLUMN surface_regime, DROP COLUMN age_gy")
        conn.execute("DELETE FROM schema_migrations WHERE version IN (48, 49, 50, 51, 52, 53)")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (47)")
        conn.commit()
    finally:
        conn.close()

    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION

    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT * FROM rogue_planets WHERE id = ?", (planet_id,)).fetchone()
    finally:
        conn.close()
    expected = rs.rogue_surface_conditions(planet.mass_kg, planet.radius_km, planet.planet_type, "sub-neptune",
                                           planet.has_moons,
                                           random.Random(name))
    assert row["surface_regime"] == expected["surface_regime"]
    assert row["age_gy"] == pytest.approx(expected["age_gy"])
    assert row["has_internal_heat"] == int(expected["has_internal_heat"])
