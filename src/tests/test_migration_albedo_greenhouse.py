# tests/test_migration_albedo_greenhouse.py

"""
GEN.412: a saved system keeps the albedo and greenhouse quantity its bodies drew; migration 0142 gives a galaxy made at
v141 both, worked out exactly from the stored scale height (a seeded class draw where that has no solution), and the stored
temperature never moves.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import climate_draws, orbit_elements
from planetgen.util import draw

pytestmark = pytest.mark.fresh_database

COLUMNS = ("bond_albedo", "greenhouse_factor", "albedo_basis")


def _system_with_moons():
    for seed in range(60):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = False
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if sum(len(p.moons) for p in system.planets if p.body_type != "a") >= 2:
            return system, cfg
    pytest.fail("no system with two moons in 60 seeds")


def _rows(conn, system_id):
    out = {}
    for table in ("planets", "moons"):
        for row in conn.execute(
                f"SELECT id, surface_temperature_k, {', '.join(COLUMNS)} FROM {table} WHERE star_system_id = ? AND body_type <> 'a'",
                (system_id,)).fetchall():
            out[(table, row["id"])] = dict(row)
    return out


def test_a_saved_system_keeps_its_draws_and_the_temperature_follows_them(mysql_config):
    system, cfg = _system_with_moons()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        rows = _rows(conn, system_id)
        assert rows
        luminosity = conn.execute("SELECT luminosity_w FROM stars WHERE star_system_id = ?", (system_id,)).fetchone()[
            "luminosity_w"]
        orbits = {row["id"]: (row["distance_km"], row["eccentricity"]) for row in conn.execute(
            "SELECT id, distance_km, eccentricity FROM planets WHERE star_system_id = ?", (system_id,)).fetchall()}
        for (table, body_id), row in rows.items():
            assert row["albedo_basis"] == climate_draws.BASIS_GENERATED
            assert 0.0 <= row["bond_albedo"] < 1.0 and row["greenhouse_factor"] >= 0.0
            if table == "planets":
                distance_km, eccentricity = orbits[body_id]
                # the climate goes by the mean light of the orbit (GEN.727)
                flux = climate_draws.starlight_w_m2(luminosity, distance_km) * orbit_elements.mean_flux_factor(eccentricity)
                # no body is held at the ionization ceiling of its air (GEN.1580: the suite makes one an error), so the draws give the temperature
                assert row["surface_temperature_k"] == pytest.approx(
                    climate_draws.surface_temperature_k(row["bond_albedo"], row["greenhouse_factor"], flux), rel=1e-9)
        loaded = store.load_star_system(conn, system_id)
        for planet, original in zip((p for p in loaded.planets if p.body_type != "a"),
                                    (p for p in system.planets if p.body_type != "a")):
            assert planet.bond_albedo == pytest.approx(original.bond_albedo)
            assert planet.greenhouse_factor == pytest.approx(original.greenhouse_factor)
            assert planet.albedo_basis == original.albedo_basis
    finally:
        conn.rollback()
        conn.close()
