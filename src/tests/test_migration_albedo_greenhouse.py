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
from planetgen.physics import climate_draws, constants
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
        distances = {row["id"]: row["distance_km"] for row in conn.execute(
            "SELECT id, distance_km FROM planets WHERE star_system_id = ?", (system_id,)).fetchall()}
        for (table, body_id), row in rows.items():
            assert row["albedo_basis"] == climate_draws.BASIS_GENERATED
            assert 0.0 <= row["bond_albedo"] < 1.0 and row["greenhouse_factor"] >= 0.0
            if table == "planets":
                flux = climate_draws.starlight_w_m2(luminosity, distances[body_id])
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


def test_a_v141_galaxy_gets_the_draws_back_and_keeps_its_temperatures(mysql_config):
    system, cfg = _system_with_moons()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        recorded = _rows(conn, system_id)
    finally:
        conn.close()
    # make it a v141 galaxy again: the columns gone, the rows kept
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            conn.execute(f"ALTER TABLE {table} ADD COLUMN evolutionary_speed VARCHAR(64) NULL")  # v143 (GEN.666)
            conn.execute(f"ALTER TABLE {table} DROP COLUMN bond_albedo, DROP COLUMN greenhouse_factor,"
                         " DROP COLUMN albedo_basis, ALGORITHM=COPY")
        conn.execute("UPDATE alembic_version SET version_num = '0141'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (141)")
        conn.commit()
        assert store._stored_version(conn, "schema_migrations") == 141
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        rebuilt = _rows(conn, system_id)
    finally:
        conn.close()
    assert set(rebuilt) == set(recorded)
    inverted = 0
    for key, row in rebuilt.items():
        before = recorded[key]
        assert row["surface_temperature_k"] == before["surface_temperature_k"]
        assert row["albedo_basis"] in (climate_draws.BASIS_INVERSE, climate_draws.BASIS_CLASS_DRAW)
        assert 0.0 <= row["bond_albedo"] < 1.0 and row["greenhouse_factor"] >= -1e-9
        if row["albedo_basis"] == climate_draws.BASIS_INVERSE:
            inverted += 1
            assert row["bond_albedo"] == pytest.approx(before["bond_albedo"], rel=1e-6, abs=1e-9)
            assert row["greenhouse_factor"] == pytest.approx(before["greenhouse_factor"], rel=1e-6, abs=1e-6)
    assert inverted, "no body could be inverted exactly"


def test_a_row_with_no_solution_takes_a_seeded_class_draw(mysql_config):
    system, cfg = _system_with_moons()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
            # an edited row: a pressure far above the barometric law cannot come from this code
            conn.execute("UPDATE planets SET atmospheric_pressure_pa = 1e12 WHERE star_system_id = ? AND body_type <> 'a'"
                         " AND atmosphere <> 'None'", (system_id,))
        edited = conn.execute("SELECT id FROM planets WHERE star_system_id = ? AND atmospheric_pressure_pa = 1e12",
                              (system_id,)).fetchall()
    finally:
        conn.close()
    if not edited:
        pytest.skip("the system has no planet with air")
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            conn.execute(f"ALTER TABLE {table} ADD COLUMN evolutionary_speed VARCHAR(64) NULL")  # v143 (GEN.666)
            conn.execute(f"ALTER TABLE {table} DROP COLUMN bond_albedo, DROP COLUMN greenhouse_factor,"
                         " DROP COLUMN albedo_basis, ALGORITHM=COPY")
        conn.execute("UPDATE alembic_version SET version_num = '0141'")
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (141)")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for row in edited:
            drawn = conn.execute("SELECT bond_albedo, greenhouse_factor, albedo_basis, surface_temperature_k FROM planets"
                                 " WHERE id = ?", (row["id"],)).fetchone()
            assert drawn["albedo_basis"] == climate_draws.BASIS_CLASS_DRAW
            assert drawn["bond_albedo"] is not None and drawn["greenhouse_factor"] >= -1e-9
    finally:
        conn.close()
