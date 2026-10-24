# tests/test_migration_block_lights.py

"""
MAP.259: migration 0147 creates `galaxy_block_lights` and fills it from the stored scatter: a row for every block of the
eleven levels above a sector that holds a bright star nothing was built from or a phenomenon not yet built.
"""

import pytest

from planetgen.db import store
from planetgen.galaxy.drill import DRILL_LEVELS, drill_chain_of
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

STARS = [(30, 3, 5), (30, 3, 6), (400, 40, 11)]
POINTS = [(30, 3, 5), (31, 3, 5)]
BUILT = (31, 3, 5)


def _old_galaxy(mysql_config):
    load_old_schema(mysql_config, 146)
    assert store.schema_status(mysql_config)[0] == 146
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for shell, band, wedge in STARS:
            conn.execute(
                "INSERT INTO bright_stars (shell_index, band_index, wedge_index, position_x_mpc, position_y_mpc,"
                " position_z_mpc, population, star_type, yerkes_class, mass_kg, radius_km, temperature_k, luminosity_w,"
                " age_gy, initial_mass_sol, seed) VALUES (?, ?, ?, 1000, 2000, 3000, 'young', 'O5V', 'V', 1, 1, 30000,"
                " 3.828e27, 0.001, 40, 1)", (shell, band, wedge))
        for shell, band, wedge in POINTS:
            conn.execute(
                "INSERT INTO phenomenon_scatter (shell_index, band_index, wedge_index, kind, position_x_mpc,"
                " position_y_mpc, position_z_mpc, seed, map_luminosity_sol, built_at) VALUES (?, ?, ?, 'black-hole', 500,"
                " 500, 500, 1, 2.0, " + ("CURRENT_TIMESTAMP" if (shell, band, wedge) == BUILT else "NULL") + ")",
                (shell, band, wedge))
        conn.commit()
    finally:
        conn.close()


def test_the_migration_fills_every_level_from_the_unbuilt_objects(mysql_config):
    _old_galaxy(mysql_config)
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        rows = conn.execute("SELECT m, shell_index, band_index, wedge_index, star_count, point_count, star_luminosity_sol,"
                            " star_lpx, star_lt FROM galaxy_block_lights").fetchall()
        assert {r["m"] for r in rows} == set(DRILL_LEVELS[:-1])
        stored = {(r["m"], r["shell_index"], r["band_index"], r["wedge_index"]): r for r in rows}
        expected = {}
        for address in STARS:
            for block in drill_chain_of(*address)[:-1]:
                expected.setdefault(tuple(block), [0, 0])[0] += 1
        for address in POINTS:
            if address != BUILT:
                for block in drill_chain_of(*address)[:-1]:
                    expected.setdefault(tuple(block), [0, 0])[1] += 1
        assert {key: [row["star_count"], row["point_count"]] for key, row in stored.items()} == {
            key: value for key, value in expected.items()}
        root = [r for key, r in stored.items() if key[0] == DRILL_LEVELS[0]]
        assert sum(r["star_luminosity_sol"] for r in root) == pytest.approx(len(STARS) * 10.0, rel=1e-6)    # 3.828e27 W
        assert sum(r["star_lpx"] for r in root) == pytest.approx(len(STARS) * 10.0 * 1.0, rel=1e-6)    # 1000 mpc = 1 pc
        assert sum(r["star_lt"] for r in root) == pytest.approx(len(STARS) * 10.0 * 30000, rel=1e-6)
        assert store.schema_status(mysql_config)[0] == store.SCHEMA_VERSION
    finally:
        conn.close()
