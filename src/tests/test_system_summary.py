# tests/test_system_summary.py

"""
GEN.1478: what each star system adds up to is kept in `system_summary` by the system save and the editor, equals the sums the readers used to make from
the rows, and is what the gravity sources and the sector paths now read.
"""

import importlib

import pytest

from planetgen.db import edits, gravity_sources, sector_paths, store, system_summary
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests.fuzz_support import deterministic_entropy

pytestmark = pytest.mark.db


def _save(mysql_config, seed):
    with deterministic_entropy(seed):
        draw.set_run_seed(seed)
        cfg = SystemConfig()
        system = StarSystem(cfg, galactic_center_dist_ly=26000.0)
    return store.save_system(system, cfg, config=mysql_config)


def _by_hand(conn, system_id):
    """The figures worked out from the rows, the way the readers did before the columns."""
    def sums(table, column):
        row = conn.execute(f"SELECT SUM({column}) AS s, COUNT(*) AS n, MAX({column}) AS m FROM {table} WHERE star_system_id = ?", (system_id,)).fetchone()
        return row["s"], row["n"], row["m"]
    stars, star_count, _ = sums("stars", "mass_kg")
    planets, planet_count, _ = sums("planets", "mass_kg")
    moons, moon_count, _ = sums("moons", "mass_kg")
    light, _, brightest = sums("stars", "luminosity_w")
    return {"stars_mass_kg": stars, "system_mass_kg": (stars or 0) + (planets or 0) + (moons or 0), "brightest_star_luminosity_w": brightest,
            "total_luminosity_w": light, "star_count": star_count, "planet_count": planet_count, "moon_count": moon_count,
            "belt_count": sums("asteroid_belts", "mass_kg")[1], "comet_count": sums("comets", "id")[1]}


def _stored(conn, system_id):
    row = conn.execute("SELECT * FROM system_summary WHERE star_system_id = ?", (system_id,)).fetchone()
    return None if row is None else {name: row[name] for name in system_summary.COLUMNS}


def test_a_saved_system_holds_the_figures_of_its_rows(mysql_config):
    ids = [_save(mysql_config, seed) for seed in range(6)]
    conn = store.get_connection(mysql_config)
    try:
        for system_id in ids:
            stored = _stored(conn, system_id)
            assert stored == pytest.approx(_by_hand(conn, system_id), rel=1e-12)
            assert stored["star_count"] >= 1 and stored["stars_mass_kg"] > 0 and stored["system_mass_kg"] >= stored["stars_mass_kg"]
            assert stored["brightest_star_luminosity_w"] <= stored["total_luminosity_w"]
    finally:
        conn.close()


def test_an_edit_that_removes_a_body_changes_the_figures(mysql_config):
    for seed in range(60):
        system_id = _save(mysql_config, seed)
        conn = store.get_connection(mysql_config)
        try:
            system = store.load_star_system(conn, system_id)
            planets = [p for p in system.planets if p.body_type != "a"]
            if len(planets) < 2:
                continue
            before = _stored(conn, system_id)
            gone = planets[-1]
            system.planets.remove(gone)
            with conn:
                edits.save_system_edits(conn, system_id, system)
            after = _stored(conn, system_id)
            assert after == pytest.approx(_by_hand(conn, system_id), rel=1e-12)
            assert after["planet_count"] == before["planet_count"] - 1
            assert after["system_mass_kg"] != before["system_mass_kg"]
            return
        finally:
            conn.close()
    pytest.fail("no system with two planets in 60 seeds")


def test_the_readers_see_the_stored_figures(mysql_config):
    system_id = _save(mysql_config, 3)
    conn = store.get_connection(mysql_config)
    try:
        sector_id = conn.execute("SELECT sector_id FROM star_systems WHERE id = ?", (system_id,)).fetchone()["sector_id"]
        stored = _stored(conn, system_id)
        if sector_id is None:
            pytest.skip("a system saved alone has no sector")
        masses = gravity_sources.system_point_masses(conn, [sector_id])
        assert float(masses.masses_solar[masses.keys.index(("star_systems", system_id))]) == pytest.approx(stored["system_mass_kg"] / 1.98892e30, rel=1e-3)
    finally:
        conn.close()


def test_the_migrations_statement_gives_the_figures_the_save_does(mysql_config):
    system_id = _save(mysql_config, 4)
    migration = importlib.import_module("planetgen.db.migrations.versions.0237_system_summary_figures")
    conn = store.get_connection(mysql_config)
    try:
        row = conn.execute(migration.SELECT.replace(":low", "?").replace(":high", "?"), (system_id, system_id)).fetchone()
        assert list(row.values())[1:] == pytest.approx([_stored(conn, system_id)[name] for name in system_summary.COLUMNS], rel=1e-12)
    finally:
        conn.close()
