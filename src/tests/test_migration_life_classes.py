# tests/test_migration_life_classes.py

"""
GEN.889: migration 0175 swaps `life_chemical` for the life class columns and adds the class score tables. A body from before v175 keeps the
class its chemistry stood for and nothing else; a saved system keeps its class, grade, complexity and the classes it supports.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import life_classes
from planetgen.util import draw
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

NEW_COLUMNS = ("life_class", "life_class_variant", "life_class_best", "life_class_best_score", "life_animal_grade", "life_complexity")


def _types(conn, table):
    return {row["COLUMN_NAME"]: row["DATA_TYPE"] for row in conn.execute(
        "SELECT COLUMN_NAME, DATA_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", (table,)).fetchall()}


def test_every_chemistry_folds_into_a_class_that_exists():
    from planetgen.db.migrations.versions import __name__ as package
    module = __import__(f"{package}.0175_life_classes", fromlist=["CLASS_OF_CHEMISTRY"])
    assert set(module.CLASS_OF_CHEMISTRY) == {"Retinal", "Melanin", "Chlorophyll a", "Blue-Optimized Porphyrins", "Bacteriochlorophylls",
                                              "Zinc-Bacteriochlorophyll"}
    for number, letter in module.CLASS_OF_CHEMISTRY.values():
        assert f"{number}{letter}" in life_classes.BY_ID


def test_a_v174_galaxy_keeps_the_class_of_its_chemistry(mysql_config):
    load_old_schema(mysql_config, 174)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("SET SESSION sql_mode = ''")
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        for table in ("planets", "moons"):
            assert "life_chemical" in _types(conn, table)
            assert not set(NEW_COLUMNS) & set(_types(conn, table))
        for name, chemistry in (("Chlorophyll world", "Chlorophyll a"), ("Blue world", "Blue-Optimized Porphyrins"), ("Dark world", "Melanin"),
                                ("Bare world", None)):
            conn.execute("INSERT INTO planets (name, body_type, life_chemical) VALUES (?, 't', ?)", (name, chemistry))
        conn.execute("INSERT INTO moons (name, body_type, life_chemical) VALUES (?, 't', ?)", ("Mat moon", "Bacteriochlorophylls"))
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            types = _types(conn, table)
            assert "life_chemical" not in types
            assert types["life_class"] == "tinyint" and types["life_class_variant"] == "char" and types["life_complexity"] == "double", table
        found = {row["name"]: (row["life_class"], row["life_class_variant"], row["life_complexity"]) for row in conn.execute(
            "SELECT name, life_class, life_class_variant, life_complexity FROM planets").fetchall()}
        assert found == {"Chlorophyll world": (1, "a", None), "Blue world": (1, "b", None), "Dark world": (3, "b", None),
                         "Bare world": (None, None, None)}
        moon = conn.execute("SELECT life_class, life_class_variant FROM moons").fetchone()
        assert (moon["life_class"], moon["life_class_variant"]) == (2, None)
        for table in ("planet_life_class_scores", "moon_life_class_scores"):
            assert conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"] == 0
    finally:
        conn.close()


def test_a_saved_system_keeps_the_class_and_the_classes_it_supports(mysql_config):
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    system = None
    for seed in range(200):
        with draw.bound(seed):
            candidate = StarSystem(system_config=cfg)
        if any(p.life_class is not None for p in candidate.all_planets):
            system = candidate
            break
    assert system is not None, "no system with a life class in 200 seeds"
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        made = {p.name: p for p in system.all_planets if p.body_type != "a"}
        checked = 0
        for planet in loaded.all_planets:
            if planet.body_type == "a":
                continue
            original = made[planet.name]
            for name in life_classes.BODY_FIELDS:
                assert getattr(planet, name) == pytest.approx(getattr(original, name)), (planet.name, name)
            assert (getattr(planet, life_classes.SCORES_FIELD) or {}) == {
                key: value for key, value in (getattr(original, life_classes.SCORES_FIELD) or {}).items() if value > 0}, planet.name
            checked += planet.life_class is not None
        assert checked > 0
    finally:
        conn.rollback()
        conn.close()
