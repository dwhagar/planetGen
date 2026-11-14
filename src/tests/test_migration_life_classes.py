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

pytestmark = pytest.mark.fresh_database


def test_every_chemistry_folds_into_a_class_that_exists():
    from planetgen.db.migrations.versions import __name__ as package
    module = __import__(f"{package}.0175_life_classes", fromlist=["CLASS_OF_CHEMISTRY"])
    assert set(module.CLASS_OF_CHEMISTRY) == {"Retinal", "Melanin", "Chlorophyll a", "Blue-Optimized Porphyrins", "Bacteriochlorophylls",
                                              "Zinc-Bacteriochlorophyll"}
    for number, letter in module.CLASS_OF_CHEMISTRY.values():
        assert f"{number}{letter}" in life_classes.BY_ID


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
