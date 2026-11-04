# tests/test_life_coexist_reads.py

"""GEN.1209 part 2: what reads the several classes of a body: the system detail's `life_classes`, the any-class filter, the class that built a
civilization, and the page's line about the other classes."""

import pytest

from planetgen.db import query, store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.population import model as population
from planetgen.util import draw
from planetgen.web.lib import planet_layout

SEEDS = range(5000, 5016)


def _config():
    config = SystemConfig()
    config.PLANETS = True
    return config


def _mixed_system():
    for seed in SEEDS:
        config = _config()
        with draw.bound(seed):
            system = StarSystem(system_config=config)
        if any(len(planet.life_living or ()) > 1 for planet in system.all_planets if planet.body_type != "a"):
            return system, config
    pytest.skip("no planet with several classes in these seeds")


def test_a_body_with_no_life_has_no_classes_and_one_stored_before_has_its_own():
    assert query._life_classes_of({"life_class": None, "life_class_variant": None, "life_stage": None}, None) == []
    old = {"life_class": 3, "life_class_variant": "a", "life_stage": "complex_cells", "life_animal_grade": 0, "life_microbial": 4, "life_flora": 0,
           "life_fauna": 0, "life_sentient": 0}
    assert query._life_classes_of(old, None) == [{"life_class_id": "3a", "share": 1.0, "life_stage": "complex_cells", "animal_grade": 0,
                                                   "kinds": {"microbial": 4, "flora": 0, "fauna": 0, "sentient": 0}}]


def test_the_dominant_class_comes_first_in_the_list():
    body = {"life_class": 4, "life_class_variant": "b", "life_stage": "multicellularity"}
    living = [{"life_class": "3a", "share": 0.4, "stage": "complex_cells", "animal_grade": 0, "kinds": [3, 0, 0, 0]},
              {"life_class": "4b", "share": 0.6, "stage": "multicellularity", "animal_grade": 1, "kinds": [5, 0, 1, 0]},
              {"life_class": "2", "share": 0.1, "stage": None, "animal_grade": 0, "kinds": [0, 0, 0, 0]}]
    listed = query._life_classes_of(body, living)
    assert [entry["life_class_id"] for entry in listed] == ["4b", "3a"], "the dominant class first; a class with no stage is not listed"
    assert listed[0]["kinds"] == {"microbial": 5, "flora": 0, "fauna": 1, "sentient": 0}


def test_the_filter_matches_a_body_when_any_class_that_lives_on_it_is_asked_for():
    clauses, params = [], []
    query._life_class_clause("p", ["3a", "4b"], clauses, params)
    sql = clauses[0]
    assert "planet_life" in sql and "EXISTS" in sql and "p.life_class" in sql
    assert params == [3, "a", 4, "b", "3a", "4b"]
    clauses, params = [], []
    query._life_class_clause("m", ["6"], clauses, params)
    assert "moon_life" in clauses[0] and params == [6, "", "6"]


def test_the_page_names_the_other_classes_that_live_beside_the_dominant_one():
    row = {"life_class_id": "4b", "life_stage": "multicellularity", "life_animal_grade": 1, "life_class_best": "4b", "life_class_best_score": 7,
           "life_class_scores": {"4b": 7, "3a": 5},
           "life_classes": [{"life_class_id": "4b", "share": 0.7, "life_stage": "multicellularity", "animal_grade": 1, "kinds": {}},
                            {"life_class_id": "3a", "share": 0.3, "life_stage": "complex_cells", "animal_grade": 0, "kinds": {}}]}
    view = planet_layout.life_class_view(row)
    assert view["alongside"].startswith("Also living here: Class 3a (30 percent of the food production, at the stage of ")
    assert view["alongside"] in view["paragraph"]
    row["life_classes"] = row["life_classes"][:1]
    assert planet_layout.life_class_view(row)["alongside"] is None


def test_the_system_detail_search_facet_and_species_read_every_class(mysql_config, monkeypatch):
    from planetgen.db import species_fit
    system, config = _mixed_system()
    system_id = store.save_system(system, config, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        detail = query.system_detail(conn, system_id)
        mixed = [planet for planet in detail["planets"] if len(planet["life_classes"]) > 1]
        assert mixed and mixed[0]["life_classes"][0]["life_class_id"] == mixed[0]["life_class_id"]
        minor = mixed[0]["life_classes"][1]["life_class_id"]
        assert minor in {entry["value"] for entry in query._life_class_facet(conn, "planets")}
        # the civilization is built by the class that reached the milestone, not by the dominant one
        planet_id = mixed[0]["id"]
        with conn:
            conn.execute("UPDATE planet_life SET life_stage = 'technological_civilization', life_milestone_age_gyr = 1.0, life_system_age_gyr = 4.0 "
                         "WHERE planet_id = ? AND life_class = ?", (planet_id, minor))
            conn.execute("UPDATE system_configs SET intelligent_life = 1 WHERE id = (SELECT system_config_id FROM star_systems WHERE id = ?)", (system_id,))
            population._set_watermark(conn, 0)
        monkeypatch.setattr(species_fit, "check_homeworld", lambda conn, planet_id: None)
        monkeypatch.setattr(population, "has_civilization", lambda *args, **kwargs: True)
        assert population.scan_life_worlds(conn) >= 1
        found = conn.execute("SELECT life_class FROM species WHERE homeworld_planet_id = ?", (planet_id,)).fetchone()
        assert found["life_class"] == minor
    finally:
        conn.close()
