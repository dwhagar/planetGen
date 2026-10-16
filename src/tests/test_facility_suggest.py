# tests/test_facility_suggest.py

"""
Suggested facility ratings and the example facility type (POP.10): the
scoring rules in `planetgen.population.facility_suggest` and the routes that
fill unset ratings, against a real throwaway database.
"""

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.population import facility_suggest as suggest
from tests.test_facility_attrs import _base, _create, world  # noqa: F401 (fixtures and helpers)
from tests.test_facility_types import admin, anonymous  # noqa: F401 (fixtures)

DIMENSIONS = ("crime", "housing", "resources", "maintenance", "health")


def _inputs(**changes):
    found = {"habitability": 0.5, "richness": 0.5, "size": 0.5, "nearness": 0.5, "affiliated": True,
             "energy": 0.5, "materials": 0.5, "information": 0.5, "medical": 0.5, "propulsion": 0.5,
             "defense": 0.5}
    found.update(changes)
    return found


@pytest.mark.parametrize("score, color", [
    (0.0, "red"), (0.32, "red"), (0.33, "yellow"), (0.65, "yellow"), (0.66, "green"), (1.0, "green"),
])
def test_scores_become_colours_at_the_thresholds(score, color):
    assert suggest.color_for(score) == color


def test_the_same_facility_always_gets_the_same_answer_and_the_noise_is_small():
    inputs = _inputs()
    first = suggest.scores(inputs, 17)
    assert suggest.scores(inputs, 17) == first
    assert suggest.scores(inputs, 18) != first
    for score in first.values():
        assert abs(score - 0.5) <= tuning.FACILITY_RATING_NOISE + 1e-9      # every weight set adds up to 1


def test_a_better_input_never_lowers_a_score_beyond_the_noise():
    low = suggest.scores(_inputs(habitability=0.0, medical=0.0, materials=0.0, energy=0.0, information=0.0,
                                 nearness=0.0, richness=0.0, size=0.0), 3)
    high = suggest.scores(_inputs(habitability=1.0, medical=1.0, materials=1.0, energy=1.0, information=1.0,
                                  nearness=1.0, richness=1.0, size=1.0), 3)
    for dimension in DIMENSIONS:
        assert high[dimension] > low[dimension]
    assert suggest.suggest(_inputs(habitability=0.0, medical=0.0, materials=0.0, energy=0.0, information=0.0,
                                   nearness=0.0, richness=0.0, size=0.0), 3) == {d: "red" for d in DIMENSIONS}
    assert set(suggest.suggest(_inputs(habitability=1.0, medical=1.0, materials=1.0, energy=1.0, information=1.0,
                                       nearness=1.0, richness=1.0, size=1.0), 3).values()) == {"green"}


def test_nobody_policing_a_facility_costs_it_in_crime_only():
    policed = suggest.scores(_inputs(), 9)
    alone = suggest.scores(_inputs(affiliated=False), 9)
    assert alone["crime"] == pytest.approx(policed["crime"] - tuning.FACILITY_RATING_UNAFFILIATED_PENALTY)
    assert {k: v for k, v in alone.items() if k != "crime"} == {k: v for k, v in policed.items() if k != "crime"}


def test_suggesting_fills_only_what_is_unset(admin, anonymous, world):
    facility_id = _create(admin, world, name="Blank")
    url = f"/api/facilities/{facility_id}/suggest-ratings"

    filled = admin.post(url).get_json()
    assert sorted(filled["filled"]) == sorted(DIMENSIONS)
    found = filled["facility"]
    assert found["ratings_source"] == "generated" and found["ratings_set_at"]
    assert all(found["ratings"][name] in ("green", "yellow", "red") for name in DIMENSIONS)
    assert admin.post(url).get_json()["filled"] == []                      # nothing left to fill

    mixed = _create(admin, world, name="Mixed", ratings={"crime": "red"})
    result = admin.post(f"/api/facilities/{mixed}/suggest-ratings").get_json()
    assert sorted(result["filled"]) == sorted(set(DIMENSIONS) - {"crime"})
    assert result["facility"]["ratings"]["crime"] == "red"                 # the admin's word stands
    assert result["facility"]["ratings_source"] == "admin"

    assert admin.post("/api/facilities/FFFFFFFFFFFFFF-FFFFFFF-FFF/suggest-ratings").status_code == 404
    assert anonymous.post(url).status_code == 401


def test_the_suggestion_reads_the_host_the_species_and_the_distance(admin, world, mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("UPDATE species SET medical_index = 7, materials_index = 7, energy_index = 7, "
                         "information_index = 7")
            conn.execute("INSERT INTO system_owners (star_system_id, polity_id, distance_ly) VALUES (?, ?, 0)",
                         (world["system_id"], world["polity_id"]))
    finally:
        conn.close()
    strong_public = _create(admin, world, name="Strong", affiliation_polity_id=world["polity_id"])
    conn = store.get_connection(mysql_config)
    try:
        strong = conn.execute("SELECT id FROM facilities WHERE name = 'Strong'").fetchone()["id"]
        inputs = suggest.gather(conn, strong)
    finally:
        conn.close()
    assert inputs["affiliated"] is True and inputs["nearness"] == 1.0
    assert (inputs["medical"], inputs["information"]) == (1.0, 1.0) and inputs["propulsion"] == 0.5
    assert inputs["habitability"] == tuning.FACILITY_RATING_ORBITAL_HABITABILITY and inputs["richness"] == 0.5
    assert inputs["size"] == tuning.FACILITY_KIND_SIZE["outpost"]
    # With the best species, the nearest system and a sealed orbit, nothing rates Red.
    assert "red" not in suggest.suggest(inputs, strong).values()


def test_the_batch_fills_every_facility_with_an_unset_rating(admin, world, mysql_config):
    first = _create(admin, world, name="One")
    second = _create(admin, world, name="Two", ratings={d: "green" for d in DIMENSIONS})
    third = _create(admin, world, name="Three", ratings={"health": "red"})
    done = admin.post("/api/facilities/suggest-ratings").get_json()
    assert done == {"facilities": 2, "ratings": 5 + 4}
    assert admin.post("/api/facilities/suggest-ratings").get_json() == {"facilities": 0, "ratings": 0}
    by_id = {item["id"]: item for item in admin.get("/api/facilities").get_json()["items"]}
    assert by_id[first]["ratings_source"] == "generated"
    assert by_id[second]["ratings_source"] == "admin"
    assert by_id[third]["ratings"]["health"] == "red" and by_id[third]["ratings"]["crime"] is not None


def test_the_example_type_installs_once_and_works(admin, anonymous, world):
    first = admin.post("/api/facility-types/examples")
    assert first.status_code == 201 and first.get_json()["added"] == ["Research outpost"]
    second = admin.post("/api/facility-types/examples")
    assert second.status_code == 200 and second.get_json()["added"] == []
    [example] = anonymous.get("/api/facility-types").get_json()["items"]
    assert example["kind"] == "outpost" and [f["key"] for f in example["fields"]] == ["staff", "focus", "armed"]

    facility_id = _create(admin, world, type_id=example["id"])
    found = anonymous.get(f"/api/facilities/{facility_id}").get_json()
    assert found["fields"] == {"staff": 12, "focus": "Survey"}              # the defaults; `armed` has none
    assert anonymous.post("/api/facility-types/examples").status_code == 401
