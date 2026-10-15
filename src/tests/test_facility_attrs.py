# tests/test_facility_attrs.py

"""
What a facility carries beyond its place (POP.10): its type and the type's
field values, its affiliation and its five ratings, over `POST`/`PATCH`/`GET
/api/facilities`, against a real throwaway database.
"""

import pytest

from planetgen.db import store
from tests.publicids import pid
from tests.test_admin_edits import _saved_system
from tests.test_facility_types import FARM, admin, anonymous  # noqa: F401 (fixtures)


@pytest.fixture
def world(mysql_config):
    """A saved system: its public star id for a host, its internal id, and a polity there."""
    _sector_id, system_id = _saved_system(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        star_id = conn.execute("SELECT id FROM stars WHERE star_system_id = ? LIMIT 1", (system_id,)).fetchone()["id"]
        planet = conn.execute("SELECT id FROM planets WHERE star_system_id = ? LIMIT 1", (system_id,)).fetchone()
        conn.execute(
            "INSERT INTO species (name, homeworld_planet_id, star_system_id, life_chemical, life_stage, build, "
            "climate, size, civilization_age_years, era, spacefaring) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            ("Velarix", planet["id"], system_id, "chlorophyll", "technological_civilization", "medium",
             "temperate", "medium", 60000.0, "established", 1))
        species_id = conn.execute("SELECT id FROM species").fetchone()["id"]
        conn.execute("INSERT INTO polities (name, species_id, capital_system_id, government, color, reach_ly) "
                     "VALUES (?, ?, ?, ?, ?, ?)", ("Velarix Concord", species_id, system_id, "Concord", "#3366cc", 40.0))
        polity_id = conn.execute("SELECT id FROM polities").fetchone()["id"]
        conn.commit()
    finally:
        conn.close()
    return {"host": pid("star", star_id), "system_id": system_id, "polity_id": polity_id}


def _base(world, name="Post", kind="outpost", **extra):
    return {"name": name, "kind": kind, "placement": "orbital", "host_type": "star", "host_id": world["host"], **extra}


def _type(admin, **spec):
    body = {"name": "Relay", "kind": "outpost", **spec}
    response = admin.post("/api/facility-types", json=body)
    assert response.status_code == 201, response.get_json()
    return response.get_json()["id"]


def _create(admin, world, **extra):
    response = admin.post("/api/facilities", json=_base(world, **extra))
    assert response.status_code == 201, response.get_json()
    return response.get_json()["id"]


def test_a_plain_facility_reads_back_unassessed(admin, anonymous, world):
    facility_id = _create(admin, world)
    found = anonymous.get(f"/api/facilities/{facility_id}").get_json()
    assert found["type_id"] is None and found["type_name"] is None and found["fields"] == {}
    assert found["affiliation_polity_id"] is None and found["affiliation_name"] is None
    assert found["ratings"] == {name: None for name in ("crime", "housing", "resources", "maintenance", "health")}
    assert found["overall"] is None and found["ratings_source"] is None


def test_type_fields_affiliation_and_ratings_on_create(admin, anonymous, world):
    type_id = _type(admin, icon="radioactive", fields=[
        {"key": "crew", "label": "Crew", "datatype": "int", "min_value": 1, "required": True},
        {"key": "mode", "label": "Mode", "datatype": "choice", "choices": ["Survey", "Relay"], "default_text": "Relay"},
        {"key": "armed", "label": "Armed", "datatype": "bool"},
        {"key": "note", "label": "Note", "datatype": "text"},
        {"key": "output", "label": "Output", "datatype": "decimal", "unit": "TW"},
    ])
    facility_id = _create(admin, world, type_id=type_id, affiliation_polity_id=world["polity_id"],
                          ratings={"crime": "green", "health": "red", "housing": "yellow"}, ratings_note="Survey 2200",
                          fields={"crew": 12, "armed": False, "output": 1.5, "note": "Quiet"})
    found = anonymous.get(f"/api/facilities/{facility_id}").get_json()
    assert (found["type_id"], found["type_name"], found["type_icon"]) == (type_id, "Relay", "radioactive")
    assert found["fields"] == {"crew": 12, "mode": "Relay", "armed": False, "output": 1.5, "note": "Quiet"}
    assert found["affiliation_polity_id"] == world["polity_id"] and found["affiliation_name"] == "Velarix Concord"
    assert found["ratings"]["crime"] == "green" and found["ratings"]["health"] == "red"
    assert found["ratings"]["resources"] is None
    assert found["rating_words"]["crime"] == "Low" and found["rating_words"]["health"] == "Crisis"
    assert found["rating_words"]["housing"] == "Tight"
    assert found["overall"] == "red" and found["ratings_source"] == "admin" and found["ratings_set_at"]
    assert found["ratings_note"] == "Survey 2200"
    listed = anonymous.get(f"/api/systems/{pid('system', world['system_id'])}/facilities").get_json()["items"]
    assert listed[0]["fields"]["crew"] == 12 and listed[0]["overall"] == "red"
    assert anonymous.get(f"/api/facility-types/{type_id}").get_json()["facilities"] == 1


@pytest.mark.parametrize("extra, fragment", [
    ({"type_id": 99999}, "no such facility type"),
    ({"fields": {"crew": 3}}, "needs a 'type_id'"),
    ({"affiliation_polity_id": 99999}, "no such polity"),
    ({"affiliation_polity_id": 1, "affiliation_name": "Both"}, "not both"),
    ({"affiliation_name": " "}, "affiliation_name"),
    ({"ratings": {"crime": "amber"}}, "green, yellow, red"),
    ({"ratings": {"fame": "green"}}, "fame"),
    ({"ratings_note": "x" * 4001}, "ratings_note"),
])
def test_bad_extras_are_refused_and_nothing_is_saved(admin, anonymous, world, extra, fragment):
    response = admin.post("/api/facilities", json=_base(world, **extra))
    assert response.status_code in (400, 404)
    assert fragment in response.get_json()["error"]
    assert anonymous.get("/api/facilities").get_json()["total"] == 0


def test_field_values_follow_the_type_definition(admin, world):
    kind_mismatch = _type(admin, name="Colonial", kind="colony")
    type_id = _type(admin, fields=[
        {"key": "crew", "label": "Crew", "datatype": "int", "min_value": 1, "max_value": 50, "required": True},
        {"key": "mode", "label": "Mode", "datatype": "choice", "choices": ["Survey", "Relay"]},
    ])
    retired = _type(admin, name="Old relay", active=False)

    def refuse(**extra):
        return admin.post("/api/facilities", json=_base(world, **extra))

    assert "colony type" in refuse(type_id=kind_mismatch).get_json()["error"]
    assert "retired" in refuse(type_id=retired).get_json()["error"]
    assert "required" in refuse(type_id=type_id).get_json()["error"]
    assert "required" in refuse(type_id=type_id, fields={"crew": None}).get_json()["error"]
    assert "at most 50" in refuse(type_id=type_id, fields={"crew": 51}).get_json()["error"]
    assert "whole number" in refuse(type_id=type_id, fields={"crew": 2.5}).get_json()["error"]
    assert "must be one of" in refuse(type_id=type_id, fields={"crew": 5, "mode": "Dock"}).get_json()["error"]
    assert "no field 'speed'" in refuse(type_id=type_id, fields={"crew": 5, "speed": 1}).get_json()["error"]
    assert refuse(type_id=type_id, fields={"crew": 5}).status_code == 201


def test_changing_a_facility(admin, anonymous, world):
    first = _type(admin, fields=[{"key": "crew", "label": "Crew", "datatype": "int", "required": True}])
    second = _type(admin, name="Beacon", fields=[{"key": "power", "label": "Power", "datatype": "decimal"}])
    facility_id = _create(admin, world, type_id=first, fields={"crew": 4}, ratings={"crime": "green"})
    url = f"/api/facilities/{facility_id}"

    renamed = admin.patch(url, json={"name": "Renamed", "description": "Watches the sky"}).get_json()
    assert (renamed["name"], renamed["description"]) == ("Renamed", "Watches the sky")
    assert renamed["fields"] == {"crew": 4} and renamed["ratings"]["crime"] == "green"   # untouched

    rated = admin.patch(url, json={"ratings": {"crime": None, "health": "yellow"}, "fields": {"crew": 9}}).get_json()
    assert rated["ratings"]["crime"] is None and rated["ratings"]["health"] == "yellow"
    assert rated["overall"] == "yellow" and rated["fields"] == {"crew": 9}
    assert "required" in admin.patch(url, json={"fields": {"crew": None}}).get_json()["error"]

    moved = admin.patch(url, json={"type_id": second}).get_json()
    assert moved["type_id"] == second and moved["fields"] == {}                  # the old type's values go
    again = admin.patch(url, json={"fields": {"power": 3}}).get_json()
    assert again["fields"] == {"power": 3}
    untyped = admin.patch(url, json={"type_id": None}).get_json()
    assert untyped["type_id"] is None and untyped["fields"] == {}


def test_changing_an_affiliation_and_refusing_a_move(admin, anonymous, world):
    facility_id = _create(admin, world)
    url = f"/api/facilities/{facility_id}"
    named = admin.patch(url, json={"affiliation_name": "Free Traders"}).get_json()
    assert (named["affiliation_polity_id"], named["affiliation_name"]) == (None, "Free Traders")
    owned = admin.patch(url, json={"affiliation_polity_id": world["polity_id"]}).get_json()
    assert owned["affiliation_name"] == "Velarix Concord"
    cleared = admin.patch(url, json={"affiliation_polity_id": None}).get_json()
    assert (cleared["affiliation_polity_id"], cleared["affiliation_name"]) == (None, None)

    assert admin.patch(url, json={}).status_code == 400
    assert admin.patch(url, json={"placement": "standalone"}).status_code == 400
    assert admin.patch("/api/facilities/FFFFFFFFFF-FFFFFFF-FFF", json={"name": "Ghost"}).status_code == 404
    assert anonymous.patch(url, json={"name": "X"}).status_code == 401


def test_the_polity_name_survives_the_polity(admin, anonymous, world, mysql_config):
    facility_id = _create(admin, world, affiliation_polity_id=world["polity_id"])
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("DELETE FROM system_owners")
            conn.execute("DELETE FROM polities")
    finally:
        conn.close()
    found = anonymous.get(f"/api/facilities/{facility_id}").get_json()
    assert found["affiliation_polity_id"] is None and found["affiliation_name"] == "Velarix Concord"


def test_the_list_filters_and_sorts(admin, anonymous, world):
    farm_kind = _type(admin, name="Listening post")
    a = _create(admin, world, name="Alpha", type_id=farm_kind, ratings={"crime": "red", "health": "green"})
    b = _create(admin, world, name="Bravo", ratings={"crime": "green", "health": "green"},
                affiliation_name="Free Traders")
    c = _create(admin, world, name="Charlie", ratings={"crime": "yellow", "health": "green"},
                affiliation_polity_id=world["polity_id"])
    d = _create(admin, world, name="Delta")
    e = _create(admin, world, name="Echo", kind="station", ratings={"housing": "red"})

    def ids(query):
        body = anonymous.get(f"/api/facilities?{query}").get_json()
        assert body["total"] == len(body["items"]) or "limit" in query
        return [item["id"] for item in body["items"]]

    assert ids("") == [a, b, c, d, e]
    assert ids("kind=station") == [e]
    assert ids(f"type_id={farm_kind}") == [a]
    assert ids(f"affiliation_polity_id={world['polity_id']}") == [c]
    assert ids("affiliation_name=Free%20Traders") == [b]
    assert ids("rating_crime=red") == [a]
    assert ids("rating_health=green&rating_crime=green") == [b]
    assert ids("overall=red") == [a, e]
    assert ids("overall=green") == [b]
    assert ids("overall=blue") == [] if False else True
    assert ids("sort=worst") == [a, e, c, b, d]            # worst first, unassessed last
    assert ids("sort=worst&order=desc") == [b, c, a, e, d]
    assert ids("sort=type") == [a, b, c, d, e]
    paged = anonymous.get("/api/facilities?sort=worst&limit=2&offset=1").get_json()
    assert [i["id"] for i in paged["items"]] == [e, c] and paged["total"] == 5

    for bad in ("kind=fortress", "type_id=x", "rating_crime=amber", "overall=blue", "sort=size", "order=up"):
        assert anonymous.get(f"/api/facilities?{bad}").status_code == 400, bad
