# tests/test_facility_types.py

"""
Facility types and their typed fields (POP.10): the rules in
`planetgen.population.facility_types` and the `/api/facility-types` routes,
against a real throwaway database. Anyone reads; only an admin writes.
"""

import pytest

from planetgen.admin import auth as adminAuth
from planetgen.api.config import Config
from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.population import facility_types as types
from planetgen.population.facility_icons import ICON_NAMES
from planetgen.web.app import create_app


@pytest.fixture
def anonymous(mysql_config):
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


@pytest.fixture
def admin(mysql_config):
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    client = app.test_client()
    store.get_connection(mysql_config).close()
    _username, password = adminAuth.bootstrap_control_schema(mysql_config)
    assert client.post("/api/auth/login", json={"username": "admin", "password": password}).status_code == 200
    assert client.post("/api/auth/change-credentials", json={
        "current_password": password, "new_username": "boss", "new_password": "a-long-new-password-1",
    }).status_code == 200
    return client


def _facility_of_type(mysql_config, type_id):
    """An orbital outpost around the star of a new system, typed `type_id`; returns its id."""
    sector = SpaceSector("Typed", edge_ly=10.0)
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    config.PLANETS = False
    config.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=config), position=(1.0, 1.0, 1.0), system_config=config)
    store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            star_id = conn.execute("SELECT id FROM stars LIMIT 1").fetchone()["id"]
            facility_id = store.add_facility(conn, "Typed post", "outpost", "orbital", "star", star_id)
            conn.execute("UPDATE facilities SET type_id = ? WHERE id = ?", (type_id, facility_id))
    finally:
        conn.close()
    return facility_id


def _put_value(mysql_config, facility_id, field_id, column, value):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute(f"INSERT INTO facility_field_values (facility_id, field_id, {column}) VALUES (?, ?, ?)",
                         (facility_id, field_id, value))
    finally:
        conn.close()


FARM = {
    "name": "Agricultural colony", "kind": "colony", "icon": "tree", "description": "Grows food.", "sort_order": 2,
    "fields": [
        {"key": "hectares", "label": "Farmland", "datatype": "decimal", "unit": "ha", "min_value": 0, "required": True},
        {"key": "crop", "label": "Main crop", "datatype": "choice", "choices": ["Grain", "Fruit", "Fungus"],
         "default_text": "Grain"},
        {"key": "irrigated", "label": "Irrigated", "datatype": "bool"},
    ],
}


def test_the_whole_life_of_a_type_over_the_api(admin, anonymous):
    created = admin.post("/api/facility-types", json=FARM)
    assert created.status_code == 201, created.get_json()
    type_id = created.get_json()["id"]

    # Anyone reads it, fields in order, with the choice options as a list.
    read = anonymous.get(f"/api/facility-types/{type_id}")
    assert read.status_code == 200
    found = read.get_json()
    assert (found["name"], found["kind"], found["icon"], found["active"], found["facilities"]) == (
        "Agricultural colony", "colony", "tree", True, 0)
    assert [f["key"] for f in found["fields"]] == ["hectares", "crop", "irrigated"]
    crop = next(f for f in found["fields"] if f["key"] == "crop")
    assert crop["choices"] == ["Grain", "Fruit", "Fungus"] and crop["default_text"] == "Grain"
    hectares = next(f for f in found["fields"] if f["key"] == "hectares")
    assert hectares["required"] is True and hectares["unit"] == "ha" and hectares["min_value"] == 0

    listing = anonymous.get("/api/facility-types").get_json()
    assert [t["id"] for t in listing["items"]] == [type_id]
    assert listing["icons"] == list(ICON_NAMES)
    assert anonymous.get("/api/facility-types?kind=starbase").get_json()["items"] == []
    assert anonymous.get("/api/facility-types?kind=bogus").status_code == 400

    # Changing it.
    changed = admin.patch(f"/api/facility-types/{type_id}", json={"name": "Farm world", "active": False,
                                                                  "icon": "leaf"})
    assert changed.status_code == 200
    assert (changed.get_json()["name"], changed.get_json()["active"]) == ("Farm world", False)
    assert anonymous.get("/api/facility-types?active=1").get_json()["items"] == []
    assert len(anonymous.get("/api/facility-types").get_json()["items"]) == 1

    # A field added, changed and deleted.
    added = admin.post(f"/api/facility-types/{type_id}/fields", json={"key": "yield", "label": "Yield",
                                                                       "datatype": "int", "max_value": 100})
    assert added.status_code == 201
    field_id = added.get_json()["id"]
    edited = admin.patch(f"/api/facility-types/{type_id}/fields/{field_id}", json={"label": "Yield per ha",
                                                                                   "max_value": 200})
    assert edited.status_code == 200
    field = next(f for f in edited.get_json()["fields"] if f["id"] == field_id)
    assert (field["label"], field["max_value"], field["key"]) == ("Yield per ha", 200, "yield")
    assert admin.delete(f"/api/facility-types/{type_id}/fields/{field_id}").status_code == 200
    assert admin.delete(f"/api/facility-types/{type_id}/fields/{field_id}").status_code == 404

    assert admin.delete(f"/api/facility-types/{type_id}").status_code == 200
    assert anonymous.get(f"/api/facility-types/{type_id}").status_code == 404
    assert anonymous.get("/api/facility-types").get_json()["items"] == []


def test_only_an_admin_writes(anonymous):
    assert anonymous.post("/api/facility-types", json=FARM).status_code == 401
    assert anonymous.patch("/api/facility-types/1", json={"name": "X"}).status_code == 401
    assert anonymous.delete("/api/facility-types/1").status_code == 401
    assert anonymous.post("/api/facility-types/1/fields", json=FARM["fields"][0]).status_code == 401
    assert anonymous.patch("/api/facility-types/1/fields/1", json={"label": "X"}).status_code == 401
    assert anonymous.delete("/api/facility-types/1/fields/1").status_code == 401


@pytest.mark.parametrize("change, fragment", [
    ({"name": "  "}, "name"),
    ({"kind": "fortress"}, "kind"),
    ({"icon": "not-an-icon"}, "icon"),
    ({"sort_order": "1"}, "sort_order"),
    ({"bogus": 1}, "bogus"),
    ({"fields": [{"key": "Bad Key", "label": "x", "datatype": "text"}]}, "key"),
    ({"fields": [{"key": "a", "label": "x", "datatype": "float"}]}, "datatype"),
    ({"fields": [{"key": "a", "label": "x", "datatype": "choice", "choices": ["only"]}]}, "2 to"),
    ({"fields": [{"key": "a", "label": "x", "datatype": "choice", "choices": ["same", "same"]}]}, "differ"),
    ({"fields": [{"key": "a", "label": "x", "datatype": "text", "choices": ["a", "b"]}]}, "choices"),
    ({"fields": [{"key": "a", "label": "x", "datatype": "int", "min_value": 5, "max_value": 1}]}, "min_value"),
    ({"fields": [{"key": "a", "label": "x", "datatype": "text", "unit": "kg"}]}, "unit"),
    ({"fields": [{"key": "a", "label": "x", "datatype": "int", "min_value": 0, "default_text": "-1"}]}, "at least"),
    ({"fields": [{"key": "a", "label": "x", "datatype": "int"}, {"key": "a", "label": "y", "datatype": "int"}]},
     "own 'key'"),
    ({"fields": [{"key": f"f{n}", "label": "x", "datatype": "text"} for n in range(types.MAX_FIELDS + 1)]},
     "at most"),
])
def test_a_bad_type_is_a_400_naming_the_problem(admin, change, fragment):
    body = {"name": "Probe", "kind": "outpost", **change}
    response = admin.post("/api/facility-types", json=body)
    assert response.status_code == 400, response.get_json()
    assert fragment in response.get_json()["error"]
    assert admin.get("/api/facility-types").get_json()["items"] == []


def test_a_name_is_unique_and_a_type_can_be_empty(admin):
    assert admin.post("/api/facility-types", json={"name": "Post", "kind": "outpost"}).status_code == 201
    clash = admin.post("/api/facility-types", json={"name": "Post", "kind": "station"})
    assert clash.status_code == 409
    other = admin.post("/api/facility-types", json={"name": "Dock", "kind": "station"}).get_json()["id"]
    assert admin.patch(f"/api/facility-types/{other}", json={"name": "Post"}).status_code == 409
    assert admin.patch(f"/api/facility-types/{other}", json={}).status_code == 400
    assert admin.patch("/api/facility-types/99999", json={"name": "Nope"}).status_code == 404


def test_a_type_in_use_cannot_be_deleted_or_change_its_kind(admin, mysql_config):
    type_id = admin.post("/api/facility-types", json=FARM).get_json()["id"]
    _facility_of_type(mysql_config, type_id)

    assert admin.get(f"/api/facility-types/{type_id}").get_json()["facilities"] == 1
    refused = admin.delete(f"/api/facility-types/{type_id}")
    assert refused.status_code == 409 and "inactive" in refused.get_json()["error"]
    assert admin.patch(f"/api/facility-types/{type_id}", json={"kind": "outpost"}).status_code == 409
    assert admin.patch(f"/api/facility-types/{type_id}", json={"kind": "colony", "active": False}).status_code == 200


def test_a_field_with_values_keeps_the_values_valid(admin, mysql_config):
    type_id = admin.post("/api/facility-types", json=FARM).get_json()["id"]
    fields = {f["key"]: f["id"] for f in admin.get(f"/api/facility-types/{type_id}").get_json()["fields"]}
    facility_id = _facility_of_type(mysql_config, type_id)
    _put_value(mysql_config, facility_id, fields["hectares"], "value_dec", 40.0)
    _put_value(mysql_config, facility_id, fields["crop"], "value_text", "Fruit")
    base = f"/api/facility-types/{type_id}/fields"

    assert admin.delete(f"{base}/{fields['hectares']}").status_code == 409
    assert admin.patch(f"{base}/{fields['hectares']}", json={"datatype": "int"}).status_code == 409
    assert admin.patch(f"{base}/{fields['hectares']}", json={"max_value": 10}).status_code == 409
    assert admin.patch(f"{base}/{fields['hectares']}", json={"min_value": 50}).status_code == 409
    assert admin.patch(f"{base}/{fields['hectares']}", json={"max_value": 500, "label": "Land"}).status_code == 200
    assert admin.patch(f"{base}/{fields['crop']}", json={"choices": ["Grain", "Fungus"]}).status_code == 409
    assert admin.patch(f"{base}/{fields['crop']}", json={"choices": ["Fruit", "Grain", "Tea"]}).status_code == 200

    # A field added to a type with facilities: a required one needs a default, and the default is filled in.
    needs_default = admin.post(base, json={"key": "crew", "label": "Crew", "datatype": "int", "required": True})
    assert needs_default.status_code == 409
    added = admin.post(base, json={"key": "crew", "label": "Crew", "datatype": "int", "required": True,
                                   "default_text": "12"})
    assert added.status_code == 201
    conn = store.get_connection(mysql_config)
    try:
        row = conn.execute("SELECT value_int FROM facility_field_values WHERE facility_id = ? AND field_id = ?",
                           (facility_id, added.get_json()["id"])).fetchone()
    finally:
        conn.close()
    assert row["value_int"] == 12


def test_a_field_belongs_to_its_type_and_the_key_is_fixed(admin):
    first = admin.post("/api/facility-types", json=FARM).get_json()["id"]
    second = admin.post("/api/facility-types", json={"name": "Other", "kind": "outpost"}).get_json()["id"]
    field_id = admin.get(f"/api/facility-types/{first}").get_json()["fields"][0]["id"]

    assert admin.patch(f"/api/facility-types/{second}/fields/{field_id}", json={"label": "X"}).status_code == 404
    assert admin.delete(f"/api/facility-types/{second}/fields/{field_id}").status_code == 404
    assert admin.patch(f"/api/facility-types/{first}/fields/{field_id}", json={"key": "renamed"}).status_code == 400
    assert admin.post("/api/facility-types/99999/fields", json=FARM["fields"][0]).status_code == 404
    dup = admin.post(f"/api/facility-types/{first}/fields", json={"key": "crop", "label": "Again", "datatype": "text"})
    assert dup.status_code == 409


@pytest.mark.parametrize("field, raw, expected", [
    ({"key": "n", "datatype": "int", "min_value": 0, "max_value": 10}, 7, 7),
    ({"key": "n", "datatype": "int"}, "12", 12),
    ({"key": "n", "datatype": "int"}, 3.0, 3),
    ({"key": "n", "datatype": "decimal"}, 2, 2.0),
    ({"key": "n", "datatype": "decimal"}, "2.5", 2.5),
    ({"key": "n", "datatype": "bool"}, True, True),
    ({"key": "n", "datatype": "bool"}, "false", False),
    ({"key": "n", "datatype": "text"}, "words", "words"),
    ({"key": "n", "datatype": "choice", "choices": ["a", "b"]}, "b", "b"),
    ({"key": "n", "datatype": "choice", "choices": "a\nb"}, "a", "a"),
])
def test_parse_value_accepts_what_fits(field, raw, expected):
    value = types.parse_value(field, raw)
    assert value == expected and type(value) is type(expected)


@pytest.mark.parametrize("field, raw", [
    ({"key": "n", "datatype": "int"}, 1.5),
    ({"key": "n", "datatype": "int"}, True),
    ({"key": "n", "datatype": "int"}, "many"),
    ({"key": "n", "datatype": "int", "max_value": 5}, 6),
    ({"key": "n", "datatype": "decimal", "min_value": 0}, -0.1),
    ({"key": "n", "datatype": "decimal"}, float("inf")),
    ({"key": "n", "datatype": "decimal"}, None),
    ({"key": "n", "datatype": "bool"}, 1),
    ({"key": "n", "datatype": "text"}, 5),
    ({"key": "n", "datatype": "text"}, "x" * 256),
    ({"key": "n", "datatype": "choice", "choices": ["a", "b"]}, "c"),
])
def test_parse_value_refuses_what_does_not(field, raw):
    with pytest.raises(types.FacilityTypeError) as caught:
        types.parse_value(field, raw)
    assert caught.value.status == 400 and "'n'" in str(caught.value)


def test_every_allowed_icon_is_a_vendored_file_and_the_licence_is_kept():
    import pathlib
    folder = pathlib.Path(__file__).resolve().parents[1] / "html" / "static" / "vendor" / "shoelace" / "assets" / "icons"
    assert len(set(ICON_NAMES)) == len(ICON_NAMES)
    assert [name for name in ICON_NAMES if not (folder / f"{name}.svg").is_file()] == []
    assert "MIT License" in (folder / "LICENSE-bootstrap-icons.txt").read_text(encoding="utf-8")
    # The chips of a rating use these three shapes and five glyphs too.
    from planetgen.web.lib import facility_display
    for icon in (*facility_display.COLOR_ICONS.values(), *facility_display.DIMENSION_ICONS.values()):
        assert icon in ICON_NAMES
