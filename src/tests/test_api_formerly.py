# tests/test_api_formerly.py

"""GEN.1140: the designations an object gave up (`designations`, GEN.1138) are shown. The JSON of an object gains `former_designations`
(newest first), its page a "Formerly" line, a search finds an object by a name it used to have and says so, and `GET /api/lookup` answers what
an Object ID or a name is now."""

import pytest

from planetgen.db import designation_rows as rows
from planetgen.db import store
from planetgen.galaxy import object_uid as ou
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.system import StarSystem
from tests.test_api import admin_client, client, default_admin_client, first_admin_password  # noqa: F401

pytestmark = pytest.mark.fresh_database

NEW_NAME = "Zandrel Prime"


@pytest.fixture
def galaxy(mysql_config, admin_client):  # noqa: F811
    sector = SpaceSector("Formerly Sector", edge_ly=20.0)
    config = SystemConfig()
    config.PLANETS = True
    config.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=config), position=(1.0, 1.0, 1.0), system_config=config)
    sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    store.save_sector(sector, config=mysql_config)
    return mysql_config, admin_client


def _first(config, table):
    conn = store.get_connection(config)
    try:
        return conn.execute(f"SELECT id, uid, name FROM {table} ORDER BY id LIMIT 1").fetchone()
    finally:
        conn.close()


def _printed(row):
    return ou.format_id(ou.from_bytes(row["uid"]))


def _rename(config, table, row, new_name=NEW_NAME):
    conn = store.get_connection(config)
    try:
        with conn:
            rows.set_name(conn, table, row["id"], new_name, "rename")
    finally:
        conn.close()


def _rename_system(config, system, new_name=NEW_NAME):
    """A person's rename of a system: the stars named after it follow (`store.rename_star_system`), as in the app."""
    conn = store.get_connection(config)
    try:
        with conn:
            store.rename_star_system(conn, system["id"], new_name)
    finally:
        conn.close()


def test_the_json_of_a_renamed_system_lists_its_old_name_newest_first(galaxy):
    config, client = galaxy
    system = _first(config, "star_systems")
    _rename(config, "star_systems", system)
    _rename(config, "star_systems", system, "Zandrel Secundus")
    for url in (f"/api/systems/{_printed(system)}", f"/api/objects/system:{_printed(system)}"):
        body = client.get(url).get_json()
        assert [(f["designation"], f["kind"], f["reason"]) for f in body["former_designations"]] == [
            (NEW_NAME, "name", "rename"), (system["name"], "name", "rename")], url
        assert all(isinstance(f["retired_epoch"], int) and len(f["retired_on"]) == 10 for f in body["former_designations"])


def test_an_object_that_never_changed_has_an_empty_list_on_every_kind_of_answer(galaxy):
    config, client = galaxy
    system, planet, rogue = _first(config, "star_systems"), _first(config, "planets"), _first(config, "rogue_planets")
    sector_id = client.get("/api/sectors").get_json()["items"][0]["id"]
    for url in (f"/api/systems/{_printed(system)}", f"/api/objects/planet:{_printed(planet)}", f"/api/objects/star:{_printed(_first(config, 'stars'))}",
                f"/api/phenomena/rogue_planet/{_printed(rogue)}", f"/api/sectors/{sector_id}"):
        response = client.get(url)
        assert response.status_code == 200 and response.get_json()["former_designations"] == [], url


def test_a_renamed_sector_and_phenomenon_carry_the_field_too(galaxy):
    config, client = galaxy
    sector, rogue = _first(config, "sectors"), _first(config, "rogue_planets")
    _rename(config, "sectors", sector)
    _rename(config, "rogue_planets", rogue)
    sector_id = client.get("/api/sectors").get_json()["items"][0]["id"]
    assert [f["designation"] for f in client.get(f"/api/sectors/{sector_id}").get_json()["former_designations"]] == [sector["name"]]
    assert [f["designation"] for f in client.get(f"/api/phenomena/rogue_planet/{_printed(rogue)}").get_json()["former_designations"]] == [rogue["name"]]


def test_the_pages_say_formerly_with_the_old_name_the_reason_and_the_epoch(galaxy):
    config, client = galaxy
    system, rogue, planet = _first(config, "star_systems"), _first(config, "rogue_planets"), _first(config, "planets")
    _rename(config, "star_systems", system)
    _rename(config, "rogue_planets", rogue)
    for url, old in ((f"/system/{_printed(system)}", system["name"]), (f"/phenomenon/rogue_planet/{_printed(rogue)}", rogue["name"])):
        html = client.get(url).get_data(as_text=True)
        assert 'class="formerly"' in html and old in html and "renamed" in html and "epoch " in html, url
    plain = client.get(f"/system/{_printed(_first(config, 'star_systems'))}").get_data(as_text=True)
    assert plain.count('class="formerly"') == 1
    untouched = client.get(f"/object/planet:{_printed(planet)}")
    assert untouched.status_code == 200 and 'class="formerly"' not in untouched.get_data(as_text=True)


def test_a_planet_page_whose_system_was_renamed_lists_the_name_it_had(galaxy):
    config, client = galaxy
    planet = _first(config, "planets")
    conn = store.get_connection(config)
    try:
        system = conn.execute("SELECT s.id FROM star_systems s JOIN planets p ON p.star_system_id = s.id WHERE p.id = ?", (planet["id"],)).fetchone()
        with conn:
            store.rename_star_system(conn, system["id"], NEW_NAME)
    finally:
        conn.close()
    html = client.get(f"/object/planet:{_printed(planet)}").get_data(as_text=True)
    assert 'class="formerly"' in html and planet["name"] in html


def test_a_search_for_a_past_name_finds_the_object_and_says_it_was_called_that(galaxy):
    config, client = galaxy
    system = _first(config, "star_systems")
    _rename(config, "star_systems", system)
    old_word = system["name"].split()[0]
    by_old = client.get(f"/api/search?system_q={old_word}").get_json()["results"]["systems"]
    [row] = [r for r in by_old["rows"] if r["id"] == _printed(system)]
    assert row["name"] == NEW_NAME and system["name"] in row["formerly"]
    by_new = client.get(f"/api/search?system_q={NEW_NAME.split()[0]}").get_json()["results"]["systems"]
    [row] = [r for r in by_new["rows"] if r["id"] == _printed(system)]
    assert "formerly" not in row


def test_the_search_page_marks_a_row_found_by_a_past_name_and_points_at_the_object(galaxy):
    config, client = galaxy
    system = _first(config, "star_systems")
    _rename_system(config, system)
    html = client.get("/search", query_string={"q": system["name"]}).get_data(as_text=True)
    assert "(formerly " in html and NEW_NAME in html
    assert "is now" in html


def test_lookup_answers_what_an_id_or_a_name_is_now(galaxy):
    config, client = galaxy
    system = _first(config, "star_systems")
    conn = store.get_connection(config)
    try:
        with conn:
            store.rename_star_system(conn, system["id"], NEW_NAME)
            old_id = rows.renumber(conn, "star_systems", system["id"], ou.from_bytes(system["uid"]) + (5000 << 12), "admin renumber")
    finally:
        conn.close()
    new_id = ou.format_id(ou.from_bytes(system["uid"]) + (5000 << 12))

    now = client.get("/api/lookup", query_string={"text": new_id})
    assert now.status_code == 200 and (now.get_json()["status"], now.get_json()["ref"]) == ("current", f"system:{new_id}")
    moved = client.get("/api/lookup", query_string={"text": ou.format_id(old_id).lower()}).get_json()
    assert (moved["status"], moved["id"], moved["name"], moved["retired"]) == ("redirect", new_id, NEW_NAME, ou.format_id(old_id).lower())
    renamed = client.get("/api/lookup", query_string={"text": system["name"]}).get_json()
    assert renamed["status"] == "redirect" and renamed["name"].startswith(NEW_NAME) and renamed["retired"] == system["name"]
    current = client.get("/api/lookup", query_string={"text": NEW_NAME}).get_json()
    assert current["status"] == "current" and current["name"] == NEW_NAME
    assert system["name"] in [f["designation"] for f in current["former_designations"]]


def test_lookup_of_a_deleted_one_is_a_410_of_one_nobody_had_a_404_and_of_nothing_a_400(galaxy):
    config, client = galaxy
    system = _first(config, "star_systems")
    _rename_system(config, system)
    conn = store.get_connection(config)
    try:
        with conn:
            conn.execute("UPDATE designations SET gone_epoch = 9 WHERE designation_kind = 'name'")
    finally:
        conn.close()
    gone = client.get("/api/lookup", query_string={"text": system["name"]})
    assert gone.status_code == 410 and gone.get_json()["status"] == "gone"
    assert client.get("/api/lookup", query_string={"text": "Nobody Ever Had This"}).status_code == 404
    assert client.get("/api/lookup", query_string={"text": ou.format_id(ou.from_bytes(system["uid"]) + (9000 << 12))}).status_code == 404
    assert client.get("/api/lookup").status_code == 400


def test_lookup_of_a_name_two_objects_once_had_answers_with_one_and_lists_the_other(galaxy):
    config, client = galaxy
    system, star = _first(config, "star_systems"), _first(config, "stars")
    assert star["name"] == system["name"], "a single star shares its system's name"
    _rename(config, "star_systems", system, "Zandrel Prime")
    _rename(config, "stars", star, "Zandrel A")
    answer = client.get("/api/lookup", query_string={"text": system["name"]}).get_json()
    assert answer["status"] == "redirect" and answer["name"] in ("Zandrel Prime", "Zandrel A")
    [other] = answer["also"]
    assert other["ref"] != answer["ref"] and {answer["name"], other["name"]} == {"Zandrel Prime", "Zandrel A"}
