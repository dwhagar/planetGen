# tests/test_api_thin_routes.py

"""
TEST.49: the thin API routes' edges -- unknown ids, an empty galaxy,
paging/limit parameters and ids that name the wrong kind of row -- for
`/api/galaxy/sectors`, `/api/galaxy/shape`, `/api/galaxy/phenomena`,
`/api/galaxy/bright-stars`, the star/planet/moon rename PATCHes,
`/api/facilities` (POST and DELETE; there is no PATCH),
`/api/admin/login-failures`, `/api/population`, `/api/species/<id>` and
`/api/systems/<id>/owner`, plus what deleting a sector does to the
facilities and wiki links hanging off it.

Every refusal must be a clean `{"error": ...}` JSON 4xx (400 for bad
parameters or bodies, 404 for an unknown id, 405 for a method a route
doesn't take), never a 500. The `client` fixture runs with `TESTING`
set, so an unhandled exception would surface here as a raised error, not
a quiet 500.

Reuses `test_api.py`'s fixtures; every test is skipped without a
reachable MySQL test server.
"""

import pytest

from planetgen.db import store

# Imported fixtures (see test_bughunt_api_gaps.py for why this works).
from tests.test_api import (  # noqa: F401
    _place_sector, _save_wide_binary_with_moons, admin_client, client, default_admin_client,
    first_admin_password, seeded_sector,
)
from tests.publicids import pid, pids

HUGE_ID = 10 ** 25
"""int: Past `BIGINT UNSIGNED`'s range; the URL's `<int:...>` converter
still accepts it, so the query itself must cope."""


def _schema(mysql_config):
    """Lays the content schema down on the empty throwaway database, so a
    read is a real "nothing there", not "no table yet"."""
    store.get_connection(mysql_config).close()


def _assert_json_error(response, status):
    assert response.status_code == status, (response.status_code, response.get_data(as_text=True)[:300])
    body = response.get_json()
    assert isinstance(body, dict) and isinstance(body.get("error"), str) and body["error"], body


def _one(mysql_config, sql, params=()):
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute(sql, params).fetchone()
    finally:
        conn.close()


# --- Galaxy reads on an empty galaxy --------------------------------------------------

def test_galaxy_reads_on_an_empty_galaxy(client, mysql_config):
    _schema(mysql_config)
    assert client.get("/api/galaxy/sectors").get_json() == {"items": []}
    assert client.get("/api/galaxy/phenomena").get_json() == {"items": []}
    shape = client.get("/api/galaxy/shape").get_json()
    assert shape["shape"] is None and shape["bright_stars"]["scattered"] is False
    assert client.get("/api/galaxy/bright-stars?ring=0&layer=0&slot=0").get_json() == {"items": []}


def test_galaxy_listings_ignore_paging_parameters(client, mysql_config):
    """The galaxy listings aren't paginated (see their docstrings), so
    `limit`/`offset` -- even nonsense ones that a paginated route would
    refuse -- are ignored rather than 400 or 500."""
    _place_sector(mysql_config, "Paged One", (100.0, 0.0, 0.0), address=(30, 0, 0))
    _place_sector(mysql_config, "Paged Two", (-100.0, 0.0, 0.0), address=(30, 0, 5))
    for path in ("/api/galaxy/sectors", "/api/galaxy/phenomena", "/api/galaxy/shape"):
        plain = client.get(path).get_json()
        for query in ("?limit=1", "?limit=-5&offset=abc", "?limit=99999999999999999999", "?offset=1"):
            response = client.get(path + query)
            assert response.status_code == 200, (path, query)
            assert response.get_json() == plain, (path, query)
    assert {item["name"] for item in client.get("/api/galaxy/sectors?limit=1").get_json()["items"]} == {
        "Paged One", "Paged Two"}


@pytest.mark.parametrize("query", [
    "", "ring=0", "ring=0&layer=0", "ring=1.5&layer=0&slot=0", "ring=&layer=0&slot=0",
    "ring=nan&layer=0&slot=0", "ring=inf&layer=0&slot=0", "ring=0&layer=0&slot=%E2%91%A0",
    "ring=0x10&layer=0&slot=0",
])
def test_bright_stars_bad_address_is_400(client, mysql_config, query):
    _schema(mysql_config)
    _assert_json_error(client.get(f"/api/galaxy/bright-stars?{query}"), 400)


@pytest.mark.parametrize("query", [
    "ring=-1&layer=0&slot=0", "ring=0&layer=-99&slot=0", f"ring={HUGE_ID}&layer=0&slot=0",
    f"ring=0&layer=0&slot=-{HUGE_ID}", "ring=+3&layer=0&slot=0&all=1", "ring=3&layer=0&slot=0&all=maybe",
])
def test_bright_stars_out_of_range_address_is_an_empty_list(client, mysql_config, query):
    """Integers no cell has (negative, past BIGINT) are still integers:
    an empty answer, not a 400 or a database overflow."""
    _schema(mysql_config)
    response = client.get(f"/api/galaxy/bright-stars?{query}")
    assert response.status_code == 200
    assert response.get_json() == {"items": []}


# --- Star/planet/moon rename PATCH ----------------------------------------------------

@pytest.mark.parametrize("kind", ["stars", "planets", "moons"])
def test_rename_unknown_ids_are_404(admin_client, kind):
    for body_id in (0, HUGE_ID):
        _assert_json_error(admin_client.patch(f"/api/{kind}/{body_id}", json={"name": "Nowhere"}), 404)
    # A negative or non-numeric id never matches the route at all.
    for raw in ("-1", "abc"):
        response = admin_client.patch(f"/api/{kind}/{raw}", json={"name": "Nowhere"})
        assert response.status_code == 404
        assert response.get_json() == {"error": "not found"}


@pytest.mark.parametrize("kind", ["stars", "planets", "moons"])
def test_rename_bad_bodies_are_400(admin_client, kind):
    """Checked before the row is looked up, so an id that exists or not
    makes no difference."""
    for kwargs in ({"data": "name=x"}, {"json": ["name"]}, {"json": {}}, {"json": {"name": None}},
                   {"data": "{not json", "content_type": "application/json"}):
        _assert_json_error(admin_client.patch(f"/api/{kind}/1", **kwargs), 400)


def test_rename_with_the_wrong_kind_of_id(admin_client, mysql_config):
    """A moon's id sent to the planet route renames whichever planet has
    that id, if any -- ids are per table -- and never the moon; with no
    such planet it's a clean 404."""
    system_id = _save_wide_binary_with_moons(mysql_config)
    moon = _one(mysql_config, "SELECT id, name FROM moons WHERE star_system_id = ? ORDER BY id DESC LIMIT 1",
                (system_id,))
    top_planet = _one(mysql_config, "SELECT MAX(id) AS top FROM planets")["top"]
    top_star = _one(mysql_config, "SELECT MAX(id) AS top FROM stars")["top"]

    wrong_id = max(moon["id"], top_planet, top_star) + 1000
    _assert_json_error(admin_client.patch(f"/api/planets/{wrong_id}", json={"name": "Misfiled"}), 404)
    _assert_json_error(admin_client.patch(f"/api/stars/{wrong_id}", json={"name": "Misfiled"}), 404)
    if moon["id"] > top_planet:
        _assert_json_error(admin_client.patch(f"/api/planets/{moon['id']}", json={"name": "Misfiled"}), 404)
    assert _one(mysql_config, "SELECT name FROM moons WHERE id = ?", (moon["id"],))["name"] == moon["name"]


# --- Facilities -----------------------------------------------------------------------

def _star_and_moon(mysql_config):
    system_id = _save_wide_binary_with_moons(mysql_config)
    star_id = _one(mysql_config, "SELECT id FROM stars WHERE star_system_id = ? ORDER BY id LIMIT 1",
                   (system_id,))["id"]
    moon_id = _one(mysql_config, "SELECT id FROM moons WHERE star_system_id = ? ORDER BY id LIMIT 1",
                   (system_id,))["id"]
    return system_id, star_id, moon_id


def _facility(**overrides):
    body = {"name": "Relay", "kind": "station", "placement": "orbital", "host_type": "star", "host_id": 1}
    body.update(overrides)
    return {key: value for key, value in body.items() if value is not None}


@pytest.mark.parametrize("body", [
    _facility(name=None),
    _facility(host_id=None),
    _facility(name="   "),
    _facility(name="x" * 10000),
    _facility(kind="castle"),
    _facility(placement="underground"),
    _facility(host_type="comet"),
    _facility(host_id=0),
    _facility(host_id=-3),
    _facility(host_id="1"),
    _facility(host_id=True),
    _facility(host_id=1.0),
    _facility(distance_km=0),
    _facility(distance_km=-5),
    _facility(distance_km=True),
    _facility(phase_deg="north"),
    _facility(offset_ly=[0, 0]),
    _facility(offset_ly=[0, 0, "x"]),
    _facility(offset_ly=[0, 0, True]),
    _facility(description=5),
    _facility(owner="me"),
])
def test_create_facility_bad_bodies_are_400(admin_client, body):
    _assert_json_error(admin_client.post("/api/facilities", json=body), 400)


def test_create_facility_non_object_bodies_are_400(admin_client):
    for kwargs in ({"json": []}, {"json": "x"}, {"data": "name=x"}, {}):
        _assert_json_error(admin_client.post("/api/facilities", **kwargs), 400)


@pytest.mark.parametrize("host_type, placement, kind", [
    ("star", "orbital", "station"),
    ("planet", "orbital", "station"),
    ("moon", "terrestrial", "colony"),
    ("asteroid_belt", "asteroid", "outpost"),
    ("asteroid_field", "asteroid", "outpost"),
    ("space", "standalone", "station"),
])
def test_create_facility_unknown_host_is_404(admin_client, host_type, placement, kind):
    body = _facility(host_type=host_type, placement=placement, kind=kind, host_id=999999999)
    _assert_json_error(admin_client.post("/api/facilities", json=body), 404)


def test_create_facility_wrong_placement_for_its_host(admin_client, mysql_config, seeded_sector):
    """A real host with a placement, kind or position it can't take is a
    400 from the placement rules, and nothing is stored."""
    _config, sector_id, _system_ids = seeded_sector
    _system_id, star_id, moon_id = _star_and_moon(mysql_config)
    refused = [
        _facility(placement="terrestrial", kind="colony", host_type="star", host_id=star_id),
        _facility(placement="orbital", kind="colony", host_type="moon", host_id=moon_id),
        _facility(placement="standalone", host_type="star", host_id=star_id),
        # An orbit past the host's sphere of influence.
        _facility(host_type="moon", host_id=moon_id, distance_km=1e15),
        # A stand-alone position on an orbital facility, and an orbit on a surface one.
        _facility(host_type="star", host_id=star_id, offset_ly=[0, 0, 0]),
        _facility(placement="terrestrial", kind="outpost", host_type="moon", host_id=moon_id, phase_deg=10),
        # seeded_sector's sector was never placed in the galaxy.
        _facility(placement="standalone", host_type="space", host_id=sector_id),
    ]
    for body in refused:
        _assert_json_error(admin_client.post("/api/facilities", json=body), 400)
    assert _one(mysql_config, "SELECT COUNT(*) AS n FROM facilities")["n"] == 0


def test_standalone_facility_offset_outside_its_sector_is_400(admin_client, mysql_config):
    sector_id = _place_sector(mysql_config, "Dock Sector", (40.0, 0.0, 0.0), edge_ly=10.0, address=(12, 0, 0))
    for offset in ([5.01, 0, 0], [0, 0, -1e300], [0, 1e308, 0]):
        body = _facility(placement="standalone", host_type="space", host_id=sector_id, offset_ly=offset)
        _assert_json_error(admin_client.post("/api/facilities", json=body), 400)
    body = _facility(placement="standalone", host_type="space", host_id=sector_id, offset_ly=[4.9, 0, 0])
    assert admin_client.post("/api/facilities", json=body).status_code == 201


def test_facility_has_no_patch_and_delete_is_404_once_gone(admin_client, mysql_config):
    _system_id, star_id, _moon_id = _star_and_moon(mysql_config)
    response = admin_client.post("/api/facilities", json=_facility(host_id=star_id))
    assert response.status_code == 201
    facility_id = response.get_json()["id"]

    # There is no facility PATCH: a clean 405, and the facility is untouched.
    response = admin_client.patch(f"/api/facilities/{facility_id}", json={"name": "Renamed"})
    assert response.status_code == 405
    assert response.get_json() == {"error": "method not allowed"}
    assert admin_client.get(f"/api/facilities/{facility_id}").get_json()["name"] == "Relay"
    assert admin_client.delete("/api/facilities").status_code == 405

    assert admin_client.delete(f"/api/facilities/{facility_id}").status_code == 200
    _assert_json_error(admin_client.delete(f"/api/facilities/{facility_id}"), 404)
    _assert_json_error(admin_client.delete(f"/api/facilities/{HUGE_ID}"), 404)
    _assert_json_error(admin_client.get(f"/api/facilities/{facility_id}"), 404)


def test_facility_reads_for_unknown_ids(client, mysql_config):
    _schema(mysql_config)
    for path in ("/api/facilities/0", f"/api/facilities/{HUGE_ID}", "/api/systems/999999999/facilities",
                 "/api/sectors/999999999/facilities"):
        _assert_json_error(client.get(path), 404)


def test_facility_writes_need_an_admin(client):
    assert client.post("/api/facilities", json=_facility()).status_code == 401
    assert client.delete("/api/facilities/1").status_code == 401


def test_facility_writes_need_a_fresh_admin(default_admin_client):
    # `client` and `default_admin_client` are one test client, so the
    # signed-out checks live in their own test above.
    assert default_admin_client.post("/api/facilities", json=_facility()).status_code == 403
    assert default_admin_client.delete("/api/facilities/1").status_code == 403


# --- /api/admin/login-failures ----------------------------------------------------------

@pytest.mark.parametrize("limit", ["0", "-1", "201", "abc", "1.5", "1e2", str(HUGE_ID)])
def test_login_failures_bad_limit_is_400(admin_client, limit):
    _assert_json_error(admin_client.get(f"/api/admin/login-failures?limit={limit}"), 400)


def test_login_failures_limits_and_contents(admin_client):
    for limit in ("", "1", "200"):
        response = admin_client.get(f"/api/admin/login-failures?limit={limit}")
        assert response.status_code == 200
        assert isinstance(response.get_json()["items"], list)

    # A refused sign-in (on a separate client, so admin_client's session
    # is untouched) shows up, newest first, capped by the limit.
    other = admin_client.application.test_client()
    for _ in range(2):
        assert other.post("/api/auth/login", json={"username": "nobody-here", "password": "wrong"}).status_code == 401
    items = admin_client.get("/api/admin/login-failures?limit=1").get_json()["items"]
    assert len(items) == 1
    assert items[0]["username"] == "nobody-here"
    assert items[0]["created_at"].endswith("Z")


def test_login_failures_needs_an_admin(client):
    _assert_json_error(client.get("/api/admin/login-failures"), 401)


def test_login_failures_needs_a_fresh_admin(default_admin_client):
    _assert_json_error(default_admin_client.get("/api/admin/login-failures"), 403)


# --- Population reads -----------------------------------------------------------------

def test_population_on_an_empty_database(client, mysql_config):
    _schema(mysql_config)
    response = client.get("/api/population?limit=-1&offset=x")
    assert response.status_code == 200
    assert set(response.get_json().values()) == {False}


def test_species_unknown_ids(client, mysql_config):
    _schema(mysql_config)
    for species_id in (0, 1, 999999999, HUGE_ID):
        _assert_json_error(client.get(f"/api/species/{species_id}"), 404)
    for raw in ("-1", "abc", "1.5"):
        response = client.get(f"/api/species/{raw}")
        assert response.status_code == 404
        assert response.get_json() == {"error": "not found"}


def test_species_listing_paging_limits(client, mysql_config):
    _schema(mysql_config)
    for query in ("limit=-1", "limit=abc", "offset=-1", "offset=1.5"):
        _assert_json_error(client.get(f"/api/species?{query}"), 400)
    body = client.get("/api/species?limit=999999&offset=999999").get_json()
    assert body["items"] == [] and body["total"] == 0 and body["offset"] == 999999


def test_system_owner_unknown_unowned_and_deleted(admin_client, seeded_sector):
    _config, _sector_id, system_ids = seeded_sector
    for system_id in (0, 999999999, HUGE_ID):
        _assert_json_error(admin_client.get(f"/api/systems/{pid('system', system_id)}/owner"), 404)
    assert admin_client.get(f"/api/systems/{pid('system', system_ids[0])}/owner").get_json() == {"owner": None}
    assert admin_client.delete(f"/api/systems/{pid('system', system_ids[0])}").status_code == 200
    _assert_json_error(admin_client.get(f"/api/systems/{pid('system', system_ids[0])}/owner"), 404)


# --- Deleting a sector with facilities and wiki links ---------------------------------

def test_deleting_a_sector_with_facilities_and_wiki_links(admin_client, mysql_config):
    """Pins what `DELETE /api/sectors/<id>` does to what hangs off it:

    - a stand-alone facility parked in the sector goes with it
      (`facilities.sector_id` is `ON DELETE CASCADE`);
    - its systems stay, detached (`sector_id` set to NULL), and keep
      their own facilities and wiki links;
    - the sector's own `wiki_url` goes with its row;
    - the sector's reads are then 404s, and so is a second delete.
    """
    sector_id = _place_sector(mysql_config, "Doomed Sector", (60.0, 10.0, 0.0), address=(17, 0, 3))
    system_id = _one(mysql_config, "SELECT id FROM star_systems WHERE sector_id = ?", (sector_id,))["id"]
    star_id = _one(mysql_config, "SELECT id FROM stars WHERE star_system_id = ?", (system_id,))["id"]

    standalone = admin_client.post("/api/facilities", json=_facility(
        name="Drift Dock", placement="standalone", host_type="space", host_id=sector_id))
    assert standalone.status_code == 201, standalone.get_json()
    standalone_id = standalone.get_json()["id"]
    orbital = admin_client.post("/api/facilities", json=_facility(name="Sun Watch", host_id=star_id))
    assert orbital.status_code == 201, orbital.get_json()
    orbital_id = orbital.get_json()["id"]

    assert admin_client.patch(f"/api/sectors/{pid('sector', sector_id)}",
                              json={"wiki_url": "https://wiki.example.com/Doomed"}).status_code == 200
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("UPDATE star_systems SET wikijs_url = ? WHERE id = ?",
                         ("https://wiki.example.com/Doomed/System", system_id))
    finally:
        conn.close()
    assert [f["id"] for f in admin_client.get(f"/api/sectors/{pid('sector', sector_id)}/facilities").get_json()["items"]] == [
        standalone_id]

    assert admin_client.delete(f"/api/sectors/{pid('sector', sector_id)}").status_code == 200

    _assert_json_error(admin_client.get(f"/api/sectors/{pid('sector', sector_id)}"), 404)
    _assert_json_error(admin_client.get(f"/api/sectors/{pid('sector', sector_id)}/facilities"), 404)
    _assert_json_error(admin_client.delete(f"/api/sectors/{pid('sector', sector_id)}"), 404)
    _assert_json_error(admin_client.patch(f"/api/sectors/{pid('sector', sector_id)}", json={"wiki_url": None}), 404)

    _assert_json_error(admin_client.get(f"/api/facilities/{standalone_id}"), 404)
    assert _one(mysql_config, "SELECT COUNT(*) AS n FROM facilities WHERE sector_id = ?", (sector_id,))["n"] == 0

    system = admin_client.get(f"/api/systems/{pid('system', system_id)}").get_json()
    assert system["sector_id"] is None
    assert system["wikijs_url"] == "https://wiki.example.com/Doomed/System"
    assert admin_client.get(f"/api/facilities/{orbital_id}").get_json()["name"] == "Sun Watch"
    assert [f["id"] for f in admin_client.get(f"/api/systems/{pid('system', system_id)}/facilities").get_json()["items"]] == [
        orbital_id]
    assert admin_client.get(f"/api/systems/{pid('system', system_id)}/owner").get_json() == {"owner": None}
    assert all(s["id"] != sector_id for s in admin_client.get("/api/galaxy/sectors").get_json()["items"])


def test_galaxy_made_lists_sectors_created_in_a_window(client, mysql_config):
    """ADM.31: `/api/galaxy/made` returns the placed sectors made between two times."""
    import time
    _place_sector(mysql_config, "Made One", (100.0, 0.0, 0.0), address=(30, 0, 0))
    _place_sector(mysql_config, "Made Two", (-100.0, 0.0, 0.0), address=(30, 0, 5))
    now = time.time()
    body = client.get(f"/api/galaxy/made?since={now - 3600}&until={now + 3600}").get_json()
    assert body["total"] == 2
    assert {i["name"] for i in body["items"]} == {"Made One", "Made Two"}
    assert {(i["ring_index"], i["layer_index"], i["ring_slot_index"]) for i in body["items"]} == {(30, 0, 0), (30, 0, 5)}
    assert client.get(f"/api/galaxy/made?since={now + 3600}").get_json() == {"total": 0, "items": []}
    assert client.get("/api/galaxy/made?since=soon").status_code == 400
