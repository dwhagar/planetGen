# tests/test_belt_member_page.py

"""
GEN.1070: the page of one body an asteroid belt keeps (`/object/belt:<ID>/member/<rank>`) and the API route behind it (`GET /api/objects/<ref>/belt-members/<rank>`), walked on the
small generated database `test_web_a11y.py` builds. Skipped without the MySQL test server.
"""

import math

import pytest

from planetgen.db import store
from planetgen.galaxy import object_uid
from planetgen.web.maps import bodylook
from tests.publicids import pid
from tests.test_web_a11y import site_app, site_db  # noqa: F401 -- fixtures


@pytest.fixture(scope="module")
def client(site_app):
    return site_app.test_client()


@pytest.fixture(scope="module")
def belt(site_db):
    conn = store.get_connection(site_db)
    try:
        row = conn.execute("SELECT b.id, (SELECT COUNT(*) FROM belt_members m WHERE m.belt_id = b.id) AS members FROM asteroid_belts b ORDER BY members DESC, b.id LIMIT 1").fetchone()
        found = (f"belt:{pid('belt', row['id'], conn)}", row["members"]) if row else None
    finally:
        conn.close()
    if found is None or found[1] == 0:
        pytest.skip("the sample galaxy holds no asteroid belt with stored bodies")
    return found


def _rank(n):
    return f"{n:04X}"


def test_the_api_answers_a_member_the_belt_page_lists(client, belt):
    ref, members = belt
    listed = client.get(f"/api/objects/{ref}/belt-page").get_json()["bodies"]
    data = client.get(f"/api/objects/{ref}/belt-members/{_rank(1)}")
    assert data.status_code == 200
    body = data.get_json()
    assert set(body) == {"designation", "rank", "member", "orbit", "share", "belt", "previous", "next", "planet"}
    assert body["member"] == listed[0] and body["designation"] == listed[0]["designation"]
    assert body["previous"] is None and (body["next"] is None) == (members == 1)
    assert 0.0 < body["share"] <= 1.0


def test_the_orbit_follows_from_the_stored_elements(client, belt):
    ref, _ = belt
    orbit = client.get(f"/api/objects/{ref}/belt-members/{_rank(1)}").get_json()["orbit"]
    a, e = orbit["semi_major_axis"]["value"], orbit["eccentricity"]
    assert orbit["closest_distance"]["value"] == pytest.approx(a * (1 - e)) and orbit["farthest_distance"]["value"] == pytest.approx(a * (1 + e))
    assert orbit["period"]["unit"] == "yr" and orbit["period"]["value"] > 0 and math.isfinite(orbit["period"]["value"])


def test_the_page_shows_the_body_and_links_its_neighbours_and_its_belt(client, belt):
    ref, members = belt
    page = client.get(f"/object/{ref}/member/{_rank(1)}")
    assert page.status_code == 200
    text = page.get_data(as_text=True)
    assert 'id="member-body"' in text and 'id="member-orbit"' in text
    assert f'href="/object/{ref}"' in text and f'href="/object/{ref}/belt"' in text
    if members > 1:
        assert f'href="/object/{ref}/member/{_rank(2)}"' in text


def test_the_page_draws_its_body_from_the_look_of_a_belt_member(client, belt):
    """MAP.413: the Sensor view is a canvas with the member's look record (`bodylook.member_look`, its seed from the designation), painted by bodyview.js."""
    import html
    import json

    ref, _members = belt
    text = client.get(f"/object/{ref}/member/{_rank(1)}").get_data(as_text=True)
    assert 'id="member-example"' in text and "bodyview.js" in text and 'class="body-view-canvas"' in text
    look = json.loads(html.unescape(text.split("data-look='", 1)[1].split("'", 1)[0]))
    assert look["kind"] == "belt_member" and look["body"] == "solid" and look["color"].startswith("#") and look["radiusKm"] > 0
    assert look["seed"] == bodylook.seed_for("belt_member", client.get(f"/api/objects/{ref}/belt-members/{_rank(1)}").get_json()["designation"])


def test_the_belt_page_links_every_body_to_its_page(client, belt):
    ref, members = belt
    text = client.get(f"/object/{ref}/belt").get_data(as_text=True)
    assert text.count(f'href="/object/{ref}/member/') == members


def test_a_rank_the_belt_does_not_keep_or_does_not_write_so_is_not_found(client, belt):
    ref, members = belt
    for rank in (_rank(10 ** 7), "0000", "1", "00001", "ZZZZ"):
        assert client.get(f"/api/objects/{ref}/belt-members/{rank}").status_code == 404
        assert client.get(f"/object/{ref}/member/{rank}").status_code == 404
    assert client.get("/api/objects/planet:FE81000A2B-0000005-003/belt-members/0001").status_code == 404
    assert client.get("/object/planet:FE81000A2B-0000005-003/member/0001").status_code == 404



@pytest.fixture(scope="module")
def rich_belt(site_db):
    """`(ref, stored members, bodies at the class D size)` of a belt of the sample galaxy that holds more of them than it stores, or a skip."""
    from planetgen.describe import belt_info
    from planetgen.generation import belt_ondemand
    conn = store.get_connection(site_db)
    try:
        for row in conn.execute("SELECT b.id, b.star_system_id, b.star_id FROM asteroid_belts b ORDER BY b.id LIMIT 60").fetchall():
            system = store.load_star_system(conn, row["star_system_id"])
            belt, star = belt_info.belt_of(system, row["id"]), belt_info.host_star(system, row["star_id"])
            total = belt_ondemand.individual_count(belt_ondemand.settled_of(belt, star))
            if belt.members and total > len(belt.members) + 50:
                return f"belt:{pid('belt', row['id'], conn)}", belt.members, total
    finally:
        conn.close()
    pytest.skip("the sample galaxy holds no belt with more bodies at the class D size than it stores")


def test_a_body_past_the_stored_ones_has_a_page_and_a_class_made_when_asked(client, rich_belt):
    """GEN.1674: a rank above the stored bodies but inside the belt's count at the Class D size is drawn from the belt's size law, with the same page and a class from the one planet generator."""
    ref, members, total = rich_belt
    past = max(member["rank"] for member in members) + 20
    data = client.get(f"/api/objects/{ref}/belt-members/{_rank(past)}")
    assert data.status_code == 200
    body = data.get_json()
    assert body["rank"] == past and body["planet"]["class"] and body["previous"] and body["member"]["mass"]["value"] >= 1.0e19 and body["belt"]["class_d_count"] == total
    assert body == client.get(f"/api/objects/{ref}/belt-members/{_rank(past)}").get_json(), "the same body whenever it is asked for"
    page = client.get(f"/object/{ref}/member/{_rank(past)}")
    assert page.status_code == 200 and f"Class {body['planet']['class']}" in page.get_data(as_text=True)
    assert client.get(f"/api/objects/{ref}/belt-members/{_rank(total * 4 + 10)}").status_code == 404
