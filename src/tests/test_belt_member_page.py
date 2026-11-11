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
    assert set(body) == {"designation", "rank", "member", "orbit", "share", "belt", "previous", "next"}
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
    for rank in (_rank(members + 1), "0000", "1", "00001", "ZZZZ"):
        assert client.get(f"/api/objects/{ref}/belt-members/{rank}").status_code == 404
        assert client.get(f"/object/{ref}/member/{rank}").status_code == 404
    assert client.get("/api/objects/planet:FE81000A2B-0000005-003/belt-members/0001").status_code == 404
    assert client.get("/object/planet:FE81000A2B-0000005-003/member/0001").status_code == 404
