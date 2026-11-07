# tests/test_field_member_page.py

"""
GEN.1105 part 2: the bodies an interstellar asteroid field keeps (`GET /api/objects/<ref>/field-members` and `.../field-members/<rank>`), the page of one of them
(`/object/asteroid_field:<ID>/field-member/<rank>`) and the list on the field's own page, walked on the small generated database `test_web_a11y.py` builds with one hosted field added. Skipped
without the MySQL test server.
"""

import pytest

from planetgen.db import store
from planetgen.galaxy import object_uid
from planetgen.generation import field_host
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.util import draw
from tests.publicids import pid
from tests.test_web_a11y import site_app, site_db  # noqa: F401 -- fixtures


def _hosted_field():
    for seed in range(1, 4000):
        with draw.bound(f"remnant:{seed}"):
            remnant = SupernovaRemnant(SystemConfig())
        hosted = field_host.hosted_field(remnant, SystemConfig(), f"k:{seed}")
        if hosted is not None and len(hosted.field.members) >= 2:
            return hosted.field
    pytest.fail("no hosted field with two members in 4000 remnants")


@pytest.fixture(scope="module")
def field(site_db):
    built = _hosted_field()
    conn = store.get_connection(site_db)
    try:
        with conn:
            field_id = store.insert_asteroid_field(conn, built)
        found = (f"asteroid_field:{pid('asteroid_field', field_id, conn)}", field_id, len(built.members), built)
    finally:
        conn.close()
    return found


@pytest.fixture(scope="module")
def client(site_app):
    return site_app.test_client()


def _rank(n):
    return f"{n:04X}"


def test_the_api_lists_the_bodies_heaviest_first_under_the_fields_id(client, field):
    ref, field_id, count, built = field
    data = client.get(f"/api/objects/{ref}/field-members")
    assert data.status_code == 200
    body = data.get_json()
    assert body["count"] == count == len(body["bodies"])
    assert [entry["rank"] for entry in body["bodies"]] == list(range(1, count + 1))
    masses = [entry["mass"]["value"] for entry in body["bodies"]]
    assert masses == sorted(masses, reverse=True)
    uid = body["field"]["designation"]
    assert body["bodies"][0]["designation"] == f"{uid}.{_rank(1)}"
    assert object_uid.parse_field_designation(body["bodies"][1]["designation"])[1] == 2


def test_the_api_answers_one_body_with_its_share_and_neighbours(client, field):
    ref, _, count, _ = field
    listed = client.get(f"/api/objects/{ref}/field-members").get_json()["bodies"]
    data = client.get(f"/api/objects/{ref}/field-members/{_rank(1)}")
    assert data.status_code == 200
    body = data.get_json()
    assert set(body) == {"designation", "rank", "member", "share", "field", "previous", "next"}
    assert body["member"] == listed[0] and body["previous"] is None and body["next"] == listed[1]["designation"]
    assert 0.0 < body["share"] <= 1.0
    last = client.get(f"/api/objects/{ref}/field-members/{_rank(count)}").get_json()
    assert last["next"] is None and last["previous"] == listed[count - 2]["designation"]


def test_a_rank_the_field_does_not_keep_or_a_wrong_kind_is_a_404(client, field):
    ref, _, count, _ = field
    assert client.get(f"/api/objects/{ref}/field-members/{_rank(count + 1)}").status_code == 404
    assert client.get(f"/api/objects/{ref}/field-members/1").status_code == 404   # not written as the designation writes it
    assert client.get(f"/api/objects/{ref}/field-members/{_rank(0)}").status_code == 404
    assert client.get("/api/objects/star:1/field-members").status_code == 404
    assert client.get(f"/api/objects/{ref}/belt-members/{_rank(1)}").status_code == 404


def test_the_page_shows_the_body_and_links_its_neighbours_and_its_field(client, field):
    ref, _, count, _ = field
    page = client.get(f"/object/{ref}/field-member/{_rank(1)}")
    assert page.status_code == 200
    text = page.get_data(as_text=True)
    designation = client.get(f"/api/objects/{ref}/field-members/{_rank(1)}").get_json()["designation"]
    assert designation in text and "It goes round no star" in text
    assert f"/field-member/{_rank(2)}" in text and f"/field-member/{_rank(0)}" not in text
    assert client.get(f"/object/{ref}/field-member/{_rank(count + 1)}").status_code == 404
    assert client.get(f"/object/{ref}/field-member/1").status_code == 404
    assert client.get(f"/object/belt:1/field-member/{_rank(1)}").status_code == 404


def test_the_fields_own_page_lists_its_bodies(client, field):
    ref, _, count, _ = field
    page = client.get(f"/phenomenon/asteroid_field/{ref.split(':')[1]}")
    assert page.status_code == 200
    text = page.get_data(as_text=True)
    assert "The bodies it keeps" in text and text.count("/field-member/") >= count
