# tests/test_web_object_lists.py

"""UX.141: `/lists` names every list and `/lists/<kind>` is a keyset data table over the object index
(`web/object_lists.py`): a page asks the API for its kind, the table route hands `after` on and gives the last row's
object ID as `next`, rows link to the object and show its class, mass and radius, and an unknown kind is a 404."""

import pytest

from planetgen.web import object_lists


def _item(index, kind="planet"):
    return {"kind": kind, "id": f"FE81000A2B-0000005-00{index}", "source_id": index, "name": f"Thing {index}",
            "sector_id": 4, "system_id": 9, "parent_id": None, "object_class": "gas_giant", "mass_kg": 1.9e27,
            "radius_km": 7.0e4, "placed": True}


@pytest.fixture
def site(monkeypatch):
    from planetgen.web.app import create_app
    from planetgen.web.lib import apiclient
    from tests.test_web_pages import FakeData, _FakeConfig

    seen = []

    def get_object_list(db, kind, limit=None, sort=None, descending=False, after=None, name=None, object_class=None,
                        sector_id=None, start=None):
        seen.append({"kind": kind, "limit": limit, "sort": sort, "descending": descending, "after": after, "start": start})
        items = [_item(index, kind) for index in range(1, 3)]
        return {"items": items, "next": items[-1]["id"], "total": 7, "total_capped": False}

    data = FakeData()
    monkeypatch.setattr(apiclient, "auth_me", data.auth_me)
    monkeypatch.setattr(apiclient, "get_object_list", get_object_list)
    app = create_app(_FakeConfig)
    app.testing = True
    return app.test_client(), seen


def test_the_lists_page_names_every_list_with_its_count(site):
    client, seen = site
    body = client.get("/lists").get_data(as_text=True)
    for kind, labels in object_lists.KIND_LABELS.items():
        assert f"/lists/{kind}" in body and labels[2] in body
    assert {call["kind"] for call in seen} == set(object_lists.KIND_LABELS)
    assert "/systems" in body and "/sectors" in body and "/phenomena" in body


def test_a_list_page_shows_the_first_rows_and_links_each_to_its_object(site):
    client, seen = site
    response = client.get("/lists/planet")
    body = response.get_data(as_text=True)
    assert response.status_code == 200 and "Thing 1" in body and "Thing 2" in body
    assert "/object/planet:FE81000A2B-0000005-001" in body
    assert "Gas giant" in body and "Sector 4" in body and seen[-1]["kind"] == "planet" and seen[-1]["after"] is None


def test_the_table_route_passes_the_key_on_and_gives_the_last_object_as_next(site):
    client, seen = site
    body = client.get("/table/list-moon?after=FE81000A2B-0000005-001&sort=mass&order=desc&limit=2").get_json()
    assert seen[-1] == {"kind": "moon", "limit": 2, "sort": "mass", "descending": True,
                        "after": "FE81000A2B-0000005-001", "start": None}
    assert body["next"] == "FE81000A2B-0000005-002" and len(body["rows"]) == 2 and body["total"] == 7
    assert body["rows"][0][0]["href"] == "/object/moon:FE81000A2B-0000005-001"


def test_an_unknown_kind_has_no_page_or_table(site):
    client, _seen = site
    assert client.get("/lists/galaxy").status_code == 404
    assert client.get("/lists/facility").status_code == 404
    assert client.get("/table/list-galaxy").status_code == 404


def test_the_jump_box_begins_a_name_sorted_list_at_a_name(site):
    client, seen = site
    assert "_from" in client.get("/lists/star").get_data(as_text=True)
    client.get("/table/list-star?from=Tr&limit=2")
    assert seen[-1]["start"] == "Tr"
