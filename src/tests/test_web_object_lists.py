# tests/test_web_object_lists.py

"""UX.141: `/lists` names every list and `/lists/<kind>` is a keyset data table over the object index
(`web/object_lists.py`): a page asks the API for its kind, the table route hands `after` on and gives the last row's
object ID as `next`, rows link to the object and show its class, mass and radius, and an unknown kind is a 404."""

import pytest

from planetgen.web import object_lists


def _item(index, kind="planet"):
    item = {"kind": kind, "id": f"FE81000A2B-0000005-00{index}", "source_id": index, "name": f"Thing {index}",
            "sector_id": 4, "system_id": 9, "parent_id": None, "object_class": "gas_giant", "mass_kg": 1.9e27,
            "radius_km": 7.0e4, "placed": True}
    if kind == "star":
        item.update(temperature_k=5772.0, luminosity_w=3.828e26, age_gy=4.6)
    if kind in ("planet", "moon"):
        item.update(distance_km=1.496e8 if index == 1 else None, eccentricity=0.0167 if index == 1 else None,
                    surface_temperature_k=288.0 if index == 1 else None)
    return item


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


def test_a_star_list_shows_its_temperature_luminosity_and_age(site):
    """UX.187: the stored physics values of a star, in the reader's units."""
    client, _seen = site
    body = client.get("/lists/star").get_data(as_text=True)
    for heading in ("Temperature", "Luminosity", "Age"):
        assert heading in body
    assert "5,772 K" in body and "1 L" in body and "4.6 Gyr" in body


def test_a_planet_or_moon_list_shows_its_orbit_eccentricity_and_surface_temperature_and_a_dash_where_none_is_stored(site):
    client, _seen = site
    for kind in ("planet", "moon"):
        body = client.get(f"/lists/{kind}").get_data(as_text=True)
        for heading in ("Orbit", "Eccentricity", "Surface temperature"):
            assert heading in body
        assert "1 AU" in body and "0.017" in body and "288 K" in body
        assert body.count("\u2013") >= 3     # the second row stores none of the three


def test_the_table_route_passes_the_key_on_and_gives_the_last_object_as_next(site):
    client, seen = site
    body = client.get("/table/list-moon?after=FE81000A2B-0000005-001&sort=mass&order=desc&limit=2").get_json()
    assert seen[-1] == {"kind": "moon", "limit": 2, "sort": "mass", "descending": True,
                        "after": "FE81000A2B-0000005-001", "start": None}
    assert body["next"] == "FE81000A2B-0000005-002" and len(body["rows"]) == 2 and body["total"] == 7
    assert body["rows"][0][0]["href"] == "/object/moon:FE81000A2B-0000005-001"


def test_planet_and_moon_lists_have_a_habitability_column_of_symbols(site, monkeypatch):
    """UX.100: a scored row draws the overall symbol and the six domains; an unscored one a dash; stars have no column."""
    from planetgen.web.lib import apiclient

    client, _seen = site
    scored = {**_item(1, "moon"), "phi4": 0.8, "tier_pressure": 0, "tier_temperature": 1, "tier_air": 2,
              "tier_water": 0, "tier_radiation": 3, "tier_suitability": 1}
    unscored = {**_item(2, "moon"), "phi4": None, **{f"tier_{d}": None for d in
                                                      ("pressure", "temperature", "air", "water", "radiation", "suitability")}}
    monkeypatch.setattr(apiclient, "get_object_list", lambda db, kind, **kw: {
        "items": [scored, unscored], "next": None, "total": 2, "total_capped": False})
    body = client.get("/lists/moon").get_data(as_text=True)
    assert 'data-col="habitability"' in body and body.count("hab-symbol ") == 7
    for label in ("Overall habitability: Green", "Air chemistry: Yellow", "Radiation: Red"):
        assert f'aria-label="{label}"' in body, label
    served = client.get("/table/list-moon").get_json()["rows"]
    assert [symbol["name"] for symbol in served[0][7]["symbols"]] == [
        "overall", "pressure", "temperature", "air", "water", "radiation", "suitability"]
    assert served[1][7] == {"text": "—", "muted": True}
    assert 'data-col="habitability"' not in client.get("/lists/star").get_data(as_text=True)
    # UX.181: each symbol of the page's rows, and of the rows the table serves, links to its habitability section.
    assert [symbol["anchor"] for symbol in served[0][7]["symbols"]][:3] == ["habitability"] * 3
    assert served[0][7]["page"].startswith("/object/moon:")


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


def test_a_slow_table_read_answers_202_with_the_job_and_where_to_poll_it(site, monkeypatch):
    """UX.136: the table route waits on a page job like the pages do; `X-Slow-Retry` reaches the cost rule."""
    from planetgen.web.lib import apiclient, slowpages

    job = {"id": "0123456789abcdef", "state": "running", "done": 300, "total": 1000, "eta_seconds": 5.0, "error": None}
    seen = []

    def slow(db, kind, **kwargs):
        seen.append(slowpages._retrying.get())
        raise slowpages.Pending(job, "/lists/planet")

    monkeypatch.setattr(apiclient, "get_object_list", slow)
    client, _seen = site
    response = client.get("/table/list-planet")
    assert response.status_code == 202 and response.headers["Cache-Control"] == "no-store"
    body = response.get_json()
    assert body["pending"] is True and body["job"] == job and body["percent"] == 30
    assert body["status_url"] == "/api/page-jobs/0123456789abcdef"
    client.get("/table/list-planet", headers={"X-Slow-Retry": "1"})
    assert seen == [False, True]


def test_planet_and_moon_lists_link_to_the_score_filters(site):
    """UX.214: the lists that carry habitability point at the search page, best habitability first."""
    client, _seen = site
    for kind in ("planet", "moon"):
        html = client.get(f"/lists/{kind}").get_data(as_text=True)
        assert f"type={kind}" in html and f"{kind}s_sort=phi4" in html and f"{kind}s_order=desc" in html
    assert "list-filter-link" not in client.get("/lists/star").get_data(as_text=True)
