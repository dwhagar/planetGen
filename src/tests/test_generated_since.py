# tests/test_generated_since.py

"""ADM.150: `GET /api/generated`, what a run has made since a time, and the admin Generated so far page that reads it."""

import pytest

from planetgen.db import query
from planetgen.web.lib import apiclient
from tests.test_api import client, seeded_sector  # noqa: F401
from tests.test_web_admin import _logged_in, fake  # noqa: F401
from tests.test_web_admin import app as admin_app  # noqa: F401
from tests.test_web_admin import client as admin_client  # noqa: F401


@pytest.mark.parametrize("kind", query.GENERATED_KINDS)
def test_everything_is_listed_from_time_zero_and_nothing_from_the_future(client, seeded_sector, kind):  # noqa: F811
    everything = client.get(f"/api/generated?kind={kind}&since=0&limit=500").get_json()
    assert everything["kind"] == kind and everything["total"] >= len(everything["items"])
    ids = [item["id"] for item in everything["items"]]
    assert len(ids) == len(set(ids))
    assert all(item["kind"] == kind[:-1] for item in everything["items"])
    later = client.get(f"/api/generated?kind={kind}&since=99999999999").get_json()
    assert later["total"] == 0 and later["items"] == []


def test_the_systems_come_newest_first_and_page(client, seeded_sector):  # noqa: F811
    first = client.get("/api/generated?kind=systems&since=0&limit=1").get_json()
    rest = client.get("/api/generated?kind=systems&since=0&limit=1&offset=1").get_json()
    assert first["total"] == rest["total"]
    if first["total"] > 1:
        assert first["items"][0]["created"] >= rest["items"][0]["created"]
        assert first["items"][0]["id"] != rest["items"][0]["id"]


@pytest.mark.parametrize("query_string", ["kind=stars&since=0", "since=0", "kind=systems", "kind=systems&since=soon",
                                          "kind=systems&since=nan", "kind=systems&since=inf"])
def test_a_bad_kind_or_time_is_a_400(client, seeded_sector, query_string):  # noqa: F811
    assert client.get(f"/api/generated?{query_string}").status_code == 400


def test_the_admin_page_lists_what_the_api_gave(admin_client, fake, monkeypatch):  # noqa: F811
    asked = []

    def generated(db, kind, since, limit=None, offset=None):
        asked.append((kind, since, limit, offset))
        return {"kind": kind, "since": since, "total": 1, "limit": limit, "offset": offset, "items": [
            {"kind": "planet", "id": "5000012345", "name": "Vega b", "planet_class": "M", "radius_km": 6371.0,
             "star_system_id": "5000067890", "system_name": "Vega"}]}

    monkeypatch.setattr(apiclient, "get_generated", generated)
    _logged_in(admin_client, fake)
    html = admin_client.get("/admin/generated?kind=planets&since=100").get_data(as_text=True)
    assert asked == [("planets", 100, 50, 0)]
    assert "Vega b" in html and "6,371 km" in html and "1 planets" in html
    assert 'data-admin-shell' in html
    assert admin_client.get("/admin/generated?kind=nonsense&since=100").status_code == 200
    assert asked[-1][0] == "systems"
