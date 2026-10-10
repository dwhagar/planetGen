# tests/test_web_object_redirect.py

"""
NAV.8: `/object/<ref>` sends any object reference to the page or row that shows it. A star, planet, moon, belt
or comet opens its system's page at that body's row; nothing but a star's system gets a page of its own.
"""

import pytest

from planetgen.web.lib import apiclient
from tests.test_web_pages import FakeData, _FakeConfig  # noqa: F401
from planetgen.web.app import create_app


@pytest.fixture
def client(monkeypatch):
    parents = {"planet:12": [{"ref": "galaxy", "kind": "galaxy"}, {"ref": "sector:4", "kind": "sector"},
                             {"ref": "system:9", "kind": "system"}]}

    def get_object(db, ref):
        if ref not in parents:
            raise apiclient.NotFoundError("no such object")
        return {"ref": ref, "parents": parents[ref]}

    monkeypatch.setattr(apiclient, "get_object", get_object)
    app = create_app(_FakeConfig)
    app.testing = True
    return app.test_client()


@pytest.mark.parametrize("ref,target", [
    ("planet:12", "/system/9#planet-12"),
    ("system:9", "/system/9"),
    ("9", "/system/9"),
    ("sector:4", "/sector/4"),
    ("nebula:3", "/phenomenon/nebula/3"),
])
def test_a_reference_opens_its_page_or_row(client, ref, target):
    response = client.get(f"/object/{ref}")
    assert response.status_code == 302
    assert response.headers["Location"].endswith(target)


@pytest.mark.parametrize("ref", ["planet:99", "nonsense", "planet:"])
def test_an_unknown_reference_is_not_found(client, ref):
    assert client.get(f"/object/{ref}").status_code == 404
