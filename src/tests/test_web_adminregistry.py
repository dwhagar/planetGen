# tests/test_web_adminregistry.py

"""
The admin registry (ADM.118, `planetgen.web.adminregistry`): every admin page route is in it or exempt with a
reason, every entry opens a real page, settings come from the model, and the search ranks as documented.

Run with: pytest tests/test_web_adminregistry.py
"""

import json

import pytest

from planetgen.util import logpaths, settings, settings_view
from planetgen.web import adminregistry
from planetgen.web.app import create_app


class _Config:
    pass


@pytest.fixture
def config(tmp_path, monkeypatch):
    path = tmp_path / "config.json"
    path.write_text(json.dumps({"secret_key": "csrf-secret"}) + "\n")
    monkeypatch.setattr(logpaths, "CONFIG_PATH", str(path))
    for variable in settings.env_names():
        monkeypatch.delenv(variable, raising=False)
    settings.reset_cache()
    yield path
    settings.reset_cache()
    settings_view.forget_startup()


@pytest.fixture(scope="module")
def app():
    from tests.test_web_admin import _FakeConfig

    return create_app(_FakeConfig)


def _admin_page_endpoints(app):
    """Endpoints of the GET routes under /admin or /account that draw an HTML page."""
    found = set()
    for rule in app.url_map.iter_rules():
        if not rule.rule.startswith(("/admin", "/account")) or "GET" not in rule.methods:
            continue
        found.add(rule.endpoint.removeprefix("web."))
    return found


def test_every_admin_page_is_in_the_registry_or_exempt(app):
    registered = {entry.endpoint for entry in adminregistry.PAGES}
    missing = _admin_page_endpoints(app) - registered - set(adminregistry.EXEMPT)
    assert not missing, f"add {sorted(missing)} to adminregistry.PAGES, or to EXEMPT with the reason"


def test_nothing_is_both_a_page_and_exempt_and_every_name_is_a_route(app):
    names = {rule.endpoint.removeprefix("web.") for rule in app.url_map.iter_rules()}
    registered = {entry.endpoint for entry in adminregistry.PAGES}
    assert not registered & set(adminregistry.EXEMPT)
    assert registered <= names
    assert set(adminregistry.EXEMPT) <= names, "an exempt endpoint that no longer exists"
    assert {endpoint for _key, _label, endpoint in adminregistry.DESTINATIONS} <= names


def test_entries_are_well_formed(app, config):
    found = adminregistry.entries()
    keys = [entry.key for entry in found]
    assert len(keys) == len(set(keys))
    destinations = {key for key, _label, _endpoint in adminregistry.DESTINATIONS}
    for entry in found:
        assert entry.kind in ("page", "task", "setting")
        assert entry.danger in adminregistry.DANGER
        assert entry.destination in destinations
        assert entry.label and entry.summary.endswith(".")
    with app.test_request_context("/"):
        from planetgen.web.helpers import page_url

        for entry in found:
            assert adminregistry.url(entry, page_url).startswith("/")


def test_settings_entries_are_the_options_an_admin_can_change(config):
    entries = adminregistry.setting_entries("admin")
    allowed = [o for o in settings_view.describe("admin")["options"] if o["allowed"]]
    assert [e.label for e in entries] == [o["title"] for o in allowed]
    assert all(e.kind == "setting" and e.destination == "settings" for e in entries)


def test_search_ranks_label_starts_first_and_needs_every_word(config):
    found = adminregistry.entries()
    labels = [e.label for e in adminregistry.search("settings", found)]
    assert labels[0] == "Settings"
    assert [e.key for e in adminregistry.search("backup restore", found)][0] == "settings-backups"
    assert adminregistry.search("zzzz no such thing", found) == []
    assert {e.kind for e in adminregistry.search("", found)} == {"page"}
    assert len(adminregistry.search("a", found, limit=5)) == 5


def test_shell_lights_the_destination_of_the_page_and_is_none_elsewhere():
    shell = adminregistry.shell("admin_facility_types", lambda endpoint: f"/{endpoint}")
    assert [d["key"] for d in shell["destinations"] if d["current"]] == ["data"]
    assert len(shell["destinations"]) == 6
    assert adminregistry.shell("galaxy", lambda endpoint: "/") is None
