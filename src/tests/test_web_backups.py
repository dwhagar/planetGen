# tests/test_web_backups.py

"""
The backup archives screen and the settings export and import (ADM.85): `/api/admin/backups`, `/api/admin/settings/
export|import` and the pages `/admin/settings/backups[/<name>]`, `/admin/settings/export` and `/admin/settings/import`,
through the real app, API and control database, over a temporary `config.json` and state folder (the settings
page itself is in test_web_settings.py).

Run with: pytest tests/test_web_backups.py
"""

import html as htmllib
import io
import json
import re

import pytest

from planetgen import backup
from planetgen.api import settings as settings_api
from planetgen.backup import core as backup_core
from planetgen.util import settings_writer
from planetgen.web import csrf

from tests.test_web_settings import START, _audit_rows, config, site  # noqa: F401 -- fixtures and helpers


def _form(client, path, **fields):
    data = {csrf.FIELD_NAME: client.csrf}
    data.update({key.replace("__", "."): value for key, value in fields.items()})
    return client.post(path, data=data)


def _doc(config):
    return json.loads(config.read_text())


def _archives(client):
    return client.get("/api/admin/backups").get_json()["archives"]


def _hidden(page):
    return {name: htmllib.unescape(value) for name, value in
            re.findall(r'<input type="hidden" name="([^"]+)" value="([^"]*)"', page)}


# --- Listing and Back up now -----------------------------------------------------------

def test_the_list_is_empty_until_a_backup_is_made(site, config):
    listing = site.get("/api/admin/backups").get_json()
    assert listing["archives"] == [] and listing["problem"] == "" and listing["kept"] == 20
    assert listing["parts"] == [{"name": "settings", "label": "Settings"}] and listing["can_download"] is False
    assert "No backup has been made yet." in site.get("/admin/settings/backups").get_data(as_text=True)


def test_back_up_now_makes_one_archive_and_an_audit_row(site, config):
    made = site.post("/api/admin/backups", json={}).get_json()
    assert made["reason"] == "manual" and made["complete"] and [p["name"] for p in made["parts"]] == ["settings"]
    assert [a["name"] for a in _archives(site)] == [made["name"]]
    conn = __import__("planetgen.db.store", fromlist=["x"]).get_control_connection(site.config)
    try:
        rows = conn.execute("SELECT action, target FROM admin_audit_log WHERE action = 'backup.take'").fetchall()
    finally:
        conn.close()
    assert [(r["action"], r["target"]) for r in rows] == [("backup.take", f"backup:{made['name']}")]


def test_the_page_makes_a_backup_and_lists_it(site, config):
    resp = _form(site, "/admin/settings/backups")
    assert resp.status_code == 303
    page = site.get("/admin/settings/backups").get_data(as_text=True)
    assert "Backed up Settings in " in page
    name = _archives(site)[0]["name"]
    assert name in page and "Compare and restore" in page and "Download is for the Owner" in page


def test_back_up_now_says_when_the_folder_cannot_be_used(site, config, monkeypatch):
    monkeypatch.setattr(backup, "write_problem", lambda *a, **k: "The backups folder is not writable.")
    assert site.post("/api/admin/backups", json={}).status_code == 409
    page = site.get("/admin/settings/backups").get_data(as_text=True)
    assert "not writable" in page and "disabled" in page


def test_a_settings_save_leaves_its_own_archive_in_the_list(site, config):
    settings_writer.save({"jobs.keep": 5}, actor="boss")
    archives = _archives(site)
    assert len(archives) == 1 and archives[0]["reason"] == "config.update" and archives[0]["taken_by"] == "boss"


def test_the_api_needs_an_admin(site):
    anonymous = site.application.test_client()
    for path in ("/api/admin/backups", "/api/admin/settings/export"):
        assert anonymous.get(path).status_code == 401


# --- Compare and restore ---------------------------------------------------------------

def test_an_archive_page_shows_what_a_restore_would_change(site, config):
    name = site.post("/api/admin/backups", json={}).get_json()["name"]
    settings_writer.save({"site_name": "Home", "page_cache.enabled": False}, actor="boss")
    shown = site.get(f"/api/admin/backups/{name}").get_json()
    assert shown["verified"] and shown["problems"] == []
    changes = {c["path"]: (c["old"], c["new"]) for c in shown["diff"]["changes"]}
    assert changes == {"site_name": ("Home", "Mine"), "page_cache.enabled": (False, True)}
    page = site.get(f"/admin/settings/backups/{name}").get_data(as_text=True)
    assert "What restoring the settings would change" in page and "Mine" in page and "needs a restart" in page
    assert 'name="part" value="settings" checked' in page and "Restore the ticked parts" in page
    assert "csrf-secret" not in page and "db-secret" not in page


def test_an_archive_that_matches_today_says_so(site, config):
    name = site.post("/api/admin/backups", json={}).get_json()["name"]
    assert "same as today's" in site.get(f"/admin/settings/backups/{name}").get_data(as_text=True)


def test_restoring_puts_the_settings_back_and_is_backed_up_and_audited(site, config):
    name = site.post("/api/admin/backups", json={}).get_json()["name"]
    settings_writer.save({"site_name": "Home", "jobs.keep": 3}, actor="boss")
    resp = _form(site, f"/admin/settings/backups/{name}/restore", part="settings",
                 expect_hash=site.get(f"/api/admin/backups/{name}").get_json()["file_hash"])
    assert resp.status_code == 303
    assert _doc(config)["site_name"] == "Mine" and _doc(config)["jobs"] == {"keep": 12}
    page = site.get(f"/admin/settings/backups/{name}").get_data(as_text=True)
    assert "Restored settings from" in page and "2 settings changed" in page and "What was replaced is kept in" in page
    rows = _audit_rows(site)
    assert ("config.restore", "config:site_name", {"old": "Home", "new": "Mine"}) in rows
    reasons = [a["reason"] for a in _archives(site)]
    assert any(reason.startswith("before restoring") for reason in reasons)
    conn = __import__("planetgen.db.store", fromlist=["x"]).get_control_connection(site.config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM admin_audit_log WHERE action = 'backup.restore'"
                            ).fetchone()["n"] == 1
    finally:
        conn.close()


def test_a_web_restore_cannot_change_the_database_password_or_the_secret_key(site, config):
    name = site.post("/api/admin/backups", json={}).get_json()["name"]
    settings_writer.save({"mysql.password": "rotated", "secret_key": "rotated-key", "site_name": "Home"},
                         actor="cli:boss", role="host")
    site.post(f"/api/admin/backups/{name}/restore", json={"parts": ["settings"]})
    document = _doc(config)
    assert document["mysql"]["password"] == "rotated" and document["secret_key"] == "rotated-key"
    assert document["site_name"] == "Mine"


def test_restoring_one_part_leaves_the_others_alone(site, config, monkeypatch):
    class Extra(backup.Part):
        name, label = "extra", "Extra data"
        value = "one"

        def dump(self, directory):
            with open(f"{directory}/data.txt", "w") as out:
                out.write(Extra.value)

        def restore(self, directory, context):
            Extra.restored = open(f"{directory}/data.txt").read()

    monkeypatch.setitem(backup_core._PARTS, "extra", Extra())
    name = site.post("/api/admin/backups", json={}).get_json()["name"]
    assert [p["label"] for p in _archives(site)[0]["parts"]] == ["Settings", "Extra data"]
    settings_writer.save({"site_name": "Home"}, actor="boss")
    Extra.value = "two"
    body = site.post(f"/api/admin/backups/{name}/restore", json={"parts": ["extra"]}).get_json()
    assert body["restored"] == ["extra"] and Extra.restored == "one" and body["safety"]
    assert _doc(config)["site_name"] == "Home"                  # the settings part was not ticked
    whole = site.post(f"/api/admin/backups/{name}/restore", json={"parts": None}).get_json()
    assert whole["restored"] == ["settings", "extra"] and _doc(config)["site_name"] == "Mine"


def test_a_part_this_build_does_not_know_is_skipped_and_said_so(site, config, monkeypatch):
    class Later(backup.Part):
        name, label = "later", "Later"

        def dump(self, directory):
            with open(f"{directory}/x", "w") as out:
                out.write("x")

    monkeypatch.setitem(backup_core._PARTS, "later", Later())
    name = site.post("/api/admin/backups", json={}).get_json()["name"]
    monkeypatch.delitem(backup_core._PARTS, "later")
    shown = site.get(f"/api/admin/backups/{name}").get_json()["archive"]
    assert [p["known"] for p in shown["parts"]] == [True, False]
    page = site.get(f"/admin/settings/backups/{name}").get_data(as_text=True)
    assert "can't restore this part" in page and "disabled" in page
    body = site.post(f"/api/admin/backups/{name}/restore", json={"parts": None}).get_json()
    assert body["restored"] == ["settings"] and body["skipped"] == ["later"]


def test_a_backup_that_no_longer_fits_the_model_is_explained_and_refused(site, config):
    config.write_text(json.dumps({"jobs": {"keep": 0}}))             # an archive of a file the model rejects
    name = backup.take(who="boss", reason="manual").name
    config.write_text(json.dumps(START))
    shown = site.get(f"/api/admin/backups/{name}").get_json()
    assert not shown["diff"]["ok"] and "jobs.keep" in shown["diff"]["errors"]
    assert "can not be restored" in site.get(f"/admin/settings/backups/{name}").get_data(as_text=True)
    assert site.post(f"/api/admin/backups/{name}/restore", json={"parts": ["settings"]}).status_code == 409
    assert _doc(config)["site_name"] == "Mine"


def test_a_damaged_archive_is_refused(site, config):
    name = site.post("/api/admin/backups", json={}).get_json()["name"]
    path = backup.backups_dir() + "/" + name
    data = bytearray(open(path, "rb").read())
    data[len(data) // 2] ^= 0xFF
    open(path, "wb").write(bytes(data))
    shown = site.get(f"/api/admin/backups/{name}").get_json()
    assert shown["verified"] is False and shown["diff"] is None
    assert "damaged" in site.get(f"/admin/settings/backups/{name}").get_data(as_text=True)
    assert site.post(f"/api/admin/backups/{name}/restore", json={"parts": ["settings"]}).status_code == 409


def test_a_restore_from_a_stale_page_is_refused(site, config):
    name = site.post("/api/admin/backups", json={}).get_json()["name"]
    stale = site.get(f"/api/admin/backups/{name}").get_json()["file_hash"]
    settings_writer.save({"jobs.keep": 3}, actor="someone")
    resp = site.post(f"/api/admin/backups/{name}/restore", json={"parts": ["settings"], "expect_hash": stale})
    assert resp.status_code == 409 and "changed by someone else" in resp.get_json()["error"]
    page = _form(site, f"/admin/settings/backups/{name}/restore", part="settings", expect_hash=stale)
    assert page.status_code == 303 and _doc(config)["jobs"] == {"keep": 3}


def test_nothing_ticked_and_unknown_names_are_refused(site, config):
    name = site.post("/api/admin/backups", json={}).get_json()["name"]
    assert site.post(f"/api/admin/backups/{name}/restore", json={"parts": []}).status_code == 400
    assert site.post(f"/api/admin/backups/{name}/restore", json={"parts": ["nope"]}).status_code == 409
    assert site.get("/api/admin/backups/20200101T000000Z.tar.bz2").status_code == 404
    assert site.get("/api/admin/backups/not-an-archive").status_code == 404
    resp = site.get("/admin/settings/backups/20200101T000000Z.tar.bz2")
    assert resp.status_code == 303
    assert "There is no such backup." in site.get("/admin/settings/backups").get_data(as_text=True)
    flash = _form(site, f"/admin/settings/backups/{name}/restore")
    assert flash.status_code == 303 and "Tick at least one part" in site.get(
        f"/admin/settings/backups/{name}").get_data(as_text=True)


# --- Download is the Owner's -----------------------------------------------------------

def test_an_admin_cannot_download_an_archive_and_the_owner_can(site, config, monkeypatch):
    name = site.post("/api/admin/backups", json={}).get_json()["name"]
    resp = site.get(f"/api/admin/backups/{name}/download")
    assert resp.status_code == 403 and "Owner" in resp.get_json()["error"]
    monkeypatch.setattr(settings_api, "ROLE", "owner")
    resp = site.get(f"/api/admin/backups/{name}/download")
    assert resp.status_code == 200 and resp.data == open(backup.backups_dir() + "/" + name, "rb").read()
    assert "attachment" in resp.headers["Content-Disposition"]
    assert "Download</a>" in site.get("/admin/settings/backups").get_data(as_text=True)


# --- Export and import -----------------------------------------------------------------

def test_an_export_holds_only_what_an_admin_may_change_and_no_secret(site, config):
    settings_writer.save({"login_allowlist": ["10.0.0.1"], "base_url": "http://example.org/"}, actor="boss",
                         role="host")
    exported = site.get("/api/admin/settings/export").get_json()
    assert exported == {"config_version": 1, "site_name": "Mine", "jobs": {"keep": 12}}
    text = json.dumps(exported)
    assert "db-secret" not in text and "csrf-secret" not in text and "mysql" not in exported
    resp = site.get("/admin/settings/export")
    assert resp.status_code == 200 and "attachment" in resp.headers["Content-Disposition"]
    assert json.loads(resp.data) == exported


def test_importing_an_export_sets_what_it_names_and_leaves_the_rest(site, config):
    document = {"config_version": 1, "site_name": "Imported", "tile_cache": {"max_mb": 99}}
    checked = site.post("/api/admin/settings/import", json={"document": document, "check_only": True}).get_json()
    assert checked["ok"] and not checked["saved"] and {c["path"] for c in checked["changes"]} == {
        "site_name", "tile_cache.max_mb"}
    assert _doc(config)["site_name"] == "Mine"
    saved = site.post("/api/admin/settings/import", json={"document": document, "expect_hash": checked["file_hash"]}
                      ).get_json()
    assert saved["saved"] and saved["backup"]
    assert _doc(config)["site_name"] == "Imported" and _doc(config)["tile_cache"] == {"max_mb": 99}
    assert _doc(config)["jobs"] == {"keep": 12} and _doc(config)["mysql"]["password"] == "db-secret"
    assert ("config.import", "config:site_name", {"old": "Mine", "new": "Imported"}) in _audit_rows(site)


def test_an_export_imports_back_to_the_same_settings(site, config):
    exported = site.get("/api/admin/settings/export").get_json()
    settings_writer.save({"site_name": "Home", "jobs.keep": 3}, actor="boss")
    assert site.post("/api/admin/settings/import", json={"document": exported}).get_json()["saved"]
    assert _doc(config)["site_name"] == "Mine" and _doc(config)["jobs"] == {"keep": 12}


@pytest.mark.parametrize("document, option", [
    ({"jobs": {"keep": 0}}, "jobs.keep"),
    ({"jobs": {"keep": "many"}}, "jobs.keep"),
    ({"mysql": {"password": "x"}}, "mysql.password"),
    ({"secret_key": "x"}, "secret_key"),
    ({"base_url": "http://a/"}, "base_url"),
    ({"no_such_option": 1}, "no_such_option"),
    ({"jobs": {"nope": 1}}, "jobs.nope"),
])
def test_a_bad_value_in_an_import_changes_nothing(site, config, document, option):
    before = config.read_bytes()
    document = {**document, "site_name": "Imported"}              # a good value beside it is not saved either
    body = site.post("/api/admin/settings/import", json={"document": document}).get_json()
    assert not body["ok"] and not body["saved"] and option in body["errors"]
    assert config.read_bytes() == before


def test_an_import_from_a_newer_release_or_not_an_object_is_refused(site, config):
    before = config.read_bytes()
    body = site.post("/api/admin/settings/import", json={"document": {"config_version": 99}}).get_json()
    assert not body["saved"] and "newer" in body["errors"][""]
    assert config.read_bytes() == before


def test_the_import_page_checks_reviews_and_saves(site, config):
    upload = json.dumps({"config_version": 1, "site_name": "Imported"}).encode()
    resp = site.post("/admin/settings/import", data={csrf.FIELD_NAME: site.csrf,
                                                     "file": (io.BytesIO(upload), "planetgen-settings.json")})
    page = resp.get_data(as_text=True)
    assert resp.status_code == 200 and "This import will change these settings" in page and "Imported" in page
    assert _doc(config)["site_name"] == "Mine"
    hidden = _hidden(page)
    saved = _form(site, "/admin/settings/import/save", document=hidden["document"], expect_hash=hidden["expect_hash"])
    assert saved.status_code == 303 and _doc(config)["site_name"] == "Imported"
    assert "Imported 1 setting." in site.get("/admin/settings").get_data(as_text=True)


@pytest.mark.parametrize("payload, shown", [
    (b"not json", "not valid JSON"),
    (b"[1, 2]", "JSON object"),
    (json.dumps({"jobs": {"keep": 0}}).encode(), "jobs &gt; keep"),
    (json.dumps({"mysql": {"password": "x"}}).encode(), "on the server"),
    (json.dumps({"no_such": 1}).encode(), "is not a setting"),
    (b"x" * (1024 * 1024 + 1), "too big"),
])
def test_the_import_page_names_the_problem_and_saves_nothing(site, config, payload, shown):
    before = config.read_bytes()
    resp = site.post("/admin/settings/import", data={csrf.FIELD_NAME: site.csrf,
                                                     "file": (io.BytesIO(payload), "settings.json")})
    assert resp.status_code in (400, 422) and shown in resp.get_data(as_text=True)
    assert config.read_bytes() == before


def test_the_import_page_without_a_file_or_with_nothing_new(site, config):
    assert site.post("/admin/settings/import", data={csrf.FIELD_NAME: site.csrf}).status_code == 400
    same = json.dumps({"site_name": "Mine"}).encode()
    resp = site.post("/admin/settings/import", data={csrf.FIELD_NAME: site.csrf, "file": (io.BytesIO(same), "s.json")})
    assert resp.status_code == 200 and "Nothing to import" in resp.get_data(as_text=True)


def test_an_import_reviewed_before_someone_elses_save_is_refused(site, config):
    upload = json.dumps({"site_name": "Imported"}).encode()
    page = site.post("/admin/settings/import", data={csrf.FIELD_NAME: site.csrf,
                                                     "file": (io.BytesIO(upload), "s.json")}).get_data(as_text=True)
    hidden = _hidden(page)
    settings_writer.save({"jobs.keep": 3}, actor="someone")
    _form(site, "/admin/settings/import/save", document=hidden["document"], expect_hash=hidden["expect_hash"])
    assert _doc(config)["site_name"] == "Mine"
    assert "changed by someone else" in site.get("/admin/settings/backups").get_data(as_text=True)


def test_the_new_pages_need_an_admin(site):
    anonymous = site.application.test_client()
    for path in ("/admin/settings/backups", "/admin/settings/export"):
        assert anonymous.get(path).status_code == 302
    assert "Backups, restore, export and import" in site.get("/admin/settings").get_data(as_text=True)
