# tests/test_settings_update.py

"""
"Update config.json" (ADM.140): the button on the admin Settings page that adds the settings this version has and the file
lacks, lists the ones a site should set, and flags stale keys, keeping every value the site has.

The writer (`settings_writer.plan_update`, `update`) is tested on its own with a file in a temporary folder; the API and the
page are tested with the same signed-in admin the settings page tests use.
"""

import copy
import json
import re

import pytest

from planetgen.util import logpaths, settings, settings_view, settings_writer
from tests.test_web_settings import START, _audit_rows, config, site  # noqa: F401 -- the fixtures of the settings page tests


def _leaves(node, prefix=""):
    for key, value in node.items():
        if isinstance(value, dict):
            yield from _leaves(value, f"{prefix}{key}.")
        else:
            yield f"{prefix}{key}", value


def _write(config, document):
    config.write_text(json.dumps(document, indent=2) + "\n")


def _fields_recommending(**wanted):
    """The model's options with `x-should-set` (and `x-recommended`) put on the named ones: `path=(why, value)` or `path=(why,)`."""
    fields = {path: (info, dict(meta)) for path, info, meta in settings.iter_fields()}
    for path, spec in wanted.items():
        fields[path][1]["x-should-set"] = spec[0]
        if len(spec) > 1:
            fields[path][1]["x-recommended"] = spec[1]
    return fields


# --- The writer ------------------------------------------------------------------------

def test_a_file_missing_keys_gets_them_with_their_defaults_and_keeps_every_old_value(config):
    before = json.loads(config.read_text())
    result = settings_writer.update(actor="boss")
    assert result.saved and result.backup
    after = json.loads(config.read_text())
    for path, value in _leaves(before):
        assert dict(_leaves(after))[path] == value                          # every old value, exactly
    added = {change.path for change in result.added}
    assert added == {path for path, _i, _m in settings.iter_fields()} - {path for path, _v in _leaves(before)}
    for path, info, _meta in settings.iter_fields():
        if path in added:
            assert dict(_leaves(after))[path] == settings.field_default(info)
    assert list(after)[:len(before)] == list(before)                       # the old keys keep their order, first
    assert settings.field_errors(after)[0] is not None                     # and the new file loads


def test_a_file_that_has_every_key_is_not_rewritten_and_no_backup_is_made(config):
    settings_writer.update(actor="boss")
    written = config.read_bytes()
    again = settings_writer.update(actor="boss")
    assert not again.saved and again.backup is None and again.changes == []
    assert config.read_bytes() == written


def test_a_file_missing_three_keys_gets_exactly_those_three(config):
    settings_writer.update(actor="boss")
    full = json.loads(config.read_text())
    for section, key in (("jobs", "keep"), ("cache", "max_mb")):
        del full[section][key]
    del full["site_name"]
    _write(config, full)
    result = settings_writer.update(actor="boss")
    assert sorted(change.path for change in result.added) == ["cache.max_mb", "jobs.keep", "site_name"]
    after = json.loads(config.read_text())
    assert dict(_leaves(after)) == {**dict(_leaves(full)), "site_name": "planetGen", "jobs.keep": 20, "cache.max_mb": 1024}


def test_the_options_the_web_cannot_change_are_added_with_their_defaults_and_never_altered(config):
    _write(config, {"site_name": "Mine"})
    result = settings_writer.update(actor="boss")
    after = json.loads(config.read_text())
    assert after["mysql"]["host"] == "127.0.0.1" and after["mysql"]["port"] == 3306
    assert "mysql.host" in {change.path for change in result.added}
    _write(config, {"mysql": {"host": "db.internal", "password": "db-secret"}})
    settings_writer.update(actor="boss")
    assert json.loads(config.read_text())["mysql"]["host"] == "db.internal"
    assert json.loads(config.read_text())["mysql"]["password"] == "db-secret"


def test_a_stale_key_is_flagged_and_survives_unless_ticked(config):
    _write(config, {**START, "old_option": 1, "jobs": {"keep": 12, "retired": "x"}, "mysql": {"host": "db.internal", "gone": 2}})
    planned = settings_writer.plan_update()
    assert planned.ok and sorted(planned.stale) == ["jobs.retired", "mysql.gone", "old_option"]
    settings_writer.update(actor="boss")
    kept = json.loads(config.read_text())
    assert kept["old_option"] == 1 and kept["jobs"]["retired"] == "x" and kept["mysql"]["gone"] == 2
    ticked = settings_writer.update(actor="boss", remove=["old_option", "jobs.retired"])
    assert ticked.saved and [change.path for change in ticked.removed] == ["jobs.retired", "old_option"]
    final = json.loads(config.read_text())
    assert "old_option" not in final and "retired" not in final["jobs"] and final["jobs"]["keep"] == 12
    assert final["mysql"]["gone"] == 2                                      # only what was ticked went


def test_only_a_stale_key_can_be_ticked_for_removal(config):
    before = config.read_bytes()
    result = settings_writer.update(actor="boss", remove=["jobs.keep"])
    assert not result.saved and "jobs.keep" in result.errors
    assert config.read_bytes() == before


def test_a_recommended_setting_that_is_absent_is_listed_and_written_when_ticked(config):
    fields = _fields_recommending(**{"jobs.keep": ("Finished jobs fill the disk without a limit.", 30)})
    _write(config, {"site_name": "Mine"})
    planned = settings_writer.plan_update(fields=fields)
    [item] = [r for r in planned.recommended if r.path == "jobs.keep"]
    assert item.why.startswith("Finished jobs") and item.recommended == 30 and item.has_value and not item.selected
    ticked = settings_writer.update(actor="boss", apply=["jobs.keep"], fields=fields)
    assert ticked.saved and json.loads(config.read_text())["jobs"]["keep"] == 30


def test_a_recommended_setting_left_unticked_gets_its_default_not_the_recommendation(config):
    fields = _fields_recommending(**{"jobs.keep": ("Reason.", 30)})
    _write(config, {"site_name": "Mine"})
    settings_writer.update(actor="boss", fields=fields)
    assert json.loads(config.read_text())["jobs"]["keep"] == settings.Settings().jobs.keep


def test_a_recommended_setting_the_site_already_has_a_value_for_is_left_alone(config):
    fields = _fields_recommending(**{"jobs.keep": ("Reason.", 30)})
    assert json.loads(config.read_text())["jobs"]["keep"] == 12
    planned = settings_writer.plan_update(fields=fields, apply=[])
    assert "jobs.keep" not in [r.path for r in planned.recommended]
    refused = settings_writer.update(actor="boss", apply=["jobs.keep"], fields=fields)
    assert not refused.saved and "jobs.keep" in refused.errors
    assert json.loads(config.read_text())["jobs"]["keep"] == 12


def test_a_recommended_setting_with_the_default_value_still_counts_as_set(config):
    fields = _fields_recommending(**{"jobs.keep": ("Reason.", 30)})
    default = settings.Settings().jobs.keep
    _write(config, {"jobs": {"keep": default}})
    assert "jobs.keep" not in [r.path for r in settings_writer.plan_update(fields=fields).recommended]


def test_a_setting_only_the_site_can_fill_is_listed_with_its_reason_and_no_value(config):
    _write(config, {"site_name": "Mine"})
    planned = settings_writer.plan_update()
    [item] = planned.recommended
    assert item.path == "base_url" and not item.has_value and item.why
    refused = settings_writer.update(actor="boss", apply=["base_url"])
    assert not refused.saved and "base_url" in refused.errors


def test_a_failed_write_leaves_the_old_file(config, monkeypatch):
    before = config.read_bytes()

    def boom(path, data):
        raise settings_writer.SettingsWriteError("disk full")

    monkeypatch.setattr(settings_writer, "_atomic_write", boom)
    with pytest.raises(settings_writer.SettingsWriteError):
        settings_writer.update(actor="boss")
    assert config.read_bytes() == before


def test_a_file_changed_since_the_page_was_drawn_is_refused(config):
    seen = settings_writer.plan_update().file_hash
    settings_writer.save({"jobs.keep": 3}, actor="someone")
    with pytest.raises(settings_writer.StaleSettingsError):
        settings_writer.update(actor="boss", expect_hash=seen)
    assert json.loads(config.read_text())["jobs"]["keep"] == 3


def test_a_section_that_is_not_an_object_is_an_error_not_an_overwrite(config):
    _write(config, {"jobs": 5})
    before = config.read_bytes()
    result = settings_writer.update(actor="boss")
    assert not result.saved and "jobs" in result.errors
    assert config.read_bytes() == before


def test_a_missing_or_damaged_file_is_explained(config):
    config.write_text("{not json")
    assert "not valid JSON" in settings_writer.plan_update().errors[""]
    config.unlink()
    assert "does not exist" in settings_writer.plan_update().errors[""]


def test_the_update_is_audited_and_the_backup_holds_the_old_file(config):
    rows = []
    result = settings_writer.update(actor="boss", audit=lambda action, target, detail: rows.append((action, target)))
    assert {action for action, _t in rows} == {settings_writer.UPDATE_ACTION}
    assert len(rows) == len(result.changes)
    assert result.backup and result.backup.endswith(".tar.bz2")


def test_a_secret_added_is_never_shown_and_a_stale_key_shows_only_its_name(config):
    _write(config, {"old_secret": "hunter2"})
    planned = settings_writer.plan_update()
    secret = next(change for change in planned.added if change.path == "mysql.password")
    assert secret.secret and secret.shown("new") == settings_writer.SECRET_MARK
    assert planned.stale == ["old_secret"]
    assert "hunter2" not in json.dumps([vars(change) for change in planned.removed], default=str)


# --- The API and the page --------------------------------------------------------------

def test_the_update_api_needs_an_admin(site, config):
    anonymous = site.application.test_client()
    assert anonymous.post("/api/admin/settings/update/plan", json={}).status_code in (401, 403)
    assert anonymous.post("/api/admin/settings/update", json={}).status_code in (401, 403)


def test_the_update_api_plans_without_writing_and_then_writes_with_a_backup_and_audit_rows(site, config):
    before = config.read_bytes()
    planned = site.post("/api/admin/settings/update/plan", json={}).get_json()
    assert planned["ok"] and not planned["saved"] and planned["added"] and planned["file_hash"]
    assert {r["path"] for r in planned["recommended"]} == {"base_url"}
    assert config.read_bytes() == before
    saved = site.post("/api/admin/settings/update", json={"expect_hash": planned["file_hash"]}).get_json()
    assert saved["saved"] and saved["backup"]
    document = json.loads(config.read_text())
    assert document["mysql"]["password"] == "db-secret" and document["site_name"] == "Mine"
    rows = _audit_rows(site)
    assert rows and {action for action, _t, _d in rows} == {settings_writer.UPDATE_ACTION}
    assert "db-secret" not in json.dumps(planned) + json.dumps(saved) + json.dumps(rows)


def test_the_update_api_refuses_a_stale_page_and_unknown_choices(site, config):
    seen = site.post("/api/admin/settings/update/plan", json={}).get_json()["file_hash"]
    settings_writer.save({"jobs.keep": 3}, actor="someone")
    assert site.post("/api/admin/settings/update", json={"expect_hash": seen}).status_code == 409
    bad = site.post("/api/admin/settings/update", json={"remove": ["jobs.keep"]}).get_json()
    assert not bad["saved"] and "jobs.keep" in bad["errors"]
    assert site.post("/api/admin/settings/update", json={"remove": "x"}).status_code == 400
    assert site.post("/api/admin/settings/update", json={"extra": 1}).status_code == 400


def test_the_settings_page_has_the_button_and_the_update_page_previews_then_writes(site, config):
    assert "Update config.json" in site.get("/admin/settings").get_data(as_text=True)
    _write(config, {**START, "old_option": 1})
    page = site.get("/admin/settings/update").get_data(as_text=True)
    assert "Recommended, not set" in page and "Settings to add" in page and "Stale keys" in page
    assert "old_option" in page and "db-secret" not in page and "csrf-secret" not in page
    assert json.loads(config.read_text())["old_option"] == 1                # a preview writes nothing
    fields = {name: value for name, value in re.findall(r'<input type="hidden" name="([^"]+)" value="([^"]*)"', page)}
    from planetgen.web import csrf as webcsrf
    data = {webcsrf.FIELD_NAME: site.csrf, "action": "save", "expect_hash": fields["expect_hash"], "remove": "old_option"}
    resp = site.post("/admin/settings/update", data=data)
    assert resp.status_code == 303 and resp.headers["Location"].endswith("/admin/settings")
    document = json.loads(config.read_text())
    assert "old_option" not in document and document["mysql"]["host"] == "db.internal" and "page_cache" in document
    assert "to config.json" in site.get("/admin/settings").get_data(as_text=True)


def test_the_update_page_says_so_when_there_is_nothing_to_add(site, config):
    settings_writer.update(actor="boss")
    document = json.loads(config.read_text())
    document["base_url"] = "http://example.org/"
    _write(config, document)
    assert "nothing to update" in site.get("/admin/settings/update").get_data(as_text=True).lower()


def test_the_update_page_needs_an_admin(site):
    assert site.application.test_client().get("/admin/settings/update").status_code == 302


def test_the_update_page_is_in_the_admin_search():
    from planetgen.web import adminregistry
    assert any(page.endpoint == "admin_settings_update" for page in adminregistry.PAGES)


def test_a_recommended_value_needs_the_reason_it_is_recommended():
    with pytest.raises(ValueError):
        settings.opt(1, "x", label="x", summary="x.", group="site", recommended=2)
    marked = settings.opt(1, "x", label="x", summary="x.", group="site", should_set="Why.", recommended=2)
    assert marked.json_schema_extra["x-should-set"] == "Why." and marked.json_schema_extra["x-recommended"] == 2
    plain = settings.opt(1, "x", label="x", summary="x.", group="site")
    assert "x-should-set" not in plain.json_schema_extra and "x-recommended" not in plain.json_schema_extra
