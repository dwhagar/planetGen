# tests/test_web_settings.py

"""
The Admin settings page (ADM.84): `planetgen.util.settings_view` (what the page lists), the `/api/admin/settings`
routes and the `/admin/settings` page with its review step, through the real app, API and control database,
over a temporary `config.json`.

The browser test of the page (editing, saving and restoring one option of each kind) is in test_web_a11y.py's
neighbour, test_web_settings_browser.py.

Run with: pytest tests/test_web_settings.py
"""

import json
import re

import pytest

from planetgen.admin import auth as adminAuth
from planetgen.api.authz import SESSION_COOKIE_NAME
from planetgen.api.config import Config
from planetgen.db import store
from planetgen.util import logpaths, settings, settings_view, settings_writer
from planetgen.web import csrf, settings_page
from planetgen.web.app import create_app

NONCE = "n" * 43
START = {"mysql": {"host": "db.internal", "password": "db-secret"}, "secret_key": "csrf-secret",
         "site_name": "Mine", "jobs": {"keep": 12}}


@pytest.fixture
def config(tmp_path, monkeypatch):
    path = tmp_path / "config.json"
    path.write_text(json.dumps(START, indent=2) + "\n")
    monkeypatch.setattr(logpaths, "CONFIG_PATH", str(path))
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(tmp_path / "state"))
    for variable in settings.env_names():
        monkeypatch.delenv(variable, raising=False)
    settings.reset_cache()
    yield path
    settings.reset_cache()
    settings_view.forget_startup()


def _options(role="admin"):
    return {option["path"]: option for option in settings_view.describe(role)["options"]}


# --- What the page lists ---------------------------------------------------------------

def test_every_option_of_the_model_is_listed_with_its_help(config):
    view = settings_view.describe("admin")
    listed = [option["path"] for option in view["options"]]
    assert listed == [path for path, _info, _meta in settings.iter_fields()]
    assert all(option["description"] and option["title"] and option["summary"]
               and option["kind"] in settings_view.KINDS and option["control"] in settings_view.CONTROL_NAMES
               for option in view["options"])
    assert [group["key"] for group in view["groups"]] == [key for key in settings_view.GROUPS
                                                          if any(o["group"] == key for o in view["options"])]
    assert [option["path"] for group in view["groups"] for option in group["options"]].count("site_name") == 1
    assert sorted(o["path"] for g in view["groups"] for o in g["options"]) == sorted(listed)


def test_the_kinds_limits_and_choices_come_from_the_model(config):
    options = _options("host")
    assert options["page_cache.enabled"]["kind"] == "bool"
    assert options["jobs.keep"]["kind"] == "int" and options["jobs.keep"]["limits"] == {"ge": 1}
    assert options["cache.max_mb"]["kind"] == "float" and options["cache.max_mb"]["unit"] == "MB"
    assert options["login_allowlist"]["kind"] == "list"
    assert options["log_rotation"]["kind"] == "choice" and options["log_rotation"]["choices"] == [
        "auto", "system", "app"]
    assert options["settings.backups_kept"]["limits"] == {"ge": 1, "le": 200}
    assert options["site_name"]["limits"] == {"min_length": 1}


def test_the_value_is_the_files_and_the_source_says_where_it_comes_from(config, monkeypatch):
    options = _options()
    assert options["site_name"]["value"] == "Mine" and options["site_name"]["source"] == "config.json"
    assert options["cache.max_mb"]["value"] == 1024 and options["cache.max_mb"]["source"] == "default"
    monkeypatch.setenv("PLANETGEN_CACHE_MAX_MB", "50")
    options = _options()
    assert options["cache.max_mb"]["source"] == "environment"
    assert options["cache.max_mb"]["env"] == "PLANETGEN_CACHE_MAX_MB"
    assert options["cache.max_mb"]["effective"] == 50 and options["cache.max_mb"]["value"] == 1024


def test_who_may_change_what(config):
    admin, owner, host = _options("admin"), _options("owner"), _options("host")
    assert admin["site_name"]["allowed"] and admin["jobs.keep"]["allowed"]
    assert not admin["mysql.password"]["allowed"] and "on the server" in admin["mysql.password"]["reason"]
    assert not admin["base_url"]["allowed"] and admin["base_url"]["reason"] == "can only be changed by the Owner"
    assert owner["base_url"]["allowed"] and not owner["mysql.password"]["allowed"]
    assert all(option["allowed"] for option in host.values())


def test_a_secret_is_never_in_the_view(config):
    options = _options("host")
    for option in options.values():
        if option["secret"]:
            assert option["value"] is None and option["default"] is None and option["effective"] is None
    assert options["secret_key"]["is_set"] is True and options["mysql.password"]["is_set"] is True
    text = json.dumps(settings_view.describe("host"))
    assert "csrf-secret" not in text and "db-secret" not in text


def test_a_file_that_does_not_load_shows_the_defaults_and_says_so(config):
    config.write_text(json.dumps({"jobs": {"keep": 0}}))
    view = settings_view.describe()
    assert "do not fit" in view["file_error"]
    assert _options()["jobs.keep"]["value"] == 20


def test_a_missing_file_is_a_problem_not_a_crash(config):
    config.unlink()
    view = settings_view.describe()
    assert "does not exist" in view["problem"] and view["file_hash"] == ""


# --- Options waiting for a restart -----------------------------------------------------

def test_a_saved_restart_option_is_pending_until_the_process_restarts(config):
    settings_view.remember_startup()
    assert settings_view.pending_restart() == []
    settings_writer.save({"page_cache.enabled": False, "jobs.keep": 5}, actor="boss")
    assert settings_view.pending_restart() == [("page_cache.enabled", "Page answer cache")]
    settings_writer.save({"page_cache.enabled": settings_writer.RESET}, actor="boss")
    assert settings_view.pending_restart() == []     # back to what the process started with


def test_nothing_is_pending_before_the_start_was_recorded_or_when_the_file_is_broken(config):
    settings_view.forget_startup()
    assert settings_view.pending_restart() == []
    settings_view.remember_startup()
    config.write_text(json.dumps({"page_cache": {"max_entries": 0}}))
    assert settings_view.pending_restart() == []


# --- Reading a submitted form ----------------------------------------------------------

class Form(dict):
    """`request.form` for `collect`: values are strings or lists of them."""

    def get(self, key, default=None):
        value = super().get(key, default)
        return value[-1] if isinstance(value, list) else value

    def getlist(self, key):
        value = super().get(key, [])
        return value if isinstance(value, list) else [value]


def _form(options, **fields):
    """The form the page draws for the file as it is, with `fields` (dotted names written `__`) changed."""
    form = Form()
    for path, option in options.items():
        if not option["allowed"]:
            continue
        if option["kind"] == "bool":
            form[settings_page.FIELD + path] = ["0"] + (["1"] if option["value"] else [])
        elif option["secret"]:
            form[settings_page.FIELD + path] = ""
        elif option["kind"] == "list":
            form[settings_page.FIELD + path] = "\n".join(option["value"])
        else:
            form[settings_page.FIELD + path] = settings_page._text(option["value"])
    for key, value in fields.items():
        form[key.replace("__", ".")] = value
    return form


def test_an_unchanged_form_asks_for_nothing(config):
    options = _options("host")
    changes, _submitted = settings_page.collect(list(options.values()), _form(options))
    assert changes == {}


def test_only_the_options_that_differ_are_sent(config):
    options = _options()
    form = _form(options, **{"opt__site_name": "Home", "opt__cache__max_mb": "300",
                             "opt__page_cache__enabled": ["0", "1"][:1], "opt__jobs__keep": "12"})
    changes, submitted = settings_page.collect(list(options.values()), form)
    assert changes == {"site_name": "Home", "cache.max_mb": 300.0, "page_cache.enabled": False}
    assert submitted["jobs.keep"] == "12" and submitted["site_name"] == "Home"


def test_text_that_is_not_a_number_goes_to_the_model_to_be_worded(config):
    options = _options()
    changes, _ = settings_page.collect(list(options.values()), _form(options, **{"opt__jobs__keep": "many"}))
    assert changes == {"jobs.keep": "many"}
    assert settings_writer.plan(changes).errors["jobs.keep"]


def test_use_the_default_removes_a_key_only_when_the_file_has_it(config):
    options = _options()
    form = _form(options, **{"reset__site_name": "1", "reset__cache__max_mb": "1"})
    changes, submitted = settings_page.collect(list(options.values()), form)
    assert changes == {"site_name": {"reset": True}}     # max_mb was already the default
    assert submitted["site_name"] == "planetGen"


def test_a_secret_box_left_empty_changes_nothing_and_clear_removes_it(config):
    secret = {"path": "x.token", "kind": "text", "secret": True, "allowed": True, "source": "config.json",
              "is_set": True, "value": None, "default": None}
    assert settings_page.collect([secret], Form({"opt.x.token": ""}))[0] == {}
    assert settings_page.collect([secret], Form({"opt.x.token": "new"}))[0] == {"x.token": "new"}
    assert settings_page.collect([secret], Form({"opt.x.token": "", "clear.x.token": "1"}))[0] == {
        "x.token": {"reset": True}}


def test_a_list_is_read_one_item_per_line_and_a_field_the_role_may_not_change_is_ignored(config):
    options = _options("owner")
    form = _form(options, **{"opt__login_allowlist": "10.0.0.0/8\r\n\r\n 192.168.1.5 \n",
                             "opt__mysql__password": "evil"})
    changes, _ = settings_page.collect(list(_options("admin").values()), form)
    assert "mysql.password" not in changes and "login_allowlist" not in changes   # Owner-only for an admin
    changes, _ = settings_page.collect(list(options.values()), form)
    assert changes == {"login_allowlist": ["10.0.0.0/8", "192.168.1.5"]}


# --- The routes and the page -----------------------------------------------------------

@pytest.fixture
def site(mysql_config, config):
    class RealConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "test-secret"

    adminAuth.bootstrap_control_schema(mysql_config)
    conn = store.get_control_connection(mysql_config)
    try:
        conn.execute("UPDATE admin_users SET must_change_credentials = 0")
        conn.commit()
        admin_id = conn.execute("SELECT id FROM admin_users ORDER BY id LIMIT 1").fetchone()["id"]
        session = adminAuth.create_session(conn, admin_id)
    finally:
        conn.close()
    app = create_app(RealConfig)
    app.testing = True
    client = app.test_client()
    client.set_cookie(csrf.COOKIE_NAME, NONCE)
    with app.app_context():
        token = csrf._sign(NONCE, session)
    client.set_cookie(SESSION_COOKIE_NAME, session)
    client.csrf, client.session_token, client.config = token, session, mysql_config
    return client


def _post(client, **fields):
    data = {csrf.FIELD_NAME: client.csrf}
    data.update({key.replace("__", "."): value for key, value in fields.items()})
    return client.post("/admin/settings", data=data)


def _audit_rows(client):
    conn = store.get_control_connection(client.config)
    try:
        return [(r["action"], r["target"], json.loads(r["detail"])) for r in conn.execute(
            "SELECT action, target, detail FROM admin_audit_log WHERE action LIKE %s ORDER BY id", ("config.%",)).fetchall()]
    finally:
        conn.close()


def test_the_api_needs_an_admin(site, config):
    anonymous = site.application.test_client()
    assert anonymous.get("/api/admin/settings").status_code == 401
    assert anonymous.post("/api/admin/settings", json={"changes": {}}).status_code in (401, 403)


def test_the_api_lists_plans_and_saves(site, config):
    listing = site.get("/api/admin/settings").get_json()
    assert listing["file_hash"] and listing["problem"] == "" and listing["pending_restart"] == []
    assert len(listing["options"]) == len(list(settings.iter_fields()))

    planned = site.post("/api/admin/settings/plan", json={"changes": {"jobs.keep": 9}})
    assert planned.status_code == 200
    body = planned.get_json()
    assert body["ok"] and not body["saved"]
    assert body["changes"] == [{"path": "jobs.keep", "old": 12, "new": 9, "secret": False, "restart": False,
                                "reset": False}]
    assert json.loads(config.read_text())["jobs"] == {"keep": 12}        # a plan writes nothing


def test_the_api_reports_each_bad_value_at_its_option(site):
    body = site.post("/api/admin/settings", json={"changes": {"jobs.keep": 0, "cache.max_mb": "lots"}}).get_json()
    assert not body["ok"] and not body["saved"]
    assert set(body["errors"]) == {"jobs.keep", "cache.max_mb"}


def test_the_api_refuses_what_an_admin_may_not_change(site, config):
    body = site.post("/api/admin/settings", json={"changes": {"mysql.password": "x", "base_url": "http://a/"}}).get_json()
    assert not body["saved"] and set(body["errors"]) == {"mysql.password", "base_url"}
    assert json.loads(config.read_text())["mysql"]["password"] == "db-secret"


def test_the_api_save_backs_up_audits_and_masks(site, config):
    saved = site.post("/api/admin/settings", json={"changes": {"site_name": "Home", "page_cache.enabled": False},
                                                   "expect_hash": site.get("/api/admin/settings").get_json()[
                                                       "file_hash"]}).get_json()
    assert saved["saved"] and saved["backup"] and saved["restart"] == ["page_cache.enabled"]
    document = json.loads(config.read_text())
    assert document["site_name"] == "Home" and document["page_cache"] == {"enabled": False}
    rows = _audit_rows(site)
    assert sorted((a, t) for a, t, _d in rows) == [("config.update", "config:page_cache.enabled"),
                                                   ("config.update", "config:site_name")]
    assert ("config.update", "config:site_name", {"old": "Mine", "new": "Home"}) in rows
    assert site.get("/api/admin/settings").get_json()["pending_restart"] == [
        {"path": "page_cache.enabled", "title": "Page answer cache"}]


def test_the_api_refuses_a_save_from_a_stale_form(site, config):
    stale = site.get("/api/admin/settings").get_json()["file_hash"]
    settings_writer.save({"jobs.keep": 3}, actor="someone")
    resp = site.post("/api/admin/settings", json={"changes": {"jobs.keep": 4}, "expect_hash": stale})
    assert resp.status_code == 409 and "changed by someone else" in resp.get_json()["error"]
    assert json.loads(config.read_text())["jobs"] == {"keep": 3}


def test_the_api_reports_a_file_it_cannot_write(site, config):
    config.unlink()
    resp = site.post("/api/admin/settings", json={"changes": {"jobs.keep": 4}})
    assert resp.status_code == 200 and not resp.get_json()["saved"] and "does not exist" in resp.get_json()["errors"][""]
    assert "does not exist" in site.get("/api/admin/settings").get_json()["problem"]


def test_the_api_takes_reset_as_an_object_and_refuses_other_objects(site, config):
    body = site.post("/api/admin/settings", json={"changes": {"site_name": {"reset": True}}}).get_json()
    assert body["saved"] and "site_name" not in json.loads(config.read_text())
    assert site.post("/api/admin/settings", json={"changes": {"site_name": {"x": 1}}}).status_code == 400
    assert site.post("/api/admin/settings", json={"changes": {}, "extra": 1}).status_code == 400


def test_the_page_lists_every_option_and_hides_every_secret(site, config):
    html = site.get("/admin/settings").get_data(as_text=True)
    for path, _info, _meta in settings.iter_fields():
        assert path in html
    assert "csrf-secret" not in html and "db-secret" not in html
    assert 'name="opt.site_name"' in html and 'value="Mine"' in html
    assert 'name="opt.mysql.password"' not in html and "on the server" in html      # listed read-only, with the reason
    assert 'aria-current="page"' in html


def test_the_settings_an_admin_cannot_change_are_hidden_until_asked_for(site, config):
    view = settings_view.describe("admin")
    fixed = [o["path"] for o in view["options"] if not o["allowed"]]
    assert view["fixed_count"] == len(fixed) and 0 < len(fixed) < len(view["options"])
    hidden = site.get("/admin/settings").get_data(as_text=True)
    shown = site.get("/admin/settings?show=fixed").get_data(as_text=True)
    assert 'data-show-fixed="false"' in hidden and 'data-show-fixed="true"' in shown
    assert f"Show settings that cannot be changed ({len(fixed)})" in hidden
    assert 'action="/admin/settings?show=fixed"' in shown and 'action="/admin/settings"' in hidden


def test_each_setting_gets_the_control_that_fits_it(site, config):
    html = site.get("/admin/settings").get_data(as_text=True)
    assert 'role="switch" id="opt-seo-json_ld"' in html                              # yes or no
    assert 'name="opt.seo.indexing" value="top_pages" checked' in html               # few choices: segments
    assert "Main pages only" in html and "Keep out of search" in html
    assert 'data-limit-slider' in html and 'id="opt-cache-max_mb"' in html           # a size: slider with a box
    assert 'type="color"' in html and 'name="rate.ratelimit.pages.search.count"' in html
    assert "Change it on the Images page" in html                                    # an uploaded file
    assert "Details" in html and "Advanced (" in html


def test_a_search_keeps_only_the_settings_it_matches(site, config):
    html = site.get("/admin/settings?q=sitemap").get_data(as_text=True)
    assert 'id="setting-seo-sitemap_enabled"' in html and 'id="setting-site_name"' not in html
    assert "No setting matches" in site.get("/admin/settings?q=zzzz-nothing").get_data(as_text=True)


@pytest.mark.parametrize("text, expected", [("30 per minute", ("30", "minute")), ("", ("", "minute")),
                                             ("5 per DAY", ("5", "day")), ("200 per day;50 per hour", None)])
def test_a_rate_is_read_as_a_number_and_a_period_or_left_to_text(text, expected):
    found = settings_page.parse_rate(text)
    assert (found and (found["limit"], found["period"])) == expected


def test_a_rate_control_is_composed_back_into_the_limit_text(config):
    options = _options()
    form = {settings_page.RATE_FIELD + "ratelimit.pages.search.count": "90",
            settings_page.RATE_FIELD + "ratelimit.pages.search.per": "hour"}
    changes, _ = settings_page.collect(list(options.values()), form)
    assert changes == {"ratelimit.pages.search": "90 per hour"}
    form[settings_page.RATE_FIELD + "ratelimit.pages.search.off"] = "1"
    assert settings_page.collect(list(options.values()), form)[0] == {"ratelimit.pages.search": ""}


@pytest.mark.parametrize("value, unit, words", [("21600", "seconds", "6 hours"), ("15", "seconds", "15 seconds"),
                                                 ("90000", "seconds", "1 day 1 hour"), ("24", "hours", "1 day"),
                                                 ("0", "seconds", "off"), ("soon", "seconds", "")])
def test_a_duration_is_said_in_words(value, unit, words):
    assert settings_page.duration_text(value, unit) == words


def test_the_page_needs_an_admin(site):
    assert site.application.test_client().get("/admin/settings").status_code == 302


def test_the_review_step_shows_the_diff_and_writes_nothing(site, config):
    resp = _post(site, opt__site_name="Home", opt__page_cache__enabled=["0"], expect_hash=site.get(
        "/api/admin/settings").get_json()["file_hash"], action="review")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 200 and "These changes will be saved" in html
    assert "Mine" in html and "Home" in html and "applies after a restart" in html
    assert json.loads(config.read_text())["site_name"] == "Mine"
    assert re.search(r'name="changes" value="[^"]*Home', html)


def test_confirming_saves_and_the_next_page_says_so(site, config):
    review = _post(site, opt__site_name="Home", action="review").get_data(as_text=True)
    hidden = {name: value for name, value in re.findall(r'<input type="hidden" name="([^"]+)" value="([^"]*)"', review)}
    import html as htmllib
    resp = _post(site, action="save", changes=htmllib.unescape(hidden["changes"]), expect_hash=hidden["expect_hash"])
    assert resp.status_code == 303 and resp.headers["Location"].endswith("/admin/settings")
    assert json.loads(config.read_text())["site_name"] == "Home"
    page = site.get("/admin/settings").get_data(as_text=True)
    assert "Saved 1 setting." in page and "backup" in page
    assert "Saved 1 setting." not in site.get("/admin/settings").get_data(as_text=True)     # shown once


def test_a_bad_value_shows_its_message_by_the_option_and_changes_nothing(site, config):
    resp = _post(site, opt__jobs__keep="0", opt__site_name="Home")
    html = resp.get_data(as_text=True)
    assert resp.status_code == 422
    assert "Nothing was saved" in html and 'id="err-jobs-keep"' in html and 'aria-invalid="true"' in html
    assert 'value="0"' in html and 'value="Home"' in html                          # the submitted text is kept
    assert json.loads(config.read_text())["jobs"] == {"keep": 12}


def test_an_unchanged_form_says_there_is_nothing_to_change(site, config):
    before = config.read_bytes()
    resp = _post(site, opt__site_name="Mine")
    assert resp.status_code == 200 and "Nothing to change" in resp.get_data(as_text=True)
    assert config.read_bytes() == before


def test_a_form_from_before_someone_elses_save_is_refused(site, config):
    stale = site.get("/api/admin/settings").get_json()["file_hash"]
    settings_writer.save({"jobs.keep": 3}, actor="someone")
    resp = _post(site, opt__site_name="Home", expect_hash=stale)
    assert resp.status_code == 409 and "changed by someone else" in resp.get_data(as_text=True)
    assert json.loads(config.read_text())["site_name"] == "Mine"


def test_a_confirmation_for_a_file_changed_meanwhile_is_refused(site, config):
    review = _post(site, opt__site_name="Home", action="review").get_data(as_text=True)
    import html as htmllib
    hidden = {name: htmllib.unescape(value) for name, value in
              re.findall(r'<input type="hidden" name="([^"]+)" value="([^"]*)"', review)}
    settings_writer.save({"jobs.keep": 3}, actor="someone")
    resp = _post(site, action="save", changes=hidden["changes"], expect_hash=hidden["expect_hash"])
    assert resp.status_code == 409
    assert json.loads(config.read_text())["site_name"] == "Mine"


def test_a_missing_file_is_explained_and_has_no_save_button(site, config):
    config.unlink()
    html = site.get("/admin/settings").get_data(as_text=True)
    assert "does not exist" in html and "Review the changes" not in html
    assert "site_name" in html


def test_the_restart_banner_shows_on_every_admin_page_until_the_restart(site, config):
    assert "waiting for a restart" not in site.get("/admin").get_data(as_text=True)
    settings_writer.save({"page_cache.enabled": False}, actor="boss")
    for path in ("/admin", "/admin/settings", "/admin/queue"):
        assert "waiting for a restart of the web server: Page answer cache" in site.get(path).get_data(as_text=True)
    assert "waiting for a restart" not in site.application.test_client().get("/").get_data(as_text=True)


def test_the_settings_tab_is_on_the_admin_pages(site):
    html = site.get("/admin").get_data(as_text=True)
    assert 'href="/admin/settings"' in html

