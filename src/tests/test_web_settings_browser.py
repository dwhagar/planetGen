# tests/test_web_settings_browser.py

"""
The Admin settings page (ADM.84) in a real browser: the form passes the same accessibility and phone-fit checks as
every other page (test_web_a11y.py only sees the page without a config.json, which has no form), and an admin can
edit, review, save and put back an option of each kind.

An admin may change only some options (`x-editable`, `x-min_role`), so the page is drawn twice: as the admin the
web has today (text, whole number, number and on/off options), and as the command line's role (`host`), which
adds the list, the choice and the secret. Needs Playwright with Chromium and the MySQL test server, like
test_web_a11y.py; skipped without them.

Run with: pytest tests/test_web_settings_browser.py
"""

import json

import pytest

sync_api = pytest.importorskip("playwright.sync_api")


from planetgen.admin import auth  # noqa: E402
from planetgen.api import settings as settings_api  # noqa: E402
from planetgen.api.authz import SESSION_COOKIE_NAME  # noqa: E402
from planetgen.api.config import Config  # noqa: E402
from planetgen.api.limiter import PAGE_LIMITS_OFF  # noqa: E402
from planetgen.db import store  # noqa: E402
from planetgen.util import logpaths, settings, settings_view  # noqa: E402
from planetgen.web.app import create_app  # noqa: E402

from tests.served import served  # noqa: E402
from tests.test_web_a11y import (  # noqa: E402,F401 -- the browser and axe fixtures, and the checks
    CSP_LISTENER, SCHEMES, VIEWPORTS, _axe_blocking, _overflow, axe_source, browser,
)

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


@pytest.fixture
def server(mysql_config, config):
    class SiteConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "settings-browser-secret"
        RATELIMIT_PAGES = PAGE_LIMITS_OFF

    auth.bootstrap_control_schema(mysql_config)
    conn = store.get_control_connection(mysql_config)
    try:
        conn.execute("UPDATE admin_users SET must_change_credentials = 0")
        conn.commit()
        admin_id = conn.execute("SELECT id FROM admin_users ORDER BY id LIMIT 1").fetchone()["id"]
        session = auth.create_session(conn, admin_id)
    finally:
        conn.close()
    app = create_app(SiteConfig)
    app.testing = True
    with served(app) as url:
        yield url, session


@pytest.fixture(params=["admin", "host"])
def role(request, monkeypatch):
    monkeypatch.setattr(settings_api, "ROLE", request.param)
    return request.param


def _page(browser, server, viewport="desktop", scheme="light"):
    base_url, session = server
    context = browser.new_context(viewport=VIEWPORTS[viewport], color_scheme=scheme, device_scale_factor=1,
                                  reduced_motion="reduce")
    context.add_cookies([{"name": SESSION_COOKIE_NAME, "value": session, "url": base_url}])
    context.add_init_script(CSP_LISTENER)
    page = context.new_page()
    problems = []
    # A refused form is a 422 page, which the browser logs; a request cut off by the next click is not a failure.
    page.on("console", lambda msg: msg.type == "error" and "status of 422" not in msg.text
            and problems.append(f"console error: {msg.text}"))
    page.on("pageerror", lambda exc: problems.append(f"uncaught exception: {exc}"))
    page.on("requestfailed", lambda req: "ERR_ABORTED" in str(req.failure)
            or problems.append(f"request failed: {req.url} ({req.failure})"))
    return context, page, problems, base_url


def _document(config):
    return json.loads(config.read_text())


@pytest.mark.parametrize("scheme", SCHEMES)
@pytest.mark.parametrize("viewport", sorted(VIEWPORTS))
def test_the_form_is_accessible_and_fits(browser, axe_source, server, role, config, viewport, scheme):
    context, page, problems, base_url = _page(browser, server, viewport, scheme)
    try:
        response = page.goto(base_url + "/admin/settings", wait_until="load")
        assert response.status == 200
        assert page.locator("form.settings-form").count() == 1
        assert page.locator("button:has-text('Review the changes')").count() == 1
        problems += _overflow(page) + _axe_blocking(page, axe_source)
        problems += [f"CSP violation: {v}" for v in page.evaluate("window.__cspViolations || []")]
    finally:
        context.close()
    assert not problems, f"/admin/settings as {role} at {viewport}/{scheme}:\n  " + "\n  ".join(problems)


def test_the_review_and_error_pages_are_accessible(browser, axe_source, server, role, config):
    context, page, problems, base_url = _page(browser, server, "phone", "dark")
    try:
        page.goto(base_url + "/admin/settings", wait_until="load")
        page.fill("#opt-site_name", "Home")
        page.click("button:has-text('Review the changes')")
        page.wait_for_selector("table.settings-diff")
        problems += _overflow(page) + _axe_blocking(page, axe_source)
        page.click("a:has-text('Go back without saving')")
        page.fill("#opt-jobs-keep", "0")
        page.click("button:has-text('Review the changes')")
        page.wait_for_selector("#err-jobs-keep")
        problems += _overflow(page) + _axe_blocking(page, axe_source)
    finally:
        context.close()
    assert not problems, f"as {role}:\n  " + "\n  ".join(problems)
    assert _document(config)["site_name"] == "Mine"


def _save(page, base_url, expect_rows):
    page.click("button:has-text('Review the changes')")
    page.wait_for_selector("table.settings-diff")
    assert page.locator("table.settings-diff tbody tr").count() == expect_rows
    page.click("button:has-text('Save these changes')")
    page.wait_for_selector("text=Saved ")


def test_edit_save_and_put_back_each_kind_of_option(browser, server, role, config):
    context, page, problems, base_url = _page(browser, server)
    try:
        page.goto(base_url + "/admin/settings", wait_until="load")
        page.fill("#opt-site_name", "Home")                       # text
        page.fill("#opt-jobs-keep", "7")                          # whole number
        page.fill("#opt-cache-max_mb", "300")                # number
        page.uncheck("#opt-page_cache-stored_counts")             # on/off
        rows = 4
        if role == "host":
            page.fill("#opt-login_allowlist", "10.0.0.0/8\n192.168.1.5")        # list
            page.check("input[name='opt.log_rotation'][value='app']")      # choice
            page.fill("#opt-wiki-wikijs-api_token", "wiki-token-123")           # secret
            rows += 3
        _save(page, base_url, rows)

        saved = _document(config)
        assert saved["site_name"] == "Home" and saved["jobs"] == {"keep": 7}
        assert saved["cache"] == {"max_mb": 300.0} and saved["page_cache"] == {"stored_counts": False}
        if role == "host":
            assert saved["login_allowlist"] == ["10.0.0.0/8", "192.168.1.5"] and saved["log_rotation"] == "app"
            assert saved["wiki"]["wikijs"]["api_token"] == "wiki-token-123"
        # What was entered is shown back, a secret is not, and nothing was left unsaved.
        assert page.input_value("#opt-site_name") == "Home" and page.input_value("#opt-jobs-keep") == "7"
        assert not page.is_checked("#opt-page_cache-stored_counts")
        assert "wiki-token-123" not in page.content()

        # Put every one of them back to its default (the keys leave the file).
        page.check("input[name='reset.cache.max_mb']")
        page.check("input[name='reset.page_cache.stored_counts']")
        page.check("input[name='reset.jobs.keep']")
        page.check("input[name='reset.site_name']")
        rows = 4
        if role == "host":
            page.check("input[name='reset.login_allowlist']")
            page.check("input[name='reset.log_rotation']")
            page.check("input[name='clear.wiki.wikijs.api_token']")
            rows += 3
        _save(page, base_url, rows)
        restored = _document(config)
        assert restored["mysql"] == START["mysql"] and restored["secret_key"] == "csrf-secret"
        for path in ("site_name", "jobs", "tile_cache", "page_cache", "login_allowlist", "log_rotation", "wiki"):
            assert path not in restored, path
        assert page.input_value("#opt-site_name") == "planetGen"
    finally:
        context.close()
    assert not problems, "\n".join(problems)


def test_the_unchangeable_switch_search_and_change_count(browser, server, config):
    context, page, problems, base_url = _page(browser, server)
    try:
        page.goto(base_url + "/admin/settings", wait_until="load")
        assert page.is_visible("#setting-site_name") and not page.is_visible("#setting-mysql-host")
        page.check("[data-show-fixed-toggle]")                       # no reload: the switch shows them at once
        assert page.is_visible("#setting-mysql-host") and "show=fixed" in page.url
        page.reload(wait_until="load")
        assert page.is_visible("#setting-mysql-host")                # remembered in the address
        page.uncheck("[data-show-fixed-toggle]")
        page.fill("#settings-q", "sitemap")
        assert page.is_visible("#setting-seo-sitemap_enabled") and not page.is_visible("#setting-site_name")
        page.fill("#settings-q", "")
        assert page.is_disabled("[data-settings-save]")
        page.fill("#opt-site_name", "Home")
        assert "1 change" in page.inner_text("[data-settings-changes]") and not page.is_disabled("[data-settings-save]")
    finally:
        context.close()
    assert not problems, problems


def test_a_bad_value_in_the_browser_is_marked_and_saves_nothing(browser, server, config):
    before = config.read_bytes()
    context, page, problems, base_url = _page(browser, server)
    try:
        page.goto(base_url + "/admin/settings", wait_until="load")
        page.fill("#opt-jobs-keep", "lots")
        page.fill("#opt-site_name", "Home")
        page.click("button:has-text('Review the changes')")
        page.wait_for_selector("#err-jobs-keep")
        assert page.get_attribute("#opt-jobs-keep", "aria-invalid") == "true"
        assert page.input_value("#opt-jobs-keep") == "lots" and page.input_value("#opt-site_name") == "Home"
        assert "Nothing was saved" in page.content()
    finally:
        context.close()
    assert config.read_bytes() == before
    assert not problems, problems


# --- Backups, restore, export and import (ADM.85) ----------------------------------------

@pytest.mark.parametrize("scheme", SCHEMES)
@pytest.mark.parametrize("viewport", sorted(VIEWPORTS))
def test_the_backup_pages_are_accessible_and_fit(browser, axe_source, server, config, viewport, scheme):
    from planetgen import backup
    from planetgen.util import settings_writer
    backup.take(who="boss", reason="manual")
    settings_writer.save({"site_name": "Home", "page_cache.enabled": False}, actor="boss")
    name = backup.list_archives()[-1].name
    context, page, problems, base_url = _page(browser, server, viewport, scheme)
    try:
        for path in ("/admin/settings/backups", f"/admin/settings/backups/{name}"):
            response = page.goto(base_url + path, wait_until="load")
            assert response.status == 200, path
            problems += [f"{path}: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]
        page.goto(base_url + "/admin/settings/backups", wait_until="load")
        page.set_input_files("#import-file", {"name": "s.json", "mimeType": "application/json",
                                              "buffer": json.dumps({"site_name": "Imported"}).encode()})
        page.click("button:has-text('Check the file')")
        page.wait_for_selector("table.settings-diff")
        problems += [f"import review: {p}" for p in _overflow(page) + _axe_blocking(page, axe_source)]
        problems += [f"CSP violation: {v}" for v in page.evaluate("window.__cspViolations || []")]
    finally:
        context.close()
    assert not problems, f"at {viewport}/{scheme}:\n  " + "\n  ".join(problems)


def test_back_up_change_compare_restore_and_import_in_the_browser(browser, server, config):
    context, page, problems, base_url = _page(browser, server)
    try:
        page.goto(base_url + "/admin/settings/backups", wait_until="load")
        page.click("button:has-text('Back up now')")
        page.wait_for_selector("text=Backed up Settings in ")

        page.goto(base_url + "/admin/settings", wait_until="load")
        page.fill("#opt-site_name", "Home")
        page.fill("#opt-jobs-keep", "7")
        _save(page, base_url, 2)
        assert _document(config)["site_name"] == "Home"

        page.goto(base_url + "/admin/settings/backups", wait_until="load")
        assert page.locator("table.backup-list tbody tr").count() == 2        # Back up now, and the save's own
        page.locator("a:has-text('Compare and restore')").last.click()        # the oldest: Back up now
        page.wait_for_selector("table.settings-diff")
        assert page.locator("table.settings-diff tbody tr").count() == 2
        page.click("button:has-text('Restore the ticked parts')")
        page.wait_for_selector("text=Restored settings from")
        assert _document(config)["site_name"] == "Mine" and _document(config)["jobs"] == {"keep": 12}
        assert page.locator("text=same as today's").count() == 1

        page.goto(base_url + "/admin/settings/backups", wait_until="load")
        page.set_input_files("#import-file", {"name": "s.json", "mimeType": "application/json",
                                              "buffer": json.dumps({"site_name": "Imported"}).encode()})
        page.click("button:has-text('Check the file')")
        page.wait_for_selector("table.settings-diff")
        page.click("button:has-text('Import these settings')")
        page.wait_for_selector("text=Imported 1 setting.")
        assert _document(config)["site_name"] == "Imported"

        page.goto(base_url + "/admin/settings/backups", wait_until="load")
        page.set_input_files("#import-file", {"name": "s.json", "mimeType": "application/json",
                                              "buffer": json.dumps({"jobs": {"keep": 0}}).encode()})
        page.click("button:has-text('Check the file')")
        page.wait_for_selector("text=Nothing was imported")
        assert _document(config)["jobs"] == {"keep": 12} or "jobs" not in _document(config)
    finally:
        context.close()
    assert not problems, "\n".join(problems)
