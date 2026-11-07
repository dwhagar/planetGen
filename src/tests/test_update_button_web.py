# tests/test_update_button_web.py

"""
SEC.43, the web side of the admin Update button: `planetgen/admin/webupdate.py` (the request file and the state view),
`/api/admin/update` (who may ask) and `/admin/update` (the page). The root side (the path unit) is
`test_update_button_units.py`.

The guards: no request without the Owner and a recent password; an API key never asks; the only thing a request does is
make the empty file `update.request`; a second request while one waits is refused; a press during an update waits for it.
"""

import os
import stat

import pytest

from planetgen.admin import auth as adminAuth
from planetgen.admin import updateguard, webupdate
from planetgen.db import store
from tests.test_owner_settings import PASSWORD, _become, _control
from tests.test_web_settings import config, site  # noqa: F401 -- the fixtures
from planetgen.web import csrf


@pytest.fixture
def units(tmp_path, monkeypatch):
    """A state folder, a unit folder with both units, and an update log, none of them the real ones."""
    state = tmp_path / "state"
    state.mkdir()
    unit_dir = tmp_path / "units"
    unit_dir.mkdir()
    (unit_dir / webupdate.PATH_UNIT).write_text("[Path]\n")
    (unit_dir / webupdate.SERVICE_UNIT).write_text("[Service]\n")
    log = tmp_path / "update.log"
    monkeypatch.setenv("PLANETGEN_STATE_DIR", str(state))
    monkeypatch.setenv(webupdate.UNIT_DIR_ENV, str(unit_dir))
    monkeypatch.setattr(webupdate, "log_path", lambda: str(log))
    monkeypatch.setattr(updateguard, "is_updating", lambda: False)
    monkeypatch.setattr(webupdate, "service_status", lambda: None)
    return type("Units", (), {"state": state, "units": unit_dir, "log": log})


# --- The request file ---------------------------------------------------------------------------------------------

def test_without_the_units_the_button_says_it_is_not_set_up_and_asks_for_nothing(units):
    (units.units / webupdate.PATH_UNIT).unlink()
    view = webupdate.view()
    assert view["state"] == "not_installed" and not view["can_request"] and "planetgen-update.path" in view["reason"]
    with pytest.raises(webupdate.UpdateRequestError):
        webupdate.request()
    assert not (units.state / webupdate.REQUEST_NAME).exists()


def test_a_request_is_one_empty_file_and_nothing_else(units):
    assert webupdate.view()["state"] == "idle"
    webupdate.request()
    made = units.state / webupdate.REQUEST_NAME
    assert made.name == "update.request" and made.stat().st_size == 0
    assert sorted(p.name for p in units.state.iterdir()) == ["update.request"]
    assert not stat.S_IMODE(made.stat().st_mode) & 0o007       # not readable or writable by others
    view = webupdate.view()
    assert view["state"] == "requested" and not view["can_request"]


def test_a_second_request_while_one_waits_is_refused(units):
    webupdate.request()
    with pytest.raises(webupdate.UpdateRequestError, match="requested"):
        webupdate.request()


def test_a_press_during_an_update_is_accepted_and_runs_after_it(units, monkeypatch):
    monkeypatch.setattr(updateguard, "is_updating", lambda: True)
    assert webupdate.view()["state"] == "running" and webupdate.view()["can_request"]
    webupdate.request()
    view = webupdate.view()
    assert view["state"] == "running" and "after the one that is running" in view["reason"]
    assert (units.state / webupdate.REQUEST_NAME).exists()


def test_a_symlink_in_the_way_is_not_followed(units, tmp_path):
    target = tmp_path / "elsewhere"
    os.symlink(target, units.state / webupdate.REQUEST_NAME)
    with pytest.raises(webupdate.UpdateRequestError):
        webupdate.request()
    assert not target.exists()


def test_an_old_request_is_called_stale(units):
    webupdate.request()
    path = units.state / webupdate.REQUEST_NAME
    old = path.stat().st_mtime - webupdate.STALE_REQUEST_SECONDS - 5
    os.utime(path, (old, old))
    assert webupdate.view()["stale_request"] is True


def test_the_log_tail_drops_the_cut_first_line_and_a_missing_log_is_empty(units):
    assert webupdate.view()["log"] == ""
    units.log.write_text("".join(f"line {n}\n" for n in range(2000)))
    text, modified = webupdate.log_tail(limit=200)
    assert text.endswith("line 1999\n") and text.startswith("line ") and len(text) <= 200 and modified
    assert webupdate.view()["log"].endswith("line 1999\n")


# --- The API and the page ------------------------------------------------------------------------------------------

def _become_owner(site, minutes_ago=0):
    _become(site, "owner", minutes_ago=minutes_ago)


def test_the_api_needs_a_signed_in_session(site, units):
    anonymous = site.application.test_client()
    assert anonymous.get("/api/admin/update").status_code == 401
    assert anonymous.post("/api/admin/update", json={}).status_code == 401


def test_an_admin_who_is_not_the_owner_is_refused(site, units):
    _become(site, "admin")
    body = site.get("/api/admin/update").get_json()
    assert body["can_ask"] is False and body["needs_password"] is False
    refused = site.post("/api/admin/update", json={})
    assert refused.status_code == 403
    assert not (units.state / webupdate.REQUEST_NAME).exists()


def test_an_owner_with_a_stale_password_is_refused_until_it_is_typed_again(site, units):
    _become_owner(site, minutes_ago=30)
    assert site.get("/api/admin/update").get_json()["needs_password"] is True
    refused = site.post("/api/admin/update", json={})
    assert refused.status_code == 403 and refused.get_json()["needs_password"] is True
    assert not (units.state / webupdate.REQUEST_NAME).exists()
    assert site.post("/api/auth/confirm-password", json={"current_password": PASSWORD}).status_code == 200
    accepted = site.post("/api/admin/update", json={})
    assert accepted.status_code == 202 and accepted.get_json()["state"] == "requested"
    assert (units.state / webupdate.REQUEST_NAME).stat().st_size == 0


def test_the_request_takes_no_field_so_nothing_can_be_passed_to_root(site, units):
    _become_owner(site)
    for body in ({"args": "--yes"}, {"cancel_running": True}, {"command": "id"}):
        assert site.post("/api/admin/update", json=body).status_code == 400
    assert not (units.state / webupdate.REQUEST_NAME).exists()


def test_a_second_press_is_a_conflict_and_the_press_is_audited(site, units):
    _become_owner(site)
    assert site.post("/api/admin/update", json={}).status_code == 202
    assert site.post("/api/admin/update", json={}).status_code == 409
    conn = _control(site)
    try:
        rows = conn.execute("SELECT action FROM admin_audit_log WHERE action = 'update.request'").fetchall()
    finally:
        conn.close()
    assert len(rows) == 1


def test_an_api_key_cannot_read_or_ask_even_the_owners(site, units):
    _become_owner(site)
    conn = _control(site)
    try:
        admin_id = conn.execute("SELECT id FROM admin_users").fetchone()["id"]
        _key_id, raw = adminAuth.create_api_key(conn, admin_id, "ci")
    finally:
        conn.close()
    keyed = site.application.test_client()
    headers = {"Authorization": f"Bearer {raw}"}
    assert keyed.get("/api/admin/update", headers=headers).status_code == 403
    assert keyed.post("/api/admin/update", json={}, headers=headers).status_code == 403
    assert not (units.state / webupdate.REQUEST_NAME).exists()


def _page_post(site, **fields):
    return site.post("/admin/update", data={csrf.FIELD_NAME: site.csrf, **fields})


def test_the_page_shows_the_state_and_an_owner_presses_update(site, units):
    _become_owner(site)
    page = site.get("/admin/update").get_data(as_text=True)
    assert "No update is waiting or running" in page and "Update the site" in page
    pressed = _page_post(site, action="update")
    assert pressed.status_code == 303
    assert (units.state / webupdate.REQUEST_NAME).exists()
    after = site.get("/admin/update").get_data(as_text=True)
    assert "Requested" in after and 'http-equiv="refresh"' in after


def test_the_page_asks_an_owner_with_a_stale_password_for_it_and_refuses_the_press(site, units):
    _become_owner(site, minutes_ago=30)
    page = site.get("/admin/update").get_data(as_text=True)
    assert 'name="current_password"' in page and "Update the site" not in page
    refused = _page_post(site, action="update")
    assert refused.status_code == 303
    assert "Updating the site is for the Owner" in site.get("/admin/update").get_data(as_text=True)
    assert not (units.state / webupdate.REQUEST_NAME).exists()
    assert _page_post(site, action="confirm-password", current_password=PASSWORD).status_code == 303
    assert "Update the site" in site.get("/admin/update").get_data(as_text=True)


def test_the_page_needs_an_admin_and_a_csrf_token(site, units):
    anonymous = site.application.test_client()
    assert anonymous.get("/admin/update").status_code == 302
    _become_owner(site)
    assert site.post("/admin/update", data={"action": "update"}).status_code == 400
    assert not (units.state / webupdate.REQUEST_NAME).exists()


def test_the_page_says_when_the_button_is_not_set_up(site, units):
    _become_owner(site)
    (units.units / webupdate.SERVICE_UNIT).unlink()
    page = site.get("/admin/update").get_data(as_text=True)
    assert "not set up on this server" in page and "Update the site" not in page


# --- UX.284: how the last update ended, and the log coloured by level ---------------------------------------------

STAMP = "2026-11-07 10:00:00 EST"
GUARD_STOP = (f"{STAMP} INFO     update: === Stage 1 of 9: Checking ===\n"
              f"{STAMP} ERROR    update: Generation is running: nothing was changed.\n")
FINISHED = f"{STAMP} INFO     update: === Update done after 4:07 ===\n{STAMP} INFO     update: No warnings.\n"


def test_a_log_that_ends_in_an_error_is_a_failed_update_with_the_error_in_the_detail(units):
    units.log.write_text(GUARD_STOP)
    view = webupdate.view()
    assert view["outcome"] == "failed" and "Generation is running" in view["outcome_detail"]


def test_a_finished_run_after_an_old_error_is_done_not_failed(units):
    units.log.write_text(GUARD_STOP + FINISHED)
    assert webupdate.view()["outcome"] == "done"
    units.log.write_text(FINISHED + GUARD_STOP)
    assert webupdate.view()["outcome"] == "failed"


def test_a_start_that_never_reached_the_log_is_failed_by_the_services_state(units, monkeypatch):
    """The unit failing before `update.log` begins (a refused request file, the start limit) leaves the log silent."""
    units.log.write_text(FINISHED)
    monkeypatch.setattr(webupdate, "service_status", lambda: {"active": "failed", "result": "start-limit-hit", "status": "0"})
    view = webupdate.view()
    assert view["outcome"] == "failed" and "reset-failed" in view["outcome_detail"]
    monkeypatch.setattr(webupdate, "service_status", lambda: {"active": "failed", "result": "exit-code", "status": "1"})
    assert "exit status 1" in webupdate.view()["outcome_detail"]


def test_no_log_and_no_failure_says_nothing(units):
    assert webupdate.view()["outcome"] is None


def test_the_page_says_the_update_failed_and_colours_the_levels(site, units):
    _become_owner(site)
    units.log.write_text(f"{STAMP} DEBUG    update: dbg\n{STAMP} INFO     update: inf\n{STAMP} WARNING  update: wrn\n"
                         f"{STAMP} ERROR    update: Generation is running.\nplain output\n")
    page = site.get("/admin/update").get_data(as_text=True)
    assert "The last update failed." in page and 'role="alert"' in page
    for level in ("debug", "info", "warning", "error"):
        assert f'<span class="log-{level}">{STAMP} {level.upper()}' in page
    assert "\nplain output\n" in page and 'data-follow-end' in page
