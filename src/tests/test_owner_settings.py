# tests/test_owner_settings.py

"""
USR.31 and ADM.86: the Owner account stub and the Owner-tier settings behind a recent password.

An admin, and an Owner whose password was typed more than ten minutes ago, are both refused an option marked
`x-min_role: owner` (the mail server, `base_url`, the From address); the Owner who typed it lately may change it; an API
key is never more than an admin; `POST /api/auth/confirm-password` reopens a stale Owner. `adminAuth.set_owner` and
`planetgen.cli.accounts --set-owner` move the Owner.
"""

import json

import pytest

from planetgen.admin import auth as adminAuth
from planetgen.cli import accounts
from planetgen.db import store
from tests.test_web_settings import _post, config, site  # noqa: F401 -- the fixtures

PASSWORD = "correct horse battery staple 42"
OWNER_CHANGE = {"base_url": "https://sky.example/"}


def _control(site):
    return store.get_control_connection(site.config)


def _become(site, role, minutes_ago=0):
    """Sets the account's role and when the session's password was typed."""
    conn = _control(site)
    try:
        conn.execute("UPDATE admin_users SET role = 'admin' WHERE role = 'owner'")
        conn.execute("UPDATE admin_users SET role = ?, password_hash = ?", (role, adminAuth.hash_password(PASSWORD)))
        conn.execute("UPDATE admin_sessions SET password_confirmed_at = CURRENT_TIMESTAMP - INTERVAL ? MINUTE",
                     (minutes_ago,))
        conn.commit()
    finally:
        conn.close()


def _save(site, changes=None):
    return site.post("/api/admin/settings", json={"changes": changes or OWNER_CHANGE}).get_json()


def test_an_admin_is_refused_an_owner_option(site, config):
    _become(site, "admin")
    body = _save(site)
    assert not body["saved"] and "base_url" in body["errors"]
    assert "base_url" not in json.loads(config.read_text())


def test_an_owner_whose_password_is_stale_is_refused_and_the_page_says_to_type_it(site, config):
    _become(site, "owner", minutes_ago=11)
    body = _save(site)
    assert not body["saved"] and "base_url" in body["errors"]
    listing = site.get("/api/admin/settings").get_json()
    option = next(o for o in listing["options"] if o["path"] == "base_url")
    assert listing["reauth"] is True and option["fixed_kind"] == "reauth" and "Type your password" in option["reason"]
    page = site.get("/admin/settings").get_data(as_text=True)
    assert "Owner settings are locked" in page and 'name="current_password"' in page


def test_an_owner_who_typed_the_password_lately_may_change_it(site, config):
    _become(site, "owner", minutes_ago=9)
    body = _save(site)
    assert body["saved"] and json.loads(config.read_text())["base_url"] == "https://sky.example/"
    listing = site.get("/api/admin/settings").get_json()
    assert listing["reauth"] is False and next(o for o in listing["options"] if o["path"] == "base_url")["allowed"]
    assert "Owner settings are locked" not in site.get("/admin/settings").get_data(as_text=True)


def test_an_api_key_is_never_more_than_an_admin_even_the_owners(site, config):
    _become(site, "owner")
    conn = _control(site)
    try:
        admin_id = conn.execute("SELECT id FROM admin_users").fetchone()["id"]
        _key_id, raw = adminAuth.create_api_key(conn, admin_id, "ci")
    finally:
        conn.close()
    keyed = site.application.test_client()
    body = keyed.post("/api/admin/settings", json={"changes": OWNER_CHANGE}, headers={"Authorization": f"Bearer {raw}"})
    assert not body.get_json()["saved"] and "base_url" in body.get_json()["errors"]


def test_typing_the_password_again_reopens_a_stale_owner_and_a_wrong_one_does_not(site, config):
    _become(site, "owner", minutes_ago=30)
    wrong = site.post("/api/auth/confirm-password", json={"current_password": "nope"})
    assert wrong.status_code == 400
    assert not _save(site)["saved"]
    right = site.post("/api/auth/confirm-password", json={"current_password": PASSWORD})
    assert right.status_code == 200 and right.get_json() == {"minutes": adminAuth.RECENT_PASSWORD_MINUTES}
    assert _save(site)["saved"]


def test_the_page_form_types_the_password_and_opens_the_owner_settings(site, config):
    _become(site, "owner", minutes_ago=30)
    refused = _post(site, action="confirm-password", current_password="nope")
    assert refused.status_code == 400 and "incorrect" in refused.get_data(as_text=True)
    opened = _post(site, action="confirm-password", current_password=PASSWORD)
    assert opened.status_code == 303
    assert "open for 10 minutes" in site.get("/admin/settings").get_data(as_text=True)
    assert _save(site)["saved"]


def test_a_new_session_counts_as_the_password_typed(site):
    conn = _control(site)
    try:
        admin_id = conn.execute("SELECT id FROM admin_users").fetchone()["id"]
        token = adminAuth.create_session(conn, admin_id)
        assert adminAuth.validate_session(conn, token)["password_recent"] == 1
        conn.execute("UPDATE admin_sessions SET password_confirmed_at = NULL")
        conn.commit()
        assert adminAuth.validate_session(conn, token)["password_recent"] == 0
    finally:
        conn.close()


def test_set_owner_moves_the_one_owner_and_names_the_old_one(site):
    conn = _control(site)
    try:
        conn.execute("UPDATE admin_users SET role = 'owner'")
        conn.execute("INSERT INTO admin_users (username, password_hash, must_change_credentials, role) "
                     "VALUES ('second', ?, 0, 'admin')", (adminAuth.hash_password(PASSWORD),))
        conn.commit()
        first = conn.execute("SELECT username FROM admin_users WHERE username != 'second'").fetchone()["username"]
        assert adminAuth.set_owner(conn, "second") == first
        roles = {r["username"]: r["role"] for r in conn.execute("SELECT username, role FROM admin_users").fetchall()}
        assert roles == {first: "admin", "second": "owner"}
        assert adminAuth.set_owner(conn, "second") is None
        with pytest.raises(adminAuth.AuthError):
            adminAuth.set_owner(conn, "nobody")
    finally:
        conn.execute("DELETE FROM admin_users WHERE username = 'second'")
        conn.commit()
        conn.close()


def test_the_accounts_command_sets_the_owner(site, capsys, monkeypatch):
    monkeypatch.setattr(store, "configured_control_database", lambda: site.config.database)
    conn = _control(site)
    try:
        first = conn.execute("SELECT username FROM admin_users").fetchone()["username"]
        conn.execute("UPDATE admin_users SET role = 'admin'")
        conn.commit()
    finally:
        conn.close()
    config = site.config
    argv = ["--mysql-host", config.host, "--mysql-port", str(config.port), "--mysql-user", config.user,
            "--mysql-password", config.password, "--mysql-database", config.database]
    assert accounts.main(["--set-owner", first, *argv]) == 0
    assert f"{first} is the Owner." in capsys.readouterr().out
    with pytest.raises(SystemExit):
        accounts.main(["--set-owner", "nobody", *argv])
