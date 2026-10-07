# tests/test_admin_auth.py

"""
Unit tests for `planetgen/admin/auth.py` (password hashing, session/
API-key lifecycle, credential rotation, audit log) against a real,
throwaway MySQL database -- see `conftest.py`'s `mysql_config` fixture.
Skipped, not failed, when no MySQL test server is reachable.

Every test here treats `mysql_config`'s throwaway database as the
*control* schema directly (`store.get_control_connection(mysql_config,
ensure_schema=True)`) rather than needing a second fixture -- the control
schema's tables (`admin_users`/`admin_sessions`/`admin_api_keys`/
`admin_audit_log`) never collide with a content schema's own tables, so
one throwaway database serves fine for testing either role.
"""

import pytest

from planetgen.db import store
from planetgen.admin import auth


@pytest.fixture
def control_conn(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    try:
        yield conn
    finally:
        conn.close()


@pytest.fixture
def admin_id(control_conn):
    """Inserts one admin (not the seeded default) and returns its id."""
    cur = control_conn.execute(
        "INSERT INTO admin_users (username, password_hash, must_change_credentials) VALUES (?, ?, 0)",
        ("tester", auth.hash_password("a-strong-test-password")),
    )
    control_conn.commit()
    return cur.lastrowid


def test_hash_and_verify_password_roundtrip():
    hashed = auth.hash_password("correct horse battery staple")
    assert hashed != "correct horse battery staple"
    assert auth.verify_password("correct horse battery staple", hashed)
    assert not auth.verify_password("wrong password", hashed)


def test_validate_password_policy_rejects_short_default_and_username_match():
    with pytest.raises(auth.AuthError):
        auth.validate_password_policy("short")
    with pytest.raises(auth.AuthError):
        auth.validate_password_policy("password")
    with pytest.raises(auth.AuthError):
        auth.validate_password_policy("MyUsername123", username="myusername123")
    # Long enough, not the default, not the username -- no exception.
    auth.validate_password_policy("a-perfectly-fine-password", username="someone-else")


def test_bootstrap_control_schema_seeds_default_admin_once(mysql_config):
    username, password = auth.bootstrap_control_schema(mysql_config)
    assert username == auth.DEFAULT_ADMIN_USERNAME

    conn = store.get_control_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute(
            "SELECT * FROM admin_users WHERE username = ?", (auth.DEFAULT_ADMIN_USERNAME,),
        ).fetchone()
        assert row is not None
        assert row["must_change_credentials"] == 1
        assert auth.verify_password(password, row["password_hash"])
        assert not auth.verify_password("password", row["password_hash"])

        # Idempotent: a second bootstrap doesn't insert a second row or
        # reset an already-changed admin back to the default.
        conn.execute(
            "UPDATE admin_users SET must_change_credentials = 0 WHERE id = ?", (row["id"],),
        )
        conn.commit()
    finally:
        conn.close()

    # Nothing seeded, so nothing to show (migrateDb prints only a new seed).
    assert auth.bootstrap_control_schema(mysql_config) is None
    conn = store.get_control_connection(mysql_config, ensure_schema=False)
    try:
        count = conn.execute(
            "SELECT COUNT(*) AS n FROM admin_users WHERE username = ?", (auth.DEFAULT_ADMIN_USERNAME,),
        ).fetchone()["n"]
        assert count == 1
        still_changed = conn.execute(
            "SELECT must_change_credentials FROM admin_users WHERE username = ?",
            (auth.DEFAULT_ADMIN_USERNAME,),
        ).fetchone()["must_change_credentials"]
        assert still_changed == 0
    finally:
        conn.close()


def test_bootstrap_seeds_a_random_first_password_meeting_the_policy(mysql_config):
    """Security #39: no published default. Each fresh control schema gets
    its own random password, long enough for the policy, stored only as a
    hash; the seeded admin must still change it before doing anything."""
    _username, first = auth.bootstrap_control_schema(mysql_config)
    auth.validate_password_policy(first, username=auth.DEFAULT_ADMIN_USERNAME)
    assert len(first) >= auth.MIN_PASSWORD_LENGTH
    conn = store.get_control_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute("SELECT * FROM admin_users").fetchone()
        assert first not in row["password_hash"]
        assert row["must_change_credentials"] == 1
        admin = auth.authenticate(conn, auth.DEFAULT_ADMIN_USERNAME, first)
        assert admin["must_change_credentials"] == 1
        with pytest.raises(auth.AuthError):
            auth.authenticate(conn, auth.DEFAULT_ADMIN_USERNAME, "password")
        # A second seed (as after deleting the admin rows) is a new password.
        conn.execute("DELETE FROM admin_users")
        conn.commit()
    finally:
        conn.close()
    _username, second = auth.bootstrap_control_schema(mysql_config)
    assert second != first


def test_bootstrap_rotates_an_older_install_still_on_the_published_password(mysql_config):
    """An install seeded before random first passwords whose admin never
    made the forced change is still on admin/password; the next
    migrateDb run gives it a random password and ends its sessions. A
    changed login, or one that no longer verifies, is left alone."""
    auth.bootstrap_control_schema(mysql_config)
    conn = store.get_control_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("UPDATE admin_users SET password_hash = ?", (auth.hash_password("password"),))
        conn.commit()
        admin = auth.authenticate(conn, auth.DEFAULT_ADMIN_USERNAME, "password")
        auth.create_session(conn, admin["id"])
    finally:
        conn.close()
    username, rotated = auth.bootstrap_control_schema(mysql_config)
    assert username == auth.DEFAULT_ADMIN_USERNAME and rotated != "password"
    conn = store.get_control_connection(mysql_config, ensure_schema=False)
    try:
        with pytest.raises(auth.AuthError):
            auth.authenticate(conn, username, "password")
        assert auth.authenticate(conn, username, rotated)["must_change_credentials"] == 1
        assert conn.execute("SELECT COUNT(*) AS n FROM admin_sessions").fetchone()["n"] == 0
    finally:
        conn.close()
    assert auth.bootstrap_control_schema(mysql_config) is None


def test_migrate_db_prints_the_first_password_only_when_seeded(monkeypatch, capsys):
    """Security #39: `planetgen.cli.migrate` shows the seeded login exactly once."""
    from planetgen.cli import migrate

    migrate.print_initial_admin_login("admin", "Zx9-random-first-pass")
    out = capsys.readouterr().out
    assert "username: admin" in out
    assert "password: Zx9-random-first-pass" in out
    assert "/login" in out and "change" in out

    monkeypatch.setattr(migrate, "_migrate_with_progress", lambda config: migrate.SCHEMA_VERSION)
    monkeypatch.setattr(migrate.sys, "argv", ["planetgen.cli.migrate"])
    for seeded, shown in ((("admin", "Zx9-random-first-pass"), True), (None, False)):
        monkeypatch.setattr(migrate.auth, "bootstrap_control_schema", lambda config, s=seeded: s)
        migrate.main()
        out = capsys.readouterr().out
        assert ("Zx9-random-first-pass" in out) is shown
        assert ("New admin login password" in out) is shown


def test_authenticate_unknown_username_still_checks_a_password_hash(control_conn, admin_id, monkeypatch):
    """Security #44: an unknown username runs `verify_password` (against
    a dummy hash made with `hash_password`) so it takes as long as a wrong
    password for a real one."""
    calls = []
    real_verify = auth.verify_password
    monkeypatch.setattr(auth, "verify_password",
                        lambda password, password_hash: calls.append(password_hash) or real_verify(password, password_hash))
    with pytest.raises(auth.AuthError):
        auth.authenticate(control_conn, "no-such-admin", "irrelevant")
    assert len(calls) == 1
    dummy = calls[0]
    assert dummy == auth._get_dummy_password_hash()
    # Same algorithm and cost as a real password hash.
    real = auth.hash_password("x" * 20)
    assert dummy.split("$", 1)[0] == real.split("$", 1)[0]
    # The dummy never matches a guess, including an empty one.
    assert not real_verify("", dummy) and not real_verify("irrelevant", dummy)

    calls.clear()
    with pytest.raises(auth.AuthError):
        auth.authenticate(control_conn, "tester", "not the password")
    assert len(calls) == 1 and calls[0] != dummy


def test_change_credentials_ends_every_session_but_keeps_api_keys(control_conn, admin_id):
    """Security #45: a credential change logs out every session of that
    admin (the API route then issues the caller a fresh one); API keys
    and other admins' sessions are untouched."""
    other = control_conn.execute(
        "INSERT INTO admin_users (username, password_hash, must_change_credentials) VALUES (?, ?, 0)",
        ("other-admin", auth.hash_password("another-strong-password")),
    ).lastrowid
    control_conn.commit()
    stolen = auth.create_session(control_conn, admin_id)
    current = auth.create_session(control_conn, admin_id)
    others_session = auth.create_session(control_conn, other)
    _key_id, raw_key = auth.create_api_key(control_conn, admin_id, "ci")

    # A failed change ends nothing.
    with pytest.raises(auth.AuthError):
        auth.change_credentials(control_conn, admin_id, "wrong", "tester", "a-strong-new-password")
    assert auth.validate_session(control_conn, stolen) is not None

    auth.change_credentials(control_conn, admin_id, "a-strong-test-password", "tester", "a-strong-new-password")
    assert auth.validate_session(control_conn, stolen) is None
    assert auth.validate_session(control_conn, current) is None
    assert auth.validate_session(control_conn, others_session)["id"] == other
    assert auth.validate_api_key(control_conn, raw_key)["id"] == admin_id
    fresh = auth.create_session(control_conn, admin_id)
    assert auth.validate_session(control_conn, fresh)["id"] == admin_id


def test_authenticate_success_and_generic_failure_message(control_conn, admin_id):
    admin = auth.authenticate(control_conn, "tester", "a-strong-test-password")
    assert admin["id"] == admin_id
    assert admin["last_login_at"] is not None

    with pytest.raises(auth.AuthError) as wrong_password:
        auth.authenticate(control_conn, "tester", "not the password")
    with pytest.raises(auth.AuthError) as unknown_user:
        auth.authenticate(control_conn, "no-such-admin", "irrelevant")
    # Same message either way -- never reveals whether a username exists.
    assert str(wrong_password.value) == str(unknown_user.value)


def test_session_lifecycle(control_conn, admin_id):
    assert auth.validate_session(control_conn, "not-a-real-token") is None
    assert auth.validate_session(control_conn, None) is None

    token = auth.create_session(control_conn, admin_id)
    admin = auth.validate_session(control_conn, token)
    assert admin is not None
    assert admin["id"] == admin_id

    auth.end_session(control_conn, token)
    assert auth.validate_session(control_conn, token) is None
    # Ending an already-gone session is a no-op, not an error.
    auth.end_session(control_conn, token)


def test_session_expiry_is_enforced(control_conn, admin_id):
    token = auth.create_session(control_conn, admin_id)
    # Force this session into the past, bypassing SESSION_TTL_HOURS, to
    # confirm expired sessions are actually rejected rather than merely
    # trusting the TTL value looks right.
    control_conn.execute(
        "UPDATE admin_sessions SET expires_at = DATE_SUB(CURRENT_TIMESTAMP, INTERVAL 1 SECOND) "
        "WHERE token_hash = ?",
        (auth._hash_token(token),),
    )
    control_conn.commit()
    assert auth.validate_session(control_conn, token) is None


def test_api_key_lifecycle(control_conn, admin_id):
    assert auth.validate_api_key(control_conn, "not-a-real-key") is None

    key_id, raw_key = auth.create_api_key(control_conn, admin_id, "test key")
    admin = auth.validate_api_key(control_conn, raw_key)
    assert admin is not None
    assert admin["id"] == admin_id

    keys = auth.list_api_keys(control_conn, admin_id)
    assert len(keys) == 1
    assert keys[0]["id"] == key_id
    assert keys[0]["revoked_at"] is None

    assert auth.revoke_api_key(control_conn, admin_id, key_id) is True
    assert auth.validate_api_key(control_conn, raw_key) is None
    # Revoking an already-revoked (or nonexistent) key is reported, not raised.
    assert auth.revoke_api_key(control_conn, admin_id, key_id) is False


def test_revoke_api_key_is_scoped_to_its_own_admin(control_conn, admin_id):
    other_admin_cur = control_conn.execute(
        "INSERT INTO admin_users (username, password_hash) VALUES (?, ?)",
        ("other-admin", auth.hash_password("another-strong-password")),
    )
    control_conn.commit()
    other_admin_id = other_admin_cur.lastrowid

    key_id, _raw_key = auth.create_api_key(control_conn, admin_id, "belongs to tester")
    # other_admin can't revoke tester's key.
    assert auth.revoke_api_key(control_conn, other_admin_id, key_id) is False
    assert auth.list_api_keys(control_conn, admin_id)[0]["revoked_at"] is None


def test_change_credentials_requires_current_password_and_policy(control_conn, admin_id):
    with pytest.raises(auth.AuthError):
        auth.change_credentials(control_conn, admin_id, "wrong-current-password", "new-username", "a-strong-new-password")

    with pytest.raises(auth.AuthError):
        auth.change_credentials(control_conn, admin_id, "a-strong-test-password", "new-username", "short")

    auth.change_credentials(control_conn, admin_id, "a-strong-test-password", "renamed-admin", "a-strong-new-password")
    row = control_conn.execute("SELECT * FROM admin_users WHERE id = ?", (admin_id,)).fetchone()
    assert row["username"] == "renamed-admin"
    assert row["must_change_credentials"] == 0
    assert auth.verify_password("a-strong-new-password", row["password_hash"])


def test_change_credentials_rejects_username_already_in_use(control_conn, admin_id):
    control_conn.execute(
        "INSERT INTO admin_users (username, password_hash) VALUES (?, ?)",
        ("taken-username", auth.hash_password("irrelevant-password")),
    )
    control_conn.commit()

    with pytest.raises(auth.AuthError):
        auth.change_credentials(
            control_conn, admin_id, "a-strong-test-password", "taken-username", "a-strong-new-password",
        )


def test_record_audit(control_conn, admin_id):
    auth.record_audit(control_conn, admin_id, "tester", "sector.create", target="sector:1", detail="name='Test'")
    row = control_conn.execute(
        "SELECT * FROM admin_audit_log WHERE admin_user_id = ?", (admin_id,),
    ).fetchone()
    assert row["admin_username"] == "tester"
    assert row["action"] == "sector.create"
    assert row["target"] == "sector:1"
