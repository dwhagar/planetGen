# tests/test_admin_auth.py

"""
Unit tests for `stellarObjects/adminAuth.py` (password hashing, session/
API-key lifecycle, credential rotation, audit log) against a real,
throwaway MySQL database -- see `conftest.py`'s `mysql_config` fixture.
Skipped, not failed, when no MySQL test server is reachable.

Every test here treats `mysql_config`'s throwaway database as the
*control* schema directly (`_db.get_control_connection(mysql_config,
ensure_schema=True)`) rather than needing a second fixture -- the control
schema's tables (`admin_users`/`admin_sessions`/`admin_api_keys`/
`admin_audit_log`) never collide with a content schema's own tables, so
one throwaway database serves fine for testing either role.
"""

import pytest

from stellarObjects import _db, adminAuth


@pytest.fixture
def control_conn(mysql_config):
    conn = _db.get_control_connection(mysql_config, ensure_schema=True)
    try:
        yield conn
    finally:
        conn.close()


@pytest.fixture
def admin_id(control_conn):
    """Inserts one admin (not the seeded default) and returns its id."""
    cur = control_conn.execute(
        "INSERT INTO admin_users (username, password_hash, must_change_credentials) VALUES (?, ?, 0)",
        ("tester", adminAuth.hash_password("a-strong-test-password")),
    )
    control_conn.commit()
    return cur.lastrowid


def test_hash_and_verify_password_roundtrip():
    hashed = adminAuth.hash_password("correct horse battery staple")
    assert hashed != "correct horse battery staple"
    assert adminAuth.verify_password("correct horse battery staple", hashed)
    assert not adminAuth.verify_password("wrong password", hashed)


def test_validate_password_policy_rejects_short_default_and_username_match():
    with pytest.raises(adminAuth.AuthError):
        adminAuth.validate_password_policy("short")
    with pytest.raises(adminAuth.AuthError):
        adminAuth.validate_password_policy(adminAuth.DEFAULT_ADMIN_PASSWORD)
    with pytest.raises(adminAuth.AuthError):
        adminAuth.validate_password_policy("MyUsername123", username="myusername123")
    # Long enough, not the default, not the username -- no exception.
    adminAuth.validate_password_policy("a-perfectly-fine-password", username="someone-else")


def test_bootstrap_control_schema_seeds_default_admin_once(mysql_config):
    adminAuth.bootstrap_control_schema(mysql_config)

    conn = _db.get_control_connection(mysql_config, ensure_schema=False)
    try:
        row = conn.execute(
            "SELECT * FROM admin_users WHERE username = ?", (adminAuth.DEFAULT_ADMIN_USERNAME,),
        ).fetchone()
        assert row is not None
        assert row["must_change_credentials"] == 1
        assert adminAuth.verify_password(adminAuth.DEFAULT_ADMIN_PASSWORD, row["password_hash"])

        # Idempotent: a second bootstrap doesn't insert a second row or
        # reset an already-changed admin back to the default.
        conn.execute(
            "UPDATE admin_users SET must_change_credentials = 0 WHERE id = ?", (row["id"],),
        )
        conn.commit()
    finally:
        conn.close()

    adminAuth.bootstrap_control_schema(mysql_config)
    conn = _db.get_control_connection(mysql_config, ensure_schema=False)
    try:
        count = conn.execute(
            "SELECT COUNT(*) AS n FROM admin_users WHERE username = ?", (adminAuth.DEFAULT_ADMIN_USERNAME,),
        ).fetchone()["n"]
        assert count == 1
        still_changed = conn.execute(
            "SELECT must_change_credentials FROM admin_users WHERE username = ?",
            (adminAuth.DEFAULT_ADMIN_USERNAME,),
        ).fetchone()["must_change_credentials"]
        assert still_changed == 0
    finally:
        conn.close()


def test_authenticate_success_and_generic_failure_message(control_conn, admin_id):
    admin = adminAuth.authenticate(control_conn, "tester", "a-strong-test-password")
    assert admin["id"] == admin_id
    assert admin["last_login_at"] is not None

    with pytest.raises(adminAuth.AuthError) as wrong_password:
        adminAuth.authenticate(control_conn, "tester", "not the password")
    with pytest.raises(adminAuth.AuthError) as unknown_user:
        adminAuth.authenticate(control_conn, "no-such-admin", "irrelevant")
    # Same message either way -- never reveals whether a username exists.
    assert str(wrong_password.value) == str(unknown_user.value)


def test_session_lifecycle(control_conn, admin_id):
    assert adminAuth.validate_session(control_conn, "not-a-real-token") is None
    assert adminAuth.validate_session(control_conn, None) is None

    token = adminAuth.create_session(control_conn, admin_id)
    admin = adminAuth.validate_session(control_conn, token)
    assert admin is not None
    assert admin["id"] == admin_id

    adminAuth.end_session(control_conn, token)
    assert adminAuth.validate_session(control_conn, token) is None
    # Ending an already-gone session is a no-op, not an error.
    adminAuth.end_session(control_conn, token)


def test_session_expiry_is_enforced(control_conn, admin_id):
    token = adminAuth.create_session(control_conn, admin_id)
    # Force this session into the past, bypassing SESSION_TTL_HOURS, to
    # confirm expired sessions are actually rejected rather than merely
    # trusting the TTL value looks right.
    control_conn.execute(
        "UPDATE admin_sessions SET expires_at = DATE_SUB(CURRENT_TIMESTAMP, INTERVAL 1 SECOND) "
        "WHERE token_hash = ?",
        (adminAuth._hash_token(token),),
    )
    control_conn.commit()
    assert adminAuth.validate_session(control_conn, token) is None


def test_api_key_lifecycle(control_conn, admin_id):
    assert adminAuth.validate_api_key(control_conn, "not-a-real-key") is None

    key_id, raw_key = adminAuth.create_api_key(control_conn, admin_id, "test key")
    admin = adminAuth.validate_api_key(control_conn, raw_key)
    assert admin is not None
    assert admin["id"] == admin_id

    keys = adminAuth.list_api_keys(control_conn, admin_id)
    assert len(keys) == 1
    assert keys[0]["id"] == key_id
    assert keys[0]["revoked_at"] is None

    assert adminAuth.revoke_api_key(control_conn, admin_id, key_id) is True
    assert adminAuth.validate_api_key(control_conn, raw_key) is None
    # Revoking an already-revoked (or nonexistent) key is reported, not raised.
    assert adminAuth.revoke_api_key(control_conn, admin_id, key_id) is False


def test_revoke_api_key_is_scoped_to_its_own_admin(control_conn, admin_id):
    other_admin_cur = control_conn.execute(
        "INSERT INTO admin_users (username, password_hash) VALUES (?, ?)",
        ("other-admin", adminAuth.hash_password("another-strong-password")),
    )
    control_conn.commit()
    other_admin_id = other_admin_cur.lastrowid

    key_id, _raw_key = adminAuth.create_api_key(control_conn, admin_id, "belongs to tester")
    # other_admin can't revoke tester's key.
    assert adminAuth.revoke_api_key(control_conn, other_admin_id, key_id) is False
    assert adminAuth.list_api_keys(control_conn, admin_id)[0]["revoked_at"] is None


def test_change_credentials_requires_current_password_and_policy(control_conn, admin_id):
    with pytest.raises(adminAuth.AuthError):
        adminAuth.change_credentials(control_conn, admin_id, "wrong-current-password", "new-username", "a-strong-new-password")

    with pytest.raises(adminAuth.AuthError):
        adminAuth.change_credentials(control_conn, admin_id, "a-strong-test-password", "new-username", "short")

    adminAuth.change_credentials(control_conn, admin_id, "a-strong-test-password", "renamed-admin", "a-strong-new-password")
    row = control_conn.execute("SELECT * FROM admin_users WHERE id = ?", (admin_id,)).fetchone()
    assert row["username"] == "renamed-admin"
    assert row["must_change_credentials"] == 0
    assert adminAuth.verify_password("a-strong-new-password", row["password_hash"])


def test_change_credentials_rejects_username_already_in_use(control_conn, admin_id):
    control_conn.execute(
        "INSERT INTO admin_users (username, password_hash) VALUES (?, ?)",
        ("taken-username", adminAuth.hash_password("irrelevant-password")),
    )
    control_conn.commit()

    with pytest.raises(adminAuth.AuthError):
        adminAuth.change_credentials(
            control_conn, admin_id, "a-strong-test-password", "taken-username", "a-strong-new-password",
        )


def test_record_audit(control_conn, admin_id):
    adminAuth.record_audit(control_conn, admin_id, "tester", "sector.create", target="sector:1", detail="name='Test'")
    row = control_conn.execute(
        "SELECT * FROM admin_audit_log WHERE admin_user_id = ?", (admin_id,),
    ).fetchone()
    assert row["admin_username"] == "tester"
    assert row["action"] == "sector.create"
    assert row["target"] == "sector:1"
