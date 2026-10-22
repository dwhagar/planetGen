# tests/test_control_accounts_schema.py

"""
DB.58 (control schema v15): accounts gain a role, an email and a status, exactly one of them can be the owner, and the
`account_tokens` and `invite_uses` tables exist. A control database from before v15 keeps its accounts as admins, the lowest id
becomes the owner, and a second deploy changes nothing.
"""

import pytest

from planetgen.admin import auth
from planetgen.db import store

_V14_ADMIN_USERS = """CREATE TABLE admin_users (
    id                        BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    username                  VARCHAR(64) NOT NULL,
    password_hash             VARCHAR(255) NOT NULL,
    must_change_credentials   TINYINT(1) NOT NULL DEFAULT 1 CHECK (must_change_credentials IN (0, 1)),
    created_at                TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    last_login_at             TIMESTAMP NULL,
    UNIQUE (username)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


@pytest.fixture
def control(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    try:
        yield conn
    finally:
        conn.close()


def _account(conn, name, **more):
    columns = {"username": name, "password_hash": "x", **more}
    return conn.execute(f"INSERT INTO admin_users ({', '.join(columns)}) VALUES ({', '.join(['?'] * len(columns))})",
                        tuple(columns.values())).lastrowid


def _roles(conn):
    return {row["username"]: row["role"] for row in conn.execute("SELECT username, role FROM admin_users").fetchall()}


def test_the_install_makes_exactly_one_owner(mysql_config):
    username, _password = auth.bootstrap_control_schema(mysql_config)
    conn = store.get_control_connection(mysql_config, ensure_schema=False)
    try:
        assert _roles(conn) == {username: "owner"}
        row = conn.execute("SELECT email, status, owner_flag FROM admin_users").fetchone()
        assert row == {"email": None, "status": "active", "owner_flag": 1}
    finally:
        conn.close()


def test_a_second_owner_is_refused_but_other_roles_are_not(control):
    _account(control, "boss", role="owner")
    with pytest.raises(Exception, match="uq_admin_users_owner"):
        _account(control, "usurper", role="owner")
    for role in ("user", "moderator", "admin"):
        _account(control, f"a-{role}", role=role)
    control.commit()
    assert _roles(control)["a-moderator"] == "moderator"
    # Demoting the owner frees the place.
    control.execute("UPDATE admin_users SET role = 'admin' WHERE username = 'boss'")
    _account(control, "next", role="owner")


def test_a_role_and_a_status_have_known_values_and_an_email_is_unique(control):
    with pytest.raises(Exception, match="chk_admin_users_role"):
        _account(control, "a", role="superuser")
    with pytest.raises(Exception, match="chk_admin_users_status"):
        _account(control, "b", status="banned")
    _account(control, "c", email="c@example.org")
    _account(control, "d", status="disabled")
    _account(control, "e")
    _account(control, "f")
    with pytest.raises(Exception, match="uq_admin_users_email"):
        _account(control, "g", email="c@example.org")
    row = control.execute("SELECT role, status FROM admin_users WHERE username = 'e'").fetchone()
    assert row == {"role": "admin", "status": "active"}


def test_tokens_have_a_purpose_a_unique_hash_and_go_with_their_account(control):
    boss = _account(control, "boss", role="owner")
    other = _account(control, "other")
    token = control.execute("INSERT INTO account_tokens (account_id, purpose, token_hash, meta) VALUES (?, 'invite', ?, ?)",
                            (boss, "a" * 64, '{"uses": 3}')).lastrowid
    with pytest.raises(Exception, match="chk_account_tokens_purpose"):
        control.execute("INSERT INTO account_tokens (account_id, purpose, token_hash) VALUES (?, 'sign_in', ?)", (boss, "b" * 64))
    with pytest.raises(Exception, match="uq_account_tokens_hash"):
        control.execute("INSERT INTO account_tokens (account_id, purpose, token_hash) VALUES (?, 'reset', ?)", (boss, "a" * 64))
    row = control.execute("SELECT used_at, expires_at FROM account_tokens WHERE id = ?", (token,)).fetchone()
    assert row == {"used_at": None, "expires_at": None}
    control.execute("INSERT INTO invite_uses (token_id, account_id) VALUES (?, ?)", (token, other))
    with pytest.raises(Exception, match="PRIMARY"):
        control.execute("INSERT INTO invite_uses (token_id, account_id) VALUES (?, ?)", (token, other))
    control.execute("DELETE FROM admin_users WHERE id = ?", (other,))
    assert control.execute("SELECT COUNT(*) AS n FROM invite_uses").fetchone()["n"] == 0
    control.execute("INSERT INTO invite_uses (token_id, account_id) VALUES (?, ?)", (token, boss))
    control.execute("DELETE FROM admin_users WHERE id = ?", (boss,))
    assert control.execute("SELECT COUNT(*) AS n FROM account_tokens").fetchone()["n"] == 0
    assert control.execute("SELECT COUNT(*) AS n FROM invite_uses").fetchone()["n"] == 0


def test_an_older_control_database_keeps_its_admins_and_the_lowest_becomes_the_owner(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute(_V14_ADMIN_USERS)
        conn.execute("CREATE TABLE control_schema_migrations (version INT NOT NULL PRIMARY KEY,"
                     " applied_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP)")
        conn.execute("INSERT INTO control_schema_migrations (version) VALUES (14)")
        for name in ("first", "second", "third"):
            conn.execute("INSERT INTO admin_users (username, password_hash) VALUES (?, 'x')", (name,))
        conn.commit()
        store._ensure_control_schema(conn)
        assert _roles(conn) == {"first": "owner", "second": "admin", "third": "admin"}
        assert conn.execute("SELECT MAX(version) AS v FROM control_schema_migrations").fetchone()["v"] == store.CONTROL_SCHEMA_VERSION
        # The keys and CHECKs are the new install's: a second owner is refused here too.
        with pytest.raises(Exception, match="uq_admin_users_owner"):
            conn.execute("UPDATE admin_users SET role = 'owner' WHERE username = 'second'")
        conn.rollback()
        # A second deploy, even one that finds the old version again, changes no role.
        conn.execute("UPDATE admin_users SET role = 'admin' WHERE username = 'first'")
        conn.execute("UPDATE admin_users SET role = 'owner' WHERE username = 'third'")
        conn.execute("DELETE FROM control_schema_migrations WHERE version > 14")
        conn.commit()
        store._ensure_control_schema(conn)
        assert _roles(conn) == {"first": "admin", "second": "admin", "third": "owner"}
    finally:
        conn.close()
