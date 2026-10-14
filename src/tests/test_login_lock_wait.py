# tests/test_login_lock_wait.py

"""
ADM.53: a login (or any request that only wants to note that it happened) must not fail or hang because a
long job holds a lock on the user's row. Against a real throwaway MySQL database (see `conftest.py`).
"""

import time

import pytest

from planetgen.admin import auth
from planetgen.db import store


@pytest.fixture
def control_conn(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True, lock_wait_s=1)
    try:
        yield conn
    finally:
        conn.close()


@pytest.fixture
def holder(mysql_config, control_conn):
    """A second connection with an open transaction holding the lock on the `tester` row."""
    control_conn.execute(
        "INSERT INTO admin_users (username, password_hash, must_change_credentials) VALUES (?, ?, 0)",
        ("tester", auth.hash_password("a-strong-test-password")))
    control_conn.commit()
    other = store.get_control_connection(mysql_config)
    other.execute("SELECT id FROM admin_users WHERE username = 'tester' FOR UPDATE")
    try:
        yield other
    finally:
        other.rollback()
        other.close()


def test_a_connection_can_be_given_a_short_lock_wait(control_conn):
    row = control_conn.execute("SELECT @@innodb_lock_wait_timeout AS t").fetchone()
    assert int(row["t"]) == 1


def test_login_succeeds_while_another_transaction_holds_the_row(control_conn, holder):
    started = time.monotonic()
    user = auth.authenticate(control_conn, "tester", "a-strong-test-password")
    assert user["username"] == "tester"
    assert time.monotonic() - started < 10  # waited about the 1 s lock wait, not the 50 s default


def test_a_session_can_be_validated_while_the_row_is_locked(control_conn, holder):
    user = control_conn.execute("SELECT id FROM admin_users WHERE username = 'tester'").fetchone()
    holder.rollback()
    token = auth.create_session(control_conn, user["id"])
    holder.execute("SELECT id FROM admin_users WHERE username = 'tester' FOR UPDATE")
    control_conn.execute("UPDATE admin_sessions SET last_seen_at = DATE_SUB(NOW(), INTERVAL 1 HOUR)")
    control_conn.commit()
    assert auth.validate_session(control_conn, token)["username"] == "tester"


def test_other_database_errors_still_raise(control_conn):
    import pymysql
    with pytest.raises(pymysql.err.Error):
        auth._write_if_not_locked(control_conn, "UPDATE no_such_table SET x = ? WHERE id = ?", (1, 1))
