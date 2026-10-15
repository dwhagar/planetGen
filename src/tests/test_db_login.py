# tests/test_db_login.py

"""
OPS.49: `db.login.connect` explains a failed login. pymysql 1.2.1 to 1.2.3 raised "'NoneType' object has no attribute
'is_auth_switch_request'" for a MySQL 8 (caching_sha2_password) login without TLS, and update.sh printed just that.
"""

import pymysql
import pytest

from planetgen.db import login

ARGS = dict(host="db.example", port=3306, user="planetgen", password="secret")


def _raising(exc):
    def connect(**kwargs):
        raise exc
    return connect


def test_a_pymysql_that_breaks_inside_the_handshake_becomes_a_login_error(monkeypatch):
    monkeypatch.setattr(pymysql, "connect",
                        _raising(AttributeError("'NoneType' object has no attribute 'is_auth_switch_request'")))
    with pytest.raises(login.DatabaseLoginError) as caught:
        login.connect(**ARGS)
    message = str(caught.value)
    assert "Couldn't log in to the database server at db.example:3306 as 'planetgen'" in message
    assert "AttributeError" in message and "is_auth_switch_request" in message
    assert "secret" not in message and "config.json" in message


def test_access_denied_becomes_a_login_error(monkeypatch):
    monkeypatch.setattr(pymysql, "connect", _raising(pymysql.err.OperationalError(1045, "Access denied for user")))
    with pytest.raises(login.DatabaseLoginError, match="Access denied for user"):
        login.connect(**ARGS)


def test_errors_that_already_say_what_is_wrong_pass_through(monkeypatch):
    for exc in (pymysql.err.OperationalError(2003, "Can't connect to MySQL server"),
                pymysql.err.OperationalError(1049, "Unknown database 'x'"), ConnectionRefusedError("refused")):
        monkeypatch.setattr(pymysql, "connect", _raising(exc))
        with pytest.raises(type(exc)) as caught:
            login.connect(**ARGS)
        assert caught.value is exc


def test_a_good_login_returns_the_connection(monkeypatch):
    monkeypatch.setattr(pymysql, "connect", lambda **kwargs: ("connection", kwargs))
    connection, kwargs = login.connect(**ARGS)
    assert connection == "connection"
    assert {k: v for k, v in kwargs.items() if k != "auth_plugin_map"} == ARGS


def test_caching_sha2_logins_use_pymysqls_step_and_accept_a_release_that_returns_nothing(monkeypatch):
    """pymysql 1.2.3's step ends without `return`, so the login loop gets `None`; a release with the fix returns the
    server's packet. The handler passes either through."""
    from pymysql import _auth
    seen = {}
    monkeypatch.setattr(pymysql, "connect", lambda **kwargs: seen.update(kwargs))
    login.connect(**ARGS)
    handler = seen["auth_plugin_map"]["caching_sha2_password"]("the connection")
    for result in (None, "the OK packet"):
        monkeypatch.setattr(_auth, "caching_sha2_password_auth", lambda conn, packet, result=result: (conn, packet, result))
        assert handler.authenticate("the packet") == ("the connection", "the packet", result)


def test_every_database_connection_goes_through_the_wrapper():
    import os
    import re
    root = os.path.join(os.path.dirname(__file__), "..", "planetgen")
    found = []
    for folder, _dirs, files in os.walk(root):
        for name in files:
            path = os.path.join(folder, name)
            if name.endswith(".py") and not path.endswith(os.path.join("db", "login.py")):
                with open(path, encoding="utf-8") as f:
                    found += [path for line in f if re.search(r"\bpymysql\.connect\(", line)]
    assert found == []


def test_a_login_error_is_still_an_operational_error_that_prints_its_message(monkeypatch):
    """Code that handles an unreachable database (`except pymysql.err.OperationalError`) must keep handling it."""
    monkeypatch.setattr(pymysql, "connect", _raising(pymysql.err.OperationalError(1045, "Access denied")))
    with pytest.raises(pymysql.err.OperationalError) as caught:
        login.connect(**ARGS)
    assert caught.value.args[0] == 1045 and str(caught.value).startswith("Couldn't log in to the database server")
