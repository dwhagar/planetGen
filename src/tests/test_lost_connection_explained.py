"""OPS.171: a statement that loses its MySQL connection leaves a warning that says what the server says about the loss."""

import pymysql
import pytest

from planetgen.db import connection_loss, contention


class _Rows:
    def __init__(self, rows):
        self.rows = rows

    def fetchall(self):
        return self.rows

    def fetchone(self):
        return self.rows[0]


class _Server:
    """A connection that answers `facts` the way a MySQL that restarted 40 s ago would."""

    def execute(self, sql, params=()):
        if "SHOW GLOBAL STATUS" in sql:
            return _Rows([{"Variable_name": "Uptime", "Value": "40"}, {"Variable_name": "Aborted_clients", "Value": "3"},
                          {"Variable_name": "Aborted_connects", "Value": "0"}])
        return _Rows([{"wait_timeout": 28800, "net_read_timeout": 30, "net_write_timeout": 60, "max_allowed_packet": 67108864}])

    def close(self):
        pass


def test_a_loss_is_described_from_the_servers_own_figures(monkeypatch):
    from planetgen.db import store
    warned = []
    monkeypatch.setattr(store, "get_connection", lambda config: _Server())
    monkeypatch.setattr(connection_loss.log, "warning", warned.append)
    connection_loss.explain_lost(object(), pymysql.err.OperationalError(2013, "Lost connection"))
    (line,) = warned
    assert "uptime is 40 s" in line and "wait_timeout 28800 s" in line and "2013" in line


def test_a_server_that_cannot_be_asked_still_leaves_a_warning_and_never_raises(monkeypatch):
    from planetgen.db import store

    def refuse(config):
        raise pymysql.err.OperationalError(2003, "can't connect")

    warned = []
    monkeypatch.setattr(store, "get_connection", refuse)
    monkeypatch.setattr(connection_loss.log, "warning", warned.append)
    connection_loss.explain_lost(object(), pymysql.err.OperationalError(2006, "gone away"))
    assert len(warned) == 1 and "could not be read" in warned[0]


def test_run_with_retries_explains_a_lost_connection_before_it_raises(monkeypatch):
    seen = []
    monkeypatch.setattr(connection_loss, "explain_lost", lambda config, exc: seen.append((config, exc.args[0])))

    class Conn:
        _config = "the-config"

        def rollback(self):
            pass

    class Gate:
        def is_closed(self):
            return False

    def lost(conn):
        raise pymysql.err.OperationalError(2013, "Lost connection to MySQL server during query")

    with pytest.raises(pymysql.err.OperationalError):
        contention.run_with_retries(Conn(), lost, "gate", "Unit 1", attempts=2, gate_object=Gate())
    assert seen == [("the-config", 2013)]
