# tests/test_migration_helpers.py

"""
DB.18: the helpers a revision uses on a big table (`planetgen/db/migrations/helpers.py`) write the `ALGORITHM=` out,
skip what is already done, update in committed pages and report progress.
"""

import pytest
import sqlalchemy as sa

from planetgen.db import alembic_runner
from planetgen.db.migrations import helpers

pytestmark = pytest.mark.fresh_database


@pytest.fixture
def connection(mysql_config):
    engine = alembic_runner._engine(mysql_config)
    with engine.connect() as conn:
        conn.execute(sa.text("CREATE TABLE helper_probe (id INT NOT NULL PRIMARY KEY, a INT NOT NULL, b INT NULL)"
                             " ENGINE=InnoDB"))
        conn.execute(sa.text("INSERT INTO helper_probe (id, a) VALUES " + ", ".join(f"({i}, {i * 2})" for i in range(1, 26))))
        conn.commit()
        try:
            yield conn
        finally:
            conn.rollback()  # an open transaction holds the tables' metadata locks and would block the drop
    engine.dispose()


@pytest.fixture
def progress(monkeypatch):
    calls = []
    monkeypatch.setattr(alembic_runner, "report_progress", lambda label, done, total: calls.append((label, done, total)))
    return calls


def test_exists_checks(connection):
    assert helpers.exists_table(connection, "helper_probe") and not helpers.exists_table(connection, "nope")
    assert helpers.exists_column(connection, "helper_probe", "a") and not helpers.exists_column(connection, "helper_probe", "z")
    assert helpers.exists_index(connection, "helper_probe", "PRIMARY")
    assert helpers.index_columns(connection, "helper_probe", "PRIMARY") == ["id"]
    assert helpers.column_type(connection, "helper_probe", "a") == "int"
    assert helpers.column_type(connection, "helper_probe", "z") is None


def test_an_instant_column_is_added_once_and_dropped_once(connection):
    assert helpers.add_column_instant(connection, "helper_probe", "c", "INT NOT NULL DEFAULT 7", after="a") is True
    assert helpers.add_column_instant(connection, "helper_probe", "c", "INT NOT NULL DEFAULT 7") is False
    assert connection.execute(sa.text("SELECT c FROM helper_probe WHERE id = 1")).scalar() == 7
    assert helpers.drop_column_instant(connection, "helper_probe", "c") is True
    assert helpers.drop_column_instant(connection, "helper_probe", "c") is False


def test_a_change_that_cannot_be_instant_is_refused_not_copied(connection):
    with pytest.raises(sa.exc.OperationalError, match="ALGORITHM=INSTANT is not supported|1845"):
        helpers.add_column_instant(connection, "helper_probe", "d", "INT NULL AUTO_INCREMENT UNIQUE")


def test_an_index_builds_online_once(connection):
    assert helpers.add_index_online(connection, "helper_probe", "idx_probe_a", ("a DESC", "b")) is True
    assert helpers.add_index_online(connection, "helper_probe", "idx_probe_a", ("a", "b")) is False
    assert helpers.index_columns(connection, "helper_probe", "idx_probe_a") == ["a", "b"]
    with pytest.raises(RuntimeError, match="exists with columns"):
        helpers.add_index_online(connection, "helper_probe", "idx_probe_a", ("b",))
    assert helpers.drop_index_online(connection, "helper_probe", "idx_probe_a") is True
    assert helpers.drop_index_online(connection, "helper_probe", "idx_probe_a") is False


def test_a_unique_index_can_be_asked_for(connection):
    helpers.add_index_online(connection, "helper_probe", "uq_probe_a", ("a",), unique=True)
    with pytest.raises(sa.exc.IntegrityError):
        connection.execute(sa.text("INSERT INTO helper_probe (id, a) VALUES (100, 2)"))
    connection.rollback()


@pytest.mark.parametrize("bad", ["a; DROP TABLE x", "a b c", "a SIDEWAYS", "1a", "a-b"])
def test_names_and_directions_are_not_spliced_into_sql(connection, bad):
    with pytest.raises(ValueError):
        helpers.add_index_online(connection, "helper_probe", "idx_x", (bad,))
    with pytest.raises(ValueError):
        helpers.add_column_instant(connection, "helper_probe", bad, "INT NULL")


def test_a_batched_update_commits_each_page_and_reports_progress(connection, progress):
    done = helpers.batched_update(connection, "helper_probe", "b = a * 2", "b IS NULL", batch=10, label="probe")
    assert done == 25
    assert progress == [("probe", 10, 25), ("probe", 20, 25), ("probe", 25, 25)]
    assert connection.execute(sa.text("SELECT COUNT(*) FROM helper_probe WHERE b IS NULL")).scalar() == 0
    assert connection.execute(sa.text("SELECT b FROM helper_probe WHERE id = 3")).scalar() == 12
    connection.rollback()  # nothing is left uncommitted: the rows are still changed
    assert connection.execute(sa.text("SELECT COUNT(*) FROM helper_probe WHERE b = a * 2")).scalar() == 25


def test_a_stopped_batched_update_starts_again_where_it_stopped(connection, progress):
    connection.execute(sa.text("UPDATE helper_probe SET b = 0 WHERE id <= 12"))
    connection.commit()
    assert helpers.batched_update(connection, "helper_probe", "b = 1", "b IS NULL", batch=5, label="probe") == 13
    assert progress[-1] == ("probe", 13, 13)
    assert connection.execute(sa.text("SELECT COUNT(*) FROM helper_probe WHERE b = 0")).scalar() == 12
    assert helpers.batched_update(connection, "helper_probe", "b = 1", "b IS NULL", label="probe") == 0
    assert progress[-1] == ("probe", 0, 0)


def test_a_batched_update_takes_values_as_parameters(connection):
    helpers.batched_update(connection, "helper_probe", "b = :v", "b IS NULL AND a > :floor", batch=7,
                           params={"v": 9, "floor": 40})
    assert connection.execute(sa.text("SELECT COUNT(*) FROM helper_probe WHERE b = 9")).scalar() == 5


def test_the_lock_wait_timeout_is_set_for_the_block_only(connection):
    before = connection.execute(sa.text("SELECT @@SESSION.lock_wait_timeout")).scalar()
    with helpers.lock_wait_timeout(connection, 7):
        assert connection.execute(sa.text("SELECT @@SESSION.lock_wait_timeout")).scalar() == 7
    assert connection.execute(sa.text("SELECT @@SESSION.lock_wait_timeout")).scalar() == before


def test_a_long_transaction_stops_the_update(connection, mysql_config):
    other = alembic_runner._engine(mysql_config).connect()
    try:
        other.execute(sa.text("SELECT COUNT(*) FROM helper_probe"))  # opens a transaction
        helpers.require_no_long_transactions(connection, seconds=3600)
        with pytest.raises(RuntimeError, match="open too long"):
            helpers.require_no_long_transactions(connection, seconds=0)
    finally:
        other.close()


def test_a_lock_wait_timeout_is_tried_again_then_raised(connection, mysql_config, monkeypatch):
    monkeypatch.setattr(helpers.time, "sleep", lambda seconds: None)
    holder = alembic_runner._engine(mysql_config).connect()
    try:
        holder.execute(sa.text("UPDATE helper_probe SET b = 1 WHERE id = 1"))  # a row lock to wait on
        connection.execute(sa.text("SET SESSION innodb_lock_wait_timeout = 1"))  # row locks wait on this one
        with pytest.raises(sa.exc.OperationalError, match="Lock wait timeout"):
            helpers._with_retries(connection, "UPDATE helper_probe SET b = 5 WHERE id = 1")
    finally:
        holder.close()
