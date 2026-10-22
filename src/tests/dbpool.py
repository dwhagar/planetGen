# tests/dbpool.py

"""
One database per test worker, emptied between tests (TEST.139 stage 3).

`mysql_config` used to create a database, apply `schema.sql` (122 tables) and drop it again for every
database test: 0.66 s on an idle machine, about a third of a database test's whole time, and the part that
doubles under load. Here a worker keeps one database with the schema already in it, and after each test
puts it back the way a new one is:

- every table the test put rows in is emptied (`DELETE`, foreign keys off), and an `AUTO_INCREMENT` the test
  advanced is set back to 1, so a test still sees row id 1 first;
- the rows a new database has (the baseline `schema_migrations` row) are put back;
- this process's cached facts about the database (id blocks, facet and count caches, full-text sizes, queue
  connections) are dropped, because they are keyed by the database name and the next test reuses the name;
- the connection pool is closed, so no session variable, lock or open transaction reaches the next test.

A test that changed the schema (a table, column, index, view, routine or trigger added, dropped or altered)
leaves a database that cannot be reused, so that database is dropped and the next test gets a new one. A test
that needs an empty database with no schema at all (it checks what a first connection does) is marked
`fresh_database` and gets the old create-and-drop fixture. `PLANETGEN_TEST_FRESH_DB=1` does the same for every
test, which is the way to tell a problem here from a problem in a test.

`PLANETGEN_TEST_DEEP_CHECK=1` makes every reset prove itself: it compares the full `SHOW CREATE TABLE` of every
table, the row counts and the auto-increment counters with the new database it started as, and fails the
test that left the difference.
"""

import os
import re
import sys
import threading
import traceback
import uuid
import warnings

import pymysql
import sqlalchemy.event
import sqlalchemy.pool

from planetgen.db import store as _db
from planetgen.db.store import MySQLConfig

# Module-level tables in planetgen that must survive a reset: the schema is still in place, and these hold
# facts about the server rather than about the data.
KEEP_CACHES = ("_pools", "_schema_ensured", "_fk_ranks", "_fk_same_rank_parents", "_max_packets",
               "_insert_shapes")


def _related(key, database, name):
    """Whether a cache key belongs to this database: its `_key()` tuple, a tuple starting with it, or its name."""
    if key == name or key == database:
        return True
    if isinstance(key, tuple):
        return bool(key) and (key[0] == database or key[0] == name or database in key)
    return False


def forget_database(config):
    """Drops every cache in `planetgen` that holds something about `config`'s database."""
    database, name = config._key(), config.database
    for module_name, module in list(sys.modules.items()):
        if module is None or not module_name.startswith("planetgen.") or module_name.startswith("planetgen.cli"):
            continue
        for attribute, value in list(vars(module).items()):
            if attribute in KEEP_CACHES or not attribute.startswith("_") or not isinstance(value, (dict, set)):
                continue
            try:
                stale = [key for key in list(value) if _related(key, database, name)]
            except RuntimeError:  # changed while we looked; the next reset will catch it
                continue
            for key in stale:
                if isinstance(value, dict):
                    value.pop(key, None)
                else:
                    value.discard(key)
    _db.forget_id_blocks(database)
    _forget_redis(name)


def _forget_redis(name):
    """Deletes what `planetgen` stored in Redis under this database's name (the stored counts, for one): the
    next test has the same name, and a stored answer from the last one would be read as its own."""
    url = os.environ.get("PLANETGEN_TEST_REDIS_URL")
    if not url:
        return
    try:
        import redis
        client = redis.Redis.from_url(url, socket_timeout=2, socket_connect_timeout=2)
        for key in client.scan_iter(match=f"*{name}*", count=1000):
            client.delete(key)
    except Exception:  # no Redis, or none that answers: nothing of ours is in it
        pass


_checked_out = {}
"""server connection id -> (thread name, stack) of every pooled connection a caller holds right now (TEST.193): when
a database cannot be dropped because a connection still holds it, this names the code that took it."""


def _note_checkout(dbapi_connection, _record, _proxy):
    try:
        _checked_out[dbapi_connection.thread_id()] = (
            threading.current_thread().name,
            "".join(traceback.format_list(traceback.extract_stack(limit=14)[:-2])))
    except Exception:  # a connection without an id (a fake): nothing to note
        pass


def _note_checkin(dbapi_connection, _record):
    try:
        _checked_out.pop(dbapi_connection.thread_id(), None)
    except Exception:
        pass


def watch_checkouts():
    """Starts noting which code holds each pooled connection (once per process)."""
    if not getattr(watch_checkouts, "started", False):
        sqlalchemy.event.listen(sqlalchemy.pool.Pool, "checkout", _note_checkout)
        sqlalchemy.event.listen(sqlalchemy.pool.Pool, "checkin", _note_checkin)
        watch_checkouts.started = True


LOCK_WAIT_TIMEOUT = 1205
"""int: MySQL's error for a statement that waited `lock_wait_timeout` seconds for a lock."""

DROP_WAIT_SECONDS = 10
"""How long a `DROP DATABASE` of a test database waits for a lock before it fails (TEST.193). Without a limit the
server waits as long as another connection holds a table of the database (it waited 853 s in a gate run)."""


GENTLE_WAIT_SECONDS = 3
"""How long `drop_database` waits for the connections that merely read from a database to let go before it ends them
too (TEST.262): another worker's `list_databases` counts the rows of every test database and holds each table's
metadata lock for a moment, and a drop that ended it failed that worker's test with "Lost connection to MySQL
server during query"."""


def database_holders(cur, name, owned_only=False):
    """`[(connection id, user, host, command, idle seconds, state)]` of the other connections that use or hold a lock
    on database `name`: those whose current database it is, and those holding a granted metadata lock on one of its
    tables (a connection that did `USE` something else but read from it; needs the server's performance schema,
    which some servers run without, and then only the first kind is found). With `owned_only` only the first kind."""
    found = {}
    cur.execute(
        "SELECT ID, USER, HOST, COMMAND, TIME, STATE FROM information_schema.PROCESSLIST"
        " WHERE DB = %s AND ID <> CONNECTION_ID()", (name,))
    for row in cur.fetchall():
        found[row[0]] = tuple(row)
    if owned_only:
        return sorted(found.values())
    try:
        cur.execute(
            "SELECT DISTINCT t.PROCESSLIST_ID, t.PROCESSLIST_USER, t.PROCESSLIST_HOST, t.PROCESSLIST_COMMAND,"
            " t.PROCESSLIST_TIME, t.PROCESSLIST_STATE FROM performance_schema.metadata_locks m"
            " JOIN performance_schema.threads t ON t.THREAD_ID = m.OWNER_THREAD_ID"
            " WHERE m.OBJECT_SCHEMA = %s AND m.LOCK_STATUS = 'GRANTED' AND t.PROCESSLIST_ID IS NOT NULL"
            " AND t.PROCESSLIST_ID <> CONNECTION_ID()", (name,))
        for row in cur.fetchall():
            found.setdefault(row[0], tuple(row))
    except pymysql.MySQLError:  # no performance schema: the first list is all there is
        pass
    return sorted(found.values())


def drop_database(server_kwargs, name, wait_seconds=DROP_WAIT_SECONDS):
    """Drops test database `name`, first ending the connections whose current database it is (TEST.193), and never
    waits more than `wait_seconds` for a lock. A connection left open by a thread that outlived its test (a server
    thread still writing a response the browser gave up on, a pooled connection checked out) would otherwise keep
    the drop waiting for minutes. A connection of another database that holds one of its tables (another worker
    counting the rows of every test database) is given `GENTLE_WAIT_SECONDS` to let go first and is ended only if
    it does not (TEST.262), so that worker's query is not cut off.

    Returns:
        list[tuple]: The connections it ended, as `database_holders` lists them (empty when none was open)."""
    conn = pymysql.connect(**server_kwargs)
    ended = []
    try:
        with conn.cursor() as cur:
            gentle = min(GENTLE_WAIT_SECONDS, int(wait_seconds))
            for attempt, (owned_only, wait) in enumerate(((True, gentle), (False, int(wait_seconds)))):
                cur.execute("SET SESSION lock_wait_timeout = %s", (wait,))
                for holder in database_holders(cur, name, owned_only=owned_only):
                    try:
                        cur.execute(f"KILL {int(holder[0])}")
                        ended.append(holder)
                    except pymysql.MySQLError:  # it ended on its own meanwhile
                        pass
                try:
                    cur.execute(f"DROP DATABASE IF EXISTS `{name}`")
                    break
                except pymysql.err.OperationalError as error:
                    if error.args[0] != LOCK_WAIT_TIMEOUT or attempt:
                        raise
        conn.commit()
    finally:
        conn.close()
    if ended:
        lines = []
        for h in ended:
            taken = _checked_out.get(h[0])
            where = f"\n  taken in thread {taken[0]} by:\n{taken[1]}" if taken else " (not from a pool)"
            lines.append(f"id {h[0]} {h[1]}@{h[2]} {h[3]} idle {h[4]} s state {h[5]}{where}")
        warnings.warn(f"TEST.193: {len(ended)} connection(s) still held {name} when its fixture ended: "
                      + "; ".join(lines))
    return ended


class WorkerDatabase:
    """The database one test process reuses."""

    def __init__(self, server_kwargs, deep_check=False):
        self.server = dict(server_kwargs)
        self.deep_check = deep_check
        self.config = None
        self.reference = None  # what a new database looks like once the schema is in it
        self.baseline = None  # the tables' creation times when the current test began
        self.dropped = 0  # databases thrown away because a test changed the schema
        self.resets = 0

    # ------------------------------------------------------------- the server ----

    def _admin(self, **extra):
        return pymysql.connect(**self.server, **extra)

    def _create(self):
        name = f"planetgen_test_{uuid.uuid4().hex[:16]}"
        conn = self._admin()
        try:
            with conn.cursor() as cur:
                cur.execute(f"CREATE DATABASE `{name}`")
            conn.commit()
        finally:
            conn.close()
        config = MySQLConfig(database=name, **self.server)
        connection = _db.get_connection(config)  # applies the schema and the baseline row
        connection.close()
        self.config = config
        return config

    def _drop(self):
        config, self.config = self.config, None
        if config is None:
            return
        _db.close_pool(config)
        forget_database(config)
        drop_database(self.server, config.database)

    def _inspect(self, cur, name):
        """(structure, created, increments): the database's shape without row data."""
        cur.execute("SELECT TABLE_NAME, TABLE_TYPE, ENGINE, CREATE_TIME FROM information_schema.TABLES "
                    "WHERE TABLE_SCHEMA = %s ORDER BY TABLE_NAME", (name,))
        tables = cur.fetchall()
        cur.execute("SELECT TABLE_NAME, COUNT(*), SUM(CRC32(CONCAT_WS('|', COLUMN_NAME, COLUMN_TYPE, IS_NULLABLE, "
                    "COALESCE(COLUMN_DEFAULT, '~'), EXTRA))) FROM information_schema.COLUMNS "
                    "WHERE TABLE_SCHEMA = %s GROUP BY TABLE_NAME ORDER BY TABLE_NAME", (name,))
        columns = cur.fetchall()
        cur.execute("SELECT TABLE_NAME, COUNT(*), SUM(CRC32(CONCAT_WS('|', INDEX_NAME, COLUMN_NAME, NON_UNIQUE, "
                    "SEQ_IN_INDEX))) FROM information_schema.STATISTICS WHERE TABLE_SCHEMA = %s "
                    "GROUP BY TABLE_NAME ORDER BY TABLE_NAME", (name,))
        indexes = cur.fetchall()
        extras = []
        for query in ("SELECT ROUTINE_NAME FROM information_schema.ROUTINES WHERE ROUTINE_SCHEMA = %s",
                      "SELECT TRIGGER_NAME FROM information_schema.TRIGGERS WHERE TRIGGER_SCHEMA = %s",
                      "SELECT CONSTRAINT_NAME FROM information_schema.TABLE_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = %s"):
            cur.execute(query + " ORDER BY 1", (name,))
            extras.append(cur.fetchall())
        structure = ([t[:3] for t in tables], columns, indexes, extras)
        return structure, {t[0]: t[3] for t in tables}

    def _auto_increments(self, cur, name):
        cur.execute("SET SESSION information_schema_stats_expiry = 0")  # MySQL 8 caches these otherwise
        cur.execute("SELECT TABLE_NAME FROM information_schema.TABLES WHERE TABLE_SCHEMA = %s "
                    "AND TABLE_TYPE = 'BASE TABLE' AND AUTO_INCREMENT > 1", (name,))
        return [row[0] for row in cur.fetchall()]

    def _tables_with_rows(self, cur, tables):
        if not tables:
            return []
        probe = " UNION ALL ".join(f"SELECT '{t}' FROM DUAL WHERE EXISTS (SELECT 1 FROM `{t}`)" for t in tables)
        cur.execute(probe)
        return [row[0] for row in cur.fetchall()]

    # ------------------------------------------------------------ the lifecycle ----

    def lease(self):
        """The database for the next test: empty, with the schema in place."""
        if self.config is None:
            self._create()
            conn = self._admin(database=self.config.database)
            try:
                with conn.cursor() as cur:
                    structure, created = self._inspect(cur, self.config.database)
                    tables = [t[0] for t in structure[0] if t[1] == "BASE TABLE"]
                    rows = {}
                    for table in self._tables_with_rows(cur, tables):
                        cur.execute(f"SELECT * FROM `{table}`")
                        rows[table] = cur.fetchall()
                    self.reference = {"structure": structure, "tables": tables, "rows": rows}
                    if self.deep_check:
                        self.reference["detail"] = self._detail(cur, tables)
            finally:
                conn.close()
        _db._schema_ensured.add(self.config._key())
        conn = self._admin(database=self.config.database)
        try:
            with conn.cursor() as cur:
                self.baseline = self._inspect(cur, self.config.database)[1]
        finally:
            conn.close()
        return self.config

    def _detail(self, cur, tables):
        detail = {}
        for table in tables:
            cur.execute(f"SHOW CREATE TABLE `{table}`")
            detail[table] = _without_counter(cur.fetchone()[1])
        return detail

    def release(self):
        """Puts the database back to how a new one is. Returns "" when it did, or why it threw the database
        away; with the deep check on, raises AssertionError describing a difference the reset left."""
        config = self.config
        if config is None:
            return "none"
        _db.close_pool(config)
        forget_database(config)
        conn = self._admin(database=config.database, read_timeout=120, write_timeout=120)
        try:
            with conn.cursor() as cur:
                cur.execute("SET SESSION lock_wait_timeout = 10")
                structure, created = self._inspect(cur, config.database)
                if structure != self.reference["structure"] or created != self.baseline:
                    reason = "the schema changed"
                else:
                    reason = ""
                if not reason:
                    cur.execute("SET SESSION FOREIGN_KEY_CHECKS = 0")
                    tables = self.reference["tables"]
                    for table in self._tables_with_rows(cur, tables):
                        cur.execute(f"DELETE FROM `{table}`")
                    for table, rows in self.reference["rows"].items():
                        marks = ", ".join(["%s"] * len(rows[0]))
                        cur.executemany(f"INSERT INTO `{table}` VALUES ({marks})", rows)
                    for table in self._auto_increments(cur, config.database):
                        if table not in self.reference["rows"]:
                            cur.execute(f"ALTER TABLE `{table}` AUTO_INCREMENT = 1")
                    cur.execute("SET SESSION FOREIGN_KEY_CHECKS = 1")
                    conn.commit()
                    self.resets += 1
                    if self.deep_check:
                        self._verify(cur, tables)
        except pymysql.MySQLError as exc:
            reason = f"the reset failed ({exc.args[0] if exc.args else exc})"
        finally:
            conn.close()
        if reason:
            self.dropped += 1
            self._drop()
        return reason

    def _verify(self, cur, tables):
        problems = []
        detail = self._detail(cur, tables)
        for table in tables:
            if detail[table] != self.reference["detail"][table]:
                problems.append(f"{table}: definition differs")
        have = self._tables_with_rows(cur, tables)
        if sorted(have) != sorted(self.reference["rows"]):
            problems.append(f"tables with rows: {sorted(have)} (a new database: {sorted(self.reference['rows'])})")
        left = self._auto_increments(cur, self.config.database)
        extra = [t for t in left if t not in self.reference["rows"]]
        if extra:
            problems.append(f"auto-increment not reset: {extra}")
        for table, rows in self.reference["rows"].items():
            cur.execute(f"SELECT * FROM `{table}`")
            if cur.fetchall() != rows:
                problems.append(f"{table}: rows differ from a new database")
        if problems:
            raise AssertionError("the database was not put back as a new one is: " + "; ".join(problems))

    def close(self):
        self._drop()


def _without_counter(create_statement):
    """A `SHOW CREATE TABLE` without its ` AUTO_INCREMENT=n` table option, which is checked on its own."""
    return re.sub(r"\s+AUTO_INCREMENT=\d+", "", create_statement)
