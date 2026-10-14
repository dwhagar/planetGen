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
import uuid

import pymysql

from planetgen.db import store as _db
from planetgen.db.store import MySQLConfig

# Module-level tables in planetgen that must survive a reset: the schema is still in place, and these hold
# facts about the server rather than about the data.
KEEP_CACHES = ("_pools", "_schema_ensured", "_fk_ranks", "_fk_same_rank_parents", "_max_packets",
               "_insert_shapes", "_mariadb_servers")


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
        conn = self._admin()
        try:
            with conn.cursor() as cur:
                cur.execute(f"DROP DATABASE IF EXISTS `{config.database}`")
            conn.commit()
        finally:
            conn.close()

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
        try:
            cur.execute("SET SESSION information_schema_stats_expiry = 0")  # MySQL 8 caches these otherwise
        except pymysql.MySQLError:  # MariaDB has no such variable and always reads them live
            pass
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
