# tests/db_schema_support.py

"""
Helpers for the migration tests (TEST.8, TEST.9): load an old released
schema into an empty database, and describe a database's whole shape so
two databases can be compared.
"""

import contextlib
import gzip
import os
import uuid

import pymysql

from stellarObjects import _db

OLD_SCHEMA_DIR = os.path.join(os.path.dirname(__file__), "fixtures", "old_schemas")


def old_schema_versions():
    """Every version with a checked-in schema, oldest first."""
    versions = []
    for name in os.listdir(OLD_SCHEMA_DIR):
        if name.startswith("schema_v") and name.endswith(".sql.gz"):
            versions.append(int(name[len("schema_v"):-len(".sql.gz")]))
    return sorted(versions)


def load_old_schema(config, version):
    """Creates `version`'s tables in the empty database `config` names,
    recorded as that version, the way that release's `_ensure_schema`
    left a new database."""
    with gzip.open(os.path.join(OLD_SCHEMA_DIR, f"schema_v{version}.sql.gz"), "rt", encoding="utf-8") as handle:
        script = handle.read()
    conn = _db.get_connection(config, ensure_schema=False)
    try:
        conn.executescript(script)
        conn.execute("INSERT INTO schema_migrations (version) VALUES (?)", (version,))
        conn.commit()
    finally:
        conn.close()


@contextlib.contextmanager
def scratch_database(beside):
    """A second throwaway database on the same server as the `MySQLConfig`
    `beside` (a test's `mysql_config`), dropped after."""
    server_kwargs = dict(host=beside.host, port=beside.port, user=beside.user, password=beside.password)
    name = f"planetgen_test_{uuid.uuid4().hex[:16]}"
    admin = pymysql.connect(**server_kwargs)
    try:
        with admin.cursor() as cur:
            cur.execute(f"CREATE DATABASE `{name}`")
    finally:
        admin.close()
    config = _db.MySQLConfig(database=name, **server_kwargs)
    try:
        yield config
    finally:
        _db.close_pool(config)
        admin = pymysql.connect(**server_kwargs)
        try:
            with admin.cursor() as cur:
                cur.execute(f"DROP DATABASE IF EXISTS `{name}`")
        finally:
            admin.close()


def schema_snapshot(config):
    """
    The shape of the database `config` names, from information_schema:
    tables and views, each column (type, nullability, default, extra,
    collation), each index (uniqueness, kind, columns in order), each
    foreign key (columns, target, rules) and each table's CHECK clauses
    (by clause: an unnamed CHECK's automatic name depends on the engine
    and on the order it was added in).
    """
    conn = _db.get_connection(config, ensure_schema=False)
    try:
        def rows(sql):
            return conn.execute(sql).fetchall()

        tables = {row["t"]: row["k"] for row in rows(
            "SELECT TABLE_NAME AS t, TABLE_TYPE AS k FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()")}
        columns = {(row["t"], row["c"]): (row["ty"], row["n"], row["d"], row["e"], row["coll"]) for row in rows(
            "SELECT TABLE_NAME AS t, COLUMN_NAME AS c, COLUMN_TYPE AS ty, IS_NULLABLE AS n, COLUMN_DEFAULT AS d,"
            " EXTRA AS e, COLLATION_NAME AS coll FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()")}
        indexes = {}
        for row in rows("SELECT TABLE_NAME AS t, INDEX_NAME AS i, NON_UNIQUE AS u, SEQ_IN_INDEX AS s,"
                        " COLUMN_NAME AS c, INDEX_TYPE AS ty FROM information_schema.STATISTICS"
                        " WHERE TABLE_SCHEMA = DATABASE() ORDER BY TABLE_NAME, INDEX_NAME, SEQ_IN_INDEX"):
            indexes.setdefault((row["t"], row["i"]), (int(row["u"]), row["ty"], []))[2].append(row["c"])
        foreign_keys = {}
        for row in rows("SELECT k.TABLE_NAME AS t, k.CONSTRAINT_NAME AS n, k.COLUMN_NAME AS c,"
                        " k.REFERENCED_TABLE_NAME AS rt, k.REFERENCED_COLUMN_NAME AS rc,"
                        " r.UPDATE_RULE AS ur, r.DELETE_RULE AS dr"
                        " FROM information_schema.KEY_COLUMN_USAGE k JOIN information_schema.REFERENTIAL_CONSTRAINTS r"
                        "   ON r.CONSTRAINT_SCHEMA = k.TABLE_SCHEMA AND r.CONSTRAINT_NAME = k.CONSTRAINT_NAME"
                        "  AND r.TABLE_NAME = k.TABLE_NAME"
                        " WHERE k.TABLE_SCHEMA = DATABASE() AND k.REFERENCED_TABLE_NAME IS NOT NULL"
                        " ORDER BY k.TABLE_NAME, k.CONSTRAINT_NAME, k.ORDINAL_POSITION"):
            foreign_keys.setdefault((row["t"], row["n"]), (row["rt"], row["ur"], row["dr"], []))[3].append(
                (row["c"], row["rc"]))
        checks = {}
        for row in rows("SELECT tc.TABLE_NAME AS t, cc.CHECK_CLAUSE AS cl FROM information_schema.CHECK_CONSTRAINTS cc"
                        " JOIN information_schema.TABLE_CONSTRAINTS tc"
                        "   ON tc.CONSTRAINT_SCHEMA = cc.CONSTRAINT_SCHEMA AND tc.CONSTRAINT_NAME = cc.CONSTRAINT_NAME"
                        "  AND tc.CONSTRAINT_TYPE = 'CHECK'"
                        " WHERE cc.CONSTRAINT_SCHEMA = DATABASE()"):
            checks.setdefault(row["t"], []).append(" ".join(row["cl"].split()))
    finally:
        conn.close()
    return {
        "table": tables,
        "column": columns,
        "index": {key: (unique, kind, tuple(cols)) for key, (unique, kind, cols) in indexes.items()},
        "foreign key": {key: (target, update, delete, tuple(cols))
                        for key, (target, update, delete, cols) in foreign_keys.items()},
        "checks on": {table: sorted(clauses) for table, clauses in checks.items()},
    }


def schema_differences(migrated, fresh):
    """One line per thing the two snapshots disagree on; empty when equal."""
    lines = []
    for part in fresh:
        have, want = migrated[part], fresh[part]
        for key in sorted(set(have) | set(want), key=str):
            if have.get(key) != want.get(key):
                lines.append(f"{part} {key}: migrated {have.get(key)!r}, fresh {want.get(key)!r}")
    return lines
