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

from planetgen.db import store

OLD_SCHEMA_DIR = os.path.join(os.path.dirname(__file__), "fixtures", "old_schemas")

MIGRATIONS_DIR = os.path.join(os.path.dirname(__file__), "..", "planetgen", "db", "migrations", "versions")


def migration_number(stem):
    """The number of the migration file `<number>_<stem>.py`, so a test that starts from the schema before it follows a renumbering."""
    for name in os.listdir(MIGRATIONS_DIR):
        if name[:4].isdigit() and name[5:] == stem + ".py":
            return int(name[:4])
    raise LookupError(stem)


ROWS_BEFORE = (
    ("brine_seep", ("UPDATE planets SET hydrosphere = 'surface ocean' WHERE hydrosphere = 'brine seep'",
                    "UPDATE moons SET hydrosphere = 'surface ocean' WHERE hydrosphere = 'brine seep'")),
)
"""tuple: `(migration stem, statements)`: what a row must hold to be kept in a schema older than that migration, when it widened a CHECK
(GEN.908: the `brine seep` hydrosphere is a `surface ocean` again). `downgrade_to` runs them before it puts the old CHECKs back."""


def _ddl_parts(conn):
    """The shape of the connected database as `{table: {"create": text, "columns": {name: line}, "keys": {name: line},
    "constraints": {name: line}}}` (views as `{"view": text}`), read from `SHOW CREATE TABLE`: both sides of a comparison
    come from the same server, so a line that has not changed is the same text on both."""
    found = {}
    for row in conn.execute("SHOW FULL TABLES").fetchall():
        values = list(row.values())
        name, kind = values[0], values[1]
        if kind == "VIEW":
            found[name] = {"view": conn.execute(f"SHOW CREATE VIEW `{name}`").fetchone()["Create View"]}
            continue
        create = conn.execute(f"SHOW CREATE TABLE `{name}`").fetchone()["Create Table"]
        parts = {"create": create, "columns": {}, "keys": {}, "constraints": {}}
        for line in create.splitlines()[1:-1]:
            line = line.strip().rstrip(",")
            if line.startswith("`"):
                parts["columns"][line.split("`")[1]] = line
            elif line.startswith("CONSTRAINT"):
                parts["constraints"][line.split("`")[1]] = line
            elif line.startswith("PRIMARY KEY"):
                parts["keys"]["PRIMARY"] = line
            else:
                parts["keys"][line.split("`")[1]] = line
        found[name] = parts
    return found


def downgrade_to(config, version):
    """
    Makes the database `config` names, which holds the current schema, the schema of `version`, with its rows kept (TEST.328):
    what a migration test starts from. The old schema is the checked-in one (`fixtures/old_schemas/schema_v<version>.sql.gz`)
    loaded into a scratch database, and every difference is undone in the test's own database (a table, column, key,
    CHECK or foreign key the old schema lacks is dropped, one it has and the database lacks is added back, a column of
    another definition is modified), so no test keeps its own list of what each newer migration added. The database is
    left recorded as `version`, in `alembic_version` and `schema_migrations`. A database with no schema yet is given the current one
    first, as `store.get_connection` does.
    """
    store.get_connection(config).close()
    with scratch_database(config) as scratch:
        load_old_schema(scratch, version)
        conn = store.get_connection(scratch, ensure_schema=False)
        try:
            want = _ddl_parts(conn)
        finally:
            conn.close()
    conn = store.get_connection(config, ensure_schema=False)
    try:
        have = _ddl_parts(conn)
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        for stem, statements in ROWS_BEFORE:
            if version < migration_number(stem):
                for statement in statements:
                    conn.execute(statement)
        for name in sorted(set(have) - set(want)):
            conn.execute(f"DROP {'VIEW' if 'view' in have[name] else 'TABLE'} `{name}`")
        for name in sorted(set(have) & set(want)):
            now, old = have[name], want[name]
            if "view" in now or "view" in old:
                if now.get("view") != old.get("view"):
                    conn.execute(f"DROP VIEW `{name}`")
                    conn.execute(old["view"])
                continue
            clauses = []
            for key, line in now["constraints"].items():
                if old["constraints"].get(key) != line:
                    clauses.append(f"DROP FOREIGN KEY `{key}`" if "FOREIGN KEY" in line else f"DROP CHECK `{key}`")
            for key, line in now["keys"].items():
                if old["keys"].get(key) != line:
                    clauses.append("DROP PRIMARY KEY" if key == "PRIMARY" else f"DROP INDEX `{key}`")
            clauses += [f"DROP COLUMN `{column}`" for column in now["columns"] if column not in old["columns"]]
            clauses += [f"MODIFY COLUMN {line}" for column, line in old["columns"].items()
                        if column in now["columns"] and now["columns"][column] != line]
            clauses += [f"ADD COLUMN {line}" for column, line in old["columns"].items() if column not in now["columns"]]
            clauses += [f"ADD {line}" for key, line in old["keys"].items() if now["keys"].get(key) != line]
            clauses += [f"ADD {line}" for key, line in old["constraints"].items() if now["constraints"].get(key) != line]
            if clauses:
                conn.execute(f"ALTER TABLE `{name}` {', '.join(clauses)}, ALGORITHM=COPY")
        for name in sorted(set(want) - set(have)):
            conn.execute(want[name].get("view") or want[name]["create"])
        conn.execute("SET FOREIGN_KEY_CHECKS = 1")
        conn.execute("UPDATE alembic_version SET version_num = ?", (f"{version:04d}",))
        conn.execute("DELETE FROM schema_migrations")
        conn.execute("INSERT INTO schema_migrations (version) VALUES (?)", (version,))
        conn.commit()
    finally:
        conn.close()


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
    left a new database. Tables already there (a pooled test database holds the
    current schema) are dropped first, so none keeps its newer columns."""
    with gzip.open(os.path.join(OLD_SCHEMA_DIR, f"schema_v{version}.sql.gz"), "rt", encoding="utf-8") as handle:
        script = handle.read()
    conn = store.get_connection(config, ensure_schema=False)
    try:
        tables = [next(iter(row.values())) for row in conn.execute("SHOW FULL TABLES WHERE Table_type = 'BASE TABLE'").fetchall()]
        if tables:
            conn.execute("SET FOREIGN_KEY_CHECKS = 0")
            for table in tables:
                conn.execute(f"DROP TABLE `{table}`")
            conn.execute("SET FOREIGN_KEY_CHECKS = 1")
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
    config = store.MySQLConfig(database=name, **server_kwargs)
    try:
        yield config
    finally:
        store.close_pool(config)
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
    conn = store.get_connection(config, ensure_schema=False)
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
        for row in check_constraints(conn):
            checks.setdefault(row["t"], []).append(" ".join(row["clause"].split()))
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


def check_constraints(conn):
    """Every CHECK in the connected database as `{t, n, clause}` rows. MySQL's CHECK_CONSTRAINTS has no
    TABLE_NAME (constraint names are unique per schema), so the table comes from TABLE_CONSTRAINTS."""
    return conn.execute(
        "SELECT tc.TABLE_NAME AS t, cc.CONSTRAINT_NAME AS n, cc.CHECK_CLAUSE AS clause"
        " FROM information_schema.CHECK_CONSTRAINTS cc JOIN information_schema.TABLE_CONSTRAINTS tc"
        "   ON tc.CONSTRAINT_SCHEMA = cc.CONSTRAINT_SCHEMA AND tc.CONSTRAINT_NAME = cc.CONSTRAINT_NAME"
        "  AND tc.CONSTRAINT_TYPE = 'CHECK'"
        " WHERE cc.CONSTRAINT_SCHEMA = DATABASE() ORDER BY tc.TABLE_NAME, cc.CONSTRAINT_NAME"
    ).fetchall()


PROBE_REVISIONS = [
    "op.execute('CREATE TABLE alembic_probe (id INT PRIMARY KEY)')",
    "op.execute('ALTER TABLE alembic_probe ADD COLUMN note VARCHAR(20)')",
]

_PROBE_TEMPLATE = """from alembic import op

revision = "{revision}"
down_revision = "{previous}"
branch_labels = None
depends_on = None


def upgrade():
    {body}
"""


def migrations_with_probes(tmp_path, count=1):
    """A copy of the real migrations directory plus `count` (1 or 2) test
    revisions after the real head, for exercising Alembic's upgrade path
    whatever the head is. Returns its path."""
    import shutil
    from planetgen.db import alembic_runner
    target = tmp_path / "migrations"
    shutil.copytree(alembic_runner.MIGRATIONS_DIR, target, ignore=shutil.ignore_patterns("__pycache__"))
    head = alembic_runner.head_version()
    for offset in range(count):
        version = head + 1 + offset
        (target / "versions" / f"{version:04d}_probe.py").write_text(_PROBE_TEMPLATE.format(
            revision=f"{version:04d}", previous=f"{version - 1:04d}", body=PROBE_REVISIONS[offset]))
    return str(target)


def table_columns(conn, table):
    """The names of `table`'s columns."""
    return {row["COLUMN_NAME"] for row in conn.execute(
        "SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?",
        (table,)).fetchall()}


def insert_filled_row(conn, table, **given):
    """A row of `table` with every NOT NULL column that has no default filled with a neutral value."""
    rows = conn.execute(
        "SELECT COLUMN_NAME, DATA_TYPE, IS_NULLABLE, COLUMN_DEFAULT, EXTRA FROM information_schema.COLUMNS"
        " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?", (table,)).fetchall()
    values = dict(given)
    for row in rows:
        name = row["COLUMN_NAME"]
        if name in values or row["IS_NULLABLE"] == "YES" or row["COLUMN_DEFAULT"] is not None or "auto_increment" in row["EXTRA"]:
            continue
        kind = row["DATA_TYPE"]
        values[name] = b"\x00" * 12 if kind in ("binary", "varbinary", "blob") else (
            "" if kind in ("char", "varchar", "text") else 0)
    names = ", ".join(f"`{name}`" for name in values)
    conn.execute(f"INSERT INTO {table} ({names}) VALUES ({', '.join('?' for _ in values)})", tuple(values.values()))
