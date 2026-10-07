# tests/test_db_check_constraints.py

"""
TEST.14: every CHECK constraint in the schema rejects a bad row on the
engine the suite runs on (CI runs MySQL 8.0, 8.4 and MariaDB 11.4). The
constraints are read from information_schema, so a new one is covered
without touching this file; each is broken on a real stored row of a
small generated galaxy (`rich_galaxy_support`).
"""

import re

import pymysql
import pytest

from planetgen.db import store as _db
from tests.db_schema_support import check_constraints
from tests.rich_galaxy_support import rich_galaxy  # noqa: F401  (fixture)

pytestmark = pytest.mark.db

CHECK_VIOLATED = (3819, 4025)
"""MySQL's "Check constraint ... is violated", MariaDB's "CONSTRAINT ... failed"."""

_IN_LIST = re.compile(r"^`?(\w+)`?\s+in\s*\(", re.IGNORECASE)
_EQUALS_ONE = re.compile(r"^`?(\w+)`?\s*=\s*1$")
_NULL_TOGETHER = re.compile(r"`?(\w+)`?\s+is\s+null\s*\)?\s*=", re.IGNORECASE)

_EMPTY_TABLE_ROWS = {
    # Only a habitable moon gets these; copied from a planet's.
    "moon_reflection_spectrum": (
        "INSERT INTO moon_reflection_spectrum (moon_id, spectrum_type, position, value)"
        " VALUES ((SELECT MIN(id) FROM moons), 'visible', 0, 'blue')"),
}


def _breaking_value(conn, table, clause, row):
    """`(column, value)` that the clause refuses for `row`."""
    clause = " ".join(clause.strip().strip("()").split())
    match = _IN_LIST.match(clause)
    if match:
        column = match[1]
        kind = conn.execute(
            "SELECT DATA_TYPE AS k FROM information_schema.COLUMNS"
            " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ? AND COLUMN_NAME = ?", (table, column),
        ).fetchone()["k"]
        return column, (7 if "int" in kind else "zz")
    match = _EQUALS_ONE.match(clause)
    if match:
        return match[1], 2
    match = _NULL_TOGETHER.search(clause)
    if match:
        column = match[1]
        assert column in row, f"{table}: {column} from ({clause}) not in {sorted(row)}"
        return column, (None if row[column] is not None else 0.0)
    raise AssertionError(f"no way to break {table}'s CHECK ({clause}) -- teach _breaking_value its shape")


def test_every_check_rejects_a_bad_row(rich_galaxy):
    conn = _db.get_connection(rich_galaxy)
    try:
        checks = check_constraints(conn)
        assert len({row["t"] for row in checks}) >= 25
        failures = []
        for check in checks:
            table, name = check["t"], check["n"]
            conn.execute("SET SESSION FOREIGN_KEY_CHECKS = 0")
            try:
                if table in _EMPTY_TABLE_ROWS and conn.execute(f"SELECT 1 FROM {table} LIMIT 1").fetchone() is None:
                    conn.execute(_EMPTY_TABLE_ROWS[table])
                row = conn.execute(f"SELECT * FROM {table} ORDER BY id LIMIT 1").fetchone()
                assert row is not None, f"the seeded galaxy has no {table} row to break {name} on"
                column, value = _breaking_value(conn, table, check["clause"], row)
                try:
                    conn.execute(f"UPDATE {table} SET {column} = ? WHERE id = ?", (value, row["id"]))
                except pymysql.err.MySQLError as exc:
                    if exc.args[0] not in CHECK_VIOLATED:
                        failures.append(f"{table}.{name}: {column} = {value!r} failed with {exc.args}")
                else:
                    failures.append(f"{table}.{name} ({check['clause']}) let {column} = {value!r} through")
            finally:
                conn.rollback()
                conn.execute("SET SESSION FOREIGN_KEY_CHECKS = 1")
        assert not failures, "\n".join(failures)
    finally:
        conn.close()


def test_mysql_enforces_every_check(rich_galaxy):
    """MySQL can store a CHECK as NOT ENFORCED; none of ours may be."""
    conn = _db.get_connection(rich_galaxy)
    try:
        if "mariadb" in conn.execute("SELECT VERSION() AS v").fetchone()["v"].lower():
            pytest.skip("MariaDB always enforces CHECK constraints")
        loose = conn.execute(
            "SELECT TABLE_NAME AS t, CONSTRAINT_NAME AS n FROM information_schema.TABLE_CONSTRAINTS"
            " WHERE TABLE_SCHEMA = DATABASE() AND CONSTRAINT_TYPE = 'CHECK' AND ENFORCED <> 'YES'"
        ).fetchall()
    finally:
        conn.close()
    assert loose == []
