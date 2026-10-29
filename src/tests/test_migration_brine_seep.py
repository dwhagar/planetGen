# tests/test_migration_brine_seep.py

"""
GEN.908: migration 0176 widens the `hydrosphere` CHECK of `planets` and `moons` with `brine seep` and relabels a stored liquid surface of
chloride brine, which the water state called a `surface ocean`. Nothing else about a body changes.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _clause(conn, name):
    row = conn.execute("SELECT CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS WHERE CONSTRAINT_SCHEMA = DATABASE()"
                       " AND CONSTRAINT_NAME = ?", (name,)).fetchone()
    return next(iter(row.values())) if row else None


def _stop_enforcing_other_checks(conn, table, keep):
    """The row below holds placeholders, which the table's other CHECKs would refuse; only the hydrosphere's is the test's business."""
    for row in conn.execute("SELECT CONSTRAINT_NAME FROM information_schema.TABLE_CONSTRAINTS WHERE TABLE_SCHEMA = DATABASE()"
                            " AND TABLE_NAME = ? AND CONSTRAINT_TYPE = 'CHECK'", (table,)).fetchall():
        name = row["CONSTRAINT_NAME"]
        if name != keep:
            conn.execute(f"ALTER TABLE {table} ALTER CHECK {name} NOT ENFORCED")


def _insert(conn, table, row_id, hydrosphere, ocean_class):
    """A row of `table` holding only what its NOT NULL columns need (foreign keys off for the test)."""
    columns = conn.execute(
        "SELECT COLUMN_NAME, DATA_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?"
        " AND IS_NULLABLE = 'NO' AND COLUMN_DEFAULT IS NULL AND INSTR(EXTRA, 'auto_increment') = 0"
        " AND INSTR(EXTRA, 'GENERATED') = 0", (table,)).fetchall()
    values = {"id": row_id}
    for column in columns:
        name, kind = column["COLUMN_NAME"], column["DATA_TYPE"]
        values.setdefault(name, 0 if kind in ("int", "bigint", "smallint", "tinyint", "double", "float", "decimal") else "x")
    values.update(hydrosphere=hydrosphere, ocean_class=ocean_class)
    names = ", ".join(f"`{name}`" for name in values)
    conn.execute(f"INSERT INTO {table} ({names}) VALUES ({', '.join('?' * len(values))})", tuple(values.values()))


def test_a_v175_galaxy_takes_the_new_state_and_relabels_its_brine_seeps(mysql_config):
    load_old_schema(mysql_config, 175)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for name in ("planets_chk_5", "moons_chk_5"):
            assert "brine seep" not in _clause(conn, name)
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        for table, offset in (("planets", 0), ("moons", 10)):
            _stop_enforcing_other_checks(conn, table, f"{table}_chk_5")
            _insert(conn, table, offset + 1, "surface ocean", "chloride brine")
            _insert(conn, table, offset + 2, "surface ocean", "neutral")
            _insert(conn, table, offset + 3, "ice-covered ocean", "chloride brine")
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for name in ("planets_chk_5", "moons_chk_5"):
            assert "brine seep" in _clause(conn, name)
        for table, offset in (("planets", 0), ("moons", 10)):
            found = {row["id"]: row["hydrosphere"] for row in conn.execute(f"SELECT id, hydrosphere FROM {table} ORDER BY id").fetchall()}
            assert found == {offset + 1: "brine seep", offset + 2: "surface ocean", offset + 3: "ice-covered ocean"}, table
    finally:
        conn.close()
