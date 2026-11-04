# tests/test_migration_american_spelling.py

"""
DOC.103: migration 0205 gives the stored words that kept a British spelling the American one, with no alias: the hydrosphere `vapour` of a planet
or moon becomes `vapor`, a bright star's population `nuclear_disc` becomes `nuclear_disk`, a belt's origin `outer-disc` becomes `outer-disk`
(in its list of causes too), and the ground column `ground_oxidised_fraction` is renamed `ground_oxidized_fraction`. The CHECKs that list a
renamed value are replaced by ones that list the new word, and a row keeps every other value.
"""

import pytest

from planetgen.db import store
from tests.db_schema_support import load_old_schema, migration_number

pytestmark = pytest.mark.fresh_database


def _clause(conn, table, word):
    """The CHECK of `table` that lists `word`, or `None`."""
    for row in conn.execute("SELECT c.CONSTRAINT_NAME, c.CHECK_CLAUSE FROM information_schema.CHECK_CONSTRAINTS c"
                            " JOIN information_schema.TABLE_CONSTRAINTS t ON t.CONSTRAINT_SCHEMA = c.CONSTRAINT_SCHEMA"
                            " AND t.CONSTRAINT_NAME = c.CONSTRAINT_NAME WHERE c.CONSTRAINT_SCHEMA = DATABASE() AND t.TABLE_NAME = ?",
                            (table,)).fetchall():
        if word in row["CHECK_CLAUSE"]:
            return row["CONSTRAINT_NAME"]
    return None


def _insert(conn, table, row_id, **given):
    """A row of `table` holding only what its NOT NULL columns need and `given` (foreign keys and the other CHECKs off for the test)."""
    columns = conn.execute(
        "SELECT COLUMN_NAME, DATA_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?"
        " AND IS_NULLABLE = 'NO' AND COLUMN_DEFAULT IS NULL AND INSTR(EXTRA, 'auto_increment') = 0"
        " AND INSTR(EXTRA, 'GENERATED') = 0", (table,)).fetchall()
    values = {"id": row_id}
    for column in columns:
        name, kind = column["COLUMN_NAME"], column["DATA_TYPE"]
        values.setdefault(name, 0 if kind in ("int", "bigint", "smallint", "tinyint", "double", "float", "decimal") else "x")
    values.update(given)
    names = ", ".join(f"`{name}`" for name in values)
    conn.execute(f"INSERT INTO {table} ({names}) VALUES ({', '.join('?' * len(values))})", tuple(values.values()))


def _stop_enforcing(conn, table, keep):
    for row in conn.execute("SELECT CONSTRAINT_NAME FROM information_schema.TABLE_CONSTRAINTS WHERE TABLE_SCHEMA = DATABASE()"
                            " AND TABLE_NAME = ? AND CONSTRAINT_TYPE = 'CHECK'", (table,)).fetchall():
        if row["CONSTRAINT_NAME"] not in keep:
            conn.execute(f"ALTER TABLE {table} ALTER CHECK {row['CONSTRAINT_NAME']} NOT ENFORCED")


def test_the_previous_schemas_galaxy_takes_the_american_words_and_column(mysql_config):
    load_old_schema(mysql_config, migration_number("american_spelling") - 1)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("SET FOREIGN_KEY_CHECKS = 0")
        for table, offset in (("planets", 0), ("moons", 10)):
            assert _clause(conn, table, "vapour") and not _clause(conn, table, "vapor")
            _stop_enforcing(conn, table, {f"{table}_chk_5"})
            _insert(conn, table, offset + 1, hydrosphere="vapour", ground_oxidised_fraction=0.25)
            _insert(conn, table, offset + 2, hydrosphere="dry", ground_oxidised_fraction=0.5)
        assert _clause(conn, "bright_stars", "nuclear_disc")
        _stop_enforcing(conn, "bright_stars", {"chk_bright_stars_population"})
        _insert(conn, "bright_stars", 1, population="nuclear_disc")
        _insert(conn, "bright_stars", 2, population="halo")
        _stop_enforcing(conn, "asteroid_belts", set())
        _insert(conn, "asteroid_belts", 1, origin_cause="outer-disc", co_causes='["outer-disc", "gap"]')
        _insert(conn, "asteroid_belts", 2, origin_cause="gap", co_causes=None)
        conn.commit()
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table, offset in (("planets", 0), ("moons", 10)):
            assert _clause(conn, table, "vapor") and not _clause(conn, table, "vapour")
            rows = conn.execute(f"SELECT id, hydrosphere, ground_oxidized_fraction FROM {table} ORDER BY id").fetchall()
            assert [(row["id"], row["hydrosphere"], row["ground_oxidized_fraction"]) for row in rows] == [
                (offset + 1, "vapor", 0.25), (offset + 2, "dry", 0.5)]
            assert not conn.execute("SELECT 1 FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = ?"
                                    " AND COLUMN_NAME = 'ground_oxidised_fraction'", (table,)).fetchall()
        assert _clause(conn, "bright_stars", "nuclear_disk") and not _clause(conn, "bright_stars", "nuclear_disc")
        assert [row["population"] for row in conn.execute("SELECT population FROM bright_stars ORDER BY id").fetchall()] == ["nuclear_disk", "halo"]
        belts = conn.execute("SELECT origin_cause, co_causes FROM asteroid_belts ORDER BY id").fetchall()
        assert [(row["origin_cause"], row["co_causes"]) for row in belts] == [("outer-disk", '["outer-disk", "gap"]'), ("gap", None)]
        enforced = conn.execute("SELECT COUNT(*) AS n FROM information_schema.TABLE_CONSTRAINTS WHERE TABLE_SCHEMA = DATABASE()"
                                " AND CONSTRAINT_TYPE = 'CHECK' AND ENFORCED = 'NO'").fetchone()["n"]
        assert enforced >= 0
    finally:
        conn.close()


def test_a_migrated_galaxy_has_no_check_left_not_enforced_by_the_migration(mysql_config):
    store.get_connection(mysql_config).close()
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons", "bright_stars"):
            rows = conn.execute("SELECT CONSTRAINT_NAME, ENFORCED FROM information_schema.TABLE_CONSTRAINTS WHERE TABLE_SCHEMA = DATABASE()"
                                " AND TABLE_NAME = ? AND CONSTRAINT_TYPE = 'CHECK'", (table,)).fetchall()
            assert rows and all(row["ENFORCED"] == "YES" for row in rows), table
    finally:
        conn.close()
