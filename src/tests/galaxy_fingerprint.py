# tests/galaxy_fingerprint.py

"""
A comparable picture of a galaxy database's generated content (GEN.39):
every row of every content table, without what legitimately differs
between two runs of one seed (row ids and the foreign keys holding them,
timestamps, name-collision decorations), as a sorted list per table.

Not GEN.58's canonical fingerprint (per-sector hashes with canonical
number formatting); just enough for the tests to say "same galaxy".
"""

import re

from stellarObjects import _db
from planetgen.names.wordlists import DIMINUTIVE_PREFIXES, GREEK_LETTERS

_DECORATION = re.compile(r"\b(?:%s) (?=[A-Z])" % "|".join(GREEK_LETTERS + DIMINUTIVE_PREFIXES))
"""A name-collision decoration word (`nameUniqueness`) before a name.
Which of two colliding systems keeps the plain name depends on which
saved first (GEN.57, phase 1), so the picture leaves the decoration out."""

_CONTROL_SCHEMA = _db.os.path.join(_db.os.path.dirname(_db.__file__), "control_schema.sql")

with open(_CONTROL_SCHEMA, encoding="utf-8") as _file:
    _CONTROL_TABLES = frozenset(re.findall(r"CREATE TABLE IF NOT EXISTS (\w+)", _file.read()))

SKIPPED_TABLES = frozenset({
    "schema_migrations", "nearest_systems", "name_registry", "galaxy_shape", "galaxy_layer", "galaxy_column",
    "generation_runs",
}) | _CONTROL_TABLES
"""Tables left out: bookkeeping (the run history among it), the control database's tables (the tests
keep them in the same database), links rebuilt from positions, and the
plan itself (the tests write it)."""

_SKIPPED_TYPES = frozenset({"datetime", "timestamp", "date", "time"})


def _columns(conn, table):
    rows = conn.execute(
        "SELECT column_name AS name, data_type AS type FROM information_schema.columns"
        " WHERE table_schema = DATABASE() AND table_name = ? ORDER BY ordinal_position", (table,)).fetchall()
    return [row["name"] for row in rows
            if row["name"] != "id" and not row["name"].endswith("_id") and row["type"] not in _SKIPPED_TYPES]


def _comparable(value):
    if isinstance(value, str):
        return repr(_DECORATION.sub("", value))
    return repr(value)


def galaxy_rows(config):
    """`{table: sorted list of row tuples}` for every non-empty content
    table of the database at `config`."""
    conn = _db.get_connection(config)
    try:
        tables = [row["name"] for row in conn.execute(
            "SELECT table_name AS name FROM information_schema.tables WHERE table_schema = DATABASE()"
            " AND table_type = 'BASE TABLE'").fetchall()]
        picture = {}
        for table in sorted(tables):
            if table in SKIPPED_TABLES:
                continue
            columns = _columns(conn, table)
            if not columns:
                continue
            rows = conn.execute(f"SELECT {', '.join(f'`{c}`' for c in columns)} FROM `{table}`").fetchall()
            if rows:
                picture[table] = sorted(tuple(_comparable(row[c]) for c in columns) for row in rows)
        return picture
    finally:
        conn.close()
