#!/usr/bin/env python3
"""
Writes the three body tables of `src/planetgen/db/schema.sql` from `src/planetgen/db/body_fields.py` (GEN.1658).

    python scripts/build_body_tables.py           rewrite the three statements in schema.sql
    python scripts/build_body_tables.py --check   fail (exit 1) when schema.sql differs from what the declaration writes

Each statement sits between the marker lines `-- BEGIN GENERATED body tables ...: <table>` and `-- END GENERATED body tables: <table>`; the first run, on a schema.sql that
has no markers yet, replaces the three hand-written statements by their generated ones. Nothing outside the markers is touched, so another branch's change to the rest of the
file never conflicts with this one; a conflict inside a marker is resolved by running this script (it reads whatever is between the markers, conflict lines included, and replaces it).
Add the Alembic revision of a new column by hand, then run this and `scripts/generate_db_models.py`.
"""

import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from planetgen.db import body_fields  # noqa: E402

SCHEMA = os.path.join(ROOT, "src", "planetgen", "db", "schema.sql")


def rewritten(text):
    """`text` (the contents of schema.sql) with each body table's statement replaced by the generated one, between its markers."""
    for letter, table in body_fields.TABLES.items():
        marked = re.compile(rf"^-- BEGIN GENERATED body tables[^\n]*: {table}\n.*?^-- END GENERATED body tables: {table}\n", re.S | re.M)
        if marked.search(text):
            text = marked.sub(lambda _match: body_fields.schema_block(letter), text, count=1)
            continue
        plain = re.compile(rf"^CREATE TABLE IF NOT EXISTS {table} \(\n.*?\n\) ENGINE=[^\n]*;\n", re.S | re.M)
        if not plain.search(text):
            raise SystemExit(f"schema.sql holds neither the markers of {table} nor its CREATE TABLE")
        text = plain.sub(lambda _match: body_fields.schema_block(letter), text, count=1)
    return text


def main(argv):
    with open(SCHEMA, encoding="utf-8") as handle:
        current = handle.read()
    wanted = rewritten(current)
    if "--check" in argv:
        if wanted != current:
            print("schema.sql is out of date with db/body_fields.py: run python scripts/build_body_tables.py", file=sys.stderr)
            return 1
        return 0
    if wanted != current:
        with open(SCHEMA, "w", encoding="utf-8") as handle:
            handle.write(wanted)
        print("schema.sql: the body tables rewritten")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
