# tests/test_body_fields.py

"""
GEN.1658: the three body tables (`planets`, `moons`, `rogue_planets`) are declared once in `db/body_fields.py`, and `schema.sql`, the unit groups of `db/column_units.py` and the
column lists are served from it. These tests hold the declaration to the files it writes and to the database a new install builds from them.
"""

import re
import subprocess
import sys

import pytest

from planetgen.db import body_fields, column_units, store
from tests.db_schema_support import scratch_database

ROOT = store.SCHEMA_PATH.rsplit("/src/", 1)[0]


def test_schema_sql_holds_exactly_what_the_declaration_writes():
    """The statements of the three tables in schema.sql are the generated ones, between their marker lines."""
    result = subprocess.run([sys.executable, f"{ROOT}/scripts/build_body_tables.py", "--check"], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_every_table_is_between_its_markers_once():
    text = open(store.SCHEMA_PATH, encoding="utf-8").read()
    for letter, table in body_fields.TABLES.items():
        assert text.count(body_fields.schema_block(letter)) == 1, table
        assert len(re.findall(rf"^CREATE TABLE IF NOT EXISTS {table} \(", text, re.M)) == 1, table


def test_a_new_database_has_the_columns_the_declaration_lists_in_its_order(mysql_config):
    with scratch_database(mysql_config) as config:
        conn = store.get_connection(config)
        try:
            for letter, table in body_fields.TABLES.items():
                found = [row["column_name"] for row in conn.execute(
                    "SELECT column_name AS column_name FROM information_schema.columns WHERE table_schema = DATABASE() AND table_name = ? ORDER BY ordinal_position", (table,)).fetchall()]
                assert found == [field.name for field in body_fields.fields_of(letter)], table
        finally:
            conn.close()


def test_no_column_is_declared_twice_in_a_table_and_a_dict_value_names_only_tables_that_have_it():
    for letter in body_fields.TABLES:
        names = [field.name for field in body_fields.fields_of(letter)]
        assert len(names) == len(set(names)), letter
    for field in body_fields.FIELDS:
        assert field.t and set(field.t) <= set(body_fields.TABLES), field.name
        for attribute in ("type", "nn", "extra", "unit", "c", "lead", "post"):
            value = getattr(field, attribute)
            if isinstance(value, dict):
                letters = "".join(value)
                assert len(letters) == len(set(letters)) and set(letters) == set(field.t), (field.name, attribute)


def test_a_table_with_its_own_order_lists_exactly_its_columns():
    for letter, order in body_fields.ORDERS.items():
        assert sorted(order) == sorted(field.name for field in body_fields.FIELDS if letter in field.t), letter


def test_the_unit_groups_of_the_three_tables_come_from_the_declaration():
    for table in body_fields.TABLES.values():
        assert column_units.TABLES[table] == body_fields.units_by_group(table)
        assert column_units.COLUMNS[(table, "id")] == "id"
    assert body_fields.field_units("moons")["radius_km"] == "km"


def test_sql_of_writes_the_type_then_not_null_then_the_rest():
    assert body_fields.sql_of(body_fields.Field("x", "DOUBLE", nn=True), "p") == "DOUBLE NOT NULL"
    assert body_fields.sql_of(body_fields.Field("y", "VARCHAR(16)", nn=True, extra="DEFAULT 'star'"), "p") == "VARCHAR(16) NOT NULL DEFAULT 'star'"


def test_a_value_that_differs_by_table_is_read_by_the_tables_letter():
    field = body_fields.Field("parent_type", "VARCHAR(16)", extra={"p": "DEFAULT 'star'", "mr": "DEFAULT 'x'"}, nn=True)
    assert body_fields.sql_of(field, "p") == "VARCHAR(16) NOT NULL DEFAULT 'star'"
    assert body_fields.sql_of(field, "r") == "VARCHAR(16) NOT NULL DEFAULT 'x'"
    with pytest.raises(KeyError):
        field.of("extra", "z")
