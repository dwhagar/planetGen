# tests/test_kind_tables.py -- DB.74: the tables that hold one kind of object are listed once (`db/kind_tables.py`), each
# group agrees with the schema, and every other place that names them agrees with the registry.

import pytest

from planetgen import tuning
from planetgen.api import ids
from planetgen.db import edits, fingerprint, gravity_sources, kind_tables, models, point_columns, query, sector_paths, store
from planetgen.galaxy import objectref
from planetgen.web.lib import slowpages

TABLES = models.metadata.tables
ROWS = {row.table: row for row in kind_tables.KIND_TABLES}


def _columns(table):
    return {column.name for column in TABLES[table].columns}


def test_every_row_is_a_table_and_listed_once():
    assert len(ROWS) == len(kind_tables.KIND_TABLES)
    assert not set(ROWS) - set(TABLES)


@pytest.mark.parametrize("group, columns", [
    ("uid", {"uid"}), ("timestamped", {"created_at", "modified_at"}), ("sector_filed", {"sector_id"}),
    ("system_filed", {"star_system_id"}), ("star_filed", {"star_id"}), ("containable", {"inside_nebula_id"}),
    ("phenomenon", {"center_x_pc", "center_y_pc", "center_z_pc", "quadrant", "sector_id"}),
    ("unique_name", {"name"}), ("named_phenomenon", {"name", "center_x_pc"}), ("body", {"star_system_id", "uid"}),
    ("parent_link", {"parent_type", "parent_id"}),
])
def test_a_group_names_tables_that_have_its_columns(group, columns):
    for table in kind_tables.tables(group):
        assert columns <= _columns(table), f"{table} is in {group} and lacks {columns - _columns(table)}"


def test_every_table_with_an_object_id_is_in_the_uid_group():
    with_uid = {name for name, table in TABLES.items() if "uid" in {column.name for column in table.columns}}
    left_over = with_uid - set(kind_tables.tables("uid"))
    assert left_over <= set(point_columns.NOT_BODIES) | set(kind_tables.SCATTER_UID_TABLES), f"{sorted(left_over)} carry a uid and are not in the registry"
    assert set(kind_tables.SCATTER_UID_TABLES) <= with_uid, "a scatter table listed for its uid has none"


def test_a_phenomenon_table_is_a_generator_type_and_a_reference_kind():
    for row in kind_tables.KIND_TABLES:
        assert ("phenomenon" in row.groups) == (row.generator_type is not None)
        assert ("phenomenon" in row.groups) == (row.kind in objectref.PHENOMENON_KINDS)
    assert set(kind_tables.generator_type_tables()) <= set(tuning.PHENOMENON_TYPE_CHOICES)


def test_the_sites_that_used_to_list_the_tables_read_the_registry():
    assert store.ID_BLOCK_TABLES == frozenset(kind_tables.tables("id_block"))
    assert store.CONTAINABLE_TABLES == kind_tables.tables("containable", order="containable")
    assert store.PLACED_PHENOMENON_TABLES == kind_tables.tables("phenomenon", order="placed")
    assert store.PHENOMENON_UID_TABLES == kind_tables.tables("phenomenon")
    assert store.TIMESTAMPED_TABLES == kind_tables.tables("timestamped", order="timestamped")
    assert store.UNIQUE_NAME_TABLES == kind_tables.tables("unique_name", order="unique_name")
    assert store.NAMED_PHENOMENON_TABLES == kind_tables.tables("named_phenomenon", order="unique_name")
    assert store.PARENT_LINK_TABLES == kind_tables.tables("parent_link")
    assert store.BODY_UID_TABLES == store._SYSTEM_CONTENT_TABLES == kind_tables.tables("body")
    assert store.SCATTER_MOTION_TABLES == kind_tables.SCATTER_UID_TABLES
    assert store.STAR_FILED_TABLES == gravity_sources.STAND_ALONE_TABLES == kind_tables.tables("star_filed")
    assert fingerprint.SECTOR_TABLES == kind_tables.tables("sector_filed")
    assert fingerprint.SYSTEM_TABLES == kind_tables.tables("system_filed")
    assert fingerprint.STAR_TABLES == kind_tables.tables("star_filed")
    assert sector_paths.PATH_TABLES == kind_tables.tables("path")
    assert edits.PHENOMENON_TABLES == kind_tables.kind_tables("phenomenon", order="placed")
    assert ids.OBJECT_TABLES == kind_tables.kind_tables("uid", order="ids")
    assert slowpages.PHENOMENON_TABLES == ("phenomenon_scatter", *kind_tables.tables("phenomenon", order="slow"))


def test_the_lists_kept_for_their_order_name_the_same_kinds():
    assert set(objectref.PHENOMENON_KINDS) == set(kind_tables.kind_tables("phenomenon"))
    assert set(ids.PHENOMENON_KINDS) == set(objectref.PHENOMENON_KINDS)
    assert {row[0] for row in query._PHENOMENON_TABLES} == set(kind_tables.tables("phenomenon"))
    assert set(objectref.TABLES) == {"sector", "system", *objectref.BODY_KINDS, *objectref.FACILITY_KINDS}
    assert all(ROWS[table].kind == kind for kind, table in objectref.TABLES.items() if kind != "sector")


def test_the_named_phenomena_are_the_uniquely_named_ones():
    unique = set(kind_tables.tables("unique_name"))
    assert set(kind_tables.tables("named_phenomenon")) == unique & set(kind_tables.tables("phenomenon"))
    assert set(kind_tables.tables("body")) == set(kind_tables.tables("system_filed")) - {"facilities"}


def test_an_unknown_group_is_refused():
    with pytest.raises(ValueError):
        kind_tables.tables("nonsense")


ORDER_GROUPS = {"placed": "phenomenon", "timestamped": "timestamped", "containable": "containable", "unique_name": "unique_name",
                "slow": "phenomenon", "ids": "uid", "phenomenon_kinds": "phenomenon"}


@pytest.mark.parametrize("order, group", sorted(ORDER_GROUPS.items()))
def test_a_named_order_lists_the_tables_of_its_group_once_each(order, group):
    members = set(kind_tables.tables(group))
    named = kind_tables.ORDERS[order]
    assert len(named) == len(set(named))
    assert members <= set(named), (order, members - set(named))


def test_every_named_order_is_one_a_reader_uses():
    assert set(kind_tables.ORDERS) == set(ORDER_GROUPS)


def test_a_reader_keeps_the_order_it_had():
    """The orders of the lists the registry replaced, written out (a progress counter counts the placed
    phenomena in one of them)."""
    assert store.PLACED_PHENOMENON_TABLES[-1] == "quasars" and store.PHENOMENON_UID_TABLES[4] == "quasars"
    assert store.TIMESTAMPED_TABLES[:2] == ("sectors", "star_systems") and store.TIMESTAMPED_TABLES[-1] == "quasars"
    assert list(ids.OBJECT_TABLES)[6] == "facility" and ids.PHENOMENON_KINDS[0] == "nebula"
