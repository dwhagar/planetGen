# tests/test_migration_designations.py

"""
GEN.1138: migration 0187 creates `designations`, empty. A galaxy made at v186 keeps every row it had and gains a table that takes a
retired ID once and a retired name any number of times, and refuses a row that is neither.
"""

import pytest

from planetgen.db import store
from tests import db_schema_support
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database

_UID = bytes(range(12))


def test_a_v186_galaxy_gets_an_empty_designations_table(mysql_config):
    number = db_schema_support.migration_number("designations")
    load_old_schema(mysql_config, number - 1)
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        conn.execute("INSERT INTO sectors (name, edge_mpc) VALUES ('Before', 4000)")
        conn.commit()
        assert conn.execute("SHOW TABLES LIKE 'designations'").fetchone() is None
    finally:
        conn.close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM designations").fetchone()["n"] == 0
        assert conn.execute("SELECT COUNT(*) AS n FROM sectors").fetchone()["n"] == 1
    finally:
        conn.close()


def _insert(conn, kind, uid=None, name=None, reason="rename"):
    conn.execute("INSERT INTO designations (designation_kind, uid, name, object_kind, object_row_id, retired_epoch, reason)"
                 " VALUES (?, ?, ?, 'moons', 1, 3, ?)", (kind, uid, name, reason))


def test_an_id_is_retired_once_and_a_name_as_often_as_it_changes(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        _insert(conn, "uid", uid=_UID)
        with pytest.raises(Exception, match="uq_designations_uid"):
            _insert(conn, "uid", uid=_UID)
        conn.rollback()
        _insert(conn, "name", name="Vestara II")
        _insert(conn, "name", name="Vestara II", reason="collision bump")
        with pytest.raises(Exception, match="chk_designations_kind"):
            _insert(conn, "uid", name="Vestara II")
        with pytest.raises(Exception, match="chk_designations_reason"):
            _insert(conn, "name", name="Vestara III", reason="whim")
    finally:
        conn.rollback()
        conn.close()
