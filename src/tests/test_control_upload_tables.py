# tests/test_control_upload_tables.py

"""
DB.62 (control schema v16): the bookkeeping tables for the uploads of a remote run. They are made empty, a unit has a known
state, a sector address is claimed by one run, the id ranges are ordered, and everything goes with its run, which goes with its
owner.
"""

import pytest

from planetgen.db import store

TABLES = ("upload_runs", "upload_batches", "upload_units", "upload_claims", "upload_id_ranges")


@pytest.fixture
def control(mysql_config):
    conn = store.get_control_connection(mysql_config, ensure_schema=True)
    try:
        yield conn
    finally:
        conn.close()


def _count(conn, table):
    return conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]


def _run(conn, owner, key=None, galaxy="galaxy_a"):
    return conn.execute("INSERT INTO upload_runs (owner_id, key_id, galaxy_database) VALUES (?, ?, ?)",
                        (owner, key, galaxy)).lastrowid


@pytest.fixture
def owner(control):
    return control.execute("INSERT INTO admin_users (username, password_hash) VALUES ('uploader', 'x')").lastrowid


def test_the_tables_are_made_empty_and_the_version_is_the_current_one(control):
    for table in TABLES:
        assert _count(control, table) == 0
    assert control.execute("SELECT MAX(version) AS v FROM control_schema_migrations").fetchone()["v"] == store.CONTROL_SCHEMA_VERSION


def test_a_control_database_at_v15_gains_the_tables_and_keeps_its_rows(control):
    for table in reversed(TABLES):
        control.execute(f"DROP TABLE {table}")
    control.execute("DELETE FROM control_schema_migrations WHERE version = 16")
    control.execute("INSERT INTO admin_users (username, password_hash) VALUES ('kept', 'x')")
    control.commit()
    store._ensure_control_schema(control)
    store._ensure_control_schema(control)
    for table in TABLES:
        assert _count(control, table) == 0
    assert control.execute("SELECT username FROM admin_users WHERE username = 'kept'").fetchone() is not None


def test_a_run_starts_open_and_a_unit_has_a_known_state(control, owner):
    run = _run(control, owner)
    row = control.execute("SELECT state, key_id FROM upload_runs WHERE id = ?", (run,)).fetchone()
    assert row == {"state": "open", "key_id": None}
    for state in ("open", "sealed", "verified", "finalized", "rejected"):
        control.execute("INSERT INTO upload_units (run_id, unit_key, state) VALUES (?, ?, ?)", (run, f"u-{state}", state))
    with pytest.raises(Exception, match="chk_upload_units_state"):
        control.execute("INSERT INTO upload_units (run_id, unit_key, state) VALUES (?, 'u-bad', 'lost')", (run,))
    with pytest.raises(Exception, match="PRIMARY"):
        control.execute("INSERT INTO upload_units (run_id, unit_key) VALUES (?, 'u-open')", (run,))


def test_a_batch_number_and_a_sector_address_are_held_once(control, owner):
    first, second = _run(control, owner), _run(control, owner)
    control.execute("INSERT INTO upload_batches (run_id, seq, sha256, bytes) VALUES (?, 1, ?, 100)", (first, "a" * 64))
    control.execute("INSERT INTO upload_batches (run_id, seq, sha256, bytes) VALUES (?, 1, ?, 100)", (second, "a" * 64))
    with pytest.raises(Exception, match="PRIMARY"):
        control.execute("INSERT INTO upload_batches (run_id, seq, sha256, bytes) VALUES (?, 1, ?, 5)", (first, "b" * 64))
    claim = "INSERT INTO upload_claims (galaxy_database, shell_index, band_index, wedge_index, run_id) VALUES (?, 3, -1, 7, ?)"
    control.execute(claim, ("galaxy_a", first))
    control.execute(claim, ("galaxy_b", second))
    with pytest.raises(Exception, match="PRIMARY"):
        control.execute(claim, ("galaxy_a", second))
    assert control.execute("SELECT state FROM upload_claims WHERE galaxy_database = 'galaxy_a'").fetchone() == {"state": "claimed"}


def test_an_id_range_runs_upward(control, owner):
    run = _run(control, owner)
    control.execute("INSERT INTO upload_id_ranges (run_id, table_name, first_id, last_id) VALUES (?, 'stars', 10, 10)", (run,))
    with pytest.raises(Exception, match="chk_upload_id_ranges_order"):
        control.execute("INSERT INTO upload_id_ranges (run_id, table_name, first_id, last_id) VALUES (?, 'stars', 30, 20)", (run,))


def test_everything_goes_with_its_run_and_a_removed_key_only_clears_the_link(control, owner):
    key = control.execute("INSERT INTO admin_api_keys (admin_user_id, label, key_hash) VALUES (?, 'k', ?)",
                          (owner, "c" * 64)).lastrowid
    run = _run(control, owner, key)
    control.execute("INSERT INTO upload_batches (run_id, seq, sha256, bytes) VALUES (?, 1, ?, 1)", (run, "a" * 64))
    control.execute("INSERT INTO upload_units (run_id, unit_key) VALUES (?, 'u')", (run,))
    control.execute("INSERT INTO upload_claims (galaxy_database, shell_index, band_index, wedge_index, run_id)"
                    " VALUES ('galaxy_a', 1, 0, 0, ?)", (run,))
    control.execute("INSERT INTO upload_id_ranges (run_id, table_name, first_id, last_id) VALUES (?, 'stars', 1, 9)", (run,))
    control.execute("DELETE FROM admin_api_keys WHERE id = ?", (key,))
    assert control.execute("SELECT key_id FROM upload_runs WHERE id = ?", (run,)).fetchone() == {"key_id": None}
    control.execute("DELETE FROM admin_users WHERE id = ?", (owner,))
    for table in TABLES:
        assert _count(control, table) == 0, table
