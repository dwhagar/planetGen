"""OPS.13, OPS.58: the version-key history an update records."""

import datetime
import hashlib

import pytest

from planetgen.db import store
from planetgen.galaxy import version_history, version_key
from planetgen.generation import generator_battery
from planetgen.util import generator_epoch

SEED = bytes(range(16))


def _control(mysql_config):
    return store.get_control_connection(mysql_config, ensure_schema=True)


def test_each_change_adds_a_row_and_the_last_ten_are_kept(mysql_config):
    conn = _control(mysql_config)
    try:
        for index in range(11):
            version_history.record(conn, "galaxy_a", SEED, lock_sha256=f"{index:064x}")
        version_history.record(conn, "galaxy_b", SEED, lock_sha256="0" * 64)
        rows = version_history.history(conn, "galaxy_a")
        other = version_history.history(conn, "galaxy_b")
    finally:
        conn.execute("DELETE FROM version_key_history")
        conn.commit()
        conn.close()
    assert len(rows) == version_history.KEEP == 10
    # The oldest of the eleven went; the newest is first.
    assert rows[0]["requirements_sha256"] == f"{10:064x}" and rows[-1]["requirements_sha256"] == f"{1:064x}"
    assert len(other) == 1


def test_a_row_holds_the_seed_the_key_and_the_release(mysql_config):
    conn = _control(mysql_config)
    try:
        key = version_history.record(conn, "galaxy_c", SEED, lock_sha256="a" * 64)
        row = version_history.history(conn, "galaxy_c")[0]
    finally:
        conn.execute("DELETE FROM version_key_history")
        conn.commit()
        conn.close()
    assert key == version_key.version_key() and row["version_key"] == key
    assert row["galaxy_seed"] == SEED.hex().upper()


def test_the_lock_hash_is_the_files_sha256(tmp_path):
    lock = tmp_path / "requirements.lock"
    lock.write_bytes(b"pin==1\n")
    assert version_history.requirements_sha256(str(lock)) == hashlib.sha256(b"pin==1\n").hexdigest()
    assert version_history.requirements_sha256(str(tmp_path / "missing")) is None



def test_the_command_only_records_and_has_no_list_option(monkeypatch, capsys):
    """OPS.60: `planetgen versions` lists; `cli.version_history` records."""
    import sys
    import pytest
    from planetgen.cli import version_history as command
    monkeypatch.setattr(sys, "argv", ["version_history", "--list"])
    with pytest.raises(SystemExit) as stopped:
        command.main()
    assert stopped.value.code == 2
    assert "unrecognized arguments: --list" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# OPS.58: the epoch, the battery and the environment, and the row only when something changed
# ---------------------------------------------------------------------------

@pytest.fixture
def control(mysql_config):
    conn = _control(mysql_config)
    try:
        yield conn
    finally:
        conn.execute("DELETE FROM version_key_history")
        conn.commit()
        conn.close()


def _count(conn, database):
    return conn.execute("SELECT COUNT(*) AS n FROM version_key_history WHERE database_name = ?",
                        (database,)).fetchone()["n"]


def test_a_row_holds_the_epoch_the_battery_and_the_environment(control):
    version_history.record(control, "galaxy_d", SEED, lock_sha256="b" * 64)
    row = version_history.history(control, "galaxy_d")[0]
    battery = generator_battery.running_digest()
    assert row["generator_epoch"] == generator_epoch.generator_epoch()
    assert row["battery_digest"] == battery.overall
    assert row["battery_cases"] == generator_battery.cases_of(battery)
    assert row["environment"] == version_history.environment()
    assert set(row["environment"]) >= {"python", "system", "machine", "libc", "numpy", "astropy", "scipy", "scikit-image"}
    assert row["last_seen_at"] == row["recorded_at"]


def test_an_update_that_finds_nothing_changed_only_sets_the_last_seen_time(control):
    version_history.record(control, "galaxy_e", SEED, lock_sha256="c" * 64)
    control.execute("UPDATE version_key_history SET recorded_at = '2026-01-01 00:00:00', last_seen_at = NULL")
    control.commit()
    version_history.record(control, "galaxy_e", SEED, lock_sha256="c" * 64)
    version_history.record(control, "galaxy_e", SEED, lock_sha256="c" * 64)
    rows = version_history.history(control, "galaxy_e")
    assert len(rows) == 1
    assert rows[0]["recorded_at"] == datetime.datetime(2026, 1, 1)
    assert rows[0]["last_seen_at"] > rows[0]["recorded_at"]


@pytest.mark.parametrize("changed", [
    {"lock_sha256": "d" * 64}, {"epoch": 999}, {"key": "0" * 22}, {"env": {"python": "0.0.0"}},
    {"battery": generator_battery.BatteryDigest("f" * 64, {name: generator_battery.CaseDigest("0" * 64, "", "")
                                                          for name in generator_battery.CASES})},
])
def test_a_change_of_one_recorded_value_adds_a_row(control, changed):
    version_history.record(control, "galaxy_f", SEED, lock_sha256="c" * 64)
    version_history.record(control, "galaxy_f", SEED, lock_sha256="c" * 64)
    assert _count(control, "galaxy_f") == 1
    version_history.record(control, "galaxy_f", SEED, **{"lock_sha256": "c" * 64, **changed})
    assert _count(control, "galaxy_f") == 2


def test_a_new_galaxy_seed_adds_a_row(control):
    version_history.record(control, "galaxy_g", SEED, lock_sha256="c" * 64)
    version_history.record(control, "galaxy_g", bytes(reversed(SEED)), lock_sha256="c" * 64)
    assert _count(control, "galaxy_g") == 2


def test_the_changes_between_two_rows_are_named(control):
    older = generator_battery.running_digest()
    cases = generator_battery.cases_of(older)
    newer = generator_battery.BatteryDigest("f" * 64, {
        name: generator_battery.CaseDigest("0" * 64 if name == "nebula-shape" else digest, "", "")
        for name, digest in cases.items()})
    env = version_history.environment()
    version_history.record(control, "galaxy_h", SEED, lock_sha256="1" * 64, epoch=7, env=env, battery=older)
    version_history.record(control, "galaxy_h", SEED, lock_sha256="2" * 64, epoch=8, battery=newer,
                           env=dict(env, numpy="0.0.1"))
    rows = version_history.history(control, "galaxy_h")
    found = version_history.changes(rows[0], rows[1])
    assert found == ["epoch 7 to 8", "lock file", f"numpy {env['numpy']} to 0.0.1",
                     "battery digest (first case nebula-shape)"]
    assert version_history.changes(rows[1], None) == []


def test_a_row_recorded_before_v18_has_no_epoch_and_is_replaced_by_the_next_update(control):
    control.execute(
        "INSERT INTO version_key_history (database_name, galaxy_seed, version_key, planetgen_version,"
        " requirements_sha256, recorded_at) VALUES ('galaxy_i', ?, ?, '8.0.1', ?, '2026-10-01 00:00:00')",
        (SEED.hex().upper(), version_key.version_key(), "c" * 64))
    control.commit()
    row = version_history.history(control, "galaxy_i")[0]
    assert row["generator_epoch"] is None and row["battery_cases"] is None and row["environment"] is None
    assert row["last_seen_at"] == row["recorded_at"]
    version_history.record(control, "galaxy_i", SEED, lock_sha256="c" * 64)
    rows = version_history.history(control, "galaxy_i")
    assert len(rows) == 2
    assert version_history.changes(rows[0], rows[1]) == []  # nothing to say about what was never recorded


def test_an_older_control_database_gains_the_columns(mysql_config):
    """A `version_key_history` as control schema v16 made it is given the five columns and keeps its rows."""
    conn = _control(mysql_config)
    try:
        conn.execute("DROP TABLE version_key_history")
        conn.execute(
            "CREATE TABLE version_key_history (id BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY, database_name VARCHAR(64) NOT NULL,"
            " galaxy_seed CHAR(32) NOT NULL, version_key CHAR(22) NOT NULL, planetgen_version VARCHAR(32) NOT NULL,"
            " requirements_sha256 CHAR(64) NULL, recorded_at DATETIME(6) NOT NULL,"
            " INDEX idx_version_key_history_database (database_name, id)) ENGINE=InnoDB")
        conn.execute("INSERT INTO version_key_history (database_name, galaxy_seed, version_key, planetgen_version, recorded_at)"
                     " VALUES ('old', ?, ?, '8.0.1', '2026-10-01 00:00:00')", (SEED.hex().upper(), "0" * 22))
        conn.commit()
        store._ensure_control_schema(conn)
        store._ensure_control_schema(conn)  # a second deploy changes nothing
        columns = {row["c"] for row in conn.execute(
            "SELECT COLUMN_NAME AS c FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
            " AND TABLE_NAME = 'version_key_history'").fetchall()}
        row = version_history.history(conn, "old")[0]
    finally:
        conn.execute("DELETE FROM version_key_history")
        conn.commit()
        conn.close()
    assert {"generator_epoch", "battery_digest", "battery_cases", "environment", "last_seen_at"} <= columns
    assert row["version_key"] == "0" * 22 and row["generator_epoch"] is None
    assert row["last_seen_at"] == datetime.datetime(2026, 10, 1)
