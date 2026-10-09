"""OPS.13: the version-key history an update records."""

import hashlib

from planetgen.db import store
from planetgen.galaxy import version_history, version_key

SEED = bytes(range(16))


def _control(mysql_config):
    return store.get_control_connection(mysql_config, ensure_schema=True)


def test_each_update_adds_a_row_and_keeps_the_last_ten(mysql_config):
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

