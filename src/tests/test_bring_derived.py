"""OPS.272: the derived stage of `bring_to_epoch`: the stored values that are a pure function of stored ones (the habitability scores and the effect ledger) are refilled by the jobs of their commands
for the systems whose stamp is not the running code's, no sector is generated again, and no sector's `modified_at` moves (so the Galaxy Map's tile cache is kept, MAP.431)."""

import pytest

from planetgen.db import derived_stamp, store
from planetgen.galaxy import bring_to_epoch, epoch_record
from tests import inline_queue
from tests.test_migration_stage_epochs import _system

pytestmark = pytest.mark.fresh_database


@pytest.fixture
def saved(mysql_config):
    """One system with planets saved in a sector, stamps current, and the sector's `modified_at` set a day back so a touch shows."""
    system, cfg = _system()
    conn = store.get_connection(mysql_config)
    with conn:
        sector_id = conn.execute("INSERT INTO sectors (name, edge_mpc, generator_epoch) VALUES ('Derived sector', 1000, 1)").lastrowid
        system_id = store.insert_star_system(conn, system, cfg)
        conn.execute("UPDATE star_systems SET sector_id = ? WHERE id = ?", (sector_id, system_id))
        store.refresh_habitability(conn, system_id)
        store.refresh_effect_ledger(conn, system_id)
        conn.execute("UPDATE sectors SET modified_at = TIMESTAMP('2026-01-01 00:00:00') WHERE id = ?", (sector_id,))
    yield conn, sector_id, system_id
    conn.close()


def _stamps(conn, system_id):
    return conn.execute("SELECT habitability_stamp, ledger_stamp FROM system_summary WHERE star_system_id = ?", (system_id,)).fetchone()


def test_current_stamps_leave_nothing_to_refill(saved):
    conn, _sector_id, _system_id = saved
    assert epoch_record.derived_stale(conn) == {}
    assert bring_to_epoch.survey(conn)["derived"] == {}


def test_an_old_stamp_is_found_reported_and_refilled_by_the_run_without_touching_the_sector(saved, mysql_config):
    conn, sector_id, system_id = saved
    with conn:
        conn.execute("UPDATE system_summary SET habitability_stamp = 1, ledger_stamp = 1 WHERE star_system_id = ?", (system_id,))
    found = bring_to_epoch.survey(conn)
    assert found["derived"] == {"habitability": 1, "ledger": 1} and not bring_to_epoch._nothing_to_do(found)
    lines = bring_to_epoch.report_lines(found)
    assert any("the habitability of 1 system(s) was written by older code" in line for line in lines)
    said = []
    with inline_queue.inline_queue():
        result = bring_to_epoch.run(mysql_config, say=said.append, wide_binaries=False)
    assert result["failed"] == [] and result["redone"] == 0 and result["derived"] == {"habitability": 1, "ledger": 1}
    assert any("Derived values: refilling the habitability of 1 system(s)." in line for line in said)
    conn.commit()           # the run wrote through connections of its own
    row = _stamps(conn, system_id)
    assert row["habitability_stamp"] == derived_stamp.current("habitability") and row["ledger_stamp"] == derived_stamp.current("ledger")
    assert str(conn.execute("SELECT modified_at FROM sectors WHERE id = ?", (sector_id,)).fetchone()["modified_at"]).startswith("2026-01-01")
    assert epoch_record.derived_stale(conn) == {}


def test_a_dry_run_names_the_refill_and_changes_nothing(saved, mysql_config):
    conn, _sector_id, system_id = saved
    with conn:
        conn.execute("UPDATE system_summary SET ledger_stamp = 1 WHERE star_system_id = ?", (system_id,))
    said = []
    result = bring_to_epoch.run(mysql_config, dry_run=True, say=said.append, wide_binaries=False)
    assert "derived" not in result and _stamps(conn, system_id)["ledger_stamp"] == 1
    assert any("the ledger of 1 system(s)" in line for line in said)
