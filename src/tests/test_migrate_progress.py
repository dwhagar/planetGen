"""
`migrate_database`'s per-step callback (planetgen.cli.migrate's progress bar) and
`schema_status` (planetgen.cli.migrate --status, which update.sh reads before
asking whether to migrate or delete the data). Needs a MySQL test server
like every other database-backed test (see `conftest.py`).
"""
from planetgen.db import store


def _roll_back_to_v28(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        # What v29 dropped, so its step has something to do again.
        conn.execute("ALTER TABLE star_systems ADD COLUMN wikitext_content LONGTEXT, "
                     "ADD COLUMN markdown_content LONGTEXT")
        conn.execute("DELETE FROM schema_migrations WHERE version > 28")
        conn.execute("INSERT IGNORE INTO schema_migrations (version) VALUES (28)")
        conn.commit()
    finally:
        conn.close()


def test_status_and_steps_report_each_pending_migration(mysql_config):
    assert store.schema_status(mysql_config) == (store.SCHEMA_VERSION, 0)

    _roll_back_to_v28(mysql_config)
    pending = store.SCHEMA_VERSION - 28
    assert store.schema_status(mysql_config) == (28, pending)

    calls = []
    assert store.migrate_database(mysql_config, on_step=lambda *a: calls.append(a)) == store.SCHEMA_VERSION
    assert calls == [(n, pending, 27 + n, 28 + n) for n in range(1, pending + 1)]

    # Already current: no steps, nothing reported.
    calls.clear()
    assert store.migrate_database(mysql_config, on_step=lambda *a: calls.append(a)) == store.SCHEMA_VERSION
    assert calls == []
    assert store.schema_status(mysql_config) == (store.SCHEMA_VERSION, 0)


def test_every_step_is_listed_once_in_order():
    targets = [target for target, _ in store._migration_steps()]
    assert targets == list(range(9, store.SCHEMA_VERSION + 1))
