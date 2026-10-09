"""
`migrate_database`'s per-step callback (planetgen.cli.migrate's progress bar) and
`schema_status` (planetgen.cli.migrate --status, which update.sh reads before
asking whether to migrate or delete the data). Needs a MySQL test server
like every other database-backed test (see `conftest.py`).
"""
from planetgen.db import alembic_runner, store
from tests.db_schema_support import migrations_with_probes


def test_status_and_steps_report_each_pending_migration(mysql_config, tmp_path, monkeypatch):
    assert store.schema_status(mysql_config) == (store.SCHEMA_VERSION, 0)
    store.get_connection(mysql_config).close()

    monkeypatch.setattr(alembic_runner, "MIGRATIONS_DIR", migrations_with_probes(tmp_path, up_to=63))
    monkeypatch.setattr(store, "SCHEMA_VERSION", 63)
    assert store.schema_status(mysql_config) == (61, 2)

    calls = []
    assert store.migrate_database(mysql_config, on_step=lambda *a: calls.append(a)) == 63
    assert calls == [(1, 2, 61, 62), (2, 2, 62, 63)]

    # Already current: no steps, nothing reported.
    calls.clear()
    assert store.migrate_database(mysql_config, on_step=lambda *a: calls.append(a)) == 63
    assert calls == []
    assert store.schema_status(mysql_config) == (63, 0)
