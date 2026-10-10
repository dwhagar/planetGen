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

    head = store.SCHEMA_VERSION
    monkeypatch.setattr(alembic_runner, "MIGRATIONS_DIR", migrations_with_probes(tmp_path, count=2))
    monkeypatch.setattr(store, "SCHEMA_VERSION", head + 2)
    assert store.schema_status(mysql_config) == (head, 2)

    calls = []
    assert store.migrate_database(mysql_config, on_step=lambda *a: calls.append(a)) == head + 2
    assert calls == [(1, 2, head, head + 1), (2, 2, head + 1, head + 2)]

    # Already current: no steps, nothing reported.
    calls.clear()
    assert store.migrate_database(mysql_config, on_step=lambda *a: calls.append(a)) == head + 2
    assert calls == []
    assert store.schema_status(mysql_config) == (head + 2, 0)


def _batch_probe(tmp_path, head):
    """A probe revision after the real head that reports 4 batches of work."""
    directory = migrations_with_probes(tmp_path, count=1)
    probe = next(path for path in (tmp_path / "migrations" / "versions").glob(f"{head + 1:04d}_probe.py"))
    probe.write_text(probe.read_text().replace(
        "op.execute('CREATE TABLE alembic_probe (id INT PRIMARY KEY)')",
        "from planetgen.db import alembic_runner\n"
        "    op.execute('CREATE TABLE alembic_probe (id INT PRIMARY KEY)')\n"
        "    for done in range(1, 5):\n"
        "        alembic_runner.report_progress('rows', done, 4)"))
    return directory


def test_a_revision_in_batches_reports_its_progress(mysql_config, tmp_path, monkeypatch):
    """DB.15: a long revision calls `report_progress`; the caller's `on_progress` hears each batch."""
    store.get_connection(mysql_config).close()
    head = store.SCHEMA_VERSION
    monkeypatch.setattr(alembic_runner, "MIGRATIONS_DIR", _batch_probe(tmp_path, head))
    monkeypatch.setattr(store, "SCHEMA_VERSION", head + 1)
    heard = []
    store.migrate_database(mysql_config, on_progress=lambda *a: heard.append(a))
    assert heard == [("rows", 1, 4), ("rows", 2, 4), ("rows", 3, 4), ("rows", 4, 4)]


def test_a_revision_in_batches_runs_without_a_listener(mysql_config, tmp_path, monkeypatch):
    store.get_connection(mysql_config).close()
    head = store.SCHEMA_VERSION
    monkeypatch.setattr(alembic_runner, "MIGRATIONS_DIR", _batch_probe(tmp_path, head))
    monkeypatch.setattr(store, "SCHEMA_VERSION", head + 1)
    assert store.migrate_database(mysql_config) == head + 1


def test_the_plain_line_reporter_speaks_every_thirty_seconds_and_each_tenth():
    from planetgen.cli import migrate

    lines, now = [], [0.0]
    reporter = migrate._LineReporter(write=lines.append, clock=lambda: now[0])
    reporter.on_step(1, 2, 70, 71)
    assert lines == ["Migrating v70 -> v71 (step 1 of 2)"]
    reporter.on_progress("rows", 1, 100)          # first: spoken
    reporter.on_progress("rows", 5, 100)          # same tenth, under 30 s: quiet
    now[0] = 31.0
    reporter.on_progress("rows", 6, 100)          # 30 s later: spoken
    reporter.on_progress("rows", 15, 100)         # a new tenth: spoken
    reporter.on_progress("rows", 100, 100)        # done: spoken once
    reporter.on_progress("rows", 100, 100)
    assert [line for line in lines if line.startswith("  ")] == [
        "  rows: 1 of 100 (1%)", "  rows: 6 of 100 (6%)", "  rows: 15 of 100 (15%)", "  rows: 100 of 100 (100%)"]


def test_on_a_terminal_the_cli_draws_the_steps_and_the_batches(monkeypatch):
    from planetgen.cli import migrate

    calls = []

    class FakeStage:
        def __init__(self, total, **kwargs):
            calls.append(("open", total, kwargs["kind"]))

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            calls.append(("close", exc[0]))

        def stage(self, description):
            calls.append(("stage", description))

        def detail(self, label, done, total):
            calls.append(("detail", label, done, total))

    def fake_migrate(config, on_step=None, on_progress=None):
        on_step(1, 1, 70, 71)
        on_progress("rows", 2, 4)
        return 71

    monkeypatch.setattr(migrate.stage_progress, "StageProgress", FakeStage)
    monkeypatch.setattr(migrate, "migrate_database", fake_migrate)
    monkeypatch.setattr(migrate.sys.stdout, "isatty", lambda: True, raising=False)
    assert migrate._migrate_with_progress(None) == 71
    assert calls == [("open", 1, "migrate"), ("stage", "Migrating v70 -> v71 (step 1 of 1)"),
                     ("detail", "rows", 2, 4), ("close", None)]
