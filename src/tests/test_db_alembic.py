# tests/test_db_alembic.py

"""
DB.11: Alembic carries the schema forward from the v61 baseline. Revision
ids are schema versions, `migrate_database` runs the revisions after the
legacy steps and mirrors each into `schema_migrations`, and the baseline
leaves an existing database alone.
"""

import pytest

from planetgen.db import alembic_runner, store

HEAD = store.SCHEMA_VERSION
NEXT = HEAD + 1
HEAD_ID = alembic_runner.revision_id(HEAD)
NEXT_ID = alembic_runner.revision_id(NEXT)

REVISION_NEXT = f'''
revision = "{NEXT_ID}"
down_revision = "{HEAD_ID}"
branch_labels = None
depends_on = None


def upgrade():
    from alembic import op
    op.execute("CREATE TABLE alembic_probe (id INT PRIMARY KEY)")
'''

REVISION_BAD_ID = f'''
revision = "abc"
down_revision = "{HEAD_ID}"
branch_labels = None
depends_on = None


def upgrade():
    pass
'''


@pytest.fixture
def scripts(tmp_path):
    """A copy of the real migrations with one more revision (the next schema version) added."""
    import shutil
    target = tmp_path / "migrations"
    shutil.copytree(alembic_runner.MIGRATIONS_DIR, target, ignore=shutil.ignore_patterns("__pycache__"))
    (target / "versions" / f"{NEXT_ID}_probe.py").write_text(REVISION_NEXT)
    return str(target)


def test_schema_version_is_the_head_revision():
    versions = alembic_runner.revisions()
    assert versions[0] == (61, "0061")  # the baseline
    assert store.SCHEMA_VERSION == versions[-1][0]


def test_revision_ids_must_be_consecutive_schema_versions(tmp_path):
    import shutil
    target = tmp_path / "migrations"
    shutil.copytree(alembic_runner.MIGRATIONS_DIR, target, ignore=shutil.ignore_patterns("__pycache__"))
    (target / "versions" / "bad.py").write_text(REVISION_BAD_ID)
    with pytest.raises(RuntimeError, match=f"should be '{NEXT_ID}'"):
        alembic_runner.revisions(str(target))


def test_pending_lists_only_later_revisions(scripts):
    assert alembic_runner.pending(HEAD, scripts) == [(NEXT, NEXT_ID)]
    assert alembic_runner.pending(NEXT, scripts) == []


def test_upgrade_stamps_the_baseline_then_runs_later_revisions(mysql_config, scripts):
    store.get_connection(mysql_config).close()
    recorded = []
    steps = []
    version = alembic_runner.upgrade(mysql_config, HEAD, recorded.append,
                                     on_step=lambda *args: steps.append(args), scripts_dir=scripts)
    assert version == NEXT
    assert recorded == [NEXT]
    assert steps == [(1, 1, HEAD, NEXT)]
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT version_num FROM alembic_version").fetchone()["version_num"] == NEXT_ID
        assert conn.execute("SELECT COUNT(*) AS n FROM alembic_probe").fetchone()["n"] == 0
    finally:
        conn.close()
    # Running again finds nothing to do.
    assert alembic_runner.upgrade(mysql_config, NEXT, recorded.append, scripts_dir=scripts) == NEXT
    assert recorded == [NEXT]


def test_the_baseline_alone_only_stamps(mysql_config):
    store.get_connection(mysql_config).close()
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        stamped = conn.execute("SELECT version_num FROM alembic_version").fetchone()["version_num"]
    finally:
        conn.close()
    assert stamped == alembic_runner.revision_id(store.SCHEMA_VERSION)


def test_migrate_database_mirrors_revisions_into_schema_migrations(mysql_config, scripts, monkeypatch):
    store.get_connection(mysql_config).close()
    monkeypatch.setattr(alembic_runner, "MIGRATIONS_DIR", scripts)
    monkeypatch.setattr(store, "SCHEMA_VERSION", NEXT)
    steps = []
    assert store.migrate_database(mysql_config, on_step=lambda *args: steps.append(args)) == NEXT
    assert steps == [(1, 1, HEAD, NEXT)]
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT MAX(version) AS v FROM schema_migrations").fetchone()["v"] == NEXT
    finally:
        conn.close()
    assert store.schema_status(mysql_config) == (NEXT, 0)
