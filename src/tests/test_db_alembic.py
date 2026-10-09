# tests/test_db_alembic.py

"""
DB.11: Alembic carries the schema forward from the v61 baseline. Revision
ids are schema versions, `migrate_database` runs the revisions after the
legacy steps and mirrors each into `schema_migrations`, and the baseline
leaves an existing database alone.
"""

import pytest

from planetgen.db import alembic_runner, store

REVISION_62 = '''
revision = "0062"
down_revision = "0061"
branch_labels = None
depends_on = None


def upgrade():
    from alembic import op
    op.execute("CREATE TABLE alembic_probe (id INT PRIMARY KEY)")
'''

REVISION_BAD_ID = '''
revision = "abc"
down_revision = "0061"
branch_labels = None
depends_on = None


def upgrade():
    pass
'''


@pytest.fixture
def scripts(tmp_path):
    """A copy of the real migrations with `0062` added."""
    import shutil
    target = tmp_path / "migrations"
    shutil.copytree(alembic_runner.MIGRATIONS_DIR, target, ignore=shutil.ignore_patterns("__pycache__"))
    (target / "versions" / "0062_probe.py").write_text(REVISION_62)
    return str(target)


def test_schema_version_is_the_head_revision():
    versions = alembic_runner.revisions()
    assert versions[0] == (61, "0061")
    assert store.SCHEMA_VERSION == versions[-1][0]


def test_revision_ids_must_be_consecutive_schema_versions(tmp_path):
    import shutil
    target = tmp_path / "migrations"
    shutil.copytree(alembic_runner.MIGRATIONS_DIR, target, ignore=shutil.ignore_patterns("__pycache__"))
    (target / "versions" / "bad.py").write_text(REVISION_BAD_ID)
    with pytest.raises(RuntimeError, match="should be '0062'"):
        alembic_runner.revisions(str(target))


def test_pending_lists_only_later_revisions(scripts):
    assert alembic_runner.pending(61, scripts) == [(62, "0062")]
    assert alembic_runner.pending(62, scripts) == []


def test_upgrade_stamps_the_baseline_then_runs_later_revisions(mysql_config, scripts):
    store.get_connection(mysql_config).close()
    recorded = []
    steps = []
    version = alembic_runner.upgrade(mysql_config, 61, recorded.append,
                                     on_step=lambda *args: steps.append(args), scripts_dir=scripts)
    assert version == 62
    assert recorded == [62]
    assert steps == [(1, 1, 61, 62)]
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT version_num FROM alembic_version").fetchone()["version_num"] == "0062"
        assert conn.execute("SELECT COUNT(*) AS n FROM alembic_probe").fetchone()["n"] == 0
    finally:
        conn.close()
    # Running again finds nothing to do.
    assert alembic_runner.upgrade(mysql_config, 62, recorded.append, scripts_dir=scripts) == 62
    assert recorded == [62]


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
    monkeypatch.setattr(store, "SCHEMA_VERSION", 62)
    steps = []
    assert store.migrate_database(mysql_config, on_step=lambda *args: steps.append(args)) == 62
    assert steps == [(1, 1, 61, 62)]
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT MAX(version) AS v FROM schema_migrations").fetchone()["v"] == 62
    finally:
        conn.close()
    assert store.schema_status(mysql_config) == (62, 0)
