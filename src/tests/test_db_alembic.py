# tests/test_db_alembic.py

"""
DB.11: Alembic carries the schema forward from the v61 baseline. Revision
ids are schema versions, `migrate_database` runs the revisions after the
legacy steps and mirrors each into `schema_migrations`, and the baseline
leaves an existing database alone.
"""

import pytest

from planetgen.db import alembic_runner, store

HEAD = alembic_runner.head_version()
PROBE = HEAD + 1

REVISION_PROBE = f'''
revision = "{PROBE:04d}"
down_revision = "{HEAD:04d}"
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
    """A copy of the real migrations with one probe revision added."""
    import shutil
    target = tmp_path / "migrations"
    shutil.copytree(alembic_runner.MIGRATIONS_DIR, target, ignore=shutil.ignore_patterns("__pycache__"))
    (target / "versions" / f"{PROBE:04d}_probe.py").write_text(REVISION_PROBE)
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
    with pytest.raises(RuntimeError, match="should be '{:04d}'".format(HEAD)):
        alembic_runner.revisions(str(target))


def test_pending_lists_only_later_revisions(scripts):
    assert alembic_runner.pending(HEAD, scripts) == [(PROBE, f"{PROBE:04d}")]
    assert alembic_runner.pending(PROBE, scripts) == []


def test_upgrade_stamps_the_baseline_then_runs_later_revisions(mysql_config, scripts):
    store.get_connection(mysql_config).close()
    recorded = []
    steps = []
    version = alembic_runner.upgrade(mysql_config, HEAD, recorded.append,
                                     on_step=lambda *args: steps.append(args), scripts_dir=scripts)
    assert version == PROBE
    assert recorded == [PROBE]
    assert steps == [(1, 1, HEAD, PROBE)]
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT version_num FROM alembic_version").fetchone()["version_num"] == f"{PROBE:04d}"
        assert conn.execute("SELECT COUNT(*) AS n FROM alembic_probe").fetchone()["n"] == 0
    finally:
        conn.close()
    # Running again finds nothing to do.
    assert alembic_runner.upgrade(mysql_config, PROBE, recorded.append, scripts_dir=scripts) == PROBE
    assert recorded == [PROBE]


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
    monkeypatch.setattr(store, "SCHEMA_VERSION", PROBE)
    steps = []
    assert store.migrate_database(mysql_config, on_step=lambda *args: steps.append(args)) == PROBE
    assert steps == [(1, 1, HEAD, PROBE)]
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        assert conn.execute("SELECT MAX(version) AS v FROM schema_migrations").fetchone()["v"] == PROBE
    finally:
        conn.close()
    assert store.schema_status(mysql_config) == (PROBE, 0)
