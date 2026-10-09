# planetgen/db/alembic_runner.py

"""
Runs the Alembic revisions in `planetgen/db/migrations/` against a
planetGen content database (DB.11).

Each revision's id is its schema version as four digits (`0062`), so the
head revision's number is `store.SCHEMA_VERSION`. The legacy
`_migrate_vN_to_vM` steps in `store.py` take a database up to the baseline
(v61); everything after that is a revision here. `store.migrate_database`
calls `upgrade` once those have run, and `record` writes the
`schema_migrations` row for each version reached.
"""

import os

import pymysql
import sqlalchemy
from alembic import command
from alembic.config import Config
from alembic.runtime.migration import MigrationContext
from alembic.script import ScriptDirectory

MIGRATIONS_DIR = os.path.join(os.path.dirname(__file__), "migrations")
"""str: The directory holding `env.py` and `versions/`."""

BASELINE_VERSION = 61
"""int: The first revision; databases below it are upgraded by the legacy
steps in `store.py` before Alembic sees them."""


def revision_id(version):
    """The revision id for a schema version: `62` -> `"0062"`."""
    return f"{version:04d}"


def _config(scripts_dir=None, connection=None):
    cfg = Config()
    cfg.set_main_option("script_location", scripts_dir or MIGRATIONS_DIR)
    if connection is not None:
        cfg.attributes["connection"] = connection
    return cfg


def revisions(scripts_dir=None):
    """Every revision, oldest first, as `(schema version, revision id)`.

    Raises:
        RuntimeError: The revisions are not one line of consecutively
            numbered ids starting at the baseline.
    """
    script = ScriptDirectory.from_config(_config(scripts_dir))
    ordered = list(reversed(list(script.walk_revisions())))
    found = []
    for position, rev in enumerate(ordered):
        version = BASELINE_VERSION + position
        if rev.revision != revision_id(version):
            raise RuntimeError(f"migration revision {rev.revision!r} should be {revision_id(version)!r}: "
                               f"ids are the schema version, in order, from {revision_id(BASELINE_VERSION)}")
        found.append((version, rev.revision))
    return found


def head_version(scripts_dir=None):
    """The schema version of the newest revision."""
    return revisions(scripts_dir)[-1][0]


def pending(version, scripts_dir=None):
    """The `(schema version, revision id)` pairs still to run on a
    database at `version` once the legacy steps have taken it to the
    baseline."""
    return [(v, r) for v, r in revisions(scripts_dir) if v > max(version, BASELINE_VERSION)]


def _engine(config):
    return sqlalchemy.create_engine(
        "mysql+pymysql://",
        creator=lambda: pymysql.connect(
            host=config.host, port=config.port, user=config.user, password=config.password,
            database=config.database, charset="utf8mb4", autocommit=False,
            init_command="SET time_zone = '+00:00'",
        ),
        poolclass=sqlalchemy.pool.NullPool,
    )


def upgrade(config, from_version, record, on_step=None, step_offset=0, total=None, scripts_dir=None):
    """
    Records the database as being at `from_version` in Alembic's own
    `alembic_version` table (a database that has only run the legacy steps
    has none yet) and runs every later revision, calling
    `record(version)` after each.

    Args:
        config (MySQLConfig): The database.
        from_version (int): The schema version it is at (at least the
            baseline).
        record (callable): Called as `record(version)` after each revision.
        on_step (callable, optional): Called as `on_step(number, total,
            from_version, to_version)` before each revision, `number`
            counting from `step_offset + 1`.
        total (int, optional): The step count `on_step` reports; defaults
            to the revisions pending here.

    Returns:
        int: The schema version the database is at now.
    """
    todo = pending(from_version, scripts_dir)
    if total is None:
        total = step_offset + len(todo)
    engine = _engine(config)
    try:
        with engine.connect() as connection:
            cfg = _config(scripts_dir, connection)
            current = MigrationContext.configure(connection).get_current_revision()
            if current != revision_id(from_version):
                command.stamp(cfg, revision_id(from_version), purge=current is not None)
                connection.commit()
            version = from_version
            for number, (target, rev) in enumerate(todo, start=step_offset + 1):
                if on_step is not None:
                    on_step(number, total, version, target)
                command.upgrade(cfg, rev)
                connection.commit()
                record(target)
                version = target
    finally:
        engine.dispose()
    return version
