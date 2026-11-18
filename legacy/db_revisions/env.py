# legacy/db_revisions/env.py

"""
Alembic environment for the legacy revisions v62 to v130 (DB.126): the steps
that took a database from the last release before 8.0 to the structure 8.0
and above accept. Copied unchanged from `planetgen/db/migrations/` as it stood
before the re-baseline (DB.98). Always run online, on the connection
`legacy/bump_database_to_8_0.py` hands in through `config.attributes["connection"]`;
there is no alembic.ini, no offline (SQL script) mode, and no database URL.
"""

from alembic import context

config = context.config


def run_migrations_online():
    connection = config.attributes.get("connection")
    if connection is None:
        raise RuntimeError("planetGen migrations run through planetgen.db.alembic_runner, "
                           "which supplies the database connection")
    context.configure(connection=connection, target_metadata=None)
    with context.begin_transaction():
        context.run_migrations()


if context.is_offline_mode():
    raise RuntimeError("planetGen migrations have no offline mode")
run_migrations_online()
