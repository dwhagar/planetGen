# planetgen/db/migrations/env.py

"""
Alembic environment for planetGen's content database. Always run online,
on the connection `planetgen.db.alembic_runner` hands in through
`config.attributes["connection"]` -- there is no alembic.ini, no offline
(SQL script) mode, and no database URL to configure.
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
