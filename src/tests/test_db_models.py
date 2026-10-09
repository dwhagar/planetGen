# tests/test_db_models.py

"""
DB.11: `planetgen/db/models.py` (generated from `schema.sql` by
`scripts/generate_db_models.py`) describes the same tables, columns, types,
nullability, keys and indexes as a database built from `schema.sql`.
"""

import pymysql
import sqlalchemy as sa
from alembic.autogenerate import compare_metadata
from alembic.runtime.migration import MigrationContext

from planetgen.db import models, store


def _differences(config):
    engine = sa.create_engine(
        "mysql+pymysql://",
        creator=lambda: pymysql.connect(host=config.host, port=config.port, user=config.user,
                                         password=config.password, database=config.database),
    )
    try:
        with engine.connect() as connection:
            context = MigrationContext.configure(
                connection,
                opts={
                    "compare_type": True,
                    "include_object": lambda obj, name, kind, reflected, other:
                        not (kind == "table" and name == "alembic_version"),
                },
            )
            return compare_metadata(context, models.metadata)
    finally:
        engine.dispose()


def test_models_match_a_database_built_from_schema_sql(mysql_config):
    store.get_connection(mysql_config).close()
    assert _differences(mysql_config) == []


def test_a_changed_column_is_noticed(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("ALTER TABLE sectors ADD COLUMN stray_column INT")
        conn.commit()
    finally:
        conn.close()
    assert [diff[0] for diff in _differences(mysql_config)] == ["remove_column"]
