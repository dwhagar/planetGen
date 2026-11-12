"""The code that generated each sector (DB.7).

Schema v67. `sectors` gains `version_key`, `planetgen_version`,
`python_version` and `platform` (the values `galaxy_shape` keeps for the
galaxy's plan) and an index on `version_key`. Existing sectors keep NULL:
they were generated before the release was recorded.
"""

import sqlalchemy as sa
from alembic import op

revision = "0067"
down_revision = "0066"
branch_labels = None
depends_on = None

_COLUMNS = (
    ("version_key", "CHAR(22)"),
    ("planetgen_version", "VARCHAR(32)"),
    ("python_version", "VARCHAR(32)"),
    ("platform", "VARCHAR(64)"),
)


def _count(connection, query, **params):
    return connection.execute(sa.text(query), params).scalar()


def upgrade():
    connection = op.get_bind()
    for name, definition in _COLUMNS:
        if not _count(connection,
                      "SELECT COUNT(*) FROM information_schema.columns WHERE table_schema = DATABASE()"
                      " AND table_name = 'sectors' AND column_name = :name", name=name):
            connection.execute(sa.text(f"ALTER TABLE sectors ADD COLUMN {name} {definition}"))
    if not _count(connection,
                  "SELECT COUNT(*) FROM information_schema.statistics WHERE table_schema = DATABASE()"
                  " AND table_name = 'sectors' AND index_name = 'idx_sectors_version_key'"):
        connection.execute(sa.text("ALTER TABLE sectors ADD KEY idx_sectors_version_key (version_key)"))
