"""A stored [Fe/H] on stars and pre-placed bright stars (GEN.158).

Schema v95. `stars` and `bright_stars` gain `metallicity_fe_h` (DOUBLE, nullable), each star's metallicity in dex
against the Sun. It is NULL on rows made before v95 and on a star whose population carries no value yet (a cluster
star's comes with the cluster table, GEN.159); a star built from such a row draws its own, as it always did. Each
step is skipped when it is already done, so a stopped run can be started again; the `bright_stars` column goes last
and marks the migration done. The columns go on the end of the table, which MySQL adds without copying the rows.
"""

import sqlalchemy as sa
from alembic import op

revision = "0095"
down_revision = "0094"
branch_labels = None
depends_on = None


def _has_column(connection, table, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar()


def upgrade():
    connection = op.get_bind()
    for table in ("stars", "bright_stars"):
        if not _has_column(connection, table, "metallicity_fe_h"):
            op.execute(f"ALTER TABLE {table} ADD COLUMN metallicity_fe_h DOUBLE NULL")
