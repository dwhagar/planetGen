"""The plan stores which shape levers it was pinned with (DB.81).

Schema v114. `galaxy_shape` gains `pinned_stellar_mass_solar` (the mass a pinned Scale dial gave the galaxy, NULL when the
seed drew it) and `galaxy_shape_pins` holds one row per lever the plan was told, so the galaxy page can say which values
were pinned and which were drawn. No row is touched: a galaxy planned before this reads as everything drawn, and its
mass chain is drawn from its seed as before. A step already done is left alone, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0114"
down_revision = "0113"
branch_labels = None
depends_on = None


def _has_column(connection, table, column):
    return bool(connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = :table_name AND COLUMN_NAME = :column_name"),
        {"table_name": table, "column_name": column}).scalar())


def upgrade():
    connection = op.get_bind()
    if not _has_column(connection, "galaxy_shape", "pinned_stellar_mass_solar"):
        op.execute("ALTER TABLE galaxy_shape ADD COLUMN pinned_stellar_mass_solar DOUBLE NULL AFTER halo_min_luminosity_sol")
    op.execute(
        "CREATE TABLE IF NOT EXISTS galaxy_shape_pins ("
        " name VARCHAR(40) NOT NULL PRIMARY KEY,"
        " value DOUBLE NOT NULL"
        ") ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci")
