"""The `nebula_field` table: the galaxy's field clouds of classes M and N, one row each (DB.24).

Schema v83. The migration only creates the empty table; `planetgen plan` fills
it from the galaxy seed, so a galaxy upgraded to v83 has no rows until its next
plan (or `plan --redo-scatters`). Skipped when it is already there, so a stopped
run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0083"
down_revision = "0082"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    exists = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'nebula_field'")
    ).scalar()
    if exists:
        return
    op.execute("""
        CREATE TABLE nebula_field (
            cell_i        INT NOT NULL,
            cell_j        INT NOT NULL,
            cell_k        INT NOT NULL,
            n             INT UNSIGNED NOT NULL,
            nebula_class  CHAR(1) NOT NULL,
            radius_ly     DOUBLE NOT NULL,
            center_x_pc   DOUBLE NOT NULL,
            center_y_pc   DOUBLE NOT NULL,
            center_z_pc   DOUBLE NOT NULL,
            uid           BINARY(10),
            stored        TINYINT(1) NOT NULL DEFAULT 0,
            PRIMARY KEY (cell_i, cell_j, cell_k, n),
            KEY idx_nebula_field_uid (uid)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""")
