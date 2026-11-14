"""Galaxy shape density settings (ADM.49).

Schema v81. `galaxy_shape` gains `arm_level` (1: the thin disk's level between
the arm crest and the inter-arm trough) and `core_amplitude` (0: no core), the
values the model had before, so an existing galaxy keeps its density.
"""

import sqlalchemy as sa
from alembic import op

revision = "0081"
down_revision = "0080"
branch_labels = None
depends_on = None


def _has(connection, column):
    return connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'galaxy_shape' AND COLUMN_NAME = :c"), {"c": column}).scalar()


def upgrade():
    connection = op.get_bind()
    if not _has(connection, "arm_level"):
        op.execute("ALTER TABLE galaxy_shape ADD COLUMN arm_level DOUBLE NOT NULL DEFAULT 1 AFTER arm_amplitude")
    if not _has(connection, "core_amplitude"):
        op.execute("ALTER TABLE galaxy_shape ADD COLUMN core_amplitude DOUBLE NOT NULL DEFAULT 0 AFTER arm_level")
