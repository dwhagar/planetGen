"""Block totals for the 3-ary block ladder (MAP.151b).

Schema v90. The ladder went from 243, 27, 3 to 243, 81, 27, 9, 3 above a sector, and the 27 and 3 blocks nest
differently, so every `galaxy_block_stats` row of the old ladder is wrong for the new one. This adds
`ladder_step` (the step between the ladder's levels the row was computed for, so the map reads only rows of its own
ladder) and empties the table. (It used to fill the table again; 0091 drops the table and 0092 makes it again and
fills it for the block ladder on shells, so that fill would only compute rows that are thrown away.)
"""

import sqlalchemy as sa
from alembic import op

revision = "0090"
down_revision = "0089"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    has_column = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'galaxy_block_stats' AND COLUMN_NAME = 'ladder_step'")).scalar()
    if not has_column:
        op.execute("ALTER TABLE galaxy_block_stats ADD COLUMN ladder_step TINYINT UNSIGNED NOT NULL DEFAULT 3")
    connection.execute(sa.text("DELETE FROM galaxy_block_stats"))
