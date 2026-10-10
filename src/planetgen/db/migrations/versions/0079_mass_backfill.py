"""The mass backfill (GEN.187).

Schema v79. `sector_stats` gains `bright_mass_sol`: the initial mass down to
which a backfill placed every living star of the sector (NULL: none did).
NULL here: sectors backfilled before v79 were filled by luminosity and keep
their `bright_level_sol` (see `schema.sql`'s "v79" note).
"""

import sqlalchemy as sa
from alembic import op

revision = "0079"
down_revision = "0078"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    present = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'sector_stats' AND COLUMN_NAME = 'bright_mass_sol'")).scalar()
    if not present:
        op.execute("ALTER TABLE sector_stats ADD COLUMN bright_mass_sol DOUBLE AFTER level_before_fill_sol")
