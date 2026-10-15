"""The `galaxy_block_stats` table: the drill-down stage sums stored per block of each ladder level (DB.25).

Schema v86. The migration creates the table (and, until 0091 dropped the table, filled it from the sectors and their
stats that the galaxy already holds, a ring at a time; now it leaves it empty). Skipped when the table is there, so a
stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0086"
down_revision = "0085"
branch_labels = None
depends_on = None

BATCH = 1000


def upgrade():
    connection = op.get_bind()
    exists = connection.execute(sa.text(
        "SELECT COUNT(*) FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
        " AND TABLE_NAME = 'galaxy_block_stats'")).scalar()
    if exists:
        return
    op.execute("""
        CREATE TABLE galaxy_block_stats (
            m                 SMALLINT UNSIGNED NOT NULL,
            ring_index        INT NOT NULL,
            wedge             INT NOT NULL,
            slab              INT NOT NULL,
            sector_count      INT UNSIGNED NOT NULL,
            systems           BIGINT UNSIGNED NOT NULL,
            expected_systems  DOUBLE NOT NULL,
            stars             BIGINT UNSIGNED NOT NULL,
            age_sum           DOUBLE NOT NULL,
            luminosity_sol    DOUBLE NOT NULL,
            PRIMARY KEY (m, ring_index, wedge, slab)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""")
    _fill(connection)


def _fill(connection):
    """Left empty: revision 0091 (DB.32) drops this table with the cylinder's block ladder, so the totals this once
    computed from the generated sectors (`db/block_stats.py`, deleted with it) would never be read."""
