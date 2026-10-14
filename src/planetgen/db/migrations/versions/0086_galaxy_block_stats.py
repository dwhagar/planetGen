"""The `galaxy_block_stats` table: the drill-down stage sums stored per block of each ladder level (DB.25).

Schema v86. The migration creates the table and fills it from the sectors and their stats that the galaxy
already holds, a ring at a time. Skipped when the table is there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner

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
    """Computes every block's sums from the generated sectors (`db.block_stats.accumulate`)."""
    from planetgen.db import block_stats

    rings = [row[0] for row in connection.execute(sa.text(
        "SELECT DISTINCT ring_index FROM sectors WHERE ring_index IS NOT NULL ORDER BY ring_index")).fetchall()]
    cells = {}
    for done, ring in enumerate(rings, start=1):
        rows = connection.execute(sa.text(
            "SELECT s.ring_index, s.layer_index, s.ring_slot_index, st.actual_systems, st.expected_systems,"
            " st.actual_stars, st.mean_age_gy, st.total_luminosity_sol FROM (SELECT * FROM sectors WHERE"
            " ring_index = :ring) s LEFT JOIN sector_stats st ON st.ring_index = s.ring_index"
            " AND st.layer_index = s.layer_index AND st.ring_slot_index = s.ring_slot_index"), {"ring": ring})
        for r in rows.fetchall():
            block_stats.accumulate(cells, (r[0], r[1], r[2]), *r[3:])
        alembic_runner.report_progress("Block totals", done, len(rings))
    insert = sa.text(
        "INSERT INTO galaxy_block_stats (" + ", ".join(block_stats.COLUMNS) + ") VALUES ("
        + ", ".join(":" + column for column in block_stats.COLUMNS) + ")")
    rows = [dict(zip(block_stats.COLUMNS, row)) for row in block_stats.cell_rows(cells)]
    for start in range(0, len(rows), BATCH):
        connection.execute(insert, rows[start:start + BATCH])
