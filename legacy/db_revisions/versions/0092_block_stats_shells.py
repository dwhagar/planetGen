"""The `galaxy_block_stats` table for the block ladder on shells (MAP.205).

Schema v92. 0091 dropped the table of the cylinder ladder; this makes it again, keyed by the block's shell, band and
wedge (`galaxy/drill.py`, docs/design/galaxy/spherical-block-ladder.md), and fills it from the sectors and their stats the
galaxy already holds, a shell at a time. Skipped when the table is there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db import alembic_runner

revision = "0092"
down_revision = "0091"
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
            shell_index       INT NOT NULL,
            band_index        INT NOT NULL,
            wedge_index       INT NOT NULL,
            sector_count      INT UNSIGNED NOT NULL,
            systems           BIGINT UNSIGNED NOT NULL,
            expected_systems  DOUBLE NOT NULL,
            stars             BIGINT UNSIGNED NOT NULL,
            age_sum           DOUBLE NOT NULL,
            luminosity_sol    DOUBLE NOT NULL,
            PRIMARY KEY (m, shell_index, band_index, wedge_index)
        ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci""")
    _fill(connection)


def _fill(connection):
    """Computes every block's sums from the generated sectors (`db.block_stats.accumulate`)."""
    from planetgen.db import block_stats

    shells = [row[0] for row in connection.execute(sa.text(
        "SELECT DISTINCT shell_index FROM sectors WHERE shell_index IS NOT NULL ORDER BY shell_index")).fetchall()]
    cells = {}
    for done, shell in enumerate(shells, start=1):
        rows = connection.execute(sa.text(
            "SELECT s.shell_index, s.band_index, s.wedge_index, st.actual_systems, st.expected_systems,"
            " st.actual_stars, st.mean_age_gy, st.total_luminosity_sol FROM (SELECT * FROM sectors WHERE"
            " shell_index = :shell) s LEFT JOIN sector_stats st ON st.shell_index = s.shell_index"
            " AND st.band_index = s.band_index AND st.wedge_index = s.wedge_index"), {"shell": shell})
        for r in rows.fetchall():
            block_stats.accumulate(cells, (r[0], r[1], r[2]), *r[3:])
        alembic_runner.report_progress("Block totals", done, len(shells))
    columns = block_stats.COLUMNS[:10]    # the table as it was at v92; the later columns are 0135's
    insert = sa.text(
        "INSERT INTO galaxy_block_stats (" + ", ".join(columns) + ") VALUES ("
        + ", ".join(":" + column for column in columns) + ")")
    rows = [dict(zip(columns, row)) for row in block_stats.cell_rows(cells)]
    for start in range(0, len(rows), BATCH):
        connection.execute(insert, rows[start:start + BATCH])
