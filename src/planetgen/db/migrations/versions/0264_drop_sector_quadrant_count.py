"""The Sectors list no longer keeps a Quadrant for each sector (DB.112).

Schema v264. `sector_system_counts` loses `quadrant` and its index `idx_sector_system_counts_quadrant` (both added by migration 0181, PERF.208): the Quadrants of the Sectors list are gone from the site
(the pages name the regions where sectors are charted instead, UX.272), so nothing filters or sorts by it. The table is derived data (`query.refresh_sector_system_counts` keeps it) and the column
is not part of any content stamp, so no sector, system or storage checksum changes. The index goes first, online (`LOCK=NONE`); the column goes instantly. A stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0264"
down_revision = "0263"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    helpers.drop_index_online(connection, "sector_system_counts", "idx_sector_system_counts_quadrant")
    helpers.drop_column_instant(connection, "sector_system_counts", "quadrant")
