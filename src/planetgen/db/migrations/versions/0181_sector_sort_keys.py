"""The Sectors list sorts by systems, density and position by key (PERF.208).

Schema v181. `sector_system_counts` gains `system_density` (systems per cubic milliparsec, NULL for a sector with no galaxy position) and
`quadrant` (`'I'` to `'IV'`, NULL for a sector with no galaxy position), each with an index ending in `sector_id`, and holds a row for every
sector, not only for those with systems: the migration fills it again. It is derived data (`query.refresh_sector_system_counts` keeps it).
The indexes are built online (`LOCK=NONE`); a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0181"
down_revision = "0180"
branch_labels = None
depends_on = None

COUNT_INDEXES = (
    ("idx_sector_system_counts_density", ("system_density", "sector_id")),
    ("idx_sector_system_counts_quadrant", ("quadrant", "sector_id")),
)
QUADRANT = ("CASE WHEN sec.center_y_pc >= 0 THEN (CASE WHEN sec.center_x_pc >= 0 THEN 'I' ELSE 'II' END) "
            "ELSE (CASE WHEN sec.center_x_pc < 0 THEN 'III' ELSE 'IV' END) END")

FILL = f"""
    INSERT INTO sector_system_counts (sector_id, system_count, system_density, quadrant)
    SELECT sec.id, COALESCE(n.system_count, 0),
           CASE WHEN sec.center_x_pc IS NULL THEN NULL ELSE COALESCE(n.system_count, 0) / POW(sec.edge_mpc, 3) END,
           CASE WHEN sec.center_x_pc IS NULL THEN NULL ELSE {QUADRANT} END
    FROM sectors sec
    LEFT JOIN (SELECT sector_id, COUNT(*) AS system_count FROM star_systems WHERE sector_id IS NOT NULL
               GROUP BY sector_id) n ON n.sector_id = sec.id
"""


def upgrade():
    connection = op.get_bind()
    helpers.add_column_instant(connection, "sector_system_counts", "system_density", "DOUBLE NULL", after="system_count")
    helpers.add_column_instant(connection, "sector_system_counts", "quadrant", "VARCHAR(8) NULL", after="system_density")
    for name, columns in COUNT_INDEXES:
        helpers.add_index_online(connection, "sector_system_counts", name, columns)
    connection.execute(sa.text("DELETE FROM sector_system_counts"))
    connection.execute(sa.text(FILL))
