"""The named regions of the Galaxy Map keep their sector count, share charted and extent (DB.107).

Schema v161. A new table, `region_stats`: for every 9-block of the block ladder that holds at least one generated sector (a
region, `galaxy/named_regions.py`), the sectors it holds (`sector_count`) and could hold (`capacity`), how many face-connected
clumps they make (`clump_count`) and the bounding sphere of the largest clump (`extent_x_pc`, `extent_y_pc`, `extent_z_pc`,
`extent_radius_pc`). Nothing is drawn, so every value is a function of the sectors stored. The rows are made by the rebuild of the
migration after it (0168, `db.block_stats.rebuild`), which reads the sector masses that migration adds; this one makes the table.
"""

import sqlalchemy as sa
from alembic import op

revision = "0161"
down_revision = "0160"
branch_labels = None
depends_on = None

TABLE = """
    CREATE TABLE IF NOT EXISTS region_stats (
        shell_index       INT NOT NULL,
        band_index        INT NOT NULL,
        wedge_index       INT NOT NULL,
        sector_count      INT UNSIGNED NOT NULL,
        capacity          INT UNSIGNED NOT NULL,
        clump_count       INT UNSIGNED NOT NULL,
        extent_x_pc       DOUBLE NOT NULL,
        extent_y_pc       DOUBLE NOT NULL,
        extent_z_pc       DOUBLE NOT NULL,
        extent_radius_pc  DOUBLE NOT NULL,
        PRIMARY KEY (shell_index, band_index, wedge_index)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    connection = op.get_bind()
    connection.execute(sa.text(TABLE))    # `schema.sql` is applied first, so the table is usually there already
