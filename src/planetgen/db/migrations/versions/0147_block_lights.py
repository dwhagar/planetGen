"""The light of the pre-placed objects in each block of the ladder (MAP.259).

Schema v147. A new table, `galaxy_block_lights`: for every block of the eleven levels above a sector that holds a
bright star no system was built around, or a scattered phenomenon not yet built, the count, summed luminosity and
luminosity-weighted position of each (and the luminosity-weighted temperature of the stars). The single-space map's
nodes draw a block they have not opened from these sums, so a node's glow and its centre of light are one row read.
The table is filled from `bright_stars` and `phenomenon_scatter` (`db.block_lights.rebuild`), a unit of shells at a time.
The table is emptied before the fill, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0147"
down_revision = "0146"
branch_labels = None
depends_on = None

TABLE = """
    CREATE TABLE IF NOT EXISTS galaxy_block_lights (
        m                    INT UNSIGNED NOT NULL,
        shell_index          INT NOT NULL,
        band_index           INT NOT NULL,
        wedge_index          INT NOT NULL,
        star_count           BIGINT UNSIGNED NOT NULL,
        star_luminosity_sol  DOUBLE NOT NULL,
        star_lpx              DOUBLE NOT NULL,
        star_lpy              DOUBLE NOT NULL,
        star_lpz              DOUBLE NOT NULL,
        star_lt              DOUBLE NOT NULL,
        point_count          BIGINT UNSIGNED NOT NULL,
        point_luminosity_sol DOUBLE NOT NULL,
        point_lpx             DOUBLE NOT NULL,
        point_lpy             DOUBLE NOT NULL,
        point_lpz             DOUBLE NOT NULL,
        PRIMARY KEY (m, shell_index, band_index, wedge_index)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    connection = op.get_bind()
    connection.execute(sa.text(TABLE))    # `schema.sql` is applied first, so the table is usually there already
    from planetgen.db import block_lights

    block_lights.rebuild(helpers.StoreConnection(connection))
