"""The parts of the block light totals (GEN.1476).

Schema v236. One table, `galaxy_block_light_parts`: for each sector cell that holds a bright star no system was built around or a scattered phenomenon not yet built, the same sums
`galaxy_block_lights` holds (see `schema.sql` and `db/block_lights.py`). The scatter writers add to it as they write, and the block light stage rolls it up. A galaxy scattered before v236 has
no parts; the table is not filled here (it is derived data, and the first block light stage checks it against the scatter and makes it again from the objects when they differ, `block_lights.roll`).
A table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

revision = "0236"
down_revision = "0235"
branch_labels = None
depends_on = None

PARTS = """
    CREATE TABLE IF NOT EXISTS galaxy_block_light_parts (
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
        PRIMARY KEY (shell_index, band_index, wedge_index)
    ) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci"""


def upgrade():
    op.get_bind().execute(sa.text(PARTS))
