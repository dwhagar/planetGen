"""The spheres the sky tiers have drawn (VIEW.25).

Schema v214. `sky_tier_draws` is new: one row for each tier of a system's sky (`generation/sky_tiers.py`) that charting has drawn over a sphere with the star backfill
(`generation/sky_charting.py`): whether the tier is a luminosity band or a mass ring, the band's range, the centre of the sphere (the 50 parsec cell of the system, the same cell
the glow band is kept for) and its radius, the limiting magnitude it was built for, the stars written, the seconds it took and the generator epoch. A cell inside a sphere is
complete down to the tier's floor without a `sector_stats` row of its own (a sky sphere holds billions of cells); the backfill and the fill read the rows here
(`db/sky_tier_ledger.py`). A table already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0214"
down_revision = "0213"
branch_labels = None
depends_on = None

TABLE = """
CREATE TABLE sky_tier_draws (
    id               BIGINT UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    kind             VARCHAR(12) NOT NULL,
    band_low         DOUBLE NOT NULL,
    band_high        DOUBLE NULL,
    center_x_pc      DOUBLE NOT NULL,
    center_y_pc      DOUBLE NOT NULL,
    center_z_pc      DOUBLE NOT NULL,
    radius_pc        DOUBLE NOT NULL,
    limit_mag        DOUBLE NOT NULL,
    stars            INT UNSIGNED NOT NULL DEFAULT 0,
    seconds          DOUBLE NOT NULL DEFAULT 0,
    generator_epoch  SMALLINT UNSIGNED NULL,
    batches_done     INT UNSIGNED NOT NULL DEFAULT 0,
    complete         TINYINT(1) NOT NULL DEFAULT 1,
    drawn_at         TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP(3),
    CONSTRAINT chk_sky_tier_draws_kind CHECK (kind IN ('luminosity', 'mass')),
    CONSTRAINT chk_sky_tier_draws_radius CHECK (radius_pc >= 0),
    UNIQUE KEY uq_sky_tier_draws_cell (kind, band_low, center_x_pc, center_y_pc, center_z_pc, radius_pc),
    KEY idx_sky_tier_draws_kind (kind, band_low)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci
"""


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_table(connection, "sky_tier_draws"):
        connection.execute(sa.text(TABLE))
