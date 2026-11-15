"""A sphere of sky can be drawn in parts (OPS.253).

Schema v262. `sky_tier_draws` gains `part_index` and `part_count` (SMALLINT UNSIGNED, 0 and 0 for a whole sphere): a sphere drawn by several units at once (`generation/sky_charting.py`, the
sliced draw of the backfill stage) keeps one row for each part, which says how many batches of its share of the sphere's rows are written, so a stopped part goes on from its last batch. The
parts' rows never count in the ledger (`sky_tier_ledger.read` reads whole spheres only); the stage's last step records the sphere itself and drops them. The unique key of a draw gains the two columns.
A column or key already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0262"
down_revision = "0261"
branch_labels = None
depends_on = None

OLD_KEY = ("kind", "band_low", "center_x_pc", "center_y_pc", "center_z_pc", "radius_pc")
NEW_KEY = (*OLD_KEY, "part_count", "part_index")


def upgrade():
    connection = op.get_bind()
    for column in ("part_index", "part_count"):
        if not helpers.exists_column(connection, "sky_tier_draws", column):
            connection.execute(sa.text(f"ALTER TABLE sky_tier_draws ADD COLUMN {column} SMALLINT UNSIGNED NOT NULL DEFAULT 0"))
    if helpers.index_columns(connection, "sky_tier_draws", "uq_sky_tier_draws_cell") != list(NEW_KEY):
        if helpers.exists_index(connection, "sky_tier_draws", "uq_sky_tier_draws_cell"):
            connection.execute(sa.text("ALTER TABLE sky_tier_draws DROP INDEX uq_sky_tier_draws_cell"))
        connection.execute(sa.text(f"ALTER TABLE sky_tier_draws ADD UNIQUE KEY uq_sky_tier_draws_cell ({', '.join(NEW_KEY)})"))
