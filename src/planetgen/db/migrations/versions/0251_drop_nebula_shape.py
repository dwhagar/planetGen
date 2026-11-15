"""A nebula's shape is the density contour, nothing is stored (GEN.1267).

Schema v251. `nebulae` loses the twelve `shape_*` columns of the GEN.75 metaball shape and the table `nebula_shape_balls` goes. A nebula is now the part of its sphere where the galaxy's
density, read for it, is at least the form level (`galaxy.nebula_density.inside`), the contour the map draws; it has no shape of its own to store. A table that no longer has them is left
as it is, so a stopped run can be started again. `inside_nebula_id` of the rows already stored is not changed here: a galaxy is made again at the epoch that carries this change.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0251"
down_revision = "0250"
branch_labels = None
depends_on = None

SHAPE_COLUMNS = (
    "shape_axis_x", "shape_axis_y", "shape_axis_z", "shape_angle_x", "shape_angle_y", "shape_angle_z",
    "shape_warp_amplitude", "shape_warp_frequency", "shape_warp_octaves", "shape_noise_seed", "shape_iso", "shape_scale",
)
"""tuple: The columns that go from `nebulae`."""


def upgrade():
    connection = op.get_bind()
    connection.execute(sa.text("DROP TABLE IF EXISTS nebula_shape_balls"))
    for column in SHAPE_COLUMNS:
        if helpers.exists_column(connection, "nebulae", column):
            helpers.drop_column_instant(connection, "nebulae", column)
