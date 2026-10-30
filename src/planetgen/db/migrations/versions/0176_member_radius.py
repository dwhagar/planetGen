"""The extent of a scatter member (GEN.952).

Schema v176. `scatter_members` gains `radius_au` (DOUBLE, NULL: the member's radius in AU), so a page can work out the X-ray luminosity of an
accreting pair of a stored tree (`physics/accretion.x_ray_luminosity_w` needs the radius of the accretor and of the donor). NULL for a row made
before v176; nothing is backfilled, a row gets its radius when the galaxy is regenerated. A column already there is skipped, so a stopped run
can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0176"
down_revision = "0175"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_column(connection, "scatter_members", "radius_au"):
        connection.execute(sa.text("ALTER TABLE scatter_members ADD COLUMN radius_au DOUBLE NULL"))
