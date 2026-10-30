"""The rotation period of a scattered star (GEN.846).

Schema v176. `bright_stars` gains `rotation_period_hours` (DOUBLE, NULL: the star's spin period in hours), so a star the scatter placed
and nothing has built yet turns on the map and its page gives a period like a built star's (`stars.rotation_period_hours` is stored
since v68). NULL for a row made before v176; nothing is backfilled, a row gets its period when the galaxy is regenerated. A column
already there is skipped, so a stopped run can be started again.
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
    if not helpers.exists_column(connection, "bright_stars", "rotation_period_hours"):
        connection.execute(sa.text("ALTER TABLE bright_stars ADD COLUMN rotation_period_hours DOUBLE NULL"))
