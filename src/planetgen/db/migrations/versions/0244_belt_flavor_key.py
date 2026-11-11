"""An asteroid belt shows a flavor line (GEN.1433).

Schema v244. `asteroid_belts` gains `flavor_key` (VARCHAR(32), nullable): the id of the flavor line the belt shows, a line of the `belt` scope of `physics/data/flavor_tree.json`.
NULL for a belt with no line and for a row made before v244, which gets it when the galaxy is regenerated (a migration does not draw). A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0244"
down_revision = "0243"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_column(connection, "asteroid_belts", "flavor_key"):
        connection.execute(sa.text("ALTER TABLE asteroid_belts ADD COLUMN flavor_key VARCHAR(32) NULL"))
