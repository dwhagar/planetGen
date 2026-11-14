"""A belt round an inner pair node of a system's orbit tree (GEN.332, slice C).

Schema v258. `asteroid_belts` gains `host_node_id` (SMALLINT, nullable), as `planets` has had since v104: a belt of the circumbinary disk of an inner `close` node of a tree
(`StarSystem.node_hosts`) is saved with no star row (`star_id` NULL), `parent_type` `node` and `host_node_id` the node's number. NULL for every other belt and for a row made before
v258. A migration does not draw: a galaxy made before slice C has no belt round such a node, and the regenerated galaxy gets them. A column already there is skipped, so a stopped
run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0258"
down_revision = "0257"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_column(connection, "asteroid_belts", "host_node_id"):
        connection.execute(sa.text("ALTER TABLE asteroid_belts ADD COLUMN host_node_id SMALLINT NULL"))
