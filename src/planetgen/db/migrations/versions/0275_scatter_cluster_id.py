"""The cluster a scattered object belongs to (GEN.163).

Schema v275 (a placeholder: the train that merges this renumbers it to the next free revision). `phenomenon_scatter` gains `cluster_id` (BIGINT UNSIGNED, the `galaxy_components.id` of the globular cluster the
object was drawn in, `component_pass.scatter_phenomena`): NULL for an object of the field or of a satellite and for a row made before v275, so nothing is backfilled (a migration does not draw). A pulsar that
has one may hold a captured giant (`generation/pulsar_planets.py`). A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0275"
down_revision = "0274"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_column(connection, "phenomenon_scatter", "cluster_id"):
        connection.execute(sa.text("ALTER TABLE phenomenon_scatter ADD COLUMN cluster_id BIGINT UNSIGNED NULL"))


def downgrade():
    connection = op.get_bind()
    if helpers.exists_column(connection, "phenomenon_scatter", "cluster_id"):
        connection.execute(sa.text("ALTER TABLE phenomenon_scatter DROP COLUMN cluster_id"))
