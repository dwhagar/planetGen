"""What the star's evolution did to a belt, stored for its page (GEN.1232).

Schema v200. `asteroid_belts` gains `history_text` (TEXT, NULL): the sentences of `physics/belt_history.py` for the events a belt went through (the part the giant branch engulfed, the share the mass loss threw inward),
as one paragraph. It is NULL for a belt nothing happened to and for a row made before v200, which gets it when the galaxy is regenerated; nothing is backfilled. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0200"
down_revision = "0199"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_column(connection, "asteroid_belts", "history_text"):
        connection.execute(sa.text("ALTER TABLE asteroid_belts ADD COLUMN history_text TEXT NULL"))
