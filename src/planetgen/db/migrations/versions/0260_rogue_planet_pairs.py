"""A rogue planet that left its donor bound to a companion keeps the link and the orbit (GEN.332, slice D).

Schema v260. `rogue_planets` gains `pair_event_index` (INT, nullable: the event index, of the same donor, of the companion; each of the two rows names the other), `pair_separation_au` (DOUBLE,
nullable: the distance between the two, AU) and `pair_period_yr` (DOUBLE, nullable: the period of the pair round each other, years). NULL for a rogue that left alone, which is every row made
before v260, so nothing is rewritten. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0260"
down_revision = "0259"
branch_labels = None
depends_on = None

_COLUMNS = (("pair_event_index", "INT"), ("pair_separation_au", "DOUBLE"), ("pair_period_yr", "DOUBLE"))
"""tuple: The columns `rogue_planets` gains, with their types."""


def upgrade():
    connection = op.get_bind()
    for column, kind in _COLUMNS:
        if not helpers.exists_column(connection, "rogue_planets", column):
            connection.execute(sa.text(f"ALTER TABLE rogue_planets ADD COLUMN {column} {kind} NULL"))
