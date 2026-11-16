"""The mass a belt's star threw inward, stored for the white dwarf's accretion rate (GEN.1231).

Schema v266. `asteroid_belts` gains `thrown_inward_kg` (DOUBLE, NULL): the mass of the belt the star's mass loss threw inward (`physics/belt_history.py`), the rock that
falls on a white dwarf. It is NULL for a belt of no remnant system and for a row made before v266, which gets it when the galaxy is regenerated; nothing is backfilled. A column already there is
skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0266"
down_revision = "0265"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    if not helpers.exists_column(connection, "asteroid_belts", "thrown_inward_kg"):
        connection.execute(sa.text("ALTER TABLE asteroid_belts ADD COLUMN thrown_inward_kg DOUBLE NULL"))
