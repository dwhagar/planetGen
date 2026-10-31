"""The age of a ring (GEN.1057).

Schema v192. `planets`, `moons` and `rogue_planets` gain `ring_age_yr` (DOUBLE, NULL), the years since the event that made the ring (a moon that fell inside the Roche limit, a comet that
passed inside it; `physics/ring_origin.py`), placed after `ring_ice_fraction` as `schema.sql` has it. NULL for a body with no ring, for a resonance ring of a small body and for a row made before
v192: nothing is backfilled, and the regenerated galaxy replaces those rows. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0192"
down_revision = "0191"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons", "rogue_planets"):
        if not helpers.exists_column(connection, table, "ring_age_yr"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN ring_age_yr DOUBLE NULL AFTER ring_ice_fraction"))
