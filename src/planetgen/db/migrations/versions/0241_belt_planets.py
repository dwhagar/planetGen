"""The heaviest round bodies of an asteroid belt are planets (GEN.1510).

Schema v241. `planets` gains `belt_order` and `belt_rank` (INT, nullable): a planet promoted from a belt (`generation/belt_planets.py`) keeps the belt's number among its
star's belts (counted as `irregular_moons.belt_order` counts) and the rank it had among the belt's members, so the designation `<belt ID>.<rank>` still resolves. NULL for every
other planet and for a row made before v241. A planet made before v241 has no promoted body; the galaxy is regenerated (a migration does not draw). A column already there is
skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0241"
down_revision = "0240"
branch_labels = None
depends_on = None

COLUMNS = ("belt_order", "belt_rank")
"""tuple: The columns `planets` gains."""


def upgrade():
    connection = op.get_bind()
    for column in COLUMNS:
        if not helpers.exists_column(connection, "planets", column):
            connection.execute(sa.text(f"ALTER TABLE planets ADD COLUMN {column} INT NULL"))
