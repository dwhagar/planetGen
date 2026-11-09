"""The light a circumbinary planet keeps after its two stars eclipse each other (GEN.1484).

Schema v241. `planets`, `moons` and `rogue_planets` gain `light_factor` (DOUBLE, nullable): the share of its host's light a planet of a close pair keeps (`physics/eclipses.py`), NULL for
every other body and for a row made before v241, which gets it when the galaxy is regenerated (a migration does not draw). A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0241"
down_revision = "0240"
branch_labels = None
depends_on = None

TABLES = ("planets", "moons", "rogue_planets")
"""tuple: The tables that gain the column."""


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        if not helpers.exists_column(connection, table, "light_factor"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN light_factor DOUBLE NULL"))
