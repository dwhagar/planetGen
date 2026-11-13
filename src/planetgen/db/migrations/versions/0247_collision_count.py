"""The mergers recorded in making a planet, kept (DB.122).

Schema v247. `planets`, `moons` and `rogue_planets` gain `collision_count` (SMALLINT UNSIGNED, nullable): how many mergers the planet architecture recorded in making the body. The oddity
table reads it (`physics/oddities.py`, the history rows), so a loaded planet that read 0 lost its collision oddity.

NULL here: a body made before v247 has no count to read back, and a NULL reads as 0. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0247"
down_revision = "0246"
branch_labels = None
depends_on = None

_TABLES = ("planets", "moons", "rogue_planets")
"""tuple: The tables that gain `collision_count`."""


def upgrade():
    connection = op.get_bind()
    for table in _TABLES:
        if not helpers.exists_column(connection, table, "collision_count"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN collision_count SMALLINT UNSIGNED NULL"))
