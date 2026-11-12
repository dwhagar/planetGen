"""What a belt and a comet are made of, kept (DB.124).

Schema v247. `asteroid_belts` and `comets` gain `makeup`, the drawn composition (`physics/small_body_composition.SmallBodyMakeup`) as one JSON object: the shares of the materials,
the porosity, the densities, the albedo and what else the draw stated. Before, a saved belt or comet came back with the words and the chondrite group but no makeup, so a loaded body could not
say its porosity, density or albedo.

NULL here: a body made before v247 has no makeup to read back (the generator epoch makes it again). A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0247"
down_revision = "0246"
branch_labels = None
depends_on = None

_TABLES = ("asteroid_belts", "comets")
"""tuple: The tables that gain `makeup`."""


def upgrade():
    connection = op.get_bind()
    for table in _TABLES:
        if not helpers.exists_column(connection, table, "makeup"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN makeup TEXT NULL"))
