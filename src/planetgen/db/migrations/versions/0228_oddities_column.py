"""The oddities a body has (GEN.1047).

Schema v228. `planets`, `moons` and `rogue_planets` gain `oddities` (VARCHAR(40), nullable): the ids of the rows of the oddity table
(`physics/oddities.py`) that gave the body its rare, physically plausible state, comma separated. NULL for no oddity and for a row made before
v228, which gets it when the galaxy is regenerated (a migration does not draw). A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0228"
down_revision = "0227"
branch_labels = None
depends_on = None

TABLES = ("planets", "moons", "rogue_planets")
"""tuple: The tables that gain the column."""


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        if not helpers.exists_column(connection, table, "oddities"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN oddities VARCHAR(40) NULL"))

