"""The share of the surface the ice covers (GEN.1567).

Schema v244. `planets`, `moons` and `rogue_planets` gain `ice_cover_fraction` (DOUBLE, nullable): the share of the surface the ice covers (the lid of an ice-covered ocean, or a body that is ice
through). `ocean_fraction` now counts the liquid sea alone (it counted liquid or frozen water, so every ice-covered world stored 1.0 and the ocean and desert classes read it as a sea), and
`land_fraction` is the rest of the surface after the liquid sea. A row made before v244 has NULL here and keeps the `ocean_fraction` and `land_fraction` it was made with until the galaxy is
regenerated (a migration does not draw). A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0244"
down_revision = "0243"
branch_labels = None
depends_on = None

TABLES = ("planets", "moons", "rogue_planets")
"""tuple: The tables that gain the column."""

COLUMN = "ice_cover_fraction"


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        if not helpers.exists_column(connection, table, COLUMN):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {COLUMN} DOUBLE NULL"))
