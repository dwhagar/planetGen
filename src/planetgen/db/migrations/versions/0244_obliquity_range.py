"""The range a planet's tilt wanders over (GEN.1535).

Schema v244. `planets`, `moons` and `rogue_planets` gain `obliquity_min_deg` and `obliquity_max_deg` (DOUBLE, nullable): the lowest and highest tilt of a free-spinning planet's axis over
the age of its star, worked from the spin-orbit secular resonances (`physics/obliquity_variation.py`), NULL for a body that was not worked (a moon, a rogue planet, a locked or resonant body) and for every row made before v244, which gets it
when the galaxy is regenerated (a migration does not draw). A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0244"
down_revision = "0243"
branch_labels = None
depends_on = None

TABLES = ("planets", "moons", "rogue_planets")
"""tuple: The tables that gain the columns."""

COLUMNS = ("obliquity_min_deg", "obliquity_max_deg")


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        for column in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} DOUBLE NULL"))
