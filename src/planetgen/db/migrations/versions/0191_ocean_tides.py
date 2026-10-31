"""The open-ocean equilibrium tide of a world with a sea (GEN.1034).

Schema v191. `planets` and `moons` gain three nullable columns, derived by `physics/ocean_tides.py` from the stored mass, radius, spin, orbits, moons and hydrosphere:
`tide_range_spring_m` and `tide_range_neap_m` (DOUBLE, metres: the high-to-low equilibrium range with every tide-raiser in line, and with the largest against the rest; 0 for an
ice-covered ocean, NULL for a body without a free sea) and `ocean_liquid` (VARCHAR(40): the liquid of that sea, water until the surface liquid inventory of GEN.1035 names others). A row
made before v191 keeps NULL there; nothing is backfilled. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0191"
down_revision = "0190"
branch_labels = None
depends_on = None

COLUMNS = (("tide_range_spring_m", "DOUBLE"), ("tide_range_neap_m", "DOUBLE"), ("ocean_liquid", "VARCHAR(40)"))


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for column, kind in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} {kind} NULL"))
