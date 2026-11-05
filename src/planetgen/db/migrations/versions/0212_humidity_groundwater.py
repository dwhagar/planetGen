"""The relative humidity and the groundwater of a body (GEN.1221).

Schema v212. `planets`, `moons` and `rogue_planets` gain `relative_humidity` (DOUBLE, the air's water vapour over saturation at the surface temperature, 0.05 to 0.99, the number the wet
bulb and the survival limits use) and `groundwater` (VARCHAR(10): `none`, `reachable` or `deep`, the liquid water a world with none at its surface holds underground,
`physics/habitability_world.py`). Both are NULL here: a row made before v212 gets them when the galaxy is regenerated. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0212"
down_revision = "0211"
branch_labels = None
depends_on = None

COLUMNS = (("relative_humidity", "DOUBLE NULL"), ("groundwater", "VARCHAR(10) NULL"))
"""tuple: The columns added to `planets`, `moons` and `rogue_planets`, each with its definition."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons", "rogue_planets"):
        for column, definition in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} {definition}"))
