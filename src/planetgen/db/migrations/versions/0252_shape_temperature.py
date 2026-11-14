"""The temperature a body took its shape at (GEN.1201).

Schema v252. `planets`, `moons` and `rogue_planets` gain `shape_temperature_k` (DOUBLE, K, nullable): the temperature the body took its shape at, worked out by `physics/shape_temperature.py` from its
surface temperature, the heat of its assembly, its radioactive content and its tides, and read by the roundness engine for the yield strength and the roundness ratio (`physics/body_strength.py`).
A row made before v252 has NULL here and keeps the cold strength it was made with until the galaxy is regenerated (a migration does not draw). A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0252"
down_revision = "0251"
branch_labels = None
depends_on = None

TABLES = ("planets", "moons", "rogue_planets")
"""tuple: The tables that gain the column."""

COLUMN = "shape_temperature_k"


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        if not helpers.exists_column(connection, table, COLUMN):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {COLUMN} DOUBLE NULL"))
