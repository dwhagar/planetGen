"""The state of a day and night, the hours of night and the solar days in a year (GEN.1097).

Schema v193. `planets` and `moons` gain `day_night_state` (VARCHAR(24): `day and night`, `tidally locked`, `spin-orbit resonance`, `non-rotating`, `polar day and night`
or `no solid surface`), `night_hours` (DOUBLE, hours of night in one solar day at the equator) and `days_per_year` (DOUBLE, solar days in one orbit), worked out by
`physics/day_night.py` from the stored spin, orbit and polar circle when a body is generated. All three are NULL here: a row made before v193 gets them when the
galaxy is regenerated. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0193"
down_revision = "0192"
branch_labels = None
depends_on = None

COLUMNS = (("day_night_state", "VARCHAR(24) NULL"), ("night_hours", "DOUBLE NULL"), ("days_per_year", "DOUBLE NULL"))
"""tuple: The columns added to `planets` and `moons`, each with its definition."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for column, definition in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} {definition}"))
