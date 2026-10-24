"""Regional temperatures and the day and night swing (GEN.656).

Schema v137. `planets` and `moons` gain, all DOUBLE and NULL for a gas giant and for a body saved without a score, the
temperatures of `physics/regional_temperature.py` (kelvin): `temp_polar_k`, `temp_north_k`, `temp_equator_k`,
`temp_south_k` and `temp_swing_k` for a world that turns on its own, and `temp_day_side_k`, `temp_terminator_k` and
`temp_night_side_k` for a tidally locked one.

They are a pure function of stored columns, so no generator epoch: `planetgen rescore-habitability`
(`store.refresh_habitability`) fills them for a galaxy made before this revision, and a body not re-scored shows no
regional line until then. Each table is one `ALGORITHM=COPY` statement (`planets` and `moons` have a FULLTEXT index);
a table already changed is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0137"
down_revision = "0136"
branch_labels = None
depends_on = None

COLUMNS = (
    "temp_polar_k", "temp_north_k", "temp_equator_k", "temp_south_k", "temp_swing_k", "temp_day_side_k",
    "temp_terminator_k", "temp_night_side_k",
)
"""tuple: The DOUBLE columns added to `planets` and `moons`."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        changes = [f"ADD COLUMN {column} DOUBLE NULL" for column in COLUMNS
                   if not helpers.exists_column(connection, table, column)]
        if changes:
            connection.execute(sa.text(f"ALTER TABLE {table} {', '.join(changes)}, ALGORITHM=COPY"))
