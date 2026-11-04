"""The hours a person can be outside in a day, and what limits them (GEN.768).

Schema v187. `planets` and `moons` gain `outside_hours_day` (DOUBLE, hours of one local day) and `outside_hours_limit` (VARCHAR(16): `daylight`, `dose`,
`heat`, `cold`, `air`, `uv` or `vehicle`), worked out by `physics/outside_hours.py` from the stored equipment, daylight and day and night swing when a body is
generated or scored again. Both are NULL here: a row made before v187 gets them when the galaxy is regenerated. A column already there is skipped, so a stopped run
can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0187"
down_revision = "0186"
branch_labels = None
depends_on = None

COLUMNS = (("outside_hours_day", "DOUBLE NULL"), ("outside_hours_limit", "VARCHAR(16) NULL"))
"""tuple: The columns added to `planets` and `moons`, each with its definition."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for column, definition in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} {definition}"))
