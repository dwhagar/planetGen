"""The day and night time shares and the hours of daylight of a world (GEN.767).

Schema v170. `planets` and `moons` gain six nullable DOUBLE columns, derived by `physics/daylight.py` from the stored rotation, orbit, tilt
and spin axis and the star's radius: `day_share_equator`, `day_share_mean`, `parent_shadow_share` (a moon in its planet's shadow),
`solar_day_hours`, `daylight_hours` and `polar_circle_deg`. A row made before v170 keeps NULL there; nothing is backfilled. A column already
there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0170"
down_revision = "0169"
branch_labels = None
depends_on = None

COLUMNS = ("day_share_equator", "day_share_mean", "parent_shadow_share", "solar_day_hours", "daylight_hours", "polar_circle_deg")


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for column in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} DOUBLE NULL"))
