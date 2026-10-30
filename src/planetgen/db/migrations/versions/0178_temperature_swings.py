"""A day and night swing for every temperature the sheet shows (GEN.848).

Schema v178. `planets` and `moons` gain `temp_polar_swing_k`, `temp_north_swing_k`, `temp_south_swing_k`, `temp_mean_swing_k`,
`temp_day_side_swing_k`, `temp_terminator_swing_k` and `temp_night_side_swing_k` (DOUBLE, kelvin of difference, peak to peak); the
equator's was `temp_swing_k` since v137.

All NULL here: the swings are worked out from the temperatures when a body is generated (`physics/regional_temperature.py`), so a row
made before v178 gets them when the galaxy is regenerated. Each column is skipped when it is already there, so a stopped run can be
started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0178"
down_revision = "0177"
branch_labels = None
depends_on = None

COLUMNS = (
    "temp_polar_swing_k", "temp_north_swing_k", "temp_south_swing_k", "temp_mean_swing_k", "temp_day_side_swing_k",
    "temp_terminator_swing_k", "temp_night_side_swing_k",
)
"""tuple: The DOUBLE columns added to `planets` and `moons`."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for column in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                op.execute(f"ALTER TABLE {table} ADD COLUMN {column} DOUBLE NULL")
