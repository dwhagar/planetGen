"""The climate of a body in seven bands of latitude (GEN.1108).

Schema v196. `planets` and `moons` gain the 49 band columns `band1_mean_k` to `band7_state` (`physics/regional_temperature.py`: for each band the mean temperature, the day and
night swing, the hours of day and of night, the longest polar day and night and the state word, `physics/latitude_bands.py`) and lose the four latitude bands, the three zones of a
locked world and their swings (`temp_polar_k`, `temp_north_k`, `temp_equator_k`, `temp_south_k`, `temp_swing_k`, `temp_day_side_k`, `temp_terminator_k`, `temp_night_side_k`,
`temp_polar_swing_k`, `temp_north_swing_k`, `temp_south_swing_k`, `temp_day_side_swing_k`, `temp_terminator_swing_k`, `temp_night_side_swing_k`). All new columns are NULL here: they are worked out
from the stored temperature, spin, orbit and pressure when a body is generated, so a regenerated galaxy fills them. `planets` and `moons` carry a FULLTEXT key, which refuses
ALGORITHM=INSTANT, so each table is rebuilt once, with every change in one statement. A column already added or already dropped is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0196"
down_revision = "0195"
branch_labels = None
depends_on = None

NEW_COLUMNS = (('band1_mean_k', 'DOUBLE NULL'), ('band1_swing_k', 'DOUBLE NULL'), ('band1_day_hours', 'DOUBLE NULL'), ('band1_night_hours', 'DOUBLE NULL'), ('band1_longest_day_days', 'DOUBLE NULL'), ('band1_longest_night_days', 'DOUBLE NULL'), ('band1_state', 'VARCHAR(24) NULL'), ('band2_mean_k', 'DOUBLE NULL'), ('band2_swing_k', 'DOUBLE NULL'), ('band2_day_hours', 'DOUBLE NULL'), ('band2_night_hours', 'DOUBLE NULL'), ('band2_longest_day_days', 'DOUBLE NULL'), ('band2_longest_night_days', 'DOUBLE NULL'), ('band2_state', 'VARCHAR(24) NULL'), ('band3_mean_k', 'DOUBLE NULL'), ('band3_swing_k', 'DOUBLE NULL'), ('band3_day_hours', 'DOUBLE NULL'), ('band3_night_hours', 'DOUBLE NULL'), ('band3_longest_day_days', 'DOUBLE NULL'), ('band3_longest_night_days', 'DOUBLE NULL'), ('band3_state', 'VARCHAR(24) NULL'), ('band4_mean_k', 'DOUBLE NULL'), ('band4_swing_k', 'DOUBLE NULL'), ('band4_day_hours', 'DOUBLE NULL'), ('band4_night_hours', 'DOUBLE NULL'), ('band4_longest_day_days', 'DOUBLE NULL'), ('band4_longest_night_days', 'DOUBLE NULL'), ('band4_state', 'VARCHAR(24) NULL'), ('band5_mean_k', 'DOUBLE NULL'), ('band5_swing_k', 'DOUBLE NULL'), ('band5_day_hours', 'DOUBLE NULL'), ('band5_night_hours', 'DOUBLE NULL'), ('band5_longest_day_days', 'DOUBLE NULL'), ('band5_longest_night_days', 'DOUBLE NULL'), ('band5_state', 'VARCHAR(24) NULL'), ('band6_mean_k', 'DOUBLE NULL'), ('band6_swing_k', 'DOUBLE NULL'), ('band6_day_hours', 'DOUBLE NULL'), ('band6_night_hours', 'DOUBLE NULL'), ('band6_longest_day_days', 'DOUBLE NULL'), ('band6_longest_night_days', 'DOUBLE NULL'), ('band6_state', 'VARCHAR(24) NULL'), ('band7_mean_k', 'DOUBLE NULL'), ('band7_swing_k', 'DOUBLE NULL'), ('band7_day_hours', 'DOUBLE NULL'), ('band7_night_hours', 'DOUBLE NULL'), ('band7_longest_day_days', 'DOUBLE NULL'), ('band7_longest_night_days', 'DOUBLE NULL'), ('band7_state', 'VARCHAR(24) NULL'))
"""tuple: The columns added to `planets` and `moons`, each with its definition."""

OLD_COLUMNS = ('temp_polar_k', 'temp_north_k', 'temp_equator_k', 'temp_south_k', 'temp_swing_k', 'temp_day_side_k', 'temp_terminator_k', 'temp_night_side_k', 'temp_polar_swing_k', 'temp_north_swing_k', 'temp_south_swing_k', 'temp_day_side_swing_k', 'temp_terminator_swing_k', 'temp_night_side_swing_k')
"""tuple: The columns dropped from `planets` and `moons`."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        changes = [f"ADD COLUMN {column} {definition}" for column, definition in NEW_COLUMNS
                   if not helpers.exists_column(connection, table, column)]
        changes += [f"DROP COLUMN {column}" for column in OLD_COLUMNS if helpers.exists_column(connection, table, column)]
        if changes:
            connection.execute(sa.text(f"ALTER TABLE {table} " + ", ".join(changes)))


def downgrade():
    raise NotImplementedError("the migrations of this database only go forward")
