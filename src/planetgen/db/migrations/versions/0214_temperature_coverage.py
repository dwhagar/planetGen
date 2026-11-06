"""What is read off the temperature of every place of a body (GEN.1265).

Schema v214. `planets`, `moons` and `rogue_planets` gain the 25 columns `temp_livable_share` to `temp_site_longest_out_days` (`physics/temperature_coverage.py`): the livable share L of the surface's area
and time and its two pieces, the Blue, Green, cold and hot shares, the area livable all the time and at least 90 percent of it, whether the livable ground is only the terminator belt,
the coldest and hottest sample anywhere with where and when, and the best site with its coldest, hottest and mean temperature and its longest spell out of the livable band. The Temperature
score is made from them (`habitability.temperature_score`). All new columns are NULL here: they are worked out from the stored temperature, spin, orbit, air, sea and ice when a body is
generated or brought to the epoch. `planets` and `moons` carry a FULLTEXT key, which refuses ALGORITHM=INSTANT, so each table is rebuilt once, with every change in one statement. A column
already added is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0214"
down_revision = "0213"
branch_labels = None
depends_on = None

NEW_COLUMNS = (
    ('temp_livable_share', 'DOUBLE NULL'),
    ('temp_livable_area', 'DOUBLE NULL'),
    ('temp_livable_time', 'DOUBLE NULL'),
    ('temp_blue_share', 'DOUBLE NULL'),
    ('temp_green_share', 'DOUBLE NULL'),
    ('temp_cold_share', 'DOUBLE NULL'),
    ('temp_hot_share', 'DOUBLE NULL'),
    ('temp_area_always', 'DOUBLE NULL'),
    ('temp_area_most', 'DOUBLE NULL'),
    ('temp_livable_belt_only', 'TINYINT NULL'),
    ('temp_coldest_k', 'DOUBLE NULL'),
    ('temp_coldest_lat_deg', 'DOUBLE NULL'),
    ('temp_coldest_lon_deg', 'DOUBLE NULL'),
    ('temp_coldest_day', 'DOUBLE NULL'),
    ('temp_hottest_k', 'DOUBLE NULL'),
    ('temp_hottest_lat_deg', 'DOUBLE NULL'),
    ('temp_hottest_lon_deg', 'DOUBLE NULL'),
    ('temp_hottest_day', 'DOUBLE NULL'),
    ('temp_site_lat_deg', 'DOUBLE NULL'),
    ('temp_site_lon_deg', 'DOUBLE NULL'),
    ('temp_site_livable_time', 'DOUBLE NULL'),
    ('temp_site_coldest_k', 'DOUBLE NULL'),
    ('temp_site_hottest_k', 'DOUBLE NULL'),
    ('temp_site_mean_k', 'DOUBLE NULL'),
    ('temp_site_longest_out_days', 'DOUBLE NULL'),
)
"""tuple: The columns added to `planets`, `moons` and `rogue_planets` (a rogue planet is a planet since GEN.1111), each with its definition."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons", "rogue_planets"):
        changes = [f"ADD COLUMN {column} {definition}" for column, definition in NEW_COLUMNS
                   if not helpers.exists_column(connection, table, column)]
        if changes:
            connection.execute(sa.text(f"ALTER TABLE {table} " + ", ".join(changes)))


def downgrade():
    raise NotImplementedError("the migrations of this database only go forward")
