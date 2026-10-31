"""The nuclear arrival age of a pre-main-sequence star (GEN.1174).

Schema v198. `stars` and `bright_stars` gain `nuclear_arrival_age_gy` (DOUBLE, NULL: the age in Gyr at which MIST's tracks put the star on the
main sequence). A star is pre-main-sequence until its radius has settled (BHAC15), the earlier age `bright_stars.phase_end_age_gy` keeps, and
this column is the later age hydrogen burning arrives, so both times are stored. NULL for every star that is not pre-main-sequence and for a row
made before v198; nothing is backfilled, a row gets its value when the galaxy is regenerated. A column already there is skipped, so a stopped
run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0198"
down_revision = "0198"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    for table in ("stars", "bright_stars"):
        if not helpers.exists_column(connection, table, "nuclear_arrival_age_gy"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN nuclear_arrival_age_gy DOUBLE NULL"))
