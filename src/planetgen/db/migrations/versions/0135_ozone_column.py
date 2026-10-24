"""The ozone column of a planet or moon, in Dobson units (GEN.413).

Schema v135. `planets` and `moons` gain `ozone_column_du` (DOUBLE, NULL): the ozone column drawn from the body's oxygen and its star's temperature
(`physics/ozone.py`), 0 for a world with no oxygen, NULL for a gas giant. The UV index and the ozone-loss flag read it.

A row made before this revision has no column to read and gets NULL (nothing was drawn for it); the regenerated galaxy replaces it with a drawn value.
A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0135"
down_revision = "0134"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        if not helpers.exists_column(connection, table, "ozone_column_du"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN ozone_column_du DOUBLE NULL"))
