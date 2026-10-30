"""The largest radius and mass of a scatter member (GEN.1009).

Schema v189. `scatter_members` gains `radius_max_au` and `mass_max_solar` (DOUBLE, NULL): the largest radius and mass the member has had
(a white dwarf's red-giant tip and birth mass), from the member record that the tree and the build both read, so the history rule has what
it needs. NULL for a row made before v189; nothing is backfilled, a row gets them when the galaxy is regenerated. A column already there is
skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0189"
down_revision = "0188"
branch_labels = None
depends_on = None

COLUMNS = ("radius_max_au", "mass_max_solar")


def upgrade():
    connection = op.get_bind()
    for column in COLUMNS:
        if not helpers.exists_column(connection, "scatter_members", column):
            connection.execute(sa.text(f"ALTER TABLE scatter_members ADD COLUMN {column} DOUBLE NULL"))
