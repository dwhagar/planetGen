"""The species of the Venus sulfur-cloud set a body's air holds (GEN.1425).

Schema v245. `planets`, `moons` and `rogue_planets` gain `cloud_chemistry` (VARCHAR(48), NULL): the ids of the species the chain gave the air (`V1,V2,V10`;
`physics/sulfur_cloud_set.py`), NULL for none. It is NULL here: the set is read from a finished air, light and ground, so a row made before v245 gets it when the
galaxy is regenerated (a migration does not draw). A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0245"
down_revision = "0244"
branch_labels = None
depends_on = None

TABLES = ("planets", "moons", "rogue_planets")
"""tuple: The tables that gain the column."""


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        if not helpers.exists_column(connection, table, "cloud_chemistry"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN cloud_chemistry VARCHAR(48) NULL AFTER life_traits"))
