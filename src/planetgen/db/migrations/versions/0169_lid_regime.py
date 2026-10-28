"""The tectonic regime of a rocky body, mobile or stagnant lid (GEN.721).

Schema v169. `planets` and `moons` gain `lid_regime` (VARCHAR(10), NULL: `mobile` for plate tectonics, `stagnant` for one rigid plate). It
was three separate coin flips (the magnetic field, the ground dose and, with GEN.721, the heat) and is now one draw kept with the body, which
the interior's heat flux law, the dynamo's lifetime and the ground dose's felsic crust all read. NULL for a gas giant and for a row made before
v169, whose heat flux was a stagnant lid's; nothing is backfilled. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0169"
down_revision = "0168"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        if not helpers.exists_column(connection, table, "lid_regime"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN lid_regime VARCHAR(10) NULL"))
