"""The phosphorus budget of a rocky body (GEN.1367).

Schema v230. `planets`, `moons` and `rogue_planets` gain `crust_p2o5_wt` (the crust's P2O5, weight percent), `crust_apatite_share`, `crust_merrillite_share` and `crust_schreibersite_share` (the shares of
the crust's phosphorus in each host, adding to 1) and `ocean_phosphate_umol_l` (the ocean's phosphate, umol/L), all DOUBLE and nullable. `phosphorus` stays and is now the band of the number
(`phosphorus_budget.phosphorus_class`). A body saved before v230 has them NULL and keeps the class it has; nothing is backfilled. Each column is skipped when it is already there, so a stopped run
can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0230"
down_revision = "0229"
branch_labels = None
depends_on = None

TABLES = ("planets", "moons", "rogue_planets")
COLUMNS = ("crust_p2o5_wt", "crust_apatite_share", "crust_merrillite_share", "crust_schreibersite_share", "ocean_phosphate_umol_l")


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        missing = [name for name in COLUMNS if not helpers.exists_column(connection, table, name)]
        if missing:
            # No ALGORITHM clause: `planets` and `moons` carry a FULLTEXT key, which refuses ALGORITHM=INSTANT; one rebuild adds all of the columns.
            connection.execute(sa.text(f"ALTER TABLE {table} " + ", ".join(f"ADD COLUMN {name} DOUBLE NULL" for name in missing)))


def downgrade():
    pass
