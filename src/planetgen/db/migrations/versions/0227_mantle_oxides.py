"""The seven minor oxides of a rocky body's rock (GEN.1362).

Schema v227. `planets`, `moons` and `rogue_planets` gain `rock_na2o_fraction`, `rock_k2o_fraction`, `rock_tio2_fraction`, `rock_cr2o3_fraction`, `rock_mno_fraction`, `rock_nio_fraction` and
`rock_p2o5_fraction` (DOUBLE, nullable): each oxide's mass fraction of the rock, from the host star's Na, K, Ti, Cr, Mn, Ni and P, the volatility the body's formation zone kept and the share the core took
(`mantle_oxides`). A body saved before v227 has them NULL and keeps the rock it has; nothing is backfilled. Each column is skipped when it is already there, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0227"
down_revision = "0226"
branch_labels = None
depends_on = None

TABLES = ("planets", "moons", "rogue_planets")
COLUMNS = ("rock_na2o_fraction", "rock_k2o_fraction", "rock_tio2_fraction", "rock_cr2o3_fraction", "rock_mno_fraction", "rock_nio_fraction", "rock_p2o5_fraction")


def upgrade():
    connection = op.get_bind()
    for table in TABLES:
        missing = [name for name in COLUMNS if not helpers.exists_column(connection, table, name)]
        if missing:
            # No ALGORITHM clause: `planets` and `moons` carry a FULLTEXT key, which refuses ALGORITHM=INSTANT; one rebuild adds all of the columns.
            connection.execute(sa.text(f"ALTER TABLE {table} " + ", ".join(f"ADD COLUMN {name} DOUBLE NULL" for name in missing)))
