"""The heat a rocky body's interior gives, radiogenic plus primordial (GEN.699).

Schema v153. `planets` and `moons` gain `mantle_eta0_pa_s` (DOUBLE, NULL: the viscosity prefactor of the mantle, Pa s),
`mantle_activation_kj_mol` (DOUBLE, NULL: its activation energy, kJ/mol) and `internal_heat_flux_w_m2` (DOUBLE, NULL: the heat
flowing out of the interior over the body's surface, W/m2). `urey_ratio`, stored since v141, keeps its column and is now an output
of the thermal history, not a draw. All NULL for a gas giant and for a row made before v153, which keeps the heat and the Urey ratio
it was generated with (the Urey ratio drawn between 0.1 and 0.9); nothing is backfilled. A column already there is skipped, so a
stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0153"
down_revision = "0152"
branch_labels = None
depends_on = None

COLUMNS = ("mantle_eta0_pa_s", "mantle_activation_kj_mol", "internal_heat_flux_w_m2")


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for column in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} DOUBLE NULL"))
