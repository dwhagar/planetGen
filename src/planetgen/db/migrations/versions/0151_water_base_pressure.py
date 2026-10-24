"""The pressure at the bottom of the water layer (GEN.461).

Schema v151. `planets` and `moons` gain `water_base_pressure_pa` (DOUBLE, NULL: the pressure under the water over the covered
area, Pa: liquid, ice or both, from the compressed column of `physics/water_layer.py`; NULL for a gas giant and a body with no
water layer).

NULL for a row made before v151, which was generated with the water layer at a fixed 1000 kg/m3: its `ocean_depth_km`, `ice_shell_km`
and `hp_ice_km` are water-equivalent kilometres, not the compressed ones, and the generator epoch of GEN.461 marks it as made by
the old rules. Nothing is backfilled. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0151"
down_revision = "0150"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        if not helpers.exists_column(connection, table, "water_base_pressure_pa"):
            connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN water_base_pressure_pa DOUBLE NULL"))
