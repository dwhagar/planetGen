"""The dose a moon takes from its parent's radiation belts (GEN.983).

Schema v195. `planets` and `moons` gain `dose_belt_msv_yr` (DOUBLE, NULL: the dose, mSv/yr, a body takes from the trapped particles of its parent
planet's field; `physics/radiation_belts.py`). It is part of `surface_dose_msv_yr`.

All NULL here: the dose is made with the body, so a row made before v195 gets it when the galaxy is regenerated, and meanwhile its surface dose is what
it was. Each column is skipped when it is already there, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0195"
down_revision = "0194"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        if not helpers.exists_column(connection, table, "dose_belt_msv_yr"):
            op.execute(f"ALTER TABLE {table} ADD COLUMN dose_belt_msv_yr DOUBLE NULL")
