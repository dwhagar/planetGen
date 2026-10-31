"""The liquid a world holds on its surface (GEN.1035).

Schema v188. `planets` and `moons` gain `surface_liquid` (the id of the liquid in `physics/surface_liquids.py`), `surface_liquid_fraction` (the share of
the surface it covers) and `surface_liquid_depth_km` (its mean depth over that share): `physics/surface_liquid_world.py`.

All NULL here: the liquid is found with the ground, so a row made before v188 gets it when the galaxy is regenerated. Each column is skipped when it is
already there, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0188"
down_revision = "0187"
branch_labels = None
depends_on = None

COLUMNS = (("surface_liquid", "VARCHAR(32)"), ("surface_liquid_fraction", "DOUBLE"), ("surface_liquid_depth_km", "DOUBLE"))
"""tuple: The columns added to `planets` and `moons`, each with its type."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for column, kind in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                op.execute(f"ALTER TABLE {table} ADD COLUMN {column} {kind} NULL")
