"""A planet or moon keeps the aquifer it holds underground (GEN.654, DB.92).

Schema v157. `planets` and `moons` gain `aquifer_depth_km` and `aquifer_thickness_km` (DOUBLE, NULL): the depth to the top of the
liquid-water layer in the pores of the crust of a world with no liquid at its surface (a dry, frozen or hot one), and the layer's
thickness, km (`planetgen.physics.aquifers`). A world with none, a gas giant and a body saved earlier keep NULL: the surface
reservoir the layer is filled from was a draw that is not stored, so an older row is not worked out again; only a world generated
now has the layer, and the habitability score of an older row stays as it was stored.

The columns go in with one ALTER for each table (no ALGORITHM=INSTANT: both carry a FULLTEXT index, which refuses it); a column
already there is left alone, so a stopped run starts again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers
from planetgen.physics import aquifers

revision = "0157"
down_revision = "0156"
branch_labels = None
depends_on = None


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        changes = [f"ADD COLUMN {column} DOUBLE NULL" for column in aquifers.AQUIFER_FIELDS
                   if not helpers.exists_column(connection, table, column)]
        if changes:
            connection.execute(sa.text(f"ALTER TABLE {table} {', '.join(changes)}"))
