"""What the surface liquids hold (GEN.1215).

Schema v197. `planets` and `moons` gain `surface_liquid_acid_wt` (the strength of the liquid sulfuric acid the air makes, percent of its mass),
`surface_liquid_acid_layer` (where it is, `pool` or `cloud`) and `surface_liquid_ammonia_wt` (the ammonia of a liquid water ocean, percent of the
liquid's mass): `physics/surface_liquid_world.py`, docs/design/composition/surface-liquid-composition.md.

All NULL here: they are found with the surface liquid, so a row made before v197 gets them when the galaxy is regenerated. Each column is skipped when it
is already there, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0197"
down_revision = "0196"
branch_labels = None
depends_on = None

COLUMNS = (("surface_liquid_acid_wt", "DOUBLE"), ("surface_liquid_acid_layer", "VARCHAR(8)"), ("surface_liquid_ammonia_wt", "DOUBLE"))
"""tuple: The columns added to `planets` and `moons`, each with its type."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        after = "surface_liquid_depth_km"
        for column, kind in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                op.execute(f"ALTER TABLE {table} ADD COLUMN {column} {kind} NULL AFTER {after}")
            after = column
