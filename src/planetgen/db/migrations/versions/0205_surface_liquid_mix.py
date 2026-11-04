"""The makeup, salt and second liquid of the surface liquids (GEN.1218b).

Schema v205. `planets` and `moons` gain `surface_liquid_methane_mol`, `surface_liquid_ethane_mol` and `surface_liquid_nitrogen_mol` (the mole fractions of
a methane, ethane or nitrogen sea), `surface_liquid_salt_wt` (the salt of a water sea, as sodium chloride, percent of the solution's mass),
`surface_liquid_2` (the id of a second liquid beside the main one) and `surface_liquid_2_fraction` (the share of the surface it covers):
`physics/surface_liquid_world.py`, docs/design/composition/surface-liquid-composition.md.

All NULL here: they are found with the surface liquid, so a row made before v205 gets them when the galaxy is regenerated. Each column is skipped when it
is already there, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0205"
down_revision = "0204"
branch_labels = None
depends_on = None

COLUMNS = (("surface_liquid_methane_mol", "DOUBLE"), ("surface_liquid_ethane_mol", "DOUBLE"), ("surface_liquid_nitrogen_mol", "DOUBLE"),
           ("surface_liquid_salt_wt", "DOUBLE"), ("surface_liquid_2", "VARCHAR(32)"), ("surface_liquid_2_fraction", "DOUBLE"))
"""tuple: The columns added to `planets` and `moons`, each with its type."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        after = "surface_liquid_ammonia_wt"
        for column, kind in COLUMNS:
            if not helpers.exists_column(connection, table, column):
                op.execute(f"ALTER TABLE {table} ADD COLUMN {column} {kind} NULL AFTER {after}")
            after = column
