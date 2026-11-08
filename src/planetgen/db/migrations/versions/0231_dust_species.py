"""The dust grains of a nebula by species (GEN.1376).

Schema v231. `nebulae`, `supernova_remnants` and `nebula_field` split `grain_silicate` into `grain_amorphous_olivine`, `grain_amorphous_pyroxene` and `grain_crystalline_silicate`, split `grain_carbon`
into `grain_amorphous_carbon` and `grain_graphite`, and gain `grain_sic` (silicon carbide). The model also changed what is in the dust (SiC, an ice mantle that holds CO, CO2, methanol, ammonia and methane),
so a table that still has the old columns has its whole composition written again by the model (`nebula_columns`, the fills of migrations 0214 and 0215) and then loses the two old columns. A table
that has the new columns already (a galaxy built from a version before 214 gets them from 0214 and 0215) is left as it is. A stopped run can be started again: the rows still empty are the ones filled, and
the old columns go last.
"""

import importlib

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0231"
down_revision = "0230"
branch_labels = None
depends_on = None

NEW_AFTER = (
    ("grain_amorphous_olivine", "dust_to_metals"), ("grain_amorphous_pyroxene", "grain_amorphous_olivine"), ("grain_crystalline_silicate", "grain_amorphous_pyroxene"),
    ("grain_amorphous_carbon", "grain_iron"), ("grain_graphite", "grain_amorphous_carbon"), ("grain_sic", "grain_pah"),
)
"""tuple: `(new column, the column it follows)`, so a migrated table lists its columns in the order of `schema.sql`."""
OLD = ("grain_silicate", "grain_carbon")
"""tuple: The columns that go."""
TABLES = {"nebulae": "DOUBLE NULL", "supernova_remnants": "DOUBLE NULL", "nebula_field": "FLOAT NULL"}
"""dict: Each table and the type of its grain columns."""


def upgrade():
    connection = op.get_bind()
    started = {}
    for table, definition in TABLES.items():
        started[table] = not helpers.exists_column(connection, table, "grain_sic")
        for column, after in NEW_AFTER:
            helpers.add_column_instant(connection, table, column, definition, after=after)
        if started[table] and helpers.exists_column(connection, table, "grain_silicate"):
            connection.execute(sa.text(f"UPDATE {table} SET mass_fraction_h = NULL"))     # every row is drawn again by the model
            connection.commit()
    first = importlib.import_module("planetgen.db.migrations.versions.0214_nebula_composition")
    second = importlib.import_module("planetgen.db.migrations.versions.0215_nebula_field_composition")
    from planetgen.physics import nebula_columns
    first._fill(connection)
    second._fill(connection, nebula_columns)
    for table in TABLES:
        for column in OLD:
            if helpers.exists_column(connection, table, column):
                helpers.drop_column_instant(connection, table, column)
