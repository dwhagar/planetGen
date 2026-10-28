"""The surface composition of every planet and moon (GEN.441).

Schema v159. `planets` and `moons` gain `ground_type` and `ground_grain` and the areal share of each ground type the generator now draws
from a body's bulk fractions, temperature, air and water: `ground_ocean_fraction`, `ground_ice_fraction`, `ground_ultramafic_fraction`,
`ground_anorthositic_fraction`, `ground_basaltic_fraction`, `ground_felsic_fraction`, `ground_clay_fraction`, `ground_oxidised_fraction`,
`ground_metal_fraction`, `ground_carbon_fraction`, `ground_evaporite_fraction` and `ground_sulfur_fraction`, which add up to 1.

All NULL here: a row made before v159 has no ground to read back (it needs the draws of the generator epoch, so the regenerated galaxy
replaces it). Each column is skipped when it is already there, so a stopped run can be started again.
"""

from alembic import op

from planetgen.db.migrations import helpers

revision = "0159"
down_revision = "0158"
branch_labels = None
depends_on = None

_COLUMNS = (
    ("ground_type", "VARCHAR(16) NULL"),
    ("ground_grain", "VARCHAR(8) NULL"),
    ("ground_ocean_fraction", "DOUBLE NULL"),
    ("ground_ice_fraction", "DOUBLE NULL"),
    ("ground_ultramafic_fraction", "DOUBLE NULL"),
    ("ground_anorthositic_fraction", "DOUBLE NULL"),
    ("ground_basaltic_fraction", "DOUBLE NULL"),
    ("ground_felsic_fraction", "DOUBLE NULL"),
    ("ground_clay_fraction", "DOUBLE NULL"),
    ("ground_oxidised_fraction", "DOUBLE NULL"),
    ("ground_metal_fraction", "DOUBLE NULL"),
    ("ground_carbon_fraction", "DOUBLE NULL"),
    ("ground_evaporite_fraction", "DOUBLE NULL"),
    ("ground_sulfur_fraction", "DOUBLE NULL"),
)
"""tuple: The `(column, definition)` pairs added to both tables."""


def upgrade():
    connection = op.get_bind()
    for table in ("planets", "moons"):
        for column, definition in _COLUMNS:
            if not helpers.exists_column(connection, table, column):
                op.execute(f"ALTER TABLE {table} ADD COLUMN {column} {definition}")
