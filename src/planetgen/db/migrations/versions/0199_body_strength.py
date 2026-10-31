"""The roundness ratio, the yield strength and the binding energy of a body (GEN.1159).

Schema v199. `planets` and `moons` gain `roundness_ratio` (DOUBLE, the mean radius over the round radius; 1 or more is round), `yield_strength_pa`
(DOUBLE, Pa, the whole-body strength of the material mix, drawn from the range of each material, GEN.1178) and `binding_energy_j` (DOUBLE,
J, the energy that assembling the body releases), worked out by `physics/body_strength.py` when a body is generated. A gas giant keeps only the binding
energy. `comets` gain those three with `density_g_cm3` and `ice_mass_fraction` (DOUBLE) from the makeup of the nucleus. `asteroid_belts` and `asteroid_fields` gain `density_g_cm3` and `ice_mass_fraction`, the mean of their members from their makeup. All are NULL here: a row made
before v199 gets them when the galaxy is regenerated. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0199"
down_revision = "0198"
branch_labels = None
depends_on = None

COLUMNS = (("roundness_ratio", "DOUBLE NULL"), ("yield_strength_pa", "DOUBLE NULL"), ("binding_energy_j", "DOUBLE NULL"))
"""tuple: The columns added to `planets` and `moons`, each with its definition."""

COMET_COLUMNS = (("density_g_cm3", "DOUBLE NULL"), ("ice_mass_fraction", "DOUBLE NULL"), *COLUMNS)
"""tuple: The columns added to `comets`."""

POPULATION_COLUMNS = (("density_g_cm3", "DOUBLE NULL"), ("ice_mass_fraction", "DOUBLE NULL"))
"""tuple: The columns added to `asteroid_belts` and `asteroid_fields`."""


def upgrade():
    connection = op.get_bind()
    for table, columns in (("planets", COLUMNS), ("moons", COLUMNS), ("comets", COMET_COLUMNS),
                           ("asteroid_belts", POPULATION_COLUMNS), ("asteroid_fields", POPULATION_COLUMNS)):
        for column, definition in columns:
            if not helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} {definition}"))
