"""The bulk density and binding energy of a star, a remnant, a rogue planet and an interstellar comet (GEN.1162).

Schema v198. `stars`, `black_holes` and `neutron_stars` gain `density_g_cm3` (DOUBLE, g/cm3) and `binding_energy_j` (DOUBLE, J; -1 for the infinite binding energy of a black hole), worked out by
`physics/body_strength.py` from the object's mass and radius. `rogue_planets` gain those with `roundness_ratio` and `yield_strength_pa` (a terrestrial one; NULL for a giant or a brown dwarf).
`interstellar_comets` gain what `comets` have: `density_g_cm3`, `ice_mass_fraction`, `roundness_ratio`, `yield_strength_pa` and `binding_energy_j`. All are NULL here: a row made before v198 gets
them when the galaxy is regenerated. A column already there is skipped, so a stopped run can be started again.
"""

import sqlalchemy as sa
from alembic import op

from planetgen.db.migrations import helpers

revision = "0198"
down_revision = "0197"
branch_labels = None
depends_on = None

STELLAR_COLUMNS = (("density_g_cm3", "DOUBLE NULL"), ("binding_energy_j", "DOUBLE NULL"))
"""tuple: The columns added to `stars`, `black_holes` and `neutron_stars`, each with its definition."""

ROGUE_COLUMNS = (("density_g_cm3", "DOUBLE NULL"), ("roundness_ratio", "DOUBLE NULL"), ("yield_strength_pa", "DOUBLE NULL"), ("binding_energy_j", "DOUBLE NULL"))
"""tuple: The columns added to `rogue_planets`."""

COMET_COLUMNS = (("density_g_cm3", "DOUBLE NULL"), ("ice_mass_fraction", "DOUBLE NULL"), ("roundness_ratio", "DOUBLE NULL"),
                 ("yield_strength_pa", "DOUBLE NULL"), ("binding_energy_j", "DOUBLE NULL"))
"""tuple: The columns added to `interstellar_comets`."""


def upgrade():
    connection = op.get_bind()
    for table, columns in (("stars", STELLAR_COLUMNS), ("black_holes", STELLAR_COLUMNS), ("neutron_stars", STELLAR_COLUMNS),
                           ("rogue_planets", ROGUE_COLUMNS), ("interstellar_comets", COMET_COLUMNS)):
        for column, definition in columns:
            if not helpers.exists_column(connection, table, column):
                connection.execute(sa.text(f"ALTER TABLE {table} ADD COLUMN {column} {definition}"))

